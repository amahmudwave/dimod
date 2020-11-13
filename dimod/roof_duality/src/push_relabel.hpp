/**
# MIT License
#
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
================================================================================================
*/
#ifndef PUSH_RELABEL_HPP_INCLUDED
#define PUSH_RELABEL_HPP_INCLUDED

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

#include "helper_data_structures.hpp"

// Maximum flow solver based on Push-Relabel algorithm, using the heuristics of
// Global Relabeling, Gap Relabeling, & Highest Vertex First. This is based on
// Cherkassky, B., Goldberg, A. On Implementing the Push-Relabel Method for the
// Maximum Flow Problem . Algorithmica 19, 390-410 (1997).
// https://doi.org/10.1007/PL00009180.
template <class EdgeType> class PushRelabelSolver {
public:
  using edge_iterator = typename vector<EdgeType>::iterator;
  using capacity_t = typename EdgeType::capacity_type;
  using edge_size_t = size_t;

  // We use preallocated vertex nodes for maintaining the linked list since we
  // know the total number of vertices. Vertex number is redundant but when
  // capacity_t is of a type consuming 8 bytes, structural padding will waste
  // the same amount of memory if it were removed.
  struct vertex_node_t {
    int vertex_number;
    int height;
    capacity_t excess;
    vertex_node_t *next;
    vertex_node_t *prev;
  };

  // Linked lists to keep track of active/inactive vertices at different
  // heights/distances from the sink.
  struct level_t {
    preallocated_linked_list<vertex_node_t> active_vertices;
    preallocated_linked_list<vertex_node_t> inactive_vertices;
  };

  PushRelabelSolver(std::vector<std::vector<EdgeType>> &adjacency_list,
                    int source, int sink);

  capacity_t computeMaximumPreflow();

private:
  // This will always be manually inlined, without relying on the compiler to do
  // so since this function is very small and the algorithm calls this the
  // highest number of times.
  // void push(edge_iterator eit) {
  //     int from_vertex = eit->from_vertex;
  //     int to_vertex = eit->to_vertex;
  //     capacity_t flow = std::min(eit->residual,
  //     _vertices[from_vertex].excess); eit->residual -= flow;
  //     _adjacency_list[to_vertex][eit->reverse_edge_index].residual += flow;
  //     _vertices[from_vertex].excess -= flow;
  //     _vertices[to_vertex].excess += flow;
  // }

  void relabel(int vertex);

  void discharge(int vertex);

  // Relabel the vertices with the current distance from the sink, found by
  // using reverse breadth first search.
  void globalRelabel();

  // When there is no vertex at a particular height/distance from the sink, all
  // the other vertices at higher heights are then disconnected from the sink,
  // so label them as such (set their height to the number of vertices) so they
  // are not processed during the computation of the preflow.
  void gapRelabel(int empty_level_height);

  // Get the iterator pair of edges which need to be processed for a vertex next
  // time it is encountered, we may want to skip the edges which have been
  // saturated already before relabeling of the vertex happened.
  std::pair<edge_iterator, edge_iterator> outEdges(int vertex) {
    return {_adjacency_list[vertex].begin(), _adjacency_list[vertex].end()};
  }

  void printLevels(); // For debuggin purposes.

private:
  int _sink;
  int _source;
  int _num_vertices;
  edge_size_t _num_edges;
  int _max_active_height, _min_active_height, _max_height;
  size_t _num_global_relabels, _num_gap_relabels, _num_pushes, _num_relabels;

  std::vector<level_t> _levels;
  std::vector<vertex_node_t> _vertices;
  vector_based_queue<int> _vertex_queue;
  std::vector<std::vector<EdgeType>> &_adjacency_list;
  std::vector<std::pair<edge_iterator, edge_iterator>> _pending_out_edges;
};

// We assume that the flow graph has the correct residuals assigned to the
// edges. Or namely initial value of residual of an edge is equal to its
// capacity.
template <class EdgeType>
PushRelabelSolver<EdgeType>::PushRelabelSolver(
    std::vector<std::vector<EdgeType>> &adjacency_list, int source, int sink)
    : _adjacency_list(adjacency_list), _source(source), _sink(sink),
      _vertex_queue(vector_based_queue<int>(adjacency_list.size())) {
  _num_global_relabels = 0;
  _num_gap_relabels = 0;
  _num_relabels = 0;
  _num_pushes = 0;
  _num_vertices = _adjacency_list.size();
  _vertices.resize(_num_vertices);
  _levels.resize(_num_vertices);
  _pending_out_edges.resize(_num_vertices);

  for (int v = 0; v < _num_vertices; v++) {
    _pending_out_edges[v] = {_adjacency_list[v].begin(),
                             _adjacency_list[v].end()};
    _vertices[v].vertex_number = v;
    _vertices[v].height = 1;
    _vertices[v].excess = 0;
  }
  _vertices[_source].height = _num_vertices;
  _vertices[_sink].height = 0;

  // Saturate the edges coming out of the source.`
  edge_iterator eit, eit_end;
  for (std::tie(eit, eit_end) = outEdges(_source); eit != eit_end; eit++) {
    capacity_t flow = eit->residual;
    _adjacency_list[eit->to_vertex][eit->reverse_edge_index].residual += flow;
    eit->residual = 0;
    _vertices[eit->to_vertex].excess += flow;
    _num_pushes++;
  }

  _max_height = 1;
  _max_active_height = 0;
  _min_active_height = _num_vertices;

  globalRelabel();
}

// Relabel the vertices with the current distance from the sink, found by
// using reverse breadth first search.
template <class EdgeType> void PushRelabelSolver<EdgeType>::globalRelabel() {
  _num_global_relabels++;
  for (int h = 0; h <= _max_height; h++) {
    _levels[h].active_vertices.clear();
    _levels[h].inactive_vertices.clear();
  }
  _max_height = 0;
  _max_active_height = 0;
  _min_active_height = _num_vertices;

  // The value of height being the number of vertices means it is not reachable
  // from the sink as the value of height can be at most one less than the
  // number of vertices.
  for (int i = 0; i < _num_vertices; i++) {
    _vertices[i].height = _num_vertices;
  }

  _vertices[_sink].height = 0;
  _vertex_queue.reset(); // Reset as we are reusing the queue.
  _vertex_queue.push(_sink);
  while (!_vertex_queue.empty()) {
    int v_parent = _vertex_queue.pop();
    int current_height = _vertices[v_parent].height + 1;
    edge_iterator eit, eit_end;
    for (std::tie(eit, eit_end) = outEdges(v_parent); eit != eit_end; eit++) {
      int to_vertex = eit->to_vertex;
      if (eit->getReverseEdgeResidual() &&
          _vertices[to_vertex].height == _num_vertices) {
        _vertices[to_vertex].height = current_height;
        if (_vertices[to_vertex].excess > 0) {
          _levels[current_height].active_vertices.push_front(
              &_vertices[to_vertex]);
        } else {
          _levels[current_height].inactive_vertices.push_front(
              &_vertices[to_vertex]);
        }
        _vertex_queue.push(to_vertex);
      }
    }

    if (!_levels[current_height].active_vertices.empty()) {
      _max_height = std::max(current_height, _max_height);
      _max_active_height = std::max(current_height, _max_active_height);
      _min_active_height = std::min(current_height, _min_active_height);
    } else if (!_levels[current_height].inactive_vertices.empty()) {
      _max_height = std::max(current_height, _max_height);
    }
  }
}

template <class EdgeType>
void PushRelabelSolver<EdgeType>::relabel(int vertex) {
  int min_relabel_height = _num_vertices;
  _vertices[vertex].height = min_relabel_height;
  edge_iterator eit, eit_end, eit_min_relabel;
  for (std::tie(eit, eit_end) = outEdges(vertex); eit != eit_end; eit++) {
    if (eit->residual) {
      int to_vertex = eit->to_vertex;
      if (_vertices[to_vertex].height < min_relabel_height) {
        min_relabel_height = _vertices[to_vertex].height;
        eit_min_relabel = eit;
      }
    }
  }
  min_relabel_height++;
  if (min_relabel_height < _num_vertices) {
    _vertices[vertex].height = min_relabel_height;
    _pending_out_edges[vertex].first = eit_min_relabel;
    _max_height = std::max(min_relabel_height, _max_height);
  }
}

template <class EdgeType>
void PushRelabelSolver<EdgeType>::discharge(int vertex) {
  assert(_vertices[vertex].excess > 0);
  while (1) {
    edge_iterator eit, eit_end;
    int pushable_height = _vertices[vertex].height - 1;
    for (std::tie(eit, eit_end) = _pending_out_edges[vertex]; eit != eit_end;
         eit++) {
      if (eit->residual) {
        int to_vertex = eit->to_vertex;
        if (_vertices[to_vertex].height == pushable_height) {
          if (to_vertex != _sink && _vertices[to_vertex].excess == 0) {
            _levels[pushable_height].inactive_vertices.erase(
                &_vertices[to_vertex]);
            _levels[pushable_height].active_vertices.push_front(
                &_vertices[to_vertex]);
          }
          // Push flow inlined here: push(eit);
          capacity_t flow = std::min(eit->residual, _vertices[vertex].excess);
          eit->residual -= flow;
          _adjacency_list[to_vertex][eit->reverse_edge_index].residual += flow;
          _vertices[vertex].excess -= flow;
          _vertices[to_vertex].excess += flow;
          if (_vertices[vertex].excess == 0) {
            break;
          }
        }
      }
    }

    if (!_levels[pushable_height].active_vertices.empty()) {
      _max_active_height = std::max(pushable_height, _max_active_height);
      _min_active_height = std::min(pushable_height, _min_active_height);
    }

    // The loop did not break thus the vertex still has some excess flow to
    // discharge so relabel.
    if (eit == eit_end) {
      int preRelabelHeight = _vertices[vertex].height;
      relabel(vertex);
      if (_levels[preRelabelHeight].active_vertices.empty() &&
          _levels[preRelabelHeight].inactive_vertices.empty()) {
        gapRelabel(preRelabelHeight);
      }
      if (_vertices[vertex].height == _num_vertices) {
        break;
      }
    } else {
      _pending_out_edges[vertex].first = eit;
      _levels[_vertices[vertex].height].inactive_vertices.push_front(
          &_vertices[vertex]);
      break;
    }
  }
  assert(_max_height >= 0 && _max_height < _num_vertices);
  assert(_max_active_height >= 0 && _max_active_height < _num_vertices);
}

// When there is no vertex at a particular height/distance from sink, all the
// other vertices at higher heights are then disconnected from the sink, so
// label them as such (set their height to number of vertices) so they are not
// processed during the computation of the preflow.
template <class EdgeType>
void PushRelabelSolver<EdgeType>::gapRelabel(int empty_level_height) {
  for (auto levelIt = _levels.begin() + empty_level_height + 1;
       levelIt <= _levels.begin() + _max_height; levelIt++) {
    // The last highest active vertex processed was from the empty_level_height.
    assert(levelIt->active_vertices.empty());
    int inactive_level_size = levelIt->inactive_vertices.size();
    vertex_node_t *vertex_node_ptr = levelIt->inactive_vertices.front();
    for (int i = 0; i < inactive_level_size; i++) {
      _vertices[vertex_node_ptr->vertex_number].height = _num_vertices;
      vertex_node_ptr = vertex_node_ptr->next;
    }
    levelIt->inactive_vertices.clear();
  }
  _max_height = empty_level_height - 1;
  _max_active_height = empty_level_height - 1;
  assert(_max_height >= 0 && _max_height < _num_vertices);
  assert(_max_active_height >= 0 && _max_active_height < _num_vertices);
}

template <class EdgeType>
typename EdgeType::capacity_type
PushRelabelSolver<EdgeType>::computeMaximumPreflow() {
  while (_max_active_height >= _min_active_height) {
    // std::cout<< " Max active height " << _max_active_height <<" Min active
    // height " << _min_active_height << std::endl;
    if (_levels[_max_active_height].active_vertices.empty()) {
      _max_active_height--;
    } else {
      vertex_node_t *vertex_node_ptr =
          _levels[_max_active_height].active_vertices.pop();
      // std::cout << " Going to discharge " << vertex_node_ptr->vertex_number
      // << " at height " << _vertices[vertex_node_ptr->vertex_number].height <<
      // std::endl;
      discharge(vertex_node_ptr->vertex_number);
    }
  }
  return _vertices[_sink].excess;
}

template <class EdgeType> void PushRelabelSolver<EdgeType>::printLevels() {
  std::cout << "Printing Levels. " << _levels.size() << " in total for "
            << _num_vertices << " vertices." << std::endl;
  for (int i = 0; i < _levels.size(); i++) {
    int size_active = _levels[i].active_vertices.size();
    int size_inactive = _levels[i].inactive_vertices.size();

    if ((size_active == 0) && (size_inactive == 0)) {
      continue;
    }

    std::cout << std::endl;
    std::cout << "Level " << i << std::endl << std::endl;
    std::cout << "Active list : " << size_active << " vertices" << std::endl;

    vertex_node_t *vertex_node_ptr = _levels[i].active_vertices.front();
    for (int n = 0; n < size_active; n++) {
      std::cout << vertex_node_ptr->vertex_number << " ";
      vertex_node_ptr = vertex_node_ptr->next;
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "Inactive list : " << size_inactive << " vertices"
              << std::endl;

    vertex_node_ptr = _levels[i].inactive_vertices.front();
    for (int n = 0; n < size_inactive; n++) {
      std::cout << vertex_node_ptr->vertex_number << " ";
      vertex_node_ptr = vertex_node_ptr->next;
    }
    std::cout << std::endl;
  }
}

#endif // MAX_PUSH_RELABEL_INCLUDED
