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
#ifndef MAX_FLOW_HPP_INCLUDED
#define MAX_FLOW_HPP_INCLUDED

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
  // know the total number of vertices. Vertex id is redundant but when
  // capacity_t is of a type consuming 8 bytes, structural padding will waste
  // the same amount of memory if it were removed.
  struct vertex_node_t {
    int id;
    int height;
    capacity_t excess;
    vertex_node_t *next;
    vertex_node_t *prev;
  };

  // Linked lists to keep track of active/inactive vertices at different
  // heights.
  struct level_t {
    preallocated_linked_list<vertex_node_t> active_vertices;
    preallocated_linked_list<vertex_node_t> inactive_vertices;
  };

  PushRelabelSolver(std::vector<std::vector<EdgeType>> &adjacency_list,
                    int source, int sink)
      : _adjacency_list(adjacency_list), _source(source), _sink(sink),
        _vertex_queue(vector_based_queue<int>(adjacency_list.size())) {
    _num_global_relabels = 0;
    _num_relabels = 0;
    _num_pushes = 0;
    _num_vertices = _adjacency_list.size();
    _vertices.resize(_num_vertices);
    _levels.resize(_num_vertices);
    _pending_out_edges.resize(_num_vertices);

    for (int v = 0; v < _num_vertices; v++) {
      _pending_out_edges[v] = {_adjacency_list[v].begin(),
                               _adjacency_list[v].end()};
      _vertices[v].id = v;
      _vertices[v].height = 1;
      _vertices[v].excess = 0;
    }
    _vertices[_source].height = _num_vertices;
    _vertices[_sink].height = 0;

    edge_iterator it, itEnd;
    for (std::tie(it, itEnd) = outEdges(_source); it != itEnd; it++) {
      capacity_t flow = it->residual;
      _adjacency_list[it->to_vertex][it->reverse_edge_index].residual += flow;
      it->residual = 0;
      _vertices[it->to_vertex].excess += flow;
      _num_pushes++;
    }

    // This is the maximum height where there can be vertices in the linked
    // lists, anything beyond this height is disconnected from the sink and does
    // not have any path to reach it.
    _max_height = _num_vertices - 1;
    _max_active_height = 0;
    _min_active_height = _num_vertices;

    global_relabel();
  }

  // Relabel the vertices with the current distance from the sink, found by
  // using reverse breadth first search.
  void global_relabel() {
    _num_global_relabels++;
    for (int h = 0; h <= _max_height; h++) {
      _levels[h].active_vertices.clear();
      _levels[h].inactive_vertices.clear();
    }
    _max_height = 0;
    _max_active_height = 0;
    _min_active_height = _num_vertices;

    for (int i = 0; i < _num_vertices; i++) {
      _vertices[i].height = _num_vertices;
    }

    _vertices[_sink].height = 0;
    _vertex_queue.reset();
    _vertex_queue.push(_sink);
    while (!_vertex_queue.empty()) {
      int v_parent = _vertex_queue.pop();
      int current_height = _vertices[v_parent].height + 1;
      edge_iterator it, itEnd;
      for (std::tie(it, itEnd) = outEdges(v_parent); it != itEnd; it++) {
        int to_vertex = it->to_vertex;
        // std::cout << " Parent " << v_parent << " " << current_height -1 <<
        // std::endl;
        if (it->getReverseEdgeResidual() &&
            _vertices[to_vertex].height == _num_vertices) {
          // std::cout << " Child " << to_vertex << " " << current_height  <<
          // std::endl;
          _vertices[to_vertex].height = current_height;
          _max_height = std::max(_max_height, current_height);
          if (_vertices[to_vertex].excess > 0) {
            add_to_active_list(to_vertex);
          } else {
            add_to_inactive_list(to_vertex);
          }
          _vertex_queue.push(to_vertex);
        }
      }
    }
  }

  void discharge(int vertex) {
    assert(_vertices[vertex].excess > 0);
    while (1) {
      edge_iterator eit, eitEnd;
      int vertex_height = _vertices[vertex].height;
      for (std::tie(eit, eitEnd) = _pending_out_edges[vertex]; eit != eitEnd;
           eit++) {
        if (eit->residual) {
          int to_vertex = eit->to_vertex;
          int to_vertex_height = _vertices[to_vertex].height;
          if (vertex_height == to_vertex_height + 1) {
            if (to_vertex != _sink && _vertices[to_vertex].excess == 0) {
              // remove_from_inactive_list(to_vertex);
              _levels[to_vertex_height].inactive_vertices.erase(
                  &_vertices[to_vertex]);
              // add_to_active_list(to_vertex);
              _levels[to_vertex_height].active_vertices.push_front(
                  &_vertices[to_vertex]);
              _max_active_height =
                  std::max(to_vertex_height, _max_active_height);
              _min_active_height =
                  std::min(to_vertex_height, _min_active_height);
              assert(_max_active_height >= 0 &&
                     _max_active_height < _num_vertices);
              assert(_min_active_height >= 0 &&
                     _min_active_height < _num_vertices);
            }
            // Push flow inlined here
            capacity_t flow = std::min(eit->residual, _vertices[vertex].excess);
            eit->residual -= flow;
            _adjacency_list[to_vertex][eit->reverse_edge_index].residual +=
                flow;
            _vertices[vertex].excess -= flow;
            _vertices[to_vertex].excess += flow;
            if (_vertices[vertex].excess == 0)
              break;
          }
        }
      }

      if (eit == eitEnd) {
        int preRelabelHeight = vertex_height;
        relabel(vertex);
        if (_levels[preRelabelHeight].active_vertices.empty() &&
            _levels[preRelabelHeight].inactive_vertices.empty()) {
          gap_relabel(preRelabelHeight);
        }
        if (_vertices[vertex].height == _num_vertices)
          break;
      } else {

        _pending_out_edges[vertex].first = eit;
        add_to_inactive_list(vertex);
        break;
      }
    }
  }

  void gap_relabel(int empty_level_height) {
    for (auto levelIt = _levels.begin() + empty_level_height + 1;
         levelIt < _levels.begin() + _max_height; levelIt++) {
      assert(levelIt->active_vertices.empty());
      int inactive_level_size = levelIt->inactive_vertices.size();
      vertex_node_t *vertex_node_ptr = levelIt->inactive_vertices.front();
      for (int i = 0; i < inactive_level_size; i++) {
        _vertices[vertex_node_ptr->id].height = _num_vertices;
        vertex_node_ptr = vertex_node_ptr->next;
      }
      levelIt->inactive_vertices.clear();
    }
    _max_height = empty_level_height - 1;
    _max_active_height = empty_level_height - 1;
    assert(_max_height >= 0 && _max_height < _num_vertices);
    assert(_max_active_height >= 0 && _max_active_height < _num_vertices);
  }

  void relabel(int vertex) {
    int min_relabel_height = _num_vertices;
    _vertices[vertex].height = min_relabel_height;
    edge_iterator eit, eitEnd, eit_min_relabel;
    for (std::tie(eit, eitEnd) = outEdges(vertex); eit != eitEnd; eit++) {
      int to_vertex = eit->to_vertex;
      if (eit->residual && _vertices[to_vertex].height < min_relabel_height) {
        min_relabel_height = _vertices[to_vertex].height;
        eit_min_relabel = eit;
      }
    }
    min_relabel_height++;
    if (min_relabel_height < _num_vertices) {
      _vertices[vertex].height = min_relabel_height;
      _pending_out_edges[vertex].first = eit_min_relabel;
      _max_height = std::max(_max_height, min_relabel_height);
    }
  }

  capacity_t maximum_preflow() {
    while (_max_active_height >= _min_active_height) {
      // std::cout<< " Max active height " << _max_active_height <<" Min active
      // height " << _min_active_height << std::endl;
      if (_levels[_max_active_height].active_vertices.empty()) {
        _max_active_height--;
      } else {
        vertex_node_t *vertex_node_ptr =
            _levels[_max_active_height].active_vertices.pop();
        // std::cout << " Going to discharge " << vertex_node_ptr->id << " at
        // height " << _vertices[vertex_node_ptr->id].height << std::endl;
        discharge(vertex_node_ptr->id);
      }
    }
    return _vertices[_sink].excess;
  }

  void push(int from_vertex, edge_iterator eit) {
    int to_vertex = eit->to_vertex;
    capacity_t flow = std::min(eit->residual, _vertices[from_vertex].excess);
    eit->residual -= flow;
    _adjacency_list[to_vertex][eit->reverse_edge_index].residual += flow;
    _vertices[from_vertex].excess -= flow;
    _vertices[to_vertex].excess += flow;
  }

  std::pair<edge_iterator, edge_iterator> outEdges(int vertex) {
    return {_adjacency_list[vertex].begin(), _adjacency_list[vertex].end()};
  }

  void add_to_active_list(int vertex) {
    int height = _vertices[vertex].height;
    _levels[height].active_vertices.push_front(&_vertices[vertex]);
    _max_active_height = std::max(height, _max_active_height);
    _min_active_height = std::min(height, _min_active_height);
  }

  void remove_from_active_list(int vertex) {
    _levels[_vertices[vertex].height].active_vertices.erase(&_vertices[vertex]);
  }

  void remove_from_inactive_list(int vertex) {
    _levels[_vertices[vertex].height].inactive_vertices.erase(
        &_vertices[vertex]);
  }

  void add_to_inactive_list(int vertex) {
    _levels[_vertices[vertex].height].inactive_vertices.push_front(
        &_vertices[vertex]);
  }

  void printLevels() {
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
        std::cout << vertex_node_ptr->id << " ";
        vertex_node_ptr = vertex_node_ptr->next;
      }
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << "Inactive list : " << size_inactive << " vertices"
                << std::endl;
      vertex_node_ptr = _levels[i].inactive_vertices.front();

      for (int n = 0; n < size_inactive; n++) {
        std::cout << vertex_node_ptr->id << " ";
        vertex_node_ptr = vertex_node_ptr->next;
      }
      std::cout << std::endl;
    }
  }

  int _source;
  int _sink;
  int _num_vertices;
  edge_size_t _num_edges;
  int _max_active_height, _min_active_height, _max_height;
  size_t _num_global_relabels, _num_pushes, _num_relabels;

  std::vector<level_t> _levels;
  std::vector<vertex_node_t> _vertices;
  vector_based_queue<int> _vertex_queue;
  std::vector<std::vector<EdgeType>> &_adjacency_list;
  std::vector<std::pair<edge_iterator, edge_iterator>> _pending_out_edges;
};

#endif // MAX_FLOW_HPP_INCLUDED
