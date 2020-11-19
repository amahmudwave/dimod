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

//#define COLLECT_STATISTICS

#ifdef COLLECT_STATISTICS
#define DEBUG_INCREMENT(x) ((x)++)
#else
#define DEBUG_INCREMENT(x)
#endif

// Maximum flow solver based on Push-Relabel algorithm, using the heuristics of
// Global Relabeling, Gap Relabeling, & Highest Vertex First. This is based on
// Cherkassky, B., Goldberg, A. On Implementing the Push-Relabel Method for the
// Maximum Flow Problem . Algorithmica 19, 390-410 (1997).
// https://doi.org/10.1007/PL00009180.
template <class EdgeType> class PushRelabelSolver {
public:
  using edge_iterator = typename vector<EdgeType>::iterator;
  using capacity_t = typename EdgeType::capacity_type;

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

  void convertPreflowToFlow(bool handle_self_loops = false);

  void printStatistics();

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

  void printLevels();

private:
  int _sink;
  int _source;
  int _num_vertices;
  size_t _num_edges;
  int _max_active_height, _min_active_height, _max_height;
  size_t _num_global_relabels, _num_gap_relabels, _num_gap_vertices,
      _num_pushes, _num_relabels;

  // Thresholds and counters for triggering global relabeling.
  size_t _relabel_work;
  size_t GLOBAL_RELABEL_THRESHOLD;
  static const int ALPHA = 6, BETA = 12;
  static constexpr double GLOBAL_RELABEL_FREQUENCY = 0.5;

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
  _num_gap_vertices = 0;
  _num_relabels = 0;
  _num_pushes = 0;
  _num_vertices = _adjacency_list.size();
  _vertices.resize(_num_vertices);
  _levels.resize(_num_vertices);
  _pending_out_edges.resize(_num_vertices);

  _num_edges = 0;
  for (int v = 0; v < _num_vertices; v++) {
    _pending_out_edges[v] = {_adjacency_list[v].begin(),
                             _adjacency_list[v].end()};
    _vertices[v].vertex_number = v;
    _vertices[v].height = 1;
    _vertices[v].excess = 0;
    _num_edges += _adjacency_list[v].size();
  }
  _vertices[_source].height = _num_vertices;
  _vertices[_sink].height = 0;

  GLOBAL_RELABEL_THRESHOLD = (ALPHA * _num_vertices) + (_num_edges / 2);

  // Saturate the edges coming out of the source.`
  edge_iterator eit, eit_end;
  for (std::tie(eit, eit_end) = outEdges(_source); eit != eit_end; eit++) {
    capacity_t flow = eit->residual;
    _adjacency_list[eit->to_vertex][eit->reverse_edge_index].residual += flow;
    eit->residual = 0;
    _vertices[eit->to_vertex].excess += flow;
    DEBUG_INCREMENT(_num_pushes);
  }

  _max_height = 1;
  _max_active_height = 0;
  _min_active_height = _num_vertices;

  globalRelabel();
}

// Relabel the vertices with the current distance from the sink, found by
// using reverse breadth first search. Assumes that the unreachable vertices are
// marked as such.
template <class EdgeType> void PushRelabelSolver<EdgeType>::globalRelabel() {
  DEBUG_INCREMENT(_num_global_relabels);
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
  int already_unreachable = 0;
  assert((_vertices[_sink] != _num_vertices) &&
         "Sink should have been reachable before globalRelabel");
  for (int i = 0; i < _num_vertices; i++) {
    if (_vertices[i].height == _num_vertices) {
      already_unreachable++;
    } else {
      _vertices[i].height = _num_vertices;
    }
  }

  int num_vertices_to_find = _num_vertices - already_unreachable;
  int num_vertices_found = 0;
  _vertices[_sink].height = 0;
  _vertex_queue.reset(); // Reset as we are reusing the queue.
  _vertex_queue.push(_sink);
  num_vertices_found++; // We had assumed sink is reachable, so reaching it.
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
        num_vertices_found++;
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

    if (num_vertices_found == num_vertices_to_find) {
      _vertex_queue.reset();
      break;
    }
  }
}

template <class EdgeType>
void PushRelabelSolver<EdgeType>::relabel(int vertex) {
  DEBUG_INCREMENT(_num_relabels);
  int min_relabel_height = _num_vertices;
  _vertices[vertex].height = min_relabel_height;
  edge_iterator eit, eit_end, eit_min_relabel;
  std::tie(eit, eit_end) = outEdges(vertex);
  _relabel_work += BETA + std::distance(eit, eit_end);
  for (; eit != eit_end; eit++) {
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
          DEBUG_INCREMENT(_num_pushes);
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
        // We just relabelled all the vertices at its previous height to a
        // height of the number of vertices, so this vertex does not have any
        // residual edge to a vertex that has height/distance one less than its
        // current height, since the level at that height is empty. So we set
        // this one to a height of _num_vertices as well.
        _vertices[vertex].height = _num_vertices;
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
  DEBUG_INCREMENT(_num_gap_relabels);
  assert(empty_level_height <= _max_height);
  auto level_it = _levels.begin() + empty_level_height + 1;
  auto level_it_end = _levels.begin() + _max_height + 1;
  for (; level_it != level_it_end; level_it++) {
    // The last highest active vertex processed was from the empty_level_height.
    assert(level_it->active_vertices.empty());
    int inactive_level_size = level_it->inactive_vertices.size();
    vertex_node_t *vertex_node_ptr = level_it->inactive_vertices.front();
    for (int i = 0; i < inactive_level_size; i++) {
      _vertices[vertex_node_ptr->vertex_number].height = _num_vertices;
      vertex_node_ptr = vertex_node_ptr->next;
      DEBUG_INCREMENT(_num_gap_vertices);
    }
    level_it->inactive_vertices.clear();
  }
  _max_height = empty_level_height - 1;
  _max_active_height = empty_level_height - 1;
  assert(_max_height >= 0 && _max_height < _num_vertices);
  assert(_max_active_height >= 0 && _max_active_height < _num_vertices);
}

template <class EdgeType>
typename EdgeType::capacity_type
PushRelabelSolver<EdgeType>::computeMaximumPreflow() {

  _relabel_work = 0;
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

    if (_relabel_work * GLOBAL_RELABEL_FREQUENCY > GLOBAL_RELABEL_THRESHOLD) {
      globalRelabel();
      _relabel_work = 0;
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

template <class EdgeType> void PushRelabelSolver<EdgeType>::printStatistics() {
#ifdef COLLECT_STATISTICS
  std::cout << std::endl;
  std::cout << "Printing Statistics : " << std::endl << std::endl;
  std::cout << "Number of Pushes : " << _num_pushes << std::endl;
  std::cout << "Number of Relabels : " << _num_relabels << std::endl;
  std::cout << "Number of Global Relabels : " << _num_global_relabels
            << std::endl;
  std::cout << "Number of Gap Relabels : " << _num_gap_relabels << std::endl;
  std::cout << "Number of Gap Vertices : " << _num_gap_vertices << std::endl;
  std::cout << std::endl;
#else
  std::cout << std::endl;
  std::cout << "Statistics not collected." << std::endl << std::endl;
#endif
}

// This is the algorithm used in boost library. An iterative depth first search
// (DFS) is applied to do a topological sort of the vertices with excess flow,
// so that flow can be pushed away from them in order and preflow be converted
// to flow. A cycle in the directed graph would prevent a proper topological
// sort thus in the same search if any flow cycle is detected, the cycle is
// cutoff by reducing the flow by the minimum amount of flow in the cycle, this
// causes the edges with flow equal to the minimum flow to get saturated. The
// depth first search then backs off or in other words, undoes what was done
// from the point of the saturated edge, so that it can go on as before from
// that point. Note there will be no explicit stack used for this depth first
// search, but instead the iterators pointing to the edges to be traversed for
// each vertex will be updated and an array containing parents of vertices will
// help simulate the stack. When we want to pop the stack we can basically look
// at the parent of the current vertex being processed and when we want to push
// a vertex, we can make it a parent for the next vertex while incrementing its
// edge iterator. For solving max-flows on graphs induced from posiforms, we do
// not need to consider self loops since we do induce terms terms like X_i*X_i'
// which might create a self loop for X_i or X_i'.
template <class EdgeType>
void PushRelabelSolver<EdgeType>::convertPreflowToFlow(bool handle_self_loops) {

  // Certain graphs, like the implication networks formed from posiforms are
  // guaranteed not to contain self loops thus processing them is optional. If
  // there are self loops then we can just eliminate the flow, as it does not
  // contribute to the excess of the vertex.
  if (handle_self_loops) {
    for (int from_vertex = 0; from_vertex < _num_vertices; from_vertex++) {
      edge_iterator eit, eit_end;
      for (std::tie(eit, eit_end) = outEdges(from_vertex); eit != eit_end;
           eit++) {
        if (eit->to_vertex == from_vertex) {
          eit->residual = eit->getCapacity();
        }
      }
    }
  }

  int topology_start_vertex = -1;
  bool topology_initialized = false;
  std::vector<int> parent(_num_vertices);

  // An entry of -1 for a vertex would mean there is no successor to that vertex
  // in the topological ordering we create for removing excess flow.
  std::vector<int> topology_next(_num_vertices, -1);

  // White - not processed yet.
  // Grey - under process.
  // Black - finished processing.
  enum DFS_COLOR { WHITE, GREY, BLACK };
  std::vector<int> dfs_color(_num_vertices, DFS_COLOR::WHITE);

  for (int vertex = 0; vertex < _num_vertices; vertex++) {
    // TODO : Remove ? We fillup the parents as we go , why initialize here ?
    parent[vertex] = vertex;
    _pending_out_edges[vertex] = outEdges(vertex);
  }

  // Iterative DFS.
  for (int parent_vertex = 0; parent_vertex < _num_vertices; parent_vertex++) {
    if ((dfs_color[parent_vertex] == DFS_COLOR::WHITE) &&
        (_vertices[parent_vertex].excess > 0) && (parent_vertex != _source) &&
        (parent_vertex != _sink)) {
      int root_vertex = parent_vertex; // The root of a dfs tree.
      dfs_color[root_vertex] = DFS_COLOR::GREY;
      while (true) {
        // We need to increment the actual iterator in the array. Think
        // of recursive dfs, the counter of a for loop belonging to a particular
        // stack frame is saved in the stack frame while the function calls
        // itself recursively so when it comes back it finds the counter with
        // the proper value.
        for (; _pending_out_edges[parent_vertex].first !=
               _pending_out_edges[parent_vertex].second;
             _pending_out_edges[parent_vertex].first++) {
          // We are traversing reverse edges with capacity 0, i.e the children
          // are sending flow to the parent. So we denote the main iterator by
          // eit_reverse unlike other places in the file.
          auto eit_reverse = _pending_out_edges[parent_vertex].first;

          // The edge brings flow to the parent_vertex, we want to topologically
          // sort such that the root receives flow from the source through its
          // children, so that we can push them back.
          if ((eit_reverse->getCapacity() == 0) &&
              (eit_reverse->residual > 0)) {
            int current_vertex = eit_reverse->to_vertex;
            if (dfs_color[current_vertex] == DFS_COLOR::WHITE) {
              dfs_color[current_vertex] = DFS_COLOR::GREY;

              // Equivalent to calling dfs, i.e pushing into the stack. The
              // while loop above will start the for loop all over again, but
              // this time the parent_vertex will be the current_vertex, and we
              // will be able to trace back as the parent-child relationship is
              // saved in the parent array.
              parent[current_vertex] = parent_vertex;
              parent_vertex = current_vertex;
              break;
            }

            // We have detected a cycle, now we need to eliminate the cycle and
            // then back the DFS up to the vertex from which emanates the first
            // edge that is saturated/flow nullified in order to remove the
            // cycle.
            else if (dfs_color[current_vertex] == DFS_COLOR::GREY) {
              // We chose to enter the condition only for reverse/residual edges
              // and also made sure we incremented the array of iterators i.e
              // the _pending_out_edges, thus we can cycle over the cycle found,
              // as the path formed by the iterators saved in _pending_out_edges
              // outlines the cycle we just discovered. The value of residual in
              // a reverse edge is equal to the flow in the orignal edge, we are
              // currently traversing reverse edges, so take the minimum of
              // residuals.
              capacity_t min_flow = eit_reverse->residual;
              int cycle_traversing_vertex = current_vertex;
              while (cycle_traversing_vertex != parent_vertex) {
                min_flow = std::min(
                    _pending_out_edges[cycle_traversing_vertex].first->residual,
                    min_flow);
                cycle_traversing_vertex =
                    _pending_out_edges[cycle_traversing_vertex]
                        .first->to_vertex;
              }

              // Reduce the flow in the edges with minimum flow in the cycle,
              // that is saturate the reverse edge or say zero out the residual.
              // Process the back edge first.
              eit_reverse = _pending_out_edges[parent_vertex].first;
              eit_reverse->residual -= min_flow;
              auto eit = _adjacency_list[eit_reverse->to_vertex].begin() +
                         eit_reverse->reverse_edge_index;
              eit->residual += min_flow;

              // Reduce the flow on the other edges of the cycle, but this time
              // starting from the current_vertex. The vertex from which the
              // first saturated vertex emanates is the one we should restart
              // our DFS from.
              cycle_traversing_vertex = current_vertex;
              int restart_vertex = parent_vertex;
              bool restart_found = false;
              while (cycle_traversing_vertex != parent_vertex) {
                eit_reverse = _pending_out_edges[cycle_traversing_vertex].first;
                eit_reverse->residual -= min_flow;
                auto eit = _adjacency_list[eit_reverse->to_vertex].begin() +
                           eit_reverse->reverse_edge_index;
                eit->residual += min_flow;

                if ((eit_reverse->residual == 0) && !restart_found) {
                  restart_found = true;
                  restart_vertex = cycle_traversing_vertex;
                  dfs_color[eit_reverse->to_vertex] = DFS_COLOR::WHITE;
                } else if (restart_found) {
                  dfs_color[eit_reverse->to_vertex] = DFS_COLOR::WHITE;
                }
                cycle_traversing_vertex = eit_reverse->to_vertex;
              }

              if (restart_vertex != parent_vertex) {
                parent_vertex = restart_vertex;
                _pending_out_edges[parent_vertex].first++;
                break;
              }
            } // else if (dfs_color[current_vertex] == DFS_COLOR::GREY)
          }   // if ((..->getCapacity() == 0) && (..->residual > 0))
        }     // for (; _pending_out_edges[parent_vertex]...

        // The current parent has finished all its edges, we can put it in its
        // proper place in the topological order, and then backtrack.
        if (_pending_out_edges[parent_vertex].first ==
            _pending_out_edges[parent_vertex].second) {
          dfs_color[parent_vertex] = DFS_COLOR::BLACK;
          if (parent_vertex != _source) {
            if (!topology_initialized) {
              topology_start_vertex = parent_vertex;
              topology_initialized = true;
            } else {
              topology_next[parent_vertex] = topology_start_vertex;
              topology_start_vertex = parent_vertex;
            }
          }

          // If it is a root, then we break the while loop and the for loop
          // above will automatically select the next vertex therwise we have to
          // do it here. It is similar to the situation that in a dfs a wrapper
          // function calls the first instance of dfs and before the call there
          // is no stackframe for dfs function. In the case it is not the root
          // we need to simulate the popping of the stack, this we do by setting
          // the parent_vertex to be the parent of the parent_vertex, and
          // incrementing the edge iterator, and the for loop will start again
          // for the parent. Note: in general the parent_vertex is changing
          // within the loop so the edge iterator for the original parent_vertex
          // does not get incremented until and unless edges which cannot be
          // traversed are encountered, we increment it when all children have
          // been traversed like here.
          if (parent_vertex != root_vertex) {
            parent_vertex = parent[parent_vertex];
            _pending_out_edges[parent_vertex].first++;
          } else {
            break;
          }
        }
      } // while(true)
    }   // if(dfs_color[parent_vertex] == DFS_COLOR::WHITE...
  }     // for(int parent_vertex = 0; ...

  // Remember topology_next[vertex] == -1 means there is nothing after the
  // vertex in topological ordering.
  if (topology_initialized) {
    for (int vertex = topology_start_vertex; vertex >= 0;
         vertex = topology_next[vertex]) {
      edge_iterator eit, eit_end;
      std::tie(eit, eit_end) = outEdges(vertex);
      while ((_vertices[vertex].excess > 0) && (eit != eit_end)) {
        DEBUG_INCREMENT(_num_pushes);
        int to_vertex = eit->to_vertex;
        capacity_t flow = std::min(eit->residual, _vertices[vertex].excess);
        eit->residual -= flow;
        _adjacency_list[to_vertex][eit->reverse_edge_index].residual += flow;
        _vertices[vertex].excess -= flow;
        _vertices[to_vertex].excess += flow;
        eit++;
      }
    }
  }
} // end convertPreflowToFlow

#endif // MAX_PUSH_RELABEL_INCLUDED
