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

// Check if the flow value in a given graph represented as an adjacency list is
// valid or not.
template <class EdgeType>
std::pair<typename EdgeType::capacity_type, bool>
isFlowValid(std::vector<std::vector<EdgeType>> &adjList, int source, int sink) {

  using capacity_t = typename EdgeType::capacity_type;

  bool valid_flow = true;
  std::vector<capacity_t> excess(adjList.size(), 0);

  std::cout << "Validating flow of flow network ..." << std::endl;

  // Since we are validating our algorithms, we will not retrieve the
  // value of residual/capacity of a reverse edge from its counterpart
  // which we generally do for performance reasons. Here we will actually
  // access the data and verify if the flow constraints hold or not.
  for (int i = 0; i < adjList.size(); i++) {
    for (int j = 0; j < adjList[i].size(); j++) {
      int to_vertex = adjList[i][j].to_vertex;
      int reverse_edge_index = adjList[i][j].reverse_edge_index;
      capacity_t edge_capacity = adjList[i][j].getCapacity();
      capacity_t edge_residual = adjList[i][j].residual;
      capacity_t reverse_edge_capacity =
          adjList[to_vertex][reverse_edge_index].getCapacity();
      capacity_t reverse_edge_residual =
          adjList[to_vertex][reverse_edge_index].residual;
      bool valid_edge =
          (adjList[i][j].getReverseEdgeCapacity() == reverse_edge_capacity) &&
          (adjList[i][j].getReverseEdgeResidual() == reverse_edge_residual) &&
          (edge_capacity >= 0) && (edge_residual >= 0);
      if (edge_capacity > 0) {
        // Residual edge having capacity 0 is a alid assumption for posiforms,
        // since no term with two variables appear multiple times with different
        // ordering of the variables. This assumption can be maintained with
        // other graphs too.
        valid_edge = valid_edge && (reverse_edge_capacity == 0) &&
                     (edge_residual <= edge_capacity) &&
                     ((edge_residual + reverse_edge_residual) == edge_capacity);

        capacity_t flow = (edge_capacity - edge_residual);
        excess[i] -= flow;
        excess[to_vertex] += flow;
      }
      if (!valid_edge) {
        std::cout << "Invalid Flow due to following edge pair :" << std::endl;
        adjList[i][j].print();
        adjList[to_vertex][reverse_edge_index].print();
      }
      valid_flow = valid_flow && valid_edge;
    }
  }

  for (int i = 0; i < excess.size(); i++) {
    if ((i == source) || (i == sink))
      continue;
    if (excess[i]) {
      std::cout << "Excess flow of " << excess[i] << " in vertex : " << i
                << std::endl;
      valid_flow = false;
    }
  }

  if (excess[sink] != -excess[source]) {
    std::cout << "Flow out of source is not equal to flow into sink."
              << std::endl;
    valid_flow = false;
  }

  return {excess[sink], valid_flow};
}

// Perform reverse breadth first search from a certain vertex, a depth of the
// number of vertices means that vertex could not be reached from the
// start_vertex, since the maximum depth can be equal to number of vertices -1.
template <class EdgeType>
void breadthFirstSearch(std::vector<std::vector<EdgeType>> &adjList,
                        int start_vertex, std::vector<int> &depth_values,
                        bool reverse = false, bool print_result = false) {
  using capacity_t = typename EdgeType::capacity_type;
  int num_vertices = adjList.size();
  depth_values.resize(num_vertices);
  std::fill(depth_values.begin(), depth_values.end(), num_vertices);
  depth_values[start_vertex] = 0;
  vector_based_queue<int> vertexQ(num_vertices);
  vertexQ.reset();
  vertexQ.push(start_vertex);

  // The check for whether the search should be reverse or not could be done
  // inside the innermost loop, but that would be detrimental to performance.
  if (reverse) {
    while (!vertexQ.empty()) {
      int v_parent = vertexQ.pop();
      int children_height = depth_values[v_parent] + 1;
      auto it = adjList[v_parent].begin();
      auto itEnd = adjList[v_parent].end();
      for (; it != itEnd; it++) {
        int to_vertex = it->to_vertex;
        if (it->getReverseEdgeResidual() &&
            depth_values[to_vertex] == num_vertices) {
          depth_values[to_vertex] = children_height;
          vertexQ.push(to_vertex);
        }
      }
    }
  } else {
    while (!vertexQ.empty()) {
      int v_parent = vertexQ.pop();
      int children_height = depth_values[v_parent] + 1;
      auto it = adjList[v_parent].begin();
      auto itEnd = adjList[v_parent].end();
      for (; it != itEnd; it++) {
        int to_vertex = it->to_vertex;
        if (it->residual && depth_values[to_vertex] == num_vertices) {
          depth_values[to_vertex] = children_height;
          vertexQ.push(to_vertex);
        }
      }
    }
  }

  if (print_result) {
    std::vector<int> level_sizes;
    std::vector<std::vector<int>> levels;
    levels.resize(num_vertices + 1);
    level_sizes.resize(num_vertices + 1, 0);
    for (int i = 0; i < depth_values.size(); i++) {
      level_sizes[depth_values[i]]++;
    }
    for (int i = 0; i < level_sizes.size(); i++) {
      levels[i].reserve(level_sizes[i]);
    }
    for (int i = 0; i < depth_values.size(); i++) {
      levels[depth_values[i]].push_back(i);
    }
    std::cout << endl;
    std::cout << "Printing " << (reverse ? "reverse " : "")
              << "breadth first search result starting from vertex : "
              << start_vertex << std::endl;
    std::cout << endl;
    for (int i = 0; i < levels.size(); i++) {
      if (!levels[i].size())
        continue;
      std::cout << "Level " << i << " has " << levels[i].size()
                << " vertices : " << std::endl;
      for (int j = 0; j < levels[i].size(); j++) {
        std::cout << levels[i][j] << " ";
      }
      std::cout << endl;
    }
    std::cout << endl;
  }
}

// Check if the flow value in a given graph represented as an adjacency list is
// a valid max-flow or not and also return the flow value.
template <class EdgeType>
std::pair<typename EdgeType::capacity_type, bool>
isMaximumFlow(std::vector<std::vector<EdgeType>> &adjList, int source,
              int sink) {

  // If the flow follows the constraints of network flow.
  auto validity_result = isFlowValid(adjList, source, sink);

  // If the flow is a maximum flow, the source will be unreachable from the sink
  // through a reverse breadth first search, meaning the source cannot reach the
  // sink through any augmenting path.
  std::vector<int> depth_values;
  int num_vertices = adjList.size();
  breadthFirstSearch(adjList, sink, depth_values, true, true);
  return {validity_result.first,
          (validity_result.second && (depth_values[source] == num_vertices))};
}

template <class EdgeType> class push_relabel {
public:
  using edgeIterator = typename vector<EdgeType>::iterator;
  using capacity_t = typename EdgeType::capacity_type;
  using edge_size_t = size_t;

  struct vertex_node_t {
    int id;
    int height;
    capacity_t excess;
    vertex_node_t *next;
    vertex_node_t *prev;
  };

  struct level_t {
    preallocated_linked_list<vertex_node_t> active_vertices;
    preallocated_linked_list<vertex_node_t> inactive_vertices;
  };

  push_relabel(std::vector<std::vector<EdgeType>> &adjList, int source,
               int sink)
      : adjList(adjList), source(source), sink(sink),
        vertexQ(vector_based_queue<int>(adjList.size())) {
    numGRelabels = 0;
    numRelabels = 0;
    numPushes = 0;
    numVertices = adjList.size();
    _vertices.resize(numVertices);
    levels.resize(numVertices);
    vCurrentEdges.resize(numVertices);

    for (int v = 0; v < numVertices; v++) {
      vCurrentEdges[v] = {adjList[v].begin(), adjList[v].end()};
      _vertices[v].id = v;
      _vertices[v].height = 1;
      _vertices[v].excess = 0;
      // std::cout << " V " << v << " height " << _vertices[v].height  <<  "
      // excess " << _vertices[v].excess << std::endl;
    }
    _vertices[source].height = numVertices;
    _vertices[sink].height = 0;

    edgeIterator it, itEnd;
    for (std::tie(it, itEnd) = outEdges(source); it != itEnd; it++) {
      capacity_t flow = it->residual;
      adjList[it->to_vertex][it->reverse_edge_index].residual += flow;
      it->residual = 0;
      _vertices[it->to_vertex].excess += flow;
      numPushes++;
    }

    maxHeight = numVertices - 1;
    maxActiveHeight = 0;
    minActiveHeight = numVertices;

    global_relabel();
  }

  void global_relabel() {

    numGRelabels++;
    for (int h = 0; h <= maxHeight; h++) {
      levels[h].active_vertices.clear();
      levels[h].inactive_vertices.clear();
    }
    maxHeight = 0;
    maxActiveHeight = 0;
    minActiveHeight = numVertices;

    for (int i = 0; i < numVertices; i++) {
      _vertices[i].height = numVertices;
    }

    _vertices[sink].height = 0;
    vertexQ.reset();
    vertexQ.push(sink);
    while (!vertexQ.empty()) {
      int v_parent = vertexQ.pop();
      int children_height = _vertices[v_parent].height + 1;
      edgeIterator it, itEnd;
      for (std::tie(it, itEnd) = outEdges(v_parent); it != itEnd; it++) {
        int to_vertex = it->to_vertex;
        // std::cout << " Parent " << v_parent << " " << children_height -1 <<
        // std::endl;
        if (it->getReverseEdgeResidual() &&
            _vertices[to_vertex].height == numVertices) {

          // std::cout << " Child " << to_vertex << " " << children_height  <<
          // std::endl;
          _vertices[to_vertex].height = children_height;
          maxHeight = std::max(maxHeight, children_height);
          if (_vertices[to_vertex].excess > 0) {
            add_to_active_list(to_vertex);
          } else {
            add_to_inactive_list(to_vertex);
          }
          vertexQ.push(to_vertex);
        }
      }
    }
  }

  void discharge(int vertex) {
    assert(_vertices[vertex].excess > 0);
    while (1) {
      edgeIterator eit, eitEnd;
      int vertexHeight = _vertices[vertex].height;
      for (std::tie(eit, eitEnd) = vCurrentEdges[vertex]; eit != eitEnd;
           eit++) {
        if (eit->residual) {
          int to_vertex = eit->to_vertex;
          int to_vertexHeight = _vertices[to_vertex].height;
          if (vertexHeight == to_vertexHeight + 1) {
            if (to_vertex != sink && _vertices[to_vertex].excess == 0) {
              // remove_from_inactive_list(to_vertex);
              levels[to_vertexHeight].inactive_vertices.erase(
                  &_vertices[to_vertex]);
              // add_to_active_list(to_vertex);
              levels[to_vertexHeight].active_vertices.push_front(
                  &_vertices[to_vertex]);
              maxActiveHeight = std::max(to_vertexHeight, maxActiveHeight);
              minActiveHeight = std::min(to_vertexHeight, minActiveHeight);
              assert(maxActiveHeight >= 0 && maxActiveHeight < numVertices);
              assert(minActiveHeight >= 0 && minActiveHeight < numVertices);
            }
            // Push flow inlined here
            capacity_t flow = std::min(eit->residual, _vertices[vertex].excess);
            eit->residual -= flow;
            adjList[to_vertex][eit->reverse_edge_index].residual += flow;
            _vertices[vertex].excess -= flow;
            _vertices[to_vertex].excess += flow;
            if (_vertices[vertex].excess == 0)
              break;
          }
        }
      }

      if (eit == eitEnd) {
        int preRelabelHeight = vertexHeight;
        relabel(vertex);
        if (levels[preRelabelHeight].active_vertices.empty() &&
            levels[preRelabelHeight].inactive_vertices.empty()) {
          gap_relabel(preRelabelHeight);
        }
        if (_vertices[vertex].height == numVertices)
          break;
      } else {

        vCurrentEdges[vertex].first = eit;
        add_to_inactive_list(vertex);
        break;
      }
    }
  }

  void gap_relabel(int emptyLevelHeight) {
    for (auto levelIt = levels.begin() + emptyLevelHeight + 1;
         levelIt < levels.begin() + maxHeight; levelIt++) {
      assert(levelIt->active_vertices.empty());
      int inactiveLevelSize = levelIt->inactive_vertices.size();
      vertex_node_t *pVertexNode = levelIt->inactive_vertices.front();
      for (int i = 0; i < inactiveLevelSize; i++) {
        _vertices[pVertexNode->id].height = numVertices;
        pVertexNode = pVertexNode->next;
      }
      levelIt->inactive_vertices.clear();
    }
    maxHeight = emptyLevelHeight - 1;
    maxActiveHeight = emptyLevelHeight - 1;
    assert(maxHeight >= 0 && maxHeight < numVertices);
    assert(maxActiveHeight >= 0 && maxActiveHeight < numVertices);
  }

  void relabel(int vertex) {
    int minRelabelHeight = numVertices;
    _vertices[vertex].height = minRelabelHeight;
    edgeIterator eit, eitEnd, eitMinRelabel;
    for (std::tie(eit, eitEnd) = outEdges(vertex); eit != eitEnd; eit++) {
      int to_vertex = eit->to_vertex;
      if (eit->residual && _vertices[to_vertex].height < minRelabelHeight) {
        minRelabelHeight = _vertices[to_vertex].height;
        eitMinRelabel = eit;
      }
    }
    minRelabelHeight++;
    if (minRelabelHeight < numVertices) {
      _vertices[vertex].height = minRelabelHeight;
      vCurrentEdges[vertex].first = eitMinRelabel;
      maxHeight = std::max(maxHeight, minRelabelHeight);
    }
  }

  capacity_t maximum_preflow() {
    while (maxActiveHeight >= minActiveHeight) {

      // std::cout<< " Max active height " << maxActiveHeight <<" Min active
      // height " << minActiveHeight << std::endl;
      if (levels[maxActiveHeight].active_vertices.empty()) {
        maxActiveHeight--;
      } else {
        vertex_node_t *pVertexNode =
            levels[maxActiveHeight].active_vertices.pop();
        // std::cout << " Going to discharge " << pVertexNode->id << " at height
        // " << _vertices[pVertexNode->id].height << std::endl;
        discharge(pVertexNode->id);
      }
    }
    return _vertices[sink].excess;
  }

  void push(int from_vertex, edgeIterator eit) {
    int to_vertex = eit->to_vertex;
    capacity_t flow = std::min(eit->residual, _vertices[from_vertex].excess);
    eit->residual -= flow;
    adjList[to_vertex][eit->reverse_edge_index].residual += flow;
    _vertices[from_vertex].excess -= flow;
    _vertices[to_vertex].excess += flow;
  }

  std::pair<edgeIterator, edgeIterator> outEdges(int vertex) {
    return {adjList[vertex].begin(), adjList[vertex].end()};
  }

  void add_to_active_list(int vertex) {
    int height = _vertices[vertex].height;
    levels[height].active_vertices.push_front(&_vertices[vertex]);
    maxActiveHeight = std::max(height, maxActiveHeight);
    minActiveHeight = std::min(height, minActiveHeight);
  }

  void remove_from_active_list(int vertex) {
    levels[_vertices[vertex].height].active_vertices.erase(&_vertices[vertex]);
  }

  void remove_from_inactive_list(int vertex) {
    levels[_vertices[vertex].height].inactive_vertices.erase(
        &_vertices[vertex]);
  }

  void add_to_inactive_list(int vertex) {
    levels[_vertices[vertex].height].inactive_vertices.push_front(
        &_vertices[vertex]);
  }

  void printLevels() {
    std::cout << "Printing Levels. " << levels.size() << " in total for "
              << numVertices << " vertices." << std::endl;
    for (int i = 0; i < levels.size(); i++) {
      int size_active = levels[i].active_vertices.size();
      int size_inactive = levels[i].inactive_vertices.size();
      if ((size_active == 0) && (size_inactive == 0))
        continue;

      std::cout << std::endl;
      std::cout << "Level " << i << std::endl << std::endl;
      std::cout << "Active list : " << size_active << " vertices" << std::endl;
      vertex_node_t *pVertexNode = levels[i].active_vertices.front();

      for (int n = 0; n < size_active; n++) {
        std::cout << pVertexNode->id << " ";
        pVertexNode = pVertexNode->next;
      }
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << "Inactive list : " << size_inactive << " vertices"
                << std::endl;
      pVertexNode = levels[i].inactive_vertices.front();

      for (int n = 0; n < size_inactive; n++) {
        std::cout << pVertexNode->id << " ";
        pVertexNode = pVertexNode->next;
      }
      std::cout << std::endl;
    }
  }

  std::vector<level_t> levels;
  std::vector<vertex_node_t> _vertices;
  vector_based_queue<int> vertexQ;
  std::vector<std::vector<EdgeType>> &adjList;
  std::vector<std::pair<edgeIterator, edgeIterator>> vCurrentEdges;
  edge_size_t numEdges;

  size_t numGRelabels, numPushes, numRelabels;
  int numVertices;
  int maxActiveHeight, minActiveHeight, maxHeight;
  int source;
  int sink;
};

#endif
