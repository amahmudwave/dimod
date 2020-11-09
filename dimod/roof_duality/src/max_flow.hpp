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
      adjList[it->toVertex][it->revEdgeIdx].residual += flow;
      it->residual = 0;
      _vertices[it->toVertex].excess += flow;
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
        int toVertex = it->toVertex;
        // std::cout << " Parent " << v_parent << " " << children_height -1 <<
        // std::endl;
        if (it->getReverseEdgeResidual() &&
            _vertices[toVertex].height == numVertices) {

          // std::cout << " Child " << toVertex << " " << children_height  <<
          // std::endl;
          _vertices[toVertex].height = children_height;
          maxHeight = std::max(maxHeight, children_height);
          if (_vertices[toVertex].excess > 0) {
            add_to_active_list(toVertex);
          } else {
            add_to_inactive_list(toVertex);
          }
          vertexQ.push(toVertex);
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
          int toVertex = eit->toVertex;
          int toVertexHeight = _vertices[toVertex].height;
          if (vertexHeight == toVertexHeight + 1) {
            if (toVertex != sink && _vertices[toVertex].excess == 0) {
              // remove_from_inactive_list(toVertex);
              levels[toVertexHeight].inactive_vertices.erase(
                  &_vertices[toVertex]);
              // add_to_active_list(toVertex);
              levels[toVertexHeight].active_vertices.push_front(
                  &_vertices[toVertex]);
              maxActiveHeight = std::max(toVertexHeight, maxActiveHeight);
              minActiveHeight = std::min(toVertexHeight, minActiveHeight);
              assert(maxActiveHeight >= 0 && maxActiveHeight < numVertices);
              assert(minActiveHeight >= 0 && minActiveHeight < numVertices);
            }
            // Push flow inlined here
            capacity_t flow = std::min(eit->residual, _vertices[vertex].excess);
            eit->residual -= flow;
            adjList[toVertex][eit->revEdgeIdx].residual += flow;
            _vertices[vertex].excess -= flow;
            _vertices[toVertex].excess += flow;
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
      int toVertex = eit->toVertex;
      if (eit->residual && _vertices[toVertex].height < minRelabelHeight) {
        minRelabelHeight = _vertices[toVertex].height;
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

  void push(int fromVertex, edgeIterator eit) {
    int toVertex = eit->toVertex;
    capacity_t flow = std::min(eit->residual, _vertices[fromVertex].excess);
    eit->residual -= flow;
    adjList[toVertex][eit->revEdgeIdx].residual += flow;
    _vertices[fromVertex].excess -= flow;
    _vertices[toVertex].excess += flow;
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
    std::cout << "Levels : " << std::endl;
    for (int i = 0; i < levels.size(); i++) {
      std::cout << "Level " << i << std::endl;

      int size = levels[i].active_vertices.size();
      std::cout << "Active list :" << size << " elements" << std::endl;
      vertex_node_t *pVertexNode = levels[i].active_vertices.front();

      for (int n = 0; n < size; n++) {
        std::cout << pVertexNode->id << " ";
        pVertexNode = pVertexNode->next;
      }
      std::cout << std::endl;

      size = levels[i].inactive_vertices.size();
      std::cout << "Inactive list :" << size << " elements" << std::endl;
      pVertexNode = levels[i].inactive_vertices.front();

      for (int n = 0; n < size; n++) {
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
