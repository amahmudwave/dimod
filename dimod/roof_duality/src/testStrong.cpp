#include <iostream>
#include <stdio.h>
#include <vector>

#include "helper_data_structures.hpp"

// Tarzan's strongly connected component algoirhtm using iterative depth first
// search. For more detail on the iterative version of the depth first search,
// look at the comments in the iterative depth first search in the
// convertPreflowToFlow function in push_relabel.hpp.
int stronglyConnectedComponents(std::vector<std::vector<int>> &adjacency_list,
                                std::vector<int> &components) {
  int num_vertices = adjacency_list.size();
  using edge_iterator = typename std::vector<int>::iterator;
  components.resize(num_vertices);
  vector_based_stack<int> component_stack(num_vertices);
  std::vector<bool> in_stack(num_vertices, false);
  int UNVISITED = num_vertices;
  std::vector<int> low_link_id(num_vertices, UNVISITED);
  std::vector<int> vertex_visit_id(num_vertices, UNVISITED);
  std::vector<int> parent(num_vertices);
  int visit_id = 0;
  int num_strong_components = 0;
  std::vector<std::pair<edge_iterator, edge_iterator>> pending_out_edges(
      num_vertices);

  for (int vertex = 0; vertex < num_vertices; vertex++) {
    pending_out_edges[vertex] = {adjacency_list[vertex].begin(),
                                 adjacency_list[vertex].end()};
  }

  // Iterative DFS.
  for (int current_vertex = 0; current_vertex < num_vertices;
       current_vertex++) {
    if (vertex_visit_id[current_vertex] == UNVISITED) {
      // The root of a DFS tree.
      int root_vertex = current_vertex;
      vertex_visit_id[current_vertex] = visit_id;
      low_link_id[current_vertex] = visit_id;
      visit_id++;
      component_stack.push(current_vertex);
      in_stack[current_vertex] = true;
      while (true) {
        for (; pending_out_edges[current_vertex].first !=
               pending_out_edges[current_vertex].second;
             pending_out_edges[current_vertex].first++) {
          int child_vertex = *(pending_out_edges[current_vertex].first);
          if (vertex_visit_id[child_vertex] == UNVISITED) {
            vertex_visit_id[child_vertex] = visit_id;
            low_link_id[child_vertex] = visit_id;
            visit_id++;
            component_stack.push(child_vertex);
            in_stack[child_vertex] = true;
            parent[child_vertex] = current_vertex;
            current_vertex = child_vertex;
            break;
          } else if (in_stack[child_vertex]) {
            low_link_id[current_vertex] = std::min(low_link_id[current_vertex],
                                                   low_link_id[child_vertex]);
          }
        }
        if (pending_out_edges[current_vertex].first ==
            pending_out_edges[current_vertex].second) {
          if (low_link_id[current_vertex] == vertex_visit_id[current_vertex]) {
            int popped_vertex = -1;
            do {
              popped_vertex = component_stack.pop();
              in_stack[popped_vertex] = false;
              components[popped_vertex] = num_strong_components;
            } while (popped_vertex != current_vertex);
            num_strong_components++;
          }
          if (current_vertex != root_vertex) {
            int completed_vertex = current_vertex;
            current_vertex = parent[current_vertex];
            low_link_id[current_vertex] = std::min(
                low_link_id[current_vertex], low_link_id[completed_vertex]);
            pending_out_edges[current_vertex].first++;
          } else {
            break;
          }
        }
      }
    }
  }
  return num_strong_components;
}

// Tarzan's strongly connected component algoirhtm using iterative depth first
// search. For more detail on the iterative version of the depth first search,
// look at the comments in the iterative depth first search in the
// convertPreflowToFlow function in push_relabel.hpp.
int DstronglyConnectedComponents(std::vector<std::vector<int>> &adjacency_list,
                                 std::vector<int> &components) {
  int num_vertices = adjacency_list.size();
  using edge_iterator = typename std::vector<int>::iterator;
  components.resize(num_vertices);
  vector_based_stack<int> component_stack(num_vertices);
  std::vector<bool> in_stack(num_vertices, false);
  int UNVISITED = num_vertices;
  std::vector<int> low_link_id(num_vertices, UNVISITED);
  std::vector<int> vertex_visit_id(num_vertices, UNVISITED);
  std::vector<int> parent(num_vertices);
  int visit_id = 0;
  int num_strong_components = 0;
  std::vector<std::pair<edge_iterator, edge_iterator>> pending_out_edges(
      num_vertices);

  for (int vertex = 0; vertex < num_vertices; vertex++) {
    pending_out_edges[vertex] = {adjacency_list[vertex].begin(),
                                 adjacency_list[vertex].end()};
  }

  // Iterative DFS.
  for (int current_vertex = 0; current_vertex < num_vertices;
       current_vertex++) {
    if (vertex_visit_id[current_vertex] == UNVISITED) {
      // The root of a DFS tree.
      int root_vertex = current_vertex;
      vertex_visit_id[current_vertex] = visit_id;
      low_link_id[current_vertex] = visit_id;
      visit_id++;
      component_stack.push(current_vertex);
      in_stack[current_vertex] = true;
      std::cout << "Pushed Root " << current_vertex << " lowlink/visit id "
                << visit_id - 1 << std::endl;
      while (true) {
        for (; pending_out_edges[current_vertex].first !=
               pending_out_edges[current_vertex].second;
             pending_out_edges[current_vertex].first++) {
          int child_vertex = *(pending_out_edges[current_vertex].first);
          if (vertex_visit_id[child_vertex] == UNVISITED) {
            vertex_visit_id[child_vertex] = visit_id;
            low_link_id[child_vertex] = visit_id;
            visit_id++;
            component_stack.push(child_vertex);
            in_stack[child_vertex] = true;
            parent[child_vertex] = current_vertex;
            std::cout << "Pushed " << child_vertex << " lowlink/visit id "
                      << visit_id - 1 << " of parent " << current_vertex
                      << " lowlink id " << low_link_id[current_vertex]
                      << " visit id " << vertex_visit_id[current_vertex]
                      << std::endl;
            current_vertex = child_vertex;
            break;
          } else if (in_stack[child_vertex]) {
            std::cout << "Met in stack  " << child_vertex << " , "
                      << low_link_id[child_vertex] << " from  "
                      << current_vertex << " ,  " << low_link_id[current_vertex]
                      << std::endl;
            low_link_id[current_vertex] = std::min(low_link_id[current_vertex],
                                                   low_link_id[child_vertex]);
          }
        }
        if (pending_out_edges[current_vertex].first ==
            pending_out_edges[current_vertex].second) {
          std::cout << "Getting out of " << current_vertex << "  owner "
                    << parent[current_vertex] << std::endl;
          if (low_link_id[current_vertex] == vertex_visit_id[current_vertex]) {
            std::cout << "Met scc owner " << current_vertex << std::endl;
            int popped_vertex = -1;
            do {
              popped_vertex = component_stack.pop();
              in_stack[popped_vertex] = false;
              components[popped_vertex] = num_strong_components;
              std::cout << "Popped " << popped_vertex << std::endl;
            } while (popped_vertex != current_vertex);
            num_strong_components++;
          }
          if (current_vertex != root_vertex) {
            int completed_vertex = current_vertex;
            current_vertex = parent[current_vertex];
            low_link_id[current_vertex] = std::min(
                low_link_id[current_vertex], low_link_id[completed_vertex]);
            pending_out_edges[current_vertex].first++;
          } else {
            break;
          }
        }
      }
    }
  }
  return num_strong_components;
}

void test(std::vector<std::vector<int>> &adj) {
  std::vector<int> components;
  int numcomp = stronglyConnectedComponents(adj, components);

  std::cout << "Numcomponents " << numcomp << std::endl << std::endl;
  for (int i = 0; i < components.size(); i++) {
    std::cout << i << " in " << components[i] << std::endl;
  }
}

int main() {

  std::vector<std::vector<int>> adj;
  adj.clear();
  adj.resize(6);
  adj[1].push_back(0);
  adj[0].push_back(2);
  adj[2].push_back(1);
  adj[0].push_back(3);
  adj[3].push_back(4);
  adj[0].push_back(5);
  adj[5].push_back(2);
  test(adj);

  adj.clear();
  adj.resize(4);
  adj[0].push_back(1);
  adj[1].push_back(2);
  adj[2].push_back(3);
  test(adj);

  adj.clear();
  adj.resize(7);
  adj[0].push_back(1);
  adj[1].push_back(2);
  adj[2].push_back(0);
  adj[1].push_back(3);
  adj[1].push_back(4);
  adj[1].push_back(6);
  adj[3].push_back(5);
  adj[4].push_back(5);
  test(adj);

  adj.clear();
  adj.resize(11);
  adj[0].push_back(1);
  adj[0].push_back(3);
  adj[1].push_back(2);
  adj[1].push_back(4);
  adj[2].push_back(0);
  adj[2].push_back(6);
  adj[3].push_back(2);
  adj[4].push_back(5);
  adj[4].push_back(6);
  adj[5].push_back(6);
  adj[5].push_back(7);
  adj[5].push_back(8);
  adj[5].push_back(9);
  adj[6].push_back(4);
  adj[7].push_back(9);
  adj[8].push_back(9);
  adj[9].push_back(8);
  test(adj);

  adj.clear();
  adj.resize(5);
  adj[0].push_back(1);
  adj[1].push_back(2);
  adj[2].push_back(3);
  adj[2].push_back(4);
  adj[3].push_back(0);
  adj[4].push_back(2);
  test(adj);

  return 0;
}
