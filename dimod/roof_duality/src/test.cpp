#include <stdio.h>
#include <vector>
#include <iostream>

// Tarzan's strongly connected component algoirhtm using iterative depth first
// search. For more detail on the iterative version of the depth first search,
// look at the comments in the iterative depth first search in the
// convertPreflowToFlow function in push_relabel.hpp.
int stronglyConnectedComponents(std::vector<std::vector<int>> &adjacency_list,
                                std::vector<int> &components) {
  int num_vertices = adjacency_list.size();
  components.resize(num_vertices);
  std::vector<bool> in_stack(num_vertices, false);
  std::vector<int> seen_stack(num_vertices);
  int seen_stack_num_elements = 0;
  std::vector<int> vertex_visit_id(num_vertices, -1);
  std::vector<int> parent(num_vertices);
  std::vector<int> low_link_id(num_vertices);
  int current_id = 0;
  int num_components = 0;
  using edge_iterator = typename std::vector<int>::iterator;
  std::vector<std::pair<edge_iterator, edge_iterator>> pending_out_edges(num_vertices);

  for (int vertex = 0; vertex < num_vertices; vertex++) {
    pending_out_edges[vertex] = {
      adjacency_list[vertex].begin(),
      adjacency_list[vertex].end()
    };
  }

  // Iterative DFS.
  for (int current_vertex = 0; current_vertex < num_vertices; current_vertex++) {
    if (vertex_visit_id[current_vertex] < 0) {
      // The root of a DFS tree.
      int root_vertex = current_vertex;
      vertex_visit_id[current_vertex] = current_id;
      low_link_id[current_vertex] = current_id;
      current_id++;
      in_stack[current_vertex] = true;
      seen_stack[seen_stack_num_elements++] = current_vertex;
      while (true) {
        for (; pending_out_edges[current_vertex].first !=
               pending_out_edges[current_vertex].second;
             pending_out_edges[current_vertex].first++) {
          int child_vertex = *(pending_out_edges[current_vertex].first);
          if (vertex_visit_id[child_vertex] < 0) {
            vertex_visit_id[child_vertex] = current_id;
            low_link_id[child_vertex] = current_id;
            current_id++;
            in_stack[child_vertex] = true;
            seen_stack[seen_stack_num_elements++] = child_vertex;
            parent[child_vertex] = current_vertex;
            current_vertex = child_vertex;
            break;
          } else if (in_stack[child_vertex]) {
            low_link_id[current_vertex] = std::min(low_link_id[current_vertex], low_link_id[child_vertex]);
            }
        }
        if (pending_out_edges[current_vertex].first ==
            pending_out_edges[current_vertex].second) {
          if (current_vertex != root_vertex) {
            int completed_vertex = current_vertex;
            current_vertex = parent[current_vertex];
            low_link_id[current_vertex] = std::min(low_link_id[current_vertex], low_link_id[completed_vertex]);
            if(low_link_id[current_vertex] == vertex_visit_id[current_vertex]) {
              int popped_vertex = -1;
              do{
                  popped_vertex = seen_stack[--seen_stack_num_elements];
		  in_stack[popped_vertex] = false;
		  components[popped_vertex] = num_components;
              } while(popped_vertex != current_vertex);
	      num_components++;
            }
            pending_out_edges[current_vertex].first++;
          } else {
            break;
          }
        }
      }
    }
  }
  return num_components;
}

int main () {
 
  std::vector<std::vector<int>> adj(5);

  adj[1].push_back(0);
  adj[0].push_back(2);
  adj[2].push_back(1);
  adj[0].push_back(3);
  adj[3].push_back(4);

  std::vector<int> components;
  int numcomp = stronglyConnectedComponents(adj, components);

  std::cout << numcomp << std::endl << std::endl;
  for(int i = 0; i < components.size(); i++) {
    std::cout << i << " in " << components[i] << std::endl;
  }
 return 0;
}
