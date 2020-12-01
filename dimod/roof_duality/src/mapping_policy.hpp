/**
# Copyright 2018 D-Wave Systems Inc.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
#
================================================================================================
*/

#ifndef IMPLICATION_NETWORK_MAPPING_POLICY_HPP_INCLUDED
#define IMPLICATION_NETWORK_MAPPING_POLICY_HPP_INCLUDED

// PLEASE DO NOT MAKE AN INTERFACE FOR MAPPING POLICY. VIRTUAL DISPATCH WILL ADD
// TO THE ALREADY EXISTING PENALTY OF USING THIS MAPPING CLASS. WE CREATED THIS
// CLASS SINCE WE CAN HAVE PERFORMANCE ADVANTAGE WITH ALTERNATE MAPPING WHILE IT
// IS EASIER TO DEBUG WITH SEQUENTIAL MAPPING, BUT IT IS HARD TO SHIFT BETWEEN
// THEM WITHOUT A DEDICATED CLASS.

// If there are n variables in the Quobo, the posiform will have 2 * ( n + 1 )
// variables, each variable with its complement and a root, X_0/X_source and its
// complement. Here we treat variable 0-n-1 as the original n Qubo variables and
// variable n as X_0/X_source. Variable n+1 to 2n+1 are their complements, thus
// variable v has its complement as (v + n + 1). This method of mapping is easy
// for debugging purpose but does not keep the ordering that may be in the qubo
// as regards to the biases and the corresponding variables.
class sequentialMapper {
public:
  sequentialMapper(int num_variables) : _num_variables(num_variables) {}

  sequentialMapper() {}
  // These functions control variable to vertex mapping. Do not manually map
  // them even though it may seem easy, since we may change the mapping by these
  // functions later. The assertions should be changed if the mapping is
  // changed.

  inline int source() { return _num_variables; }

  inline int sink() { return 2 * _num_variables + 1; }

  inline int complement(int vertex) {
    assert((vertex <= 2 * _num_variables + 1) && (vertex >= 0));
    return (vertex <= _num_variables) ? (vertex + _num_variables + 1)
                                      : (vertex - _num_variables - 1);
  }

  inline int non_complemented_vertex(int vertex) {
    assert((vertex <= 2 * _num_variables + 1) && (vertex >= 0));
    return (vertex <= _num_variables) ? vertex : (vertex - _num_variables - 1);
  }

  inline int vertex_to_variable(int vertex) {
    assert((vertex < _num_variables) && (vertex >= 0));
    return vertex;
  }

  inline int variable_to_vertex(int variable) {
    assert((variable < _num_variables) && (variable >= 0));
    return variable;
  }

  inline bool complement_maintains_order() { return false; }

private:
  int _num_variables;
};

#endif
