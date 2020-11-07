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
#ifndef POSIFORM_INFORMATION_HPP_INCLUDED
#define POSIFORM_INFORMATION_HPP_INCLUDED

#include <iostream>
#include <unordered_map>
#include <vector>

// This class should contain all the information needed to recreate
// a posiform corresponding to a BQM. The intention is to reduce
// the memory footprint as much as possible, thus requires some
// processing before the stored data can be used.
//
// If linear bias for variable X_i is negative with value -L
// it means in the posiform we have the term  L * X_0 * X_i'
// (X_i' being the complement of X_i)
//
// If a quadratic bias between X_i, X_j is negative with value -L
// with i < j it means in the posiform we have the term L*X_i*X_j'
// (X_i' being the complement of X_i)
//
// For each variable X_i the posiform keeps the iterators in the
// provided bqm for biases starting from X_j, (j>i) to X_n.
//
// The linear biases are stored in integer format and can be used
// directly. But the quadratic biases are not stored, instead
// the iterators in the qubo are stored, thus the convertTo**
// function must be applied first on the biases before use.
//
// Note the number of variables and biases exposed should correspond
// to the integer version of the qubo matrix and may be different
// to the numbers corresponding to the floating point based numbers
// as many biases may be flushed to zeroes.
template <class BQM> class PosiformInfo {
public:
  using bias_type = typename BQM::bias_type;
  using quadratic_iterator_type = typename BQM::const_outvars_iterator;
  using variable_type = typename BQM::variable_type;

  PosiformInfo(const BQM &bqm) {

    _constant_posiform = 0;
    _max_absolute_value = 0;
    _num_linear_integer_biases = 0;
    _num_bqm_variables = bqm.num_variables();
    _quadratic_iterators.resize(_num_bqm_variables);
    _linear_double_biases.resize(_num_bqm_variables, 0);
    _linear_integer_biases.resize(_num_bqm_variables, 0);
    _num_quadratic_integer_biases.resize(_num_bqm_variables, 0);

    // Apart from finding the maximum absolute value in the qubo, we must
    // consider the sum of the absolute values of linear values found by using
    // non integer bqm-biases for calculating the conversion ratio. As that
    // value converted to integer is the upper bound for max flow. This is
    // because the linear values correspond to the capacities of the edges in
    // the implication network connecting to source/sink. Thus we calculate the
    // linear values of posiform in double format. Note in the posiform all
    // linear values are positive, here we keep the sign to indicate whether X_0
    // is connected to X_i or X_i'.
    for (int i = 0; i < _num_bqm_variables; i++) {
      auto bqm_linear = bqm.linear(i);
      auto bqm_linear_abs = std::fabs(bqm_linear);
      _linear_double_biases[i] = bqm_linear;
      if (_max_absolute_value < bqm_linear_abs) {
        _max_absolute_value = bqm_linear_abs;
      }
      auto span = bqm.neighborhood(i);
      auto it =
          std::lower_bound(span.first, span.second, i + 1,
                           dimod::utils::comp_v<variable_type, bias_type>);
      _quadratic_iterators[i] = {it, span.second};
      if (it != span.second) {
        for (auto itEnd = span.second; it != itEnd; it++) {
          auto bqm_quadratic = it->second;
          auto bqm_quadratic_abs = std::fabs(bqm_quadratic);
          if (bqm_quadratic < 0) {
            _linear_double_biases[i] += bqm_quadratic;
          }
          if (_max_absolute_value < bqm_quadratic_abs) {
            _max_absolute_value = bqm_quadratic_abs;
          }
        }
      }
    }

    // See comment above regarding calculating conversion rasio.
    _posiform_linear_sum_non_integer = 0;
    for (int i = 0; i < _linear_double_biases.size(); i++) {
      _posiform_linear_sum_non_integer += std::fabs(_linear_double_biases[i]);
    }

    // TODO: Uncomment this code, see comment above.
    // Consider the upper limit of max-flow in implication graph for calculating
    // ratio conversion ratio.

    // if( _max_absolute_value < posiform_linear_sum_non_integer) {
    //	_max_absolute_value = posiform_linear_sum_non_integer;
    // }

    assert(_max_absolute_value != 0);
    _bias_conversion_ratio =
        static_cast<double>(std::numeric_limits<long long int>::max()) /
        _max_absolute_value;

    // TODO: We should be removing this if we consider the sum of the linear
    // biases as mentioned above, but we may divide by 2 to be safe. We must not
    // clamp the ratio to 1 as it is done now.
    _bias_conversion_ratio /= static_cast<double>(1LL << 10);
    if (_bias_conversion_ratio < 1)
      _bias_conversion_ratio = 1;

    _bias_conversion_ratio_limit =
        static_cast<double>(std::numeric_limits<long long int>::max()) /
        std::fabs(_posiform_linear_sum_non_integer);

    for (int i = 0; i < _num_bqm_variables; i++) {
      _linear_integer_biases[i] = convertToLL(bqm.linear(i));
      int num_nonZero_quadratic_biases_in_row = 0;
      auto it = _quadratic_iterators[i].first;
      auto itEnd = _quadratic_iterators[i].second;
      for (; it != itEnd; it++) {
        auto biasQuadLL = convertToLL(it->second);
        if (biasQuadLL) {
          _num_quadratic_integer_biases[it->first]++;
          if (biasQuadLL < 0) {
            _linear_integer_biases[i] += biasQuadLL;
          }
          num_nonZero_quadratic_biases_in_row++;
        }
      }
      _num_quadratic_integer_biases[i] += num_nonZero_quadratic_biases_in_row;
    }

    _posiform_linear_sum_integer = 0; // For debugging purposes.
    for (int i = 0; i < _linear_integer_biases.size(); i++) {
      if (_linear_integer_biases[i]) {
        _num_linear_integer_biases++;
        _posiform_linear_sum_integer += std::abs(_linear_integer_biases[i]);
      }
      if (_linear_integer_biases[i] < 0) {
        _constant_posiform += _linear_integer_biases[i];
      }
    }

    _posiform_to_bqm_variable_map.reserve(_num_bqm_variables);
    _num_posiform_variables = 0;
    for (int i = 0; i < _num_bqm_variables; i++) {
      if (_linear_integer_biases[i] || _num_quadratic_integer_biases[i]) {
        _posiform_to_bqm_variable_map.push_back(i);
        _bqm_to_posiform_variable_map.insert({i, _num_posiform_variables});
        _num_posiform_variables++;
      }
    }
  }

  void print() {
    std::cout << std::endl;
    std::cout << "Posiform Information : " << std::endl << std::endl;
    std::cout << "Number of BQM Variables : " << _num_bqm_variables
              << std::endl;
    std::cout << "Constant : " << _constant_posiform << std::endl;
    std::cout << "Maximum Absolute Value : " << _max_absolute_value
              << std::endl;
    std::cout << "Numeric Limit of long long int "
              << std::numeric_limits<long long int>::max() << std::endl;
    std::cout << "Linear Sum Unconverted : " << _posiform_linear_sum_non_integer
              << std::endl;
    std::cout << "Linear Sum Converted : "
              << convertToLL(_posiform_linear_sum_non_integer) << std::endl;
    std::cout << "Linear Sum After Summing : " << _posiform_linear_sum_integer
              << std::endl;
    std::cout << "Chosen Ratio : " << _bias_conversion_ratio << std::endl;
    std::cout << "Ratio Considering Linear Sum : "
              << _bias_conversion_ratio_limit << std::endl;
    std::cout << std::endl;

    // Printing out in convluted manner, to verify the mappings.
    std::cout << "Used Variables (bqm --> posiform) : "
              << _posiform_to_bqm_variable_map.size() << std::endl;
    for (int i = 0; i < _posiform_to_bqm_variable_map.size(); i++) {
      std::cout
          << _posiform_to_bqm_variable_map[i] << " --> "
          << _bqm_to_posiform_variable_map[_posiform_to_bqm_variable_map[i]]
          << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Linear (posiform, bqm, value) : " << std::endl;
    for (int i = 0; i < _num_bqm_variables; i++) {
      if (_linear_integer_biases[i]) {
        std::cout << _bqm_to_posiform_variable_map[i] << ", " << i << ", "
                  << _linear_integer_biases[i] << std::endl;
      }
    }

    std::cout << std::endl;
    std::cout << "Quadratic (posiform-posiform, bqm-bqm, value): " << std::endl;
    for (int i = 0; i < _num_bqm_variables; i++) {
      auto it = _quadratic_iterators[i].first;
      auto itEnd = _quadratic_iterators[i].second;
      for (; it != itEnd; it++) {
        std::cout << _bqm_to_posiform_variable_map[i] << " "
                  << _bqm_to_posiform_variable_map[it->first] << ", ";
        std::cout << i << " " << it->first << ",  " << convertToLL(it->second)
                  << std::endl;
      }
    }

    std::cout << std::endl;
  }

  inline int getNumVariables() { return _num_posiform_variables; }

  inline int getNumLinear() { return _num_linear_integer_biases; }

  inline long long int getLinear(int i) {
    return _linear_integer_biases[_posiform_to_bqm_variable_map[i]];
  }

  inline int getNumQuadratic(int i) {
    return _num_quadratic_integer_biases[_posiform_to_bqm_variable_map[i]];
  }

  // These functions do not expose the bqm, but the posiform, but for iterating
  // ver the quadratic biases, we need the convertToLL and gttMappedVariable
  // functions.

  inline std::pair<quadratic_iterator_type, quadratic_iterator_type>
  getQuadratic(int i) {
    return _quadratic_iterators[_posiform_to_bqm_variable_map[i]];
  }

  // Convert bqm variable to posiform variable.
  inline int getMappedVariable(int i) {
    return _bqm_to_posiform_variable_map[i];
  }

  inline long long int convertToLL(bias_type bias) {
    return static_cast<long long int>(bias * _bias_conversion_ratio);
  }

private:
  double _max_absolute_value;
  double _bias_conversion_ratio;
  double _bias_conversion_ratio_limit;
  double _posiform_linear_sum_non_integer;
  long long int _posiform_linear_sum_integer;
  int _num_bqm_variables;
  int _num_posiform_variables;
  int _num_linear_integer_biases;
  long long int _constant_posiform;
  std::vector<int> _num_quadratic_integer_biases;
  std::vector<int> _posiform_to_bqm_variable_map;
  std::unordered_map<int, int> _bqm_to_posiform_variable_map;
  std::vector<double> _linear_double_biases;
  std::vector<long long int> _linear_integer_biases;
  std::vector<std::pair<quadratic_iterator_type, quadratic_iterator_type>>
      _quadratic_iterators;
};

#endif // POSIFORM_INFORMATION_INCLUDED
