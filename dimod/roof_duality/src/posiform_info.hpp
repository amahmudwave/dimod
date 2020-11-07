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
# ================================================================================================
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
template <class BQM>
class PosiformInfo {
	public: 
		using bias_type = typename BQM::bias_type;
		using quadratic_iterator_type = typename BQM::const_outvars_iterator;
		using variable_type = typename BQM::variable_type;

		PosiformInfo(const BQM& bqm) {
			_num_bqm_variables = bqm.num_variables();
			_linear_integer_biases.resize(_num_bqm_variables, 0);
			quadraticIterators.resize(_num_bqm_variables);
			numQuadratic.resize(_num_bqm_variables,0);
			std::vector<bool> variableUsed(_num_bqm_variables, false);
			
			_num_linear_integers = 0;
			maxAbsValue = 0;
			cnst = 0;
			for(int i = 0; i < _num_bqm_variables; i++) {
				if(maxAbsValue < std::fabs(bqm.linear(i))) {
					maxAbsValue = std::fabs(bqm.linear(i));
				}
				auto span = bqm.neighborhood(i);
				auto it = std::lower_bound(span.first, span.second, i+1,
						dimod::utils::comp_v<variable_type, bias_type>);
				quadraticIterators[i] = {it, span.second};
				if(it != span.second) {
					for(auto itEnd = span.second; it != itEnd; it++) {
						if(maxAbsValue < std::fabs(it->second)) {
							maxAbsValue = std::fabs(it->second);
						}
					}
				}
			}


			if (maxAbsValue != 0)
				ratio = static_cast<double>(std::numeric_limits<long long int>::max()) / maxAbsValue;
			ratio /= static_cast<double>(1LL << 10);
			if (ratio < 1)
				ratio = 1;

			// TODO Ideally we will take a second pass to find maxAbsValue
			// considering the array of linears. And then recalculate linears
			// and the constant. 
			for(int i = 0; i < _num_bqm_variables; i++){
				_linear_integer_biases[i] = convertToLL(bqm.linear(i));
			        auto it = quadraticIterators[i].first;
				auto itEnd = quadraticIterators[i].second;
				int numBiases = 0;
				for(; it != itEnd; it++) {
					auto biasQuadLL = convertToLL(it->second);
					if( biasQuadLL ) {
						variableUsed[i] = true;
						variableUsed[it->first] = true;
       						numQuadratic[it->first]++;
       						if( biasQuadLL < 0) {
       							_linear_integer_biases[i]+= biasQuadLL;
       						}
						numBiases++;
					}
				}
				numQuadratic[i]+= numBiases;
			}
			
			for(int i = 0; i < _linear_integer_biases.size(); i++){
			       if(_linear_integer_biases[i]) {
				       variableUsed[i] = true;
				       _num_linear_integers++;
			       }
			       if(_linear_integer_biases[i] < 0) cnst += _linear_integer_biases[i];
			}
			
			usedVariables.reserve(_num_bqm_variables);
			numUsedVars = 0;
			for(int i = 0; i < _num_bqm_variables; i++) {
				if(variableUsed[i]) {
					usedVariables.push_back(i);
					variableMap.insert({i, numUsedVars});
					numUsedVars++;
				}
			}
			assert(numUsedVars != usedVariables.size());

		}

		void print() {
			std::cout <<"Posiform Information : " << _num_bqm_variables << endl;
			std::cout <<"Number of BQM Variables : " << _num_bqm_variables << endl;
			std::cout <<"Constant : " << cnst << std::endl;
			std::cout <<"Maximum Absolute Value : " << maxAbsValue << std::endl;
			std::cout <<"Linear : " << std::endl;
			for(int i =0; i < _num_bqm_variables; i++) {
				if(_linear_integer_biases[i])
					std::cout << variableMap[i] <<" " <<  i <<" " << _linear_integer_biases[i] << std::endl;
		        }

			std::cout << "Quadratic : " << std::endl;
			for(int i =0; i < _num_bqm_variables; i++) {
	         		auto it = quadraticIterators[i].first;
				auto itEnd = quadraticIterators[i].second;
				for(; it != itEnd; it++){
					std::cout << variableMap[i] <<" " << variableMap[it->first] << " "; 
					std::cout << i << " " << it->first <<" " << convertToLL(it->second) << std::endl;
				}
			}  
		
			std::cout << "Used Variables :" << usedVariables.size() <<  std::endl;
			for(int i = 0; i < usedVariables.size(); i++) {
				std::cout << usedVariables[i] << " --> " << variableMap[usedVariables[i]] << std::endl;
			}

		}

		inline int getNumVariables() { return numUsedVars; }
		inline int getNumLinear() { return _num_linear_integers; }
		inline long long int getLinear(int i) { return _linear_integer_biases[usedVariables[i]]; }
		inline int getNumQuadratic(int i) { return numQuadratic[usedVariables[i]]; }
		inline int getMappedVariable(int i) { return variableMap[i]; }

		inline std::pair<quadratic_iterator_type, quadratic_iterator_type> getQuadratic(int i) { 
			return quadraticIterators[usedVariables[i]];
		}

		inline long long int convertToLL(bias_type bias) {
	           return static_cast<long long int>(bias * ratio);
		}

	private:
		std::vector<std::pair<quadratic_iterator_type,quadratic_iterator_type>> quadraticIterators;
		std::vector<long long int> _linear_integer_biases;
		std::vector<int> numQuadratic;
		std::vector<int> usedVariables;
		std::unordered_map<int, int> variableMap;
		long long int cnst;
		double maxAbsValue;
		double ratio;
		int _num_bqm_variables;
		int _num_linear_integers;
		int numUsedVars;
};

#endif // POSIFORM_INFORMATION_INCLUDED
