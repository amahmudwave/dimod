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

// This class should contain all the information to recreate
// a posiform corresponding to a BQM. The intention is to reduce
// the memory footprint as much as possible. All bias values
// must be multiplied by the ratio before using. Not storing
// the results to avoid extra O(N^2) storage.
template <class BQM>
class PosiformInfo {
	public: 
		using bias_type = typename BQM::bias_type;
		using variable_type = typename BQM::variable_type;
		using t_quadIter = typename BQM::const_outvars_iterator;

		PosiformInfo(const BQM& bqm) {
			numVars = bqm.num_variables();
			linear.resize(numVars, 0);
			quadraticIterators.resize(numVars);
			numQuadratic.resize(numVars,0);
			std::vector<bool> variableUsed(numVars, false);
			
			numLinear = 0;
			maxAbsValue = 0;
			cnst = 0;
			for(int i = 0; i < numVars; i++) {
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
			for(int i = 0; i < numVars; i++){
				linear[i] = convertToLL(bqm.linear(i));
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
       							linear[i]+= biasQuadLL;
       						}
						numBiases++;
					}
				}
				numQuadratic[i]+= numBiases;
			}
			
			for(int i = 0; i < linear.size(); i++){
			       if(linear[i]) {
				       variableUsed[i] = true;
				       numLinear++;
			       }
			       if(linear[i] < 0) cnst += linear[i];
			}
			
			usedVariables.reserve(numVars);
			numUsedVars = 0;
			for(int i = 0; i < numVars; i++) {
				if(variableUsed[i]) {
					usedVariables.push_back(i);
					variableMap.insert({i, numUsedVars});
					numUsedVars++;
				}
			}
			assert(numUsedVars != usedVariables.size());

		}

		void print() {
			std::cout <<"Posiform Information : " << numVars << endl;
			std::cout <<"Num Variables : " << numVars << endl;
			std::cout <<"Constant : " << cnst << std::endl;
			std::cout <<"maxAbsValue : " << maxAbsValue << std::endl;
			std::cout << "Linear : " << std::endl;
			for(int i =0; i < numVars; i++) {
				if(linear[i])
					std::cout << variableMap[i] <<" " <<  i <<" " << linear[i] << std::endl;
		        }

			std::cout << "Quadratic : " << std::endl;
			for(int i =0; i < numVars; i++) {
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
		inline int getNumLinear() { return numLinear; }
		inline long long int getLinear(int i) { return linear[usedVariables[i]]; }
		inline int getNumQuadratic(int i) { return numQuadratic[usedVariables[i]]; }
		inline int getMappedVariable(int i) { return variableMap[i]; }

		inline std::pair<t_quadIter, t_quadIter> getQuadratic(int i) { 
			return quadraticIterators[usedVariables[i]];
		}

		inline long long int convertToLL(bias_type bias) {
	           return static_cast<long long int>(bias * ratio);
		}

	private:
		std::vector<std::pair<t_quadIter,t_quadIter>> quadraticIterators;
		std::vector<long long int> linear;
		std::vector<int> numQuadratic;
		std::vector<int> usedVariables;
		std::unordered_map<int, int> variableMap;
		long long int cnst;
		double maxAbsValue;
		double ratio;
		int numVars;
		int numLinear;
		int numUsedVars;
};

#endif // POSIFORM_INFORMATION_INCLUDED
