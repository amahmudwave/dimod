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
#ifndef FIX_VARIABLES_HPP_INCLUDED
#define FIX_VARIABLES_HPP_INCLUDED

#include <vector>
#include <utility>

#include "compressed_matrix.hpp"
#include "dimod/adjarraybqm.h"
#include "dimod/adjmapbqm.h"
#include "dimod/adjvectorbqm.h"

#include <iostream>
using namespace std;

namespace fix_variables_
{

struct FixVariablesResult
{
	std::vector<std::pair<int,  int> > fixedVars; //1-based
	compressed_matrix::CompressedMatrix<double> newQ; //0-based
	double offset;
};

FixVariablesResult fixQuboVariables(const compressed_matrix::CompressedMatrix<double>& Q, int method);

std::vector<std::pair<int,  int> > fixQuboVariablesMap(std::map<std::pair<int, int>, double> QMap, int QSize, int method);

template<class V, class B>
std::vector<std::pair<int,  int>> fixQuboVariables(dimod::AdjVectorBQM<V,B>& refBQM, int method)
{
      // Temporary code to maintain compatibility with legacy code
      std::map<std::pair<int, int>, double> QMap;
      int QSize;
      for(int i = 0; i < refBQM.adj.size(); i++){
	auto& iNeighbors = refBQM.adj[i].first;
	auto linear = refBQM.adj[i].second;
        QMap[{i,i}] = linear;
        auto it = std::lower_bound(iNeighbors.begin(), iNeighbors.end(), i+1, dimod::utils::comp_v<V,B>);
        if(it != iNeighbors.end()) {
          for(auto itEnd = iNeighbors.end(); it != itEnd; it++) {
             QMap[{i,it->first}] = it->second;
	  }
	}
      }
     cout << " Calling processed map based function " << endl; 
     return fixQuboVariablesMap(QMap, refBQM.adj.size(), method); 
}

template<class V, class B>
std::vector<std::pair<int,  int>> fixQuboVariables(dimod::AdjArrayBQM<V,B>& arrayBQM, int method) 
{

}

template<class V, class B>
std::vector<std::pair<int,  int>> fixQuboVariables(dimod::AdjMapBQM<V,B>& mapBQM, int method) 
{

}


} // namespace fix_variables_

#endif // FIX_VARIABLES_HPP_INCLUDED



