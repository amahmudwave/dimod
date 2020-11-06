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
#include "posiform_info.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/strong_components.hpp>

//for debugging only
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>


#include <iostream>
using namespace std;

namespace
{

class compClass
{
public:
	bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b)
	{
		if (a.second != b.second)
			return !(a.second < b.second);
		else
			return a.first < b.first;
	}
};

bool compareAbs(double a, double b)
{
	return std::abs(a) < std::abs(b);
}

//the index is 1-based, according to the paper
struct Posiform
{
	std::vector<std::pair<std::pair<int, int>, long long int> > quadratic;
	std::vector<std::pair<int, long long int> > linear;
	long long int cst;
	int numVars;
};

void printPosiform(Posiform& pos) {
    
     std::cout <<"Num Variables : " << pos.numVars << endl;
     std::cout <<"Constant : " << pos.cst << std::endl;
     std::cout <<"Linear : " << std::endl;
     for(size_t i = 0 ; i < pos.linear.size(); i++) {
	     std::cout << pos.linear[i].first <<" " 
                       << pos.linear[i].second << std::endl;
     }
     std::cout <<"Quadratic : " << std::endl;
     for(size_t i = 0 ; i < pos.quadratic.size(); i++) {
	     std::cout << pos.quadratic[i].first.first <<" " 
		       << pos.quadratic[i].first.second <<" "
		       << pos.quadratic[i].second << std::endl;
     }
}

template <typename capacity_t>
class  ImplicationEdge {
public:
    typedef capacity_t capacity_type;
    ImplicationEdge(int toVertex, capacity_t capacity): toVertex(toVertex), residual(capacity) { }
    int toVertex;
    int revEdgeIdx;
    int symmEdgeIdx;
    capacity_t residual;
};


template <class capacity_t>
class ImplicationNetwork {
	public:
	template <class BQM>
        ImplicationNetwork(BQM& bqm) {

		numVariables = bqm.getNumVariables();
		numVertices = 2 * numVariables + 2;
		source = numVariables;
		sink = 2 * numVariables + 1;
		adjList.resize(2 * numVariables + 2);
		sizeEstimates.resize(numVertices, 0); 

		// After this point complement function 
		// should work.
		assert(sink == complement(source));
		assert(source == complement(sink));
 
		int numLinear = bqm.getNumLinear();
		adjList[source].reserve(numLinear);
		adjList[sink].reserve(numLinear);
	        sizeEstimates[source] = numLinear;
		sizeEstimates[sink] = numLinear;

	        // There are reverse edges for each edge created
	        // in the implication graph. Depending on the sign
	        // of the bias, an edge may start from v or v' but
	        // reverse edges makes the edges coming out of v and
	        // v' both equal to the number of quadratic biases
	        // in which v contributes + 1 due to the linear term.
                for(int u = 0; u <numVariables; u++) {
		  int uComp = complement(u);
		  int numEntries = bqm.getNumQuadratic(u);
		  auto linear = bqm.getLinear(u);
		  if(linear) numEntries++;
		  adjList[u].reserve(numEntries);
		  adjList[uComp].reserve(numEntries);
		  sizeEstimates[u] = numEntries;
		  sizeEstimates[uComp] = numEntries;

		  if(linear > 0) {
			createImplicationNetworkEdges(source, uComp, linear);
		  } else if ( linear < 0) {
		 	createImplicationNetworkEdges(source, u, -linear);
		  } 

		  auto quadraticSpan = bqm.getQuadratic(u);
		  auto it = quadraticSpan.first;
                  auto itEnd = quadraticSpan.second;
                  for(; it != itEnd; it++){
                     auto bias =  bqm.convertToLL(it->second);
		     int v = bqm.getMappedVariable(it->first);
		     if(bias > 0) {
                        createImplicationNetworkEdges(u, complement(v), bias);
		     } else if ( bias < 0) {
                        createImplicationNetworkEdges(u, v, -bias);
		     }
                  }
		}
	}
	
	inline int complement( int v ) { 
		if ( v <= numVariables)
			return (v + numVariables + 1);
		else 
			return (v - numVariables - 1);
	}

	void fillLastOutEdgeReferences(int from, int to) {
		auto& edge = adjList[from].back();
	        edge.revEdgeIdx = adjList[to].size() - 1;
		int symmetricFrom = complement(to);
		edge.symmEdgeIdx = adjList[symmetricFrom].size() - 1;
	}

	void createImplicationNetworkEdges(int from, int to, capacity_t capacity) {
		int fromComp = complement(from);	
		int toComp  = complement(to);
		adjList[from].emplace_back(ImplicationEdge<capacity_t>(to, capacity));
		adjList[to].emplace_back(ImplicationEdge<capacity_t>(from, 0));
		adjList[toComp].emplace_back(ImplicationEdge<capacity_t>(fromComp, capacity));
		adjList[fromComp].emplace_back(ImplicationEdge<capacity_t>(toComp, 0));
		fillLastOutEdgeReferences(from, to);
		fillLastOutEdgeReferences(to, from);
		fillLastOutEdgeReferences(toComp, fromComp);
		fillLastOutEdgeReferences(fromComp, toComp);
	}

	void print() {
		std::cout <<"Implication Graph Information : " << std::endl;
		std::cout <<"Num Variables : " << numVariables << std::endl;
		std::cout <<"Num Vertices : " << numVertices << std::endl;
	        for(int i = 0; i < adjList.size(); i++) {
		  if(adjList[i].size() != sizeEstimates[i]) {
 		    std::cout <<"Inaccurate size estimate for out edges for " << i <<" " << sizeEstimates[i] << std::endl;
		  }

		  for(int j = 0; j < adjList[i].size(); j++) {
		    auto& node = adjList[i][j];
		    std::cout << "{ " << i << " --> " << node.toVertex << " " << node.residual << " ";
		    std::cout <<  node.revEdgeIdx << " " << node.symmEdgeIdx << " } " << std::endl;	    
                  }
		  std::cout << endl;
		  assert(adjList[i].size() == sizeEstimates[i]);
		}	
	}

	vector<vector<ImplicationEdge<capacity_t>>> adjList;
	// TODO : Verify size estimates are correct, for debugging, remove it. 
	vector<int> sizeEstimates;
	int numVariables;
	int numVertices;
	int source;
	int sink;
};

template <typename node>
    class linked_list
    {
        node _head = node { }, _tail = node { };
        std::size_t _size { 0 };
    public:
        linked_list ( )
        {
            _head . next = &_tail;
            _tail . prev = &_head;
	    _size = 0;
        }

        node * pop ( ) noexcept
        {
            auto * ret = _head . next;
            _head . next = _head . next -> next;
            _head . next -> prev = &_head;
            --_size;
            return ret;
        }

        void push ( node * n ) noexcept
        {
            n -> next = &_tail;
            n -> prev = _tail . prev;
            _tail . prev -> next = n;
            _tail . prev = n;
            ++_size;
        }

        void push_front ( node * n ) noexcept
        {
            n -> next = _head . next;
            n -> prev = &_head;
            _head . next -> prev = n;
            _head . next = n;
            ++_size;
        }

        void erase ( node * n ) noexcept
        {
            n -> prev -> next = n -> next;
            n -> next -> prev = n -> prev;
            --_size;
        }

        node * front ( ) const noexcept
        {
            return _head . next;
        }

        node * back ( ) const noexcept
        {
            return _tail . prev;
        }

        bool empty ( ) const noexcept
        {
            return _size == 0;
        }

        void clear ( ) noexcept
        {
            _head . next = &_tail;
            _tail . prev = &_head;
            _size = 0;
        }

        std::size_t size ( ) const noexcept
        {
            return _size;
        }

        void append_list ( linked_list & other ) noexcept
        {
            if ( other . empty () )
                return;
            auto other_head = other . front ();
            auto other_tail = other . back ();
            this -> back () -> next = other_head;
            other_head -> prev = this -> back ();
            _tail . prev = other_tail;
            other_tail -> next = &_tail;
            _size += other . size ();
            other . clear ();
        }
    };

template <class T>
class vecQueue {
	public:
   vecQueue(int size) { front = 0; back = 0; data.resize(size); }
   void reset () { front = 0; back = 0; }
   void push(T val) { data[back++] = val;}
   T pop() { return data[front++]; }
   bool empty() { return (front == back); }
   int front;
   int back;
   std::vector<T> data;
};

template <class EdgeType>
class push_relabel {
	public:
	using edgeIterator = typename vector<EdgeType>::iterator;
	using capacity_t = typename EdgeType::capacity_type; 
	using edge_size_t = size_t;

  	struct vertex_node_t {
	  int id;  
	  int height;
	  capacity_t excess;
	  vertex_node_t * next;
	  vertex_node_t * prev;
	};

	push_relabel(std::vector<std::vector<EdgeType>>& adjList, int source, int sink) : adjList(adjList), source(source), sink(sink), vertexQ(vecQueue<int>(adjList.size())) 
	{
		numGRelabels = 0;
		numRelabels = 0;
		numPushes = 0;
		numVertices = adjList.size();
		_vertices.resize(numVertices);
                levels.resize(numVertices);
		vCurrentEdges.resize(numVertices);


		for(int v = 0 ; v < numVertices; v++) {
			vCurrentEdges[v] = {adjList[v].begin(), adjList[v].end()};
			_vertices[v].id = v;
			_vertices[v].height = 1;
			_vertices[v].excess = 0;
			//std::cout << " V " << v << " height " << _vertices[v].height  <<  " excess " << _vertices[v].excess << std::endl;
		}
		_vertices[source].height = numVertices;
		_vertices[sink].height = 0;
		
		edgeIterator it, itEnd;
                for(std::tie(it, itEnd) = outEdges(source); it != itEnd; it++) {
	            capacity_t flow  = it->residual;
		    adjList[it->toVertex][it->revEdgeIdx].residual+= flow;
		    it->residual = 0;
                    _vertices[it->toVertex].excess+= flow;
		    numPushes++;
	        }		
		
		maxHeight = numVertices -1;
		maxActiveHeight = 0;
		minActiveHeight = numVertices;
				
		global_relabel();
	}
	
	void global_relabel() {

		numGRelabels++;
		for(int h = 0; h <= maxHeight; h++)
		{
			levels[h].active_vertices.clear();
			levels[h].inactive_vertices.clear();
		}	
                maxHeight = 0;
		maxActiveHeight =0;
		minActiveHeight = numVertices;

		for(int i = 0; i < numVertices; i++) {
			_vertices[i].height = numVertices;
		}

		_vertices[sink].height = 0;
		vertexQ.reset();
		vertexQ.push(sink);
 		while(!vertexQ.empty()) {
		   int v_parent = vertexQ.pop();
		   int children_height = _vertices[v_parent].height + 1;
		   edgeIterator it, itEnd;
		   for(std::tie(it, itEnd) = outEdges(v_parent); it != itEnd; it++) {
			int toVertex = it->toVertex;
			std::cout << " Parent " << v_parent << " " << children_height -1 << std::endl;
			if(adjList[toVertex][it->revEdgeIdx].residual && _vertices[toVertex].height == numVertices)
			{

				std::cout << " Child " << toVertex << " " << children_height  << std::endl;
				_vertices[toVertex].height = children_height;
				maxHeight = std::max(maxHeight, children_height);
				if(_vertices[toVertex].excess > 0){
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
		while(1) {
			edgeIterator eit, eitEnd;
			int vertexHeight = _vertices[vertex].height;
			for(std::tie(eit, eitEnd) = vCurrentEdges[vertex]; eit != eitEnd; eit++) {
				if(eit->residual) {
					int toVertex = eit->toVertex;
					int toVertexHeight = _vertices[toVertex].height;
					if(vertexHeight == toVertexHeight + 1) {
						if(toVertex != sink && _vertices[toVertex].excess == 0){
							// remove_from_inactive_list(toVertex);
							levels[toVertexHeight].inactive_vertices.erase(&_vertices[toVertex]);
							// add_to_active_list(toVertex);
							levels[toVertexHeight].active_vertices.push_front(&_vertices[toVertex]);
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
						if(_vertices[vertex].excess  == 0) break;
					}
				}
			}

			if( eit == eitEnd ) {
				int preRelabelHeight = vertexHeight;
				relabel(vertex);
				if(levels[preRelabelHeight].active_vertices.empty() && 
				   levels[preRelabelHeight].inactive_vertices.empty()) {
					gap_relabel(preRelabelHeight);
				}
				if(_vertices[vertex].height == numVertices) break;
			} else {
					
				vCurrentEdges[vertex].first = eit; 
				add_to_inactive_list(vertex);
				break;
			}
		}
	}


	void gap_relabel(int emptyLevelHeight) {
		for(auto levelIt = levels.begin() + emptyLevelHeight + 1; levelIt < levels.begin() + maxHeight; levelIt++) {
			assert(levelIt->active_vertices.empty());
			int inactiveLevelSize = levelIt->inactive_vertices.size();
			vertex_node_t* pVertexNode = levelIt->inactive_vertices.front();
		        for(int i = 0; i < inactiveLevelSize; i++) {
				_vertices[pVertexNode->id].height = numVertices;
				pVertexNode = pVertexNode->next;
			}	
			levelIt->inactive_vertices.clear();
		}
		maxHeight = emptyLevelHeight -1;
		maxActiveHeight = emptyLevelHeight -1;
		assert(maxHeight >= 0 && maxHeight < numVertices);
		assert(maxActiveHeight >= 0 && maxActiveHeight < numVertices);
	}



        void relabel(int vertex) {
		int minRelabelHeight = numVertices;
		_vertices[vertex].height = minRelabelHeight;
		edgeIterator eit, eitEnd, eitMinRelabel;
		for(std::tie(eit, eitEnd) = outEdges(vertex); eit != eitEnd; eit++) {
			int toVertex = eit->toVertex;
			if(eit->residual && _vertices[toVertex].height < minRelabelHeight) {
				minRelabelHeight = _vertices[toVertex].height;
				eitMinRelabel = eit;	
			}
		}
		minRelabelHeight++;
		if(minRelabelHeight < numVertices) {
			_vertices[vertex].height = minRelabelHeight;
			vCurrentEdges[vertex].first = eitMinRelabel;
			maxHeight = std::max(maxHeight, minRelabelHeight);
		}
	}

		
	capacity_t maximum_preflow() {
		while(maxActiveHeight >= minActiveHeight) {

			std::cout<< " Max active height " << maxActiveHeight <<" Min active height " << minActiveHeight << std::endl;
			if(levels[maxActiveHeight].active_vertices.empty()) {
				maxActiveHeight--;
			}
			else{
				vertex_node_t* pVertexNode = levels[maxActiveHeight].active_vertices.pop();
				std::cout << " Going to discharge " << pVertexNode->id << " at height " << _vertices[pVertexNode->id].height << std::endl;
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
	 	return {adjList[vertex].begin(), adjList[vertex].end()} ;
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
	  levels[_vertices[vertex].height].inactive_vertices.erase(&_vertices[vertex]);
	}

      	void add_to_inactive_list(int vertex) {
	   levels[_vertices[vertex].height].inactive_vertices.push_front(&_vertices[vertex]);
	}
	
	void printLevels() {
		std::cout <<"Levels : " << std::endl;
	 	for(int i = 0; i < levels.size(); i++) {
		  std::cout<< "Level " << i << std::endl;
		  
		  int size = levels[i].active_vertices.size();
		  std::cout <<"Active list :" << size << " elements" << std::endl;
		  vertex_node_t* pVertexNode  = levels[i].active_vertices.front();

  		  for(int n = 0; n < size; n++){
		     std::cout << pVertexNode->id  << " ";
		     pVertexNode = pVertexNode->next;
	          }		  
		  std::cout<< std::endl;

		  size = levels[i].inactive_vertices.size();
		  std::cout <<"Inactive list :" << size << " elements" << std::endl;
		  pVertexNode  = levels[i].inactive_vertices.front();

  		  for(int n = 0; n < size; n++){
		     std::cout << pVertexNode->id  << " ";
		     pVertexNode = pVertexNode->next;
	          }		  

		  std::cout<< std::endl;
		}	
	}


        struct level_t {
		linked_list<vertex_node_t> active_vertices;
		linked_list<vertex_node_t> inactive_vertices;
	};
	std::vector<level_t> levels;
	std::vector<vertex_node_t> _vertices;
	vecQueue<int> vertexQ;
	std::vector<std::vector<EdgeType>>&  adjList;
	std::vector<std::pair<edgeIterator, edgeIterator>> vCurrentEdges;
	edge_size_t numEdges;
	
	size_t numGRelabels, numPushes, numRelabels;
	int numVertices;
	int maxActiveHeight, minActiveHeight, maxHeight;
	int source;
	int sink;
};

struct SC
{
	std::vector<int> original;
	std::vector<int> positive;
	std::vector<int> negative;
};

struct SCRet
{
	std::vector<SC> S;
	compressed_matrix::CompressedMatrix<long long int> G;
};


Posiform BQPToPosiform(const compressed_matrix::CompressedMatrix<long long int>& Q, bool makeRandom = false)
{
	int numVariables = Q.numRows();
	Posiform ret;
	ret.numVars = numVariables;
	std::vector<long long int> q(numVariables);

	//q = diag(Q);
	for (int i = 0; i < numVariables; i++)
		q[i] = Q.get(i, i);

	//tmpQ = triu(Q,1);
	std::vector<int> tmpQRowOffsets(numVariables + 1);
	std::vector<int> tmpQColIndices; tmpQColIndices.reserve(Q.nnz());
	std::vector<long long int> tmpQValues; tmpQValues.reserve(Q.nnz());
	int index = 0;
	for (int i = 0; i < Q.numRows(); i++)
	{
		int start = Q.rowOffsets()[i];
		int end = Q.rowOffsets()[i + 1];
		tmpQRowOffsets[i] = index;
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = Q.colIndices()[j];
			if (r != c)
			{
				tmpQColIndices.push_back(c);
				tmpQValues.push_back(Q.values()[j]); //Q.values()[j] is Q(r, c)
				++index;
			}
		}
	}
	tmpQRowOffsets.back() = index; //number of non-zero elements

	compressed_matrix::CompressedMatrix<long long int> tmpQ(numVariables, numVariables, tmpQRowOffsets, tmpQColIndices, tmpQValues);


	if (!makeRandom) //make a deterministic posiform from Q
	{
		//cPrime = q + sum(tmpQ.*(tmpQ<0),2); and quadratic term
		std::vector<long long int> cPrime = q;
		for (int i = 0; i < tmpQ.numRows(); i++)
		{
			int start = tmpQ.rowOffsets()[i];
			int end = tmpQ.rowOffsets()[i + 1];
			for (int j = start; j < end; j++)
			{
				int r = i;
				int c = tmpQ.colIndices()[j];
				if (tmpQ.values()[j] < 0)
				{
					cPrime[i] += tmpQ.values()[j];
					ret.quadratic.push_back(std::make_pair(std::make_pair(r + 1, -(c + 1)), -tmpQ.values()[j])); //+1 makes it 1-based
				}
				else if (tmpQ.values()[j] > 0)
					ret.quadratic.push_back(std::make_pair(std::make_pair(r + 1, c + 1), tmpQ.values()[j])); //+1 makes it 1-based
			}
		}

		//constant term and linear term
		ret.cst = 0;
		for (int i = 0; i < numVariables; i++)
		{
			if (cPrime[i] < 0)
			{
				ret.cst += cPrime[i]; //constant term
				ret.linear.push_back(std::make_pair(-(i + 1), -cPrime[i])); //+1 makes it 1-based
			}
			else if (cPrime[i] > 0)
				ret.linear.push_back(std::make_pair(i + 1, cPrime[i])); //+1 makes it 1-based
		}
	}
	else
	{
		std::vector<long long int> linear = q;
		std::set<std::pair<int, int> > negPairs;
		for (int i = 0; i < tmpQ.numRows(); i++)
		{
			int start = tmpQ.rowOffsets()[i];
			int end = tmpQ.rowOffsets()[i + 1];
			for (int j = start; j < end; j++)
			{
				int r = i;
				int c = tmpQ.colIndices()[j];
				if (tmpQ.values()[j] < 0)
				{
					if (rand() % 2 == 0)
						negPairs.insert(std::make_pair(-(r + 1), c + 1)); //+1 makes it 1-based
					else
						negPairs.insert(std::make_pair(r + 1, -(c + 1))); //+1 makes it 1-based
				}
			}
		}

		for (int i = 0; i < tmpQ.numRows(); i++)
		{
			int start = tmpQ.rowOffsets()[i];
			int end = tmpQ.rowOffsets()[i + 1];
			for (int j = start; j < end; j++)
			{
				int r = i;
				int c = tmpQ.colIndices()[j];
				if (tmpQ.values()[j] > 0)
					ret.quadratic.push_back(std::make_pair(std::make_pair(r + 1, c + 1), tmpQ.values()[j]));
				else
				{
					if (negPairs.find(std::make_pair(r + 1, -(c + 1)))!=negPairs.end())
						ret.quadratic.push_back(std::make_pair(std::make_pair(r + 1, -(c + 1)), -tmpQ.values()[j]));
					else
						ret.quadratic.push_back(std::make_pair(std::make_pair(-(r + 1), c + 1), -tmpQ.values()[j]));
				}
			}
		}

		for (std::set<std::pair<int, int> >::iterator it = negPairs.begin(); it != negPairs.end(); ++it)
		{
			if (it->first > 0)
				linear[it->first - 1] += tmpQ(it->first - 1, -it->second - 1);
			else
				linear[it->second - 1] += tmpQ(-it->first - 1, it->second - 1);
		}

		ret.cst = 0;

		for (int i = 0; i < numVariables; i++)
		{
			if (linear[i] < 0)
			{
				ret.cst += linear[i];
				ret.linear.push_back(std::make_pair(-(i + 1), -linear[i]));
			}
			else if (linear[i] > 0)
				ret.linear.push_back(std::make_pair(i + 1, linear[i]));
		}
	}

	return ret;
}

compressed_matrix::CompressedMatrix<long long int> posiformToImplicationNetwork_1(const Posiform& p)
{
	int n = p.numVars;
	int numVertices = 2 * n + 2;
	int source = 0;
	int sink = n + 1;

	//n = p.numVars
	//0: source;
	//1 to n: x_1 to x_n
	//n+1: sink
	//n+2 to 2*n+1: \overline_x_1 to \overline_x_n
	std::map<std::pair<int, int>, long long int> m;
	for (int i = 0; i < p.linear.size(); i++)
	{
		int v = p.linear[i].first;

		long long int capacity = p.linear[i].second; // originally p.linear[i].second/2
		if (v > 0)
		{
			m[std::make_pair(source, v + n + 1)] = capacity;
			m[std::make_pair(v, sink)] = capacity;
		}
		else
		{
			v = std::abs(v);
			m[std::make_pair(source, v)] = capacity;
			m[std::make_pair(v + n + 1, sink)] = capacity;
		}
	}

	for (int i = 0; i < p.quadratic.size(); i++)
	{
		int v_1 = p.quadratic[i].first.first;
		int v_2 = p.quadratic[i].first.second;

		long long int capacity = p.quadratic[i].second; // originally p.quadratic[i].second/2
		if (v_1 < 0)
		{
			v_1 = std::abs(v_1);
			m[std::make_pair(v_1 + n + 1, v_2 + n + 1)] = capacity;
			m[std::make_pair(v_2, v_1)] = capacity;
		}
		else if (v_2 < 0)
		{
			v_2 = std::abs(v_2);
			m[std::make_pair(v_1, v_2)] = capacity;
			m[std::make_pair(v_2 + n + 1, v_1 + n + 1)] = capacity;
		}
		else
		{
			m[std::make_pair(v_1, v_2 + n + 1)] = capacity;
			m[std::make_pair(v_2, v_1 + n + 1)] = capacity;
		}
	}

	return compressed_matrix::CompressedMatrix<long long int>(numVertices, numVertices, m);
}

compressed_matrix::CompressedMatrix<long long int> posiformToImplicationNetwork_2(const Posiform& p)
{
	int n = p.numVars;
	int numVertices = 2 * n + 2;
	int source = n;
	int sink = 2 * n + 1;

	//n = p.numVars
	//0 to n-1: x_1 to x_n
	//n: source;
	//n+1 to 2*n: \overline_x_1 to \overline_x_n
	//2*n+1: sink
	std::map<std::pair<int, int>, long long int> m;
	for (int i = 0; i < p.linear.size(); i++)
	{
		int v = p.linear[i].first;
		long long int capacity = p.linear[i].second; // originally p.linear[i].second/2
		if (v > 0)
		{
			m[std::make_pair(source, v + n)] = capacity;
			m[std::make_pair(v - 1, sink)] = capacity;
		}
		else
		{
			v = std::abs(v);
			m[std::make_pair(source, v - 1)] = capacity;
			m[std::make_pair(v + n, sink)] = capacity;
		}
	}

	for (int i = 0; i < p.quadratic.size(); i++)
	{
		int v_1 = p.quadratic[i].first.first;
		int v_2 = p.quadratic[i].first.second;
		long long int capacity = p.quadratic[i].second; // originally p.quadratic[i].second/2
		if (v_1 < 0)
		{
			v_1 = std::abs(v_1);
			m[std::make_pair(v_1 + n, v_2 + n)] = capacity;
			m[std::make_pair(v_2 - 1, v_1 - 1)] = capacity;
		}
		else if (v_2 < 0)
		{
			v_2 = std::abs(v_2);
			m[std::make_pair(v_1 - 1, v_2 - 1)] = capacity;
			m[std::make_pair(v_2 + n, v_1 + n)] = capacity;
		}
		else
		{
			m[std::make_pair(v_1 - 1, v_2 + n)] = capacity;
			m[std::make_pair(v_2 - 1, v_1 + n)] = capacity;
		}
	}

	return compressed_matrix::CompressedMatrix<long long int>(numVertices, numVertices, m);
}

//this version returns R instead of F
compressed_matrix::CompressedMatrix<long long int> maxFlow(const compressed_matrix::CompressedMatrix<long long int>& A)
{
	int numVertices = A.numRows();
	int numVariables = numVertices / 2 - 1;

	clock_t curr_1 = clock();
        clock_t curr_2;

	using namespace boost;

	typedef adjacency_list_traits<vecS, vecS, directedS> Traits;
	typedef adjacency_list<vecS, vecS, directedS, property<vertex_name_t, std::string>, property<edge_capacity_t, long long int, property<edge_residual_capacity_t, long long int, property<edge_reverse_t, Traits::edge_descriptor> > > > Graph; //for edge capacity is long long int
	//typedef adjacency_list<vecS, vecS, directedS, property<vertex_name_t, std::string>, property<edge_capacity_t, double, property<edge_residual_capacity_t, double, property<edge_reverse_t, Traits::edge_descriptor> > > > Graph; //for edge capacity is double

	Graph g;

	property_map<Graph, edge_capacity_t>::type capacity = get(edge_capacity, g);
	property_map<Graph, edge_reverse_t>::type reverse_edge = get(edge_reverse, g);
	property_map<Graph, edge_residual_capacity_t>::type residual_capacity = get(edge_residual_capacity, g);

	std::vector<Traits::vertex_descriptor> verts(numVertices);
	for (int i = 0; i < numVertices; ++i)
		verts[i] = add_vertex(g);

	Traits::vertex_descriptor s = verts[0];
	Traits::vertex_descriptor t = verts[numVariables + 1];

	for (int i = 0; i < A.numRows(); i++)
	{
		int start = A.rowOffsets()[i];
		int end = A.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = A.colIndices()[j];
			long long int cap = A.values()[j];

			Traits::edge_descriptor e1, e2;
			bool in1, in2;
			boost::tie(e1, in1) = add_edge(verts[r], verts[c], g);
			boost::tie(e2, in2) = add_edge(verts[c], verts[r], g);

			capacity[e1] = cap;
			capacity[e2] = 0;
			reverse_edge[e1] = e2;
			reverse_edge[e2] = e1;
		}
	}


	//curr_2 = clock();
	//mexPrintf("inside maxFlow int version: Time elapsed_constructing_graph_for_max_flow: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

   	long long int flowValueBoost = push_relabel_max_flow(g, s, t);
	std::cout <<"Flow Value from Boost : " << flowValueBoost << std::endl;

	curr_2 = clock();
	printf("inside maxFlow int version: Time elapsed_for_boost_max_flow: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	std::vector<int> RRowOffsets(numVertices+1);
	std::vector<int> RColIndices;
	std::vector<long long int> RValues;
	int offset = 0;
	int currRow = 0;
	graph_traits<Graph>::vertex_iterator u_iter, u_end;
	graph_traits<Graph>::out_edge_iterator ei, e_end;
	for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
		for (boost::tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
		{
			if (capacity[*ei] > 0)
			{
				if (currRow <= *u_iter)
				{
					for (int i = currRow; i <= *u_iter; i++)
						RRowOffsets[i] = offset;
					currRow = static_cast<int>((*u_iter) + 1);
				}

				RColIndices.push_back(static_cast<int>(target(*ei, g)));
				RValues.push_back(residual_capacity[*ei]);
				++offset;
			}
		}
	for (int i = currRow; i < RRowOffsets.size(); i++)
		RRowOffsets[i] = offset;

	compressed_matrix::CompressedMatrix<long long int> R(numVertices, numVertices, RRowOffsets, RColIndices, RValues); //residual

	//curr_2 = clock();
	//mexPrintf("inside maxFlow int version: Time elapsed_get_R: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	std::map<std::pair<int, int>, long long int> rm = compressed_matrix::compressedMatrixToMap(R);
	std::map<std::pair<int, int>, long long int> rmMissing;

	//curr_2 = clock();
	//mexPrintf("inside maxFlow int version: Time elapsed_copy_R_to_map: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	for (int i = 0; i < R.numRows(); i++)
	{
		int start = R.rowOffsets()[i];
		int end = R.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = R.colIndices()[j];
			long long int Arc = A.get(r, c);
			long long int Rcr = R.get(c, r);
			if (Arc != 0 && R.values()[j] + Rcr - Arc != 0) //it->second is: R(r, c); here it means R(r, c)+R(c, r)!=A(r, c), so R(c, r) is missing
				rmMissing.insert(std::make_pair(std::make_pair(c, r), Arc - R.values()[j]));
		}
	}

	for (std::map<std::pair<int, int>, long long int>::iterator it = rmMissing.begin(), end = rmMissing.end(); it != end; ++it)
		rm.insert(*it);

	//curr_2 = clock();
	//mexPrintf("inside maxFlow int version: Time elapsed_get_missing_elements_in_R: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	return compressed_matrix::CompressedMatrix<long long int>(numVertices, numVertices, rm);
}

std::vector<int> bfs_for_method_2(const compressed_matrix::CompressedMatrix<long long int>& g, int s)
{
	std::queue<int> q;
	q.push(s);
	int numVertices = g.numRows();
	std::vector<int> visited(numVertices, -1);
	visited[s] = 1;

	while (!q.empty())
	{
		int curr = q.front(); q.pop();
		int start = g.rowOffsets()[curr];
		int end = g.rowOffsets()[curr + 1];
		for (int j = start; j < end; j++)
		{
			int c = g.colIndices()[j];
			if (g.values()[j] != 0 && visited[c] == -1)
			{
				visited[c] = 1;
				q.push(c);
			}
		}
	}

	return visited;
}

compressed_matrix::CompressedMatrix<long long int> makeResidualSymmetric(const compressed_matrix::CompressedMatrix<long long int>& R)
{
	//clock_t curr_1 = clock();
	//clock_t curr_2;

	int numVertices = R.numRows();
	int numVariables = numVertices / 2 - 1;
	std::map<std::pair<int, int>, long long int> MsymR;
	for (int i = 0; i < R.numRows(); i++)
	{
		int start = R.rowOffsets()[i];
		int end = R.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = R.colIndices()[j];
			if (r != c)
			{
				MsymR[std::make_pair(r, c)] += R.values()[j];
				int compR;
				if (r <= numVariables)
					compR = r + (numVariables + 1);
				else
					compR = r - (numVariables + 1);
				int compC;
				if (c <= numVariables)
					compC = c + (numVariables + 1);
				else
					compC = c - (numVariables + 1);
				MsymR[std::make_pair(compC, compR)] += R.values()[j];
			}
		}
	}

	//curr_2 = clock();
	//mexPrintf("inside makeResidualSym: Time elapsed_get_MsymR: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	//remove extra 0s in MsymR, required because there will possibly be 0 in MsymR and the 0 will add non-exist edges in stronglyConnectedComponents()
	//cause the connected components results incorrct !!!
	std::map<std::pair<int, int>, long long int> MsymRWithoutZero;
	for (std::map<std::pair<int, int>, long long int>::const_iterator it = MsymR.begin(), end = MsymR.end(); it != end; ++it)
	{
		if (it->second != 0 && it->first.first != numVariables + 1 && it->first.second != 0) //add clearing R here !!!
			MsymRWithoutZero.insert(*it);
	}

	//curr_2 = clock();
	//mexPrintf("inside makeResidualSym: Time elapsed_remove_zero_in_MsymR: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	return compressed_matrix::CompressedMatrix<long long int>(numVertices, numVertices, MsymRWithoutZero);
}

SCRet stronglyConnectedComponents(const compressed_matrix::CompressedMatrix<long long int>& R)
{
	//clock_t curr_1 = clock();
	//clock_t curr_2;

	int numVertices = R.numRows();
	int numVariables = numVertices / 2 - 1;

	using namespace boost;

	typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
	typedef adjacency_list < vecS, vecS, directedS,
		property < vertex_name_t, std::string,
		property < vertex_index_t, long,
		property < vertex_color_t, boost::default_color_type,
		property < vertex_distance_t, long,
		property < vertex_predecessor_t, Traits::edge_descriptor > > > > >
	> Graph;

	typedef graph_traits<Graph>::vertex_descriptor Vertex;

	Graph G;

	std::vector<Traits::vertex_descriptor> root(numVertices);
	for (int i = 0; i < numVertices; ++i)
		root[i] = add_vertex(G);

	std::vector<int> component(num_vertices(G)), discover_time(num_vertices(G));
	std::vector<default_color_type> color(num_vertices(G));

	for (int i = 0; i < R.numRows(); i++)
	{
		int start = R.rowOffsets()[i];
		int end = R.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = R.colIndices()[j];
			add_edge(root[r], root[c], G);
		}
	}

	//curr_2 = clock();
	//mexPrintf("Time elapsed_constructing_graph_for_strongly_connected_components_graph: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	//int num = strong_components(G, &component[0], root_map(&root[0]).color_map(&color[0]).discover_time_map(&discover_time[0]));
	int num = strong_components(G, boost::make_iterator_property_map(component.begin(), boost::get(boost::vertex_index, G)),
			                                root_map(boost::make_iterator_property_map(root.begin(), boost::get(boost::vertex_index, G))).color_map(boost::make_iterator_property_map(color.begin(), boost::get(boost::vertex_index, G))).discover_time_map(boost::make_iterator_property_map(discover_time.begin(), boost::get(boost::vertex_index, G))));

	//curr_2 = clock();
	//mexPrintf("Time elapsed_boost_strongly_connected_components_graph: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	//here original's content is from 0 to 2*n+1, 0: source, 1 to n, variables, n+1: sink, n+2 to 2*n+1: bar variables
	//positive: from 1 to n (1-based variables index)
	//negative: from 1 to n (1-based variables index)
	std::vector<SC> vsc(num);
	for (int i = 0; i < component.size(); i++)
	{
		int whichComponent = component[i];
		vsc[whichComponent].original.push_back(i);
		if (vsc[whichComponent].original.back() <= numVariables)
			vsc[whichComponent].positive.push_back(vsc[whichComponent].original.back());
		else
			vsc[whichComponent].negative.push_back(vsc[whichComponent].original.back() - (numVariables + 1));
	}

	std::map<std::pair<int, int>, long long int> M;

	//speed up the program a lot !!!
	//build the graph G for strongly connected components
	for (int i = 0; i < R.numRows(); i++)
	{
		int start = R.rowOffsets()[i];
		int end = R.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = R.colIndices()[j];
			int sc1 = component[r];
			int sc2 = component[c];
			if (sc1 != sc2 && R.values()[j] > 0 && M.find(std::make_pair(sc1, sc2)) == M.end())
				M.insert(std::make_pair(std::make_pair(sc1, sc2), 1));
		}
	}

	//curr_2 = clock();
	//mexPrintf("Time elapsed_after_boost_scc: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	SCRet ret;
	ret.S = vsc;
	ret.G = compressed_matrix::CompressedMatrix<long long int>(num, num, M);

	return ret;
}

std::vector<int> classifyStronglyConnectedComponents(const std::vector<SC>& S)
{
	std::vector<int> ret(S.size()); //0: self-complement, 1: non-self-complement

	for (int i = 0; i < S.size(); i++)
	{
		bool flag = true;
		std::set<int> pos(S[i].positive.begin(), S[i].positive.end());
		for (int j = 0; j < S[i].negative.size(); j++)
		{
			if (pos.find(S[i].negative[j]) != pos.end())
			{
				flag = false;
				break;
			}
		}

		if (flag)
			ret[i] = 1;
		else
			ret[i] = 0;
	}

	return ret;
}

void shortestPath(int src, const compressed_matrix::CompressedMatrix<long long int>& g, std::vector<int>& visited)
{
	std::queue<int> q;
	q.push(src);
	visited.resize(g.numRows(), 0);
	visited[src] = 1;

	while (!q.empty())
	{
		int curr = q.front(); q.pop();

		int start = g.rowOffsets()[curr];
		int end = g.rowOffsets()[curr + 1];
		for (int j = start; j < end; j++)
		{
			int c = g.colIndices()[j];
			if (g.values()[j] != 0 && !visited[c])
			{
				visited[c] = 2;
				q.push(c);
			}
		}
	}
}

std::vector<std::pair<int, int> > fixVariables(const std::vector<int>& classifiedSC, const SCRet& scRet, int numVariables)
{
	std::vector<std::pair<int, int> > ret;

	int x0Location;
	int x0BarLocation;
	std::vector<int> I;

	for (int i = 0; i < classifiedSC.size(); i++)
	{
		if (classifiedSC[i] == 1) //1: non self-complement
		{
			if (scRet.S[i].original[0] == 0)
				x0Location = i;
			else if (scRet.S[i].original[0] == numVariables + 1)
				x0BarLocation = i;
			else
				I.push_back(i);
		}
	}

	//calculate if there is a path from x0 to other nodes
	std::vector<int> visited;
	shortestPath(x0Location, scRet.G, visited); //speed up the program

	//find complement pairs
	std::vector<int> positive(numVariables + 1, -1);
	std::vector<int> negative(numVariables + 1, -1);
	for (int i = 0; i < I.size(); i++)
	{
		for (int k = 0; k < scRet.S[I[i]].positive.size(); k++)
			positive[scRet.S[I[i]].positive[k]] = I[i];
		for (int k = 0; k < scRet.S[I[i]].negative.size(); k++)
			negative[scRet.S[I[i]].negative[k]] = I[i];
	}

	std::map<int, int> complementPairs;
	for (int i = 0; i < positive.size(); i++)
	{
		if (positive[i] != -1 && negative[i] != -1)
		{
			complementPairs[positive[i]] = negative[i];
			complementPairs[negative[i]] = positive[i];
		}
	}

	std::queue<int> q;

	std::map<std::pair<int, int>, long long int> GTransMap;
	for (int i = 0; i < scRet.G.numRows(); i++)
	{
		int start = scRet.G.rowOffsets()[i];
		int end = scRet.G.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
			GTransMap.insert(std::make_pair(std::make_pair(scRet.G.colIndices()[j], i), scRet.G.values()[j]));
	}
	compressed_matrix::CompressedMatrix<long long int> GTrans(scRet.G.numRows(), scRet.G.numRows(), GTransMap);

	std::vector<int> outDegrees(scRet.G.numRows(), -1);
	for (int i = 0; i < classifiedSC.size(); i++)
	{
		if (classifiedSC[i] == 1) //non self-complement components
		{
			outDegrees[i] = scRet.G.rowOffsets()[i + 1] - scRet.G.rowOffsets()[i];
			if (outDegrees[i] == 0 && i != x0Location && i != x0BarLocation) //need to push the original outdegree 0 nodes into q !!!
				q.push(i);
		}
	}
	outDegrees[x0Location] = -1;
	outDegrees[x0BarLocation] = -1;

	for (int i = 0; i < visited.size(); i++)
	{
		if (visited[i] == 2) //exclude x0
		{
			for (int k = 0; k < scRet.S[i].positive.size(); k++)
				ret.push_back(std::make_pair(scRet.S[i].positive[k], 1));

			for (int k = 0; k < scRet.S[i].negative.size(); k++)
				ret.push_back(std::make_pair(scRet.S[i].negative[k], 0));

			int complement = complementPairs[i];

			outDegrees[i] = -1;
			outDegrees[complement] = -1;
		}
	}


	//decrease the outdegrees of node which has outgoing edges to i and to complement
	//push node which has 0 outdegrees into the queue
	for (int i = 0; i < visited.size(); i++)
	{
		if (visited[i] == 2) //exclude x0
		{
			int start = GTrans.rowOffsets()[i];
			int end = GTrans.rowOffsets()[i + 1];
			for (int k = start; k < end; k++)
			{
				if (outDegrees[GTrans.colIndices()[k]] > 0)
				{
					--outDegrees[GTrans.colIndices()[k]];
					if (outDegrees[GTrans.colIndices()[k]] == 0)
						q.push(GTrans.colIndices()[k]);
				}
			}

			int complement = complementPairs[i];
			start = GTrans.rowOffsets()[complement];
			end = GTrans.rowOffsets()[complement + 1];
			for (int k = start; k < end; k++)
			{
				if (outDegrees[GTrans.colIndices()[k]] > 0)
				{
					--outDegrees[GTrans.colIndices()[k]];
					if (outDegrees[GTrans.colIndices()[k]] == 0)
						q.push(GTrans.colIndices()[k]);
				}
			}
		}
	}


	while (!q.empty())
	{
		int curr = q.front(); q.pop();
		if (outDegrees[curr] == 0)
		{
			outDegrees[curr] = -1;
			int complement = complementPairs[curr];
			outDegrees[complement] = -1;

			//fixed all variables in component curr
			for (int k = 0; k < scRet.S[curr].positive.size(); k++)
				ret.push_back(std::make_pair(scRet.S[curr].positive[k], 1));

			for (int k = 0; k < scRet.S[curr].negative.size(); k++)
				ret.push_back(std::make_pair(scRet.S[curr].negative[k], 0));

			//decrease the outdegrees of node which has outgoing edges to i and to complement
			//push node which has 0 outdegrees into the queue
			int start = GTrans.rowOffsets()[curr];
			int end = GTrans.rowOffsets()[curr + 1];
			for (int k = start; k < end; k++)
			{
				if (outDegrees[GTrans.colIndices()[k]] > 0)
				{
					--outDegrees[GTrans.colIndices()[k]];
					if (outDegrees[GTrans.colIndices()[k]] == 0)
						q.push(GTrans.colIndices()[k]);
				}
			}

			start = GTrans.rowOffsets()[complement];
			end = GTrans.rowOffsets()[complement + 1];
			for (int k = start; k < end; k++)
			{
				if (outDegrees[GTrans.colIndices()[k]] > 0)
				{
					--outDegrees[GTrans.colIndices()[k]];
					if (outDegrees[GTrans.colIndices()[k]] == 0)
						q.push(GTrans.colIndices()[k]);
				}
			}
		}
	}

	std::sort(ret.begin(), ret.end(), compClass());

	return ret;
}

compressed_matrix::CompressedMatrix<double> computeNewQAndOffset(const compressed_matrix::CompressedMatrix<double>& Q, const std::vector<std::pair<int, int> >& fixed, double& offset)
{
	//fixed now is 1-based
	if (!fixed.empty())
	{
		int numVariables = Q.numRows();

		std::map<int, int> mI;
		for (int i = 0; i < fixed.size(); i++)
			mI[fixed[i].first - 1] = i; //-1 to make it 0-based

		std::vector<int> J;
		std::map<int, int> mJ;
		int cnt = 0;
		for (int i = 0; i < numVariables; i++)
		{
			if (mI.find(i) == mI.end())
			{
				J.push_back(i);
				mJ[i] = cnt++;
			}
		}

		std::map<std::pair<int, int>, double> QII;
		std::map<std::pair<int, int>, double> QJJ;
		std::map<std::pair<int, int>, double> QIJ;
		std::map<std::pair<int, int>, double> QJI;

		for (int i = 0; i < Q.numRows(); i++)
		{
			int start = Q.rowOffsets()[i];
			int end = Q.rowOffsets()[i + 1];
			for (int j = start; j < end; j++)
			{
				int r = i;
				int c = Q.colIndices()[j];
				if (mI.find(r) != mI.end() && mI.find(c) != mI.end())
					QII[std::make_pair(mI[r], mI[c])] = Q.values()[j];
				if (mJ.find(r) != mJ.end() && mJ.find(c) != mJ.end())
					QJJ[std::make_pair(mJ[r], mJ[c])] = Q.values()[j];
				if (mI.find(r) != mI.end() && mJ.find(c) != mJ.end())
					QIJ[std::make_pair(mI[r], mJ[c])] = Q.values()[j];
				if (mJ.find(r) != mJ.end() && mI.find(c) != mI.end())
					QJI[std::make_pair(mJ[r], mI[c])] = Q.values()[j];
			}
		}

		//off_set = x0'*Q_II*x0;
		offset = 0;
		for (std::map<std::pair<int, int>, double>::const_iterator it = QII.begin(), end = QII.end(); it != end; ++it)
		{
			int r = it->first.first;
			int c = it->first.second;
			offset += fixed[r].second * it->second * fixed[c].second;
		}

		//x0'*Q_IJ
		std::vector<double> tmp2(numVariables - fixed.size(), 0);
		for (std::map<std::pair<int, int>, double>::const_iterator it = QIJ.begin(), end = QIJ.end(); it != end; ++it)
		{
			int r = it->first.first;
			int c = it->first.second;
			tmp2[c] += it->second * fixed[r].second;
		}

		//Q_JI*x0
		std::vector<double> tmp3(numVariables - fixed.size(), 0);
		for (std::map<std::pair<int, int>, double>::const_iterator it = QJI.begin(), end = QJI.end(); it != end; ++it)
		{
			int r = it->first.first;
			int c = it->first.second;
			tmp3[r] += it->second * fixed[c].second;
		}

		for (int i = 0; i < (int)tmp3.size(); i++)
			tmp2[i] += tmp3[i];

		cnt = 0;
		for (int i = 0; i < (int)mJ.size() * (int)mJ.size(); i += numVariables - static_cast<int>(fixed.size()) + 1)
		{
			int r = i / static_cast<int>(mJ.size());
			int c = i % mJ.size();
			QJJ[std::make_pair(r, c)] += tmp2[cnt++];
		}

		std::map<std::pair<int, int>, double> mfixedQ;
		for (std::map<std::pair<int, int>, double>::const_iterator it = QJJ.begin(), end = QJJ.end(); it != end; ++it)
		{
			int r = it->first.first;
			int c = it->first.second;
			if (it->second != 0)
				mfixedQ[std::make_pair(J[r], J[c])] = it->second;
		}

		return compressed_matrix::CompressedMatrix<double>(numVariables, numVariables, mfixedQ);
	}
	else
	{
		offset = 0;
		return Q;
	}
}

std::vector<std::pair<int, int> > applyImplication(const compressed_matrix::CompressedMatrix<long long int>& A)
{
	int numVertices = A.numRows();
	int numVariables = numVertices / 2 - 1;

	//debuging only
	//clock_t curr_1 = clock();
	//clock_t curr_2;

	using namespace boost;

	typedef adjacency_list_traits<vecS, vecS, directedS> Traits;
	typedef adjacency_list<vecS, vecS, directedS, property<vertex_name_t, std::string>, property<edge_capacity_t, long long int, property<edge_residual_capacity_t, long long int, property<edge_reverse_t, Traits::edge_descriptor> > > > Graph; //for edge capacity is long long int
	//typedef adjacency_list<vecS, vecS, directedS, property<vertex_name_t, std::string>, property<edge_capacity_t, double, property<edge_residual_capacity_t, double, property<edge_reverse_t, Traits::edge_descriptor> > > > Graph; //for edge capacity is double

	Graph g;

	property_map<Graph, edge_capacity_t>::type capacity = get(edge_capacity, g);
	property_map<Graph, edge_reverse_t>::type reverse_edge = get(edge_reverse, g);
	property_map<Graph, edge_residual_capacity_t>::type residual_capacity = get(edge_residual_capacity, g);

	std::vector<Traits::vertex_descriptor> verts(numVertices);
	for (int i = 0; i < numVertices; ++i)
		verts[i] = add_vertex(g);

	Traits::vertex_descriptor s = verts[numVariables];
	Traits::vertex_descriptor t = verts[2 * numVariables + 1];

	for (int i = 0; i < A.numRows(); i++)
	{
		int start = A.rowOffsets()[i];
		int end = A.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = A.colIndices()[j];
			long long int cap = A.values()[j];

			Traits::edge_descriptor e1, e2;
			bool in1, in2;
			boost::tie(e1, in1) = add_edge(verts[r], verts[c], g);
			boost::tie(e2, in2) = add_edge(verts[c], verts[r], g);

			capacity[e1] = cap;
			capacity[e2] = 0;
			reverse_edge[e1] = e2;
			reverse_edge[e2] = e1;
		}
	}

	//curr_2 = clock();
	//mexPrintf("inside applyImplication int version: Time elapsed_building_graph: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	long long int flowValueBoost =  push_relabel_max_flow(g, s, t);
	std::cout <<"Flow Value from Boost : " << flowValueBoost << std::endl;

	//curr_2 = clock();
	//mexPrintf("inside applyImplication int version: Time elapsed_max_flow: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	//residual
	std::vector<int> RRowOffsets(numVertices+1);
	std::vector<int> RColIndices;
	std::vector<long long int> RValues;
	//flow
	std::vector<int> FRowOffsets(numVertices+1);
	std::vector<int> FColIndices;
	std::vector<long long int> FValues;

	int offset = 0;
	int currRow = 0;
	graph_traits<Graph>::vertex_iterator u_iter, u_end;
	graph_traits<Graph>::out_edge_iterator ei, e_end;
	for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
		for (boost::tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
		{
			if (capacity[*ei] > 0)
			{
				if (currRow <= *u_iter)
				{
					for (int i = currRow; i <= *u_iter; i++)
					{
						RRowOffsets[i] = offset;
						FRowOffsets[i] = offset;
					}
					currRow = static_cast<int>((*u_iter) + 1);
				}
				RColIndices.push_back(static_cast<int>(target(*ei, g)));
				RValues.push_back(residual_capacity[*ei]);
				FColIndices.push_back(static_cast<int>(target(*ei, g)));
				FValues.push_back(capacity[*ei] - residual_capacity[*ei]);
				++offset;
			}
		}
	for (int i = currRow; i < RRowOffsets.size(); i++)
		RRowOffsets[i] = offset;
	for (int i = currRow; i < FRowOffsets.size(); i++)
		FRowOffsets[i] = offset;

	compressed_matrix::CompressedMatrix<long long int> R(numVertices, numVertices, RRowOffsets, RColIndices, RValues); //residual
	compressed_matrix::CompressedMatrix<long long int> F(numVertices, numVertices, FRowOffsets, FColIndices, FValues); //flow

	//curr_2 = clock();
	//mexPrintf("inside applyImplication int version: Time elapsed_get_RF: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	//resid = (A>0).*R + F'.*(A'>0);
	std::map<std::pair<int, int>, long long int> residM;
	for (int i = 0; i < A.numRows(); i++)
	{
		int start = A.rowOffsets()[i];
		int end = A.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = A.colIndices()[j];
			long long int RValue = R.get(r, c);
			if (A.values()[j] > 0 && RValue != 0 && r != 2 * numVariables + 1 && c != numVariables)
				residM[std::make_pair(r, c)] += RValue;

			long long int FTValue = F.get(r, c); //not F(c, r) !!!
			if (A.values()[j] > 0 && FTValue != 0 && c != 2 * numVariables + 1 && r != numVariables)
				residM[std::make_pair(c, r)] += FTValue;
		}
	}

	compressed_matrix::CompressedMatrix<long long int> resid(numVertices, numVertices, residM);

	//curr_2 = clock();
	//mexPrintf("inside applyImplication int version: Time elapsed_finally_get_resid: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	std::vector<int> forced = bfs_for_method_2(resid, numVariables);

	//curr_2 = clock();
	//mexPrintf("inside applyImplication int version: Time elapsed_bfs_for_method_2: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	std::vector<std::pair<int, int> > fixed;
	for (int i = 0; i < forced.size(); i++)
	{
		if (forced[i] > 0)
		{
			if (i <= numVariables - 1)
				fixed.push_back(std::make_pair(i + 1, 1)); //1-based
			else if (i >= numVariables + 1 && i <= 2 * numVariables)
				fixed.push_back(std::make_pair(i - numVariables, 0)); //1-based
		}
	}

	std::sort(fixed.begin(), fixed.end(), compClass());

	//curr_2 = clock();
	//mexPrintf("inside applyImplication int version: Time elapsed_fix_vars_and_sorting: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	return fixed;
}

} // anonymous namespace



/*
template <class BQM>
void BQMToPosiform(const BQM& bqm, Posiform& posiform, std::vector<bool>& usedVars, bool randomize = false)
{
      int numVars = bqm.num_variables();
      std::vector<bool> isUsed(numVars, false);

      // Temporary code to maintain compatibility with legacy code
      std::map<std::pair<int, int>, double> QMap;
      for(int i = 0; i < numVars; i++){
	auto span = bqm.neighborhood(i);
	auto linear = bqm.linear(i);
        QMap[{i,i}] = linear;
        auto it = std::lower_bound(span.first, span.second, i+1, dimod::utils::comp_v<V,B>);
        if(it != span.second) {
          for(auto itEnd = span.second; it != itEnd; it++) {
             QMap[{i,it->first}] = it->second;
	  }
	}
      }
}
*/


namespace fix_variables_
{

struct FixVariablesResult
{
	std::vector<std::pair<int,  int> > fixedVars; //1-based
	compressed_matrix::CompressedMatrix<double> newQ; //0-based
	double offset;
};

FixVariablesResult fixQuboVariables(const compressed_matrix::CompressedMatrix<double>& Q, int method)
{
	static int di=0;
        printf("Starting function %d  ******\n", ++di);

	//Q needs to be a square matrix
	if (Q.numRows() != Q.numCols())
        throw std::invalid_argument("Q's size is not correct.");

	if (!(method == 1 || method == 2))
        throw std::invalid_argument("method must be an integer of 1 or 2.");

	FixVariablesResult ret;

	//check if Q is empty
	if (Q.numRows() == 0 || Q.numCols() == 0)
	{
		ret.offset = 0;
		return ret;
	}

	int numVariables = Q.numRows();

	clock_t curr_1 = clock();
	clock_t curr_2;

	//uTriQ = triu(Q) + tril(Q,-1)'; //make upper triangular
	std::map<std::pair<int, int>, double> uTriQMap;
	std::set<int> usedVariables;
	for (int i = 0; i < Q.numRows(); i++)
	{
		int start = Q.rowOffsets()[i];
		int end = Q.rowOffsets()[i + 1];
		for (int j = start; j < end; j++)
		{
			int r = i;
			int c = Q.colIndices()[j];
			uTriQMap[std::make_pair(std::min(r, c), std::max(r, c))] += Q.values()[j];
			usedVariables.insert(r);
			usedVariables.insert(c);
		}
	}

	compressed_matrix::CompressedMatrix<double> uTriQ(numVariables, numVariables, uTriQMap);

	double maxAbsValue = 0;

	if (!uTriQ.values().empty())
		maxAbsValue = std::fabs(*std::max_element(uTriQ.values().begin(), uTriQ.values().end(), compareAbs));

	std::cout << "Max abs value " << maxAbsValue << std::endl;

	double ratio = 1;

	if (maxAbsValue != 0)
		ratio = static_cast<double>(std::numeric_limits<long long int>::max()) / maxAbsValue;

	ratio /= static_cast<double>(1LL << 10);

	if (ratio < 1)
		ratio = 1;

	std::map<std::pair<int, int>, long long int> uTriQMapLLI;

	for (std::map<std::pair<int, int>, double>::iterator it = uTriQMap.begin(), end = uTriQMap.end(); it != end; ++it)
		uTriQMapLLI[it->first] = static_cast<long long int>(it->second * ratio);

	compressed_matrix::CompressedMatrix<long long int> uTriQLLI(numVariables, numVariables, uTriQMapLLI);
	curr_2 = clock();
	printf("Time elapsed_make upper triangular: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	curr_1 = curr_2;

	Posiform p = BQPToPosiform(uTriQLLI);
	printPosiform(p);
	curr_2 = clock();
	printf("Time elapsed_BQPToPosiform: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	curr_1 = curr_2;
	printf("Method : %d \n " , method );
	if (method == 1)
	{
		compressed_matrix::CompressedMatrix<long long int> A = posiformToImplicationNetwork_1(p);

	curr_2 = clock();
	 printf("Time elapsed_posiformToImplicationNetwork_1: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	 curr_1 = curr_2;

		compressed_matrix::CompressedMatrix<long long int> R = maxFlow(A);  //use this

		curr_2 = clock();
	    printf("Time elapsed_maxFlow: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	    curr_1 = curr_2;

		compressed_matrix::CompressedMatrix<long long int> symR = makeResidualSymmetric(R); //use this

	curr_2 = clock();
	    printf("Time elapsed_makeRSym: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	    curr_1 = curr_2;

		//add clearing R in makeResidualSym, so here just use symR directly !!!
		SCRet scRet = stronglyConnectedComponents(symR);


		curr_2 = clock();
	   printf("Time elapsed_stronglyConnectedComponents: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	    curr_1 = curr_2;

		std::vector<int> classifiedSC = classifyStronglyConnectedComponents(scRet.S);

		curr_2 = clock();
	    printf("Time elapsed_classifyStrongComponents: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	    curr_1 = curr_2;

		ret.fixedVars = fixVariables(classifiedSC, scRet, numVariables);

		curr_2 = clock();
	    printf("Time elapsed_fixVarsUsingOutDegree: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	    curr_1 = curr_2;
	}
	else if (method == 2)
	{
		compressed_matrix::CompressedMatrix<long long int> A = posiformToImplicationNetwork_2(p);

		//curr_2 = clock();
	    //mexPrintf("Time elapsed_posiformToImplicationNetwork_2: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	    //curr_1 = curr_2;

		ret.fixedVars = applyImplication(A);

		//curr_2 = clock();
	    //mexPrintf("Time elapsed_applyImplication: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	    //curr_1 = curr_2;
	}

	if (ret.fixedVars.size() > numVariables)
        throw std::logic_error("ret.fixedVars has wrong size.");

	// remove unused variables from ret.fixedVars
	std::vector<std::pair<int, int> > updatedFixedVars;
	for (int i = 0; i < ret.fixedVars.size(); ++i)
	{
		if (usedVariables.find(ret.fixedVars[i].first - 1) != usedVariables.end()) // -1 to make it 0-based since usedVariables is 0-based
			updatedFixedVars.push_back(ret.fixedVars[i]);
	}

	ret.fixedVars = updatedFixedVars;

	ret.newQ = computeNewQAndOffset(uTriQ, ret.fixedVars, ret.offset);

	//curr_2 = clock();
	//mexPrintf("Time elapsed_computeNewQAndOffset: %f\n", ((double)curr_2 - curr_1) / CLOCKS_PER_SEC);
	//curr_1 = curr_2;

	return ret;
}


std::vector<std::pair<int,  int> > fixQuboVariablesMap(std::map<std::pair<int, int>, double> QMap, int QSize, int mtd)
{
    compressed_matrix::CompressedMatrix<double> QInput(QSize, QSize, QMap);

    FixVariablesResult ret = fixQuboVariables(QInput, mtd);

    return ret.fixedVars;
}

template<class V, class B>
std::vector<std::pair<int,  int>> fixQuboVariables(dimod::AdjVectorBQM<V,B>& bqm, int method)
{
      int numVars = bqm.num_variables();
      std::vector<bool> isUsed(numVars, false);

      // Temporary code to maintain compatibility with legacy code
      std::map<std::pair<int, int>, double> QMap;
      for(int i = 0; i < numVars; i++){
	auto span = bqm.neighborhood(i);
	auto linear = bqm.linear(i);
        QMap[{i,i}] = linear;
        auto it = std::lower_bound(span.first, span.second, i+1, dimod::utils::comp_v<V,B>);
        if(it != span.second) {
          for(auto itEnd = span.second; it != itEnd; it++) {
             QMap[{i,it->first}] = it->second;
	  }
	}
      }

     PosiformInfo<dimod::AdjVectorBQM<V,B>> pi(bqm);
     pi.print();

     std::cout << "Size of Implication Node without type : " << sizeof(ImplicationEdge<long long int>) << std::endl;

     ImplicationNetwork<long long int> implicationNet(pi);
     implicationNet.print();
    
     push_relabel<ImplicationEdge<long long int>> pushRelab(implicationNet.adjList, implicationNet.source, implicationNet.sink);
     implicationNet.print();

     pushRelab.global_relabel();
     pushRelab.printLevels();
    
     long long int preflow = pushRelab.maximum_preflow();
     std::cout <<"Preflow from written code : " << preflow <<std::endl;

     printf(" Calling processed map based function \n"); 
     return fixQuboVariablesMap(QMap, numVars, method); 
}

} // namespace fix_variables_

#endif // FIX_VARIABLES_HPP_INCLUDED numVars
