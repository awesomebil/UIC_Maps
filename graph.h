// graph.h <Starter Code>
// Bilal Suleman
//
// Basic graph class using adjacency matrix representation.  Currently
// limited to a graph with at most 100 vertices.
//
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
  unordered_map<VertexT, unordered_map<VertexT, WeightT> > adjList;
  int numEdges;
  vector<VertexT> vertices;
  unordered_map<VertexT, set<VertexT>> neighborSet;

  //
  // _LookupVertex
  //
  // Finds the vertex in the Vertices vector and returns it's
  // index position if found, otherwise returns -1.
  //
  int _LookupVertex(VertexT v) const {
    if (adjList.count(v) > 0) {
    	return adjList.count(v);
    }

    // if get here, not found:
    return -1;
  }

 public:
  //
  // constructor:
  //
  // Constructs an empty graph where n is the max # of vertices
  // you expect the graph to contain.
  //
  // NOTE: the graph is implemented using an adjacency matrix.
  // If n exceeds the dimensions of this matrix, an exception
  // will be thrown to let you know that this implementation
  // will not suffice.
  //
  graph() {
  	adjList.clear();
  	numEdges = 0;
  }

  //
  // NumVertices
  //
  // Returns the # of vertices currently in the graph.
  //
  int NumVertices() const {
    return static_cast<int>(this->adjList.size());
  }

  //
  // NumEdges
  //
  // Returns the # of edges currently in the graph.
  //
  int NumEdges() const {
    return numEdges;
  }

  //
  // addVertex
  //
  // Adds the vertex v to the graph if there's room, and if so
  // returns true.  If the graph is full, or the vertex already
  // exists in the graph, then false is returned.
  //
  bool addVertex(VertexT v) {
    //
    // is the vertex already in the graph?  If so, we do not
    // insert again otherwise Vertices may fill with duplicates:
    //
    if (_LookupVertex(v) >= 0) {
      return false;
    }

    //
    // if we get here, vertex does not exist so insert.  Where
    // we insert becomes the rows and col position for this
    // vertex in the adjacency matrix.
    //
    this->adjList[v];
    vertices.push_back(v);

    return true;
  }

  //
  // addEdge
  //
  // Adds the edge (from, to, weight) to the graph, and returns
  // true.  If the vertices do not exist or for some reason the
  // graph is full, false is returned.
  //
  // NOTE: if the edge already exists, the existing edge weight
  // is overwritten with the new edge weight.
  //
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    if (_LookupVertex(from) == -1) {
    	return false;
    } else if (_LookupVertex(to) == -1) {
    	return false;
    }
    if (adjList[from].count(to) == 1) {
    	numEdges--;
    }
    adjList[from][to] = weight;
    neighborSet[from].insert(to);
    numEdges++;
    return true;
  }

  //
  // getWeight
  //
  // Returns the weight associated with a given edge.  If
  // the edge exists, the weight is returned via the reference
  // parameter and true is returned.  If the edge does not
  // exist, the weight parameter is unchanged and false is
  // returned.
  //
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
  	if (_LookupVertex(from) == -1) {
    	return false;
    } else if (_LookupVertex(to) == -1) {
    	return false;
    } else if (adjList.at(from).count(to) == 0) {
    	return false;
    }

    //
    // Okay, the edge exists, return the weight via the
    // reference parameter:
    //
    weight = adjList.at(from).at(to);

    return true;
  }

  //
  // neighbors
  //
  // Returns a set containing the neighbors of v, i.e. all
  // vertices that can be reached from v along one edge.
  // Since a set is returned, the neighbors are returned in
  // sorted order; use foreach to iterate through the set.
  //
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> S;
    if (neighborSet.count(v) == 0) {
    	return S;
    }
    return neighborSet.at(v);
  }

  //
  // getVertices
  //
  // Returns a vector containing all the vertices currently in
  // the graph.
  //
  vector<VertexT> getVertices() const {
    return vertices;
  }

  //
  // dump
  //
  // Dumps the internal state of the graph for debugging purposes.
  //
  // Example:
  //    graph<string,int>  G(26);
  //    ...
  //    G.dump(cout);  // dump to console
  //
  void dump(ostream& output) const {
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    int i = 0;
    for (auto v : adjList) {
    	++i;
    	output << " " << i << ". " << v.first << endl; 
    }
    i = 0;

    output << endl;
    output << "**Edges:" << endl;
    for (auto v : adjList) {
      ++i;
      output << " Vertex: " << v.first << ": ";

      for (auto e : v.second) {
        output << "(" << e.first << ", " << e.second << ") " << "->"; 
      }
      output << endl;
    }
    output << "**************************************************" << endl;
  }
};
