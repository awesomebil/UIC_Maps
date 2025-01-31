// testing.cpp <Starter Code>
// <Your name>
//
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
// This file is used for testing graph.h.  We encourage you to use Google'
// test framework for this project, but it is not required (because we will
// not be grading it).  
//

#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <string>
#include <fstream>

#include "graph.h"

using namespace std;


//
// buildGraph:
//
// Inputs the graph vertices and edges from the given file, building
// the graph g.  File format:
//   vertex 
//   vertex
//   ...
//   #
//   src dest weight
//   src dest weight
//   ... 
//   #
//
void buildGraph(string filename, graph<string,int>& G)
{
  ifstream file(filename);
  string   v;

  if (!file.good())
  {
    cout << endl;
    cout << "**Error: unable to open input file '" << filename << "'." << endl;
    cout << endl;
    return;
  }

  //
  // Input vertices as single uppercase letters:  A B C ... #
  //
  file >> v;

  while (v != "#")
  {
    if (!G.addVertex(v))
      cout << "**Error: unable to add vertex '" << v << "', why not?" << endl;

    file >> v;
  }

  //
  // Now input edges:  Src Dest Weight ... #
  //
  string src, dest;
  int  weight;

  file >> src;

  while (src != "#")
  {
    file >> dest;
    file >> weight;

    if (!G.addEdge(src, dest, weight))
      cout << "**Error: unable to add edge (" << src << "," << dest << "," << weight << "), why not?" << endl;

    file >> src;
  }
}

//
// outputGraph:
//
// Outputs graph g to the console.
//
void outputGraph(graph<string,int>& G)
{
  vector<string> vertices = G.getVertices();

  cout << "**Vertices: ";

  for (string v : vertices)
  {
    cout << v << " ";
  }

  cout << endl;

  cout << "**Edges: ";

  for (string v : vertices)
  {
    set<string> neighbors = G.neighbors(v);

    for (string n : neighbors)
    {
      int weight;
      
      if (G.getWeight(v, n, weight))
      {
        cout << "(" << v << "," << n << "," << weight << ") ";
      }
      else
      {
        cout << "(" << v << "," << n << "," << "???" << ") ";
      }
    }
  }

  cout << endl;
}


int main()
{
  graph<string,int> G;
  string filename;
  string startV;

  cout << "Enter filename containing graph data> ";
  cin >> filename;
  cout << endl;

  //
  // Let's input the graph, and then output to see what we have:
  //
  buildGraph(filename, G);
  
  outputGraph(G);
  
  G.dump(cout);
  int weight;
  G.getWeight("G", "A", weight);
  cout << weight << endl;
  if (G.getWeight("F", "G", weight)) {
	cout << weight << endl;	
  } else {
  	cout << "Edge does not exist." << endl;
  }
  
  if (G.getWeight("Random", "Bilal", weight)) {
	cout << weight << endl;	
  } else {
  	cout << "Edge does not exist." << endl;
  }
  G.addVertex("Bilal");
  G.addVertex("Suleman");
  G.addEdge("BIlal", "Suleman", 100);
  if (G.getWeight("Bilal", "Suleman", weight)) {
	cout << weight << endl;	
  } else {
  	cout << "Edge does not exist." << endl;
  }
  G.addEdge("G", "A", 50);
  G.dump(cout);

  //
  // done:
  //
  return 0;
}
