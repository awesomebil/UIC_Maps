// application.cpp <Starter Code>
// Bilal Suleman
//
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <queue>
#include <limits>
#include <stack>

#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
#include "graph.h"

using namespace std;
using namespace tinyxml2;

// Class defining how th priority queue will compare items.
class prioritize {
	public:
	    bool operator()(const pair<long long, double>& p1,
	    	const pair<long long, double>& p2) const {
	      return p1.second > p2.second;
	    }
};

double INF = numeric_limits<double>::max();

// Helper function for main Dijkstra algorithm function. Takes care of loop.
void DijkstraHelper(graph<long long, double>& G, long long startV,
	map<long long, double>& distances, map<long long, long long>& pred,
	vector<long long>& visited, set<long long>& visitedSet,
	priority_queue<pair<long long, double>, vector<pair<long long, double>>,
	prioritize>& unvisitedQueue) {
		while(!unvisitedQueue.empty()) {
			long long current = unvisitedQueue.top().first;
			unvisitedQueue.pop();
			if (distances[current] == INF) {
			  break;
			} else if (visitedSet.count(current) > 0) {
			  continue;
			} else {
			  visitedSet.insert(current);
			  visited.push_back(current);
			}
			for (long long adjacent : G.neighbors(current)) {
			  double edgeW = 0;
			  G.getWeight(current, adjacent, edgeW);
			  double pathDistance = distances[current] + edgeW;
			  if (pathDistance < distances[adjacent]) {
			    distances[adjacent] = pathDistance;
			    pred[adjacent] = current;
			  }
			  unvisitedQueue.push(make_pair(adjacent, pathDistance));
			}
		}		
}

// Main Dijkstra function, takes care of setting up algoritm and calling helper.
vector<long long> Dijkstra(graph<long long, double>& G, long long startV,
	  map<long long, double>& distances, map<long long, long long>& pred) {
	vector<long long>  visited;
	set<long long> visitedSet;
	priority_queue<pair<long long, double>, vector<pair<long long, double>>,
		prioritize> unvisitedQueue;
	
	for(long long vertex : G.getVertices()) {
		distances[vertex] = INF;
		unvisitedQueue.push(make_pair(vertex, INF));
		pred[vertex] = 0;
	}
	
	distances[startV] = 0;
	unvisitedQueue.push(make_pair(startV, 0));
	DijkstraHelper(G, startV, distances, pred, visited, visitedSet, 
		unvisitedQueue);
	return visited;
}

// Print map stats to screen, includes number of vertices, edges, nodes etc.
void printStats(map<long long, Coordinates>& Nodes, 
	 vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
	 graph<long long, double>& Map) {
	cout << endl;
	cout << "# of nodes: " << Nodes.size() << endl;
	cout << "# of footways: " << Footways.size() << endl;
	cout << "# of buildings: " << Buildings.size() << endl;
	cout << "# of vertices: " << Map.NumVertices() << endl;
	cout << "# of edges: " << Map.NumEdges() << endl;
	cout << endl;
}

// Constructs the graph from the given map information.
void buildGraph(map<long long, Coordinates>& Nodes,
	vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
	graph<long long, double>& Map) {
  for (auto i : Nodes) {
	Map.addVertex(i.first);
  }
  for (FootwayInfo f : Footways) {
  	for (long unsigned int i = 0; i < f.Nodes.size() - 1; ++i) {
  		long long id = f.Nodes[i];
  		long long next = f.Nodes[i + 1];
  		double lat1 = Nodes[id].Lat, lon1 = Nodes[id].Lon;
  		double lat2 = Nodes[next].Lat, lon2 = Nodes[next].Lon;
  		double distance = distBetween2Points(lat1, lon1, lat2, lon2);
  		Map.addEdge(id, next, distance);
  		Map.addEdge(next, id, distance);
  	}
  }
}

// Verifies the start and end points to be valid.
void verifyInput(vector<BuildingInfo>& Buildings, BuildingInfo& start,
	BuildingInfo& dest, bool& startFound, bool& destFound,
	string& startBuilding, string& destBuilding) {
	for (BuildingInfo b : Buildings) {
    	if (b.Abbrev == startBuilding) {
    		start = b;
    		startFound = true;
    	} 
    	if (b.Abbrev == destBuilding) {
    		dest = b;
    		destFound = true;
    	}
    	if (startFound && destFound) {
    		break;
    	}
    }
    if (!startFound || !destFound) {
    	for (BuildingInfo b : Buildings) {
    		if (b.Fullname.find(startBuilding) != std::string::npos) {
	    		start = b;
	    		startFound = true;
	    	}
	    	if (b.Fullname.find(destBuilding) != std::string::npos) {
	    		dest = b;
	    		destFound = true;
	    	}
	    	if (startFound && destFound) {
	    		break;
	    	}
    	}
    }
}

// Starting prompt for loop asking user to input start and destination.
bool prompt(string& startBuilding, string& destBuilding) {
	cout << endl;
    cout << "Enter start (partial name or abbreviation), or #> ";
    getline(cin, startBuilding);
    if (startBuilding == "#") {
    	return false;
    }
    cout << "Enter destination (partial name or abbreviation)> ";
    getline(cin, destBuilding);
    return true;
}

// Loads map info for program.
bool entry(XMLDocument& xmldoc) {
  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;

  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }
  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return false;
  }
  return true;
}

// Helper function to print out the starting and destinaton points as well
// as their coordinates.
void printPoints(BuildingInfo& start, BuildingInfo& dest, double startLat,
	 double startLon, double destLat, double destLon) {
    cout << "Starting point:" << endl;
    cout << " " << start.Fullname << endl;
    cout << " " << "(" << startLat << ", " << startLon << ")" << endl;
    cout << "Destination point:" << endl;
    cout << " " << dest.Fullname << endl;
    cout << " " << "(" << destLat << ", " << destLon << ")" << endl;
    cout << endl;
}

// Function that finds the minimum distance for start and destinaton nodes.
void findMin(long long& startNode, long long& destNode,
	 map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
	 double& startLat, double& startLon, double& destLat, double& destLon) {
	double minStart = 100000, minDest = 100000;
	for (FootwayInfo f : Footways) {
    	for (long long node : f.Nodes) {
    		Coordinates near = Nodes[node];
    		double distStart, distDest;
    		distStart = distBetween2Points(startLat, startLon, near.Lat, 
    					near.Lon);
    		distDest = distBetween2Points(destLat, destLon, near.Lat, near.Lon);
    		if (distStart < minStart) {
    			minStart = distStart;
    			startNode = node;
    		}
    		if (distDest < minDest) {
    			minDest = distDest;
    			destNode = node;
    		}
    	}
    }
}

// Helper function to print out the path.
void printPath(long long& destNode, map<long long, long long>& pred) {
	stack<long long> path;
	long long nodeId = destNode;
	while (nodeId != 0) {
		path.push(nodeId);
		nodeId = pred[nodeId];
	}
    while (!path.empty()) {
    	if (path.size() > 1) {
    		cout << path.top() << "->";
    	} else {
    		cout << path.top();
    	}
    	path.pop();
    }
    cout << endl;
}

void navigationUI(vector<FootwayInfo>& Footways, 
	vector<BuildingInfo>& Buildings, map<long long, Coordinates>& Nodes,
	graph<long long, double>& Map) {
  //
  // Navigation from building to building
  //
  string startBuilding = "not hash", destBuilding;

  while (startBuilding != "#") {
	if (!prompt(startBuilding, destBuilding)) {
		break;
	}
    bool startFound = false, destFound = false;
    BuildingInfo start, dest;
    verifyInput(Buildings, start, dest, startFound, destFound, startBuilding,
    	destBuilding);
    if (!startFound) {
    	cout << "Start building not found" << endl;
    	continue;
    }
    if (!destFound) {
    	cout << "Destination building not found" << endl;
    	continue;
    }
    double startLat = start.Coords.Lat, startLon = start.Coords.Lon;
    double destLat = dest.Coords.Lat, destLon = dest.Coords.Lon;
    printPoints(start, dest, startLat, startLon, destLat, destLon);
    long long startNode, destNode;
    findMin(startNode, destNode, Nodes, Footways, startLat, startLon,
    	destLat, destLon);
    double sNearLat = Nodes[startNode].Lat, sNearLon = Nodes[startNode].Lon;
    double dNearLat = Nodes[destNode].Lat, dNearLon = Nodes[destNode].Lon;
    
    cout << "Nearest start node: " << endl;
    cout << " " << startNode << endl;
    cout << " " << "(" << sNearLat << ", " << sNearLon << ")" << endl; 
    cout << "Nearest destination node: " << endl;
    cout << " " << destNode << endl;
    cout << " " << "(" << dNearLat << ", " << dNearLon << ")" << endl;
    cout << endl;
	map<long long, double> distances;
	map<long long, long long> pred;
	vector<long long> nodes;
	cout << "Navigating with Dijkstra..." << endl;
	nodes = Dijkstra(Map, startNode, distances, pred);
	if (distances[destNode] == INF) {
		cout << "Sorry, destination unreachable" << endl;
		continue;
	}
	cout << "Distance to dest: " << distances[destNode] << " miles" << endl;
	cout << "Path: ";
	printPath(destNode, pred);
  }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;
  if (!entry(xmldoc)) {
  	return 0;
  }
  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);
  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);
  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);
  
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());
  graph<long long, double> Map;
  
  buildGraph(Nodes, Footways, Buildings, Map);
  printStats(Nodes, Footways, Buildings, Map);
  navigationUI(Footways, Buildings, Nodes, Map);
  //
  // done:
  //
  cout << "** Done **" << endl;
  return 0;
}
