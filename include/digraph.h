/*
 * Directed Graph
 *
 * Defines a directed graph data structure with vertices indexed in range 0..n-1
 */


#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

#include <cstdint>
#include <vector>
//#include <iostream>
//#include <iterator>
#include <map>

namespace directed_graph {

  template<
    typename TI,  // type of node indices (int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class DirectedGraph
  {
  public:
    DirectedGraph() = default;

    // Number of nodes of the directed graph
    TI n;
    // Number of arcs
    TI m;
    // Mapping from edges to weight
    std::map<std::pair<TI,TI>,TV> edges;
    // Lists of neighbors with edge weights
    std::vector<std::vector<std::pair<TI,TV> > > out_neighbors;
    std::vector<std::vector<std::pair<TI,TV> > > in_neighbors;

    // Constructor
    explicit DirectedGraph(const TI number_of_nodes);

    // Destructor
    virtual ~DirectedGraph();

    // Empty and reset the directed graph
    void reset(const TI number_of_nodes);

    //void read_from_file(std::string const& filename);
    
    // Return the number of nodes
    TI order() const { return n; }

    // Return the number of edges
    TI size() const {return m; }

    // Add an edge from u to v with weight w
    void add_edge(const TI u, const TI v, const TV w);

    // Check whether the directed graph contains an edge from u to v
    bool has_edge(const TI u, const TI v) const { return (edges.count(std::make_pair(u, v))>0?true:false); }

    // Return the out degree of node u
    TI out_degree(const TI u) const { return out_neighbors[u].size(); }
  
    // Return the in degree of node u
    TI in_degree(const TI u) const { return in_neighbors[u].size(); }

    // Return the weight of edge from u to v
    TV weight(const TI u, const TI v);
  };

  // Constructor
  template<typename TI, typename TV>
  DirectedGraph<TI, TV>::DirectedGraph(const TI number_of_nodes):
    n(number_of_nodes), m(0)
  {
    reset(number_of_nodes);
  }

  // Destructor
  template<typename TI, typename TV>
  DirectedGraph<TI, TV>::~DirectedGraph()
  {
    reset(0);
  }

  // Empty the directed graph
  template<typename TI, typename TV>
  void DirectedGraph<TI, TV>::reset(const TI number_of_nodes)
  {
    n = number_of_nodes;
    m = 0;
    edges.clear();
    out_neighbors.clear();
    out_neighbors.resize(n);
    in_neighbors.clear();
    in_neighbors.resize(n);
  }

  // Add an edge from u to v with weight w
  template<typename TI, typename TV>
  void DirectedGraph<TI, TV>::add_edge(const TI u, const TI v, const TV w)
  {
    edges[std::make_pair(u, v)] = w;
    out_neighbors[u].push_back(std::make_pair(v, w));
    in_neighbors[v].push_back(std::make_pair(u, w));
    m++;
  }

  // Return the weight of edge from u to v
  template<typename TI, typename TV>
  TV DirectedGraph<TI, TV>::weight(const TI u, const TI v)
  {
    return edges[std::make_pair(u, v)];
  }


  
}

#endif
