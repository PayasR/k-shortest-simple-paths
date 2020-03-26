/*
 * Easy Directed Graph
 *
 * Defines a directed graph data structure with templated vertex and edge
 * weights types. Essentially, this data structure is an interface for
 * the DirectedGraph class that considers vertices indexed in 0..n-1.
 *
 */


#ifndef EASY_DIRECTED_GRAPH_H
#define EASY_DIRECTED_GRAPH_H

#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
//#include <iterator>
#include <map>
#include <cstring>
#include <algorithm>

#include "digraph.h"
#include "yen.h"
#include "node_classification.h"
#include "sidetrack_based.h"
#include "parsimonious_sidetrack_based.h"


#define ERROR(x) do { std::cerr << "Error: " << x << std::endl;\
                      std::exit(EXIT_FAILURE); } while (0)


namespace directed_graph {

  template<
    typename TI,    // type of node indices (int, uint32_t, etc.)
    typename TV,    // type of edge weights (int, unsigned int, float, etc.)
    typename DGTI   // type of node indices (int, uint32_t, etc.)
    >
  class EasyDirectedGraph
  {
  public:
    TI n;  // Number of nodes of the directed graph
    TI m;  // Number of arcs
    std::vector<TI> int_to_vertex;     // mapping from 0..n-1 to vertex indexes
    std::map<TI,size_t> vertex_to_int; // mapping from vertex indexes to integers in 0..n-1

    EasyDirectedGraph() = default;

    // Constructor - Read graph from file
    explicit EasyDirectedGraph(std::string const& filename);

    // Destructor
    virtual ~EasyDirectedGraph();

    // Empty and reset the directed graph
    void reset(const TI number_of_nodes);

    //void read_from_file(std::string const& filename);
    
    // Return the number of nodes
    TI order() const { return n; }

    // Return the number of edges
    TI size() const {return m; }

    // Check whether the directed graph contains an edge from u to v
    bool has_edge(const TI u, const TI v) const { return g->has_edge(vertex_to_int[u], vertex_to_int[v]); }

    // Return the out degree of node u
    TI out_degree(const TI u) const { return g->out_degree(vertex_to_int[u]); }
  
    // Return the in degree of node u
    TI in_degree(const TI u) const { return g->in_degree(vertex_to_int[u]); }

    // Return the weight of edge from u to v
    TV weight(const TI u, const TI v) const { return g->weight(vertex_to_int[u], vertex_to_int[v]); }


    // reset all algorithms
    void reset_kssp();

    // initialize algorithm
    void init_kssp(std::string const& algorithm, TI source, TI target,
		   size_t version, bool verbose);

    // Return next shortest simple path from source to target using specified algorithm
    std::pair<std::vector<TI>, TV> next_path();

    size_t used_trees();

    void report();

  private:

    DirectedGraph<DGTI,TV> *g;

    kssp::Yen<DGTI, TV> *Y;
    kssp::NodeClassification<DGTI, TV> *NC;
    kssp::SidetrackBased<DGTI, TV> *SB;
    kssp::ParsimoniousSidetrackBased<DGTI, TV> *PSB;

    void read_from_file(std::string const& filename);

    std::vector<TI> convert_path(std::vector<DGTI> path);

  };

  // Constructor
  template<typename TI, typename TV, typename DGTI>
  EasyDirectedGraph<TI, TV, DGTI>::EasyDirectedGraph(std::string const& filename)
  {
    read_from_file(filename);
    Y = nullptr;
    NC = nullptr;
    SB = nullptr;
    PSB = nullptr;
  }

  // Destructor
  template<typename TI, typename TV, typename DGTI>
  EasyDirectedGraph<TI, TV, DGTI>::~EasyDirectedGraph()
  {
    reset_kssp();
    delete g;
  }

  // Read graph from file
  template<typename TI, typename TV, typename DGTI>
  void EasyDirectedGraph<TI, TV, DGTI>::read_from_file(std::string const& filename)
  {
    std::ifstream file(filename);
    if (!file.is_open())
      ERROR("Could not open file " << filename);

    std::streamsize ignore_count = std::numeric_limits<std::streamsize>::max();
    std::map<std::pair<TI, TI>, TV> tmp_edges;
    tmp_edges.clear();
    std::string a, source_string, target_string, weight_string;
    int_to_vertex.clear();

    while (file.peek() != std::ifstream::traits_type::eof())
      {
	if ( ((file >> std::ws).peek() == 'c') || ((file >> std::ws).peek() == 'p') )
	  {
	    file.ignore(ignore_count, '\n');
	    continue;
	  }

	file >> a >> source_string >> target_string >> weight_string;
	auto source_id = std::stoi(source_string);
	auto target_id = std::stoi(target_string);
	auto weight = std::stoi(weight_string);
	assert((source_id >= 0) && (target_id >= 0));

	tmp_edges[std::make_pair(source_id, target_id)] = weight;
	int_to_vertex.push_back(source_id);
	int_to_vertex.push_back(target_id);
      }

    // remove duplicate vertex ids
    std::sort(int_to_vertex.begin(), int_to_vertex.end());
    auto new_end = std::unique(int_to_vertex.begin(), int_to_vertex.end());
    int_to_vertex.erase(new_end, int_to_vertex.end());
    n = int_to_vertex.size();

    // compute mapping
    vertex_to_int.clear();
    size_t i = 0;
    for (TI u: int_to_vertex)
      vertex_to_int[u] = i++;

    // build directed graph
    g = new DirectedGraph<DGTI,TV>(int_to_vertex.size());
    for (auto const& edge: tmp_edges)
      {
	DGTI u = vertex_to_int[edge.first.first];
	DGTI v = vertex_to_int[edge.first.second];
	g->add_edge(u, v, edge.second);
      }
    m = g->m;

    tmp_edges.clear();
  }

  // reset all algorithms
  template<typename TI, typename TV, typename DGTI>
  void EasyDirectedGraph<TI, TV, DGTI>::report()
  {
    //if (Y != nullptr) delete Y;
    //if (NC != nullptr) delete F;
    if (SB != nullptr) {SB->verbose = true; SB->report();}
    if (PSB != nullptr) {PSB->verbose = true; PSB->report();}
  }

  
  // reset all algorithms
  template<typename TI, typename TV, typename DGTI>
  void EasyDirectedGraph<TI, TV, DGTI>::reset_kssp()
  {
    if (Y != nullptr) delete Y;
    if (NC != nullptr) delete NC;
    if (SB != nullptr) delete SB;
    if (PSB != nullptr) delete PSB;
    Y = nullptr;
    NC = nullptr;
    SB = nullptr;
    PSB = nullptr;
  }

  // initialize algorithms
  template<typename TI, typename TV, typename DGTI>
  void EasyDirectedGraph<TI, TV, DGTI>::init_kssp(std::string const& algorithm,
						  TI source, TI target,
						  size_t version, bool verbose)
  {
    if ((algorithm.compare("PSB") == 0) || (algorithm.compare("Parsimonious_Sidetrack_Based") == 0))
      PSB = new kssp::ParsimoniousSidetrackBased<DGTI, TV>(g, vertex_to_int[source], vertex_to_int[target],
							   version, verbose);
    else if ((algorithm.compare("SB") == 0) || (algorithm.compare("Sidetrack_Based") == 0))
      SB = new kssp::SidetrackBased<DGTI, TV>(g, vertex_to_int[source], vertex_to_int[target], version, verbose);
    else if ((algorithm.compare("NC") == 0) || (algorithm.compare("Node_Classification") == 0))
      NC = new kssp::NodeClassification<DGTI, TV>(g, vertex_to_int[source], vertex_to_int[target]);
    else if ((algorithm.compare("Y") == 0) || (algorithm.compare("Yen") == 0))
      Y = new kssp::Yen<DGTI, TV>(g, vertex_to_int[source], vertex_to_int[target]);
    else
      ERROR("Unkown algorithm " << algorithm);
  }


  // Return next shortest simple path from source to target using specified algorithm
  template<typename TI, typename TV, typename DGTI>
  std::pair<std::vector<TI>, TV> EasyDirectedGraph<TI, TV, DGTI>::next_path()
  {
    if (PSB != nullptr)
      {
	auto res = PSB->next_path();
	return std::make_pair(convert_path(res.first), res.second);
      }
    else if (SB != nullptr)
      {
	auto res = SB->next_path();
	return std::make_pair(convert_path(res.first), res.second);
      }
    else if (NC != nullptr)
      {
	auto res = NC->next_path();
	return std::make_pair(convert_path(res.first), res.second);
      }
    else if (Y != nullptr)
      {
	auto res = Y->next_path();
	return std::make_pair(convert_path(res.first), res.second);
      }
    else
      ERROR("No initialized algorithm");
  }

  // convert a path from internal indexes (0..n-1) to external ids
  template<typename TI, typename TV, typename DGTI>
  std::vector<TI> EasyDirectedGraph<TI, TV, DGTI>::convert_path(std::vector<DGTI> path)
  {
    std::vector<TI> new_path;
    new_path.clear();
    for (DGTI u: path)
      new_path.push_back(int_to_vertex[u]);
    return new_path;
  }

  template<typename TI, typename TV, typename DGTI>
  size_t EasyDirectedGraph<TI, TV, DGTI>::used_trees()
  {
    if (PSB != nullptr)
      return PSB->cpt_used_trees;
    else if (SB != nullptr)
      return SB->cpt_used_trees;
    else if (NC != nullptr)
      return NC->cpt_used_trees;
    else if (Y != nullptr)
      return Y->cpt_used_trees;
    else
      return 0;
  }
  
}

#endif
