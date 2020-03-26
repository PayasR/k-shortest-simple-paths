/*
 * Yen's k shortest simple paths algorithm
 *
 * Implements the Yen's k shortest simple paths algorithm as described in [1,2].
 * Our implementation does not assume that k is known in advance. Instead,
 * method next_path() acts as an iterator and returns the next shortest path.
 * When no more path can be found, the method returns an empty path with null
 * weight.
 *
 * [1] Jin Y. Yen, "Finding the k Shortest Loopless Paths in a Network".
 *     Management Science. 17 (11): 712-716, 1971. doi:10.1287/mnsc.17.11.712
 * [2] https://en.wikipedia.org/wiki/Yen%27s_algorithm
 */


#ifndef YEN_H
#define YEN_H

#include <cstdint>
#include <limits>       /* used to get maximum possible value of type TV */
#include "digraph.h"
#include "tools.h"
#include "pairing_heap.h"
#include "dijkstra.h"

namespace kssp {

  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class Yen
  {
  public:
    Yen() = default;

    // Constructor
    explicit Yen(directed_graph::DirectedGraph<TI,TV> *g,
		 const TI ssource,
		 const TI ttarget);

    // Destructor
    virtual ~Yen();

    // Check whether there is a candidate path
    bool empty();

    // Get path from/to u
    std::pair<std::vector<TI>,TV> next_path();

    size_t cpt_yielded_paths;
    size_t cpt_used_trees;

  private:
    directed_graph::DirectedGraph<TI,TV> *graph;
    TI source;
    TI target;
    TV MAX_WEIGHT;

    heap::PairingHeap<tools_for_kssp::CandidatePath<TI,TV> *,TV> *heap_sorted_paths;
    std::vector<tools_for_kssp::CandidatePath<TI,TV> *> yielded_paths;
    
    std::set<TI> forbidden_vertices;
    std::vector<std::set<std::pair<TI,TI> > > forbidden_edges;

    // Compute and store deviations of this path
    void compute_deviations(tools_for_kssp::CandidatePath<TI,TV> *candidate_path);
  };


  // Constructor
  template<typename TI, typename TV>
  Yen<TI, TV>::Yen(directed_graph::DirectedGraph<TI,TV> *g,
		   const TI ssource, const TI ttarget):
    graph(g), source(ssource), target(ttarget)
  {
    heap_sorted_paths = new heap::PairingHeap<tools_for_kssp::CandidatePath<TI,TV> *, TV>();
    yielded_paths.clear();
    forbidden_vertices.clear();
    forbidden_edges.clear();
    MAX_WEIGHT = std::numeric_limits<TV>::max();
    cpt_yielded_paths = 0;
    cpt_used_trees = 0;
    
    // Compute shortest path from source to target and store it
    if (source == target)
      {
	// trivial case
	std::vector<TI> path;
	path.clear();
	path.push_back(source);
	tools_for_kssp::CandidatePath<TI,TV> *C = new tools_for_kssp::CandidatePath<TI,TV>(path, 0, 0);
	heap_sorted_paths->push(C, 0);
      }
    else
      {
	// We call Dijkstra from source to target
	dijkstra::Dijkstra<TI,TV> *D = new dijkstra::Dijkstra<TI,TV>(graph, source, false);
	D->run(target);
	cpt_used_trees++;

	// and store computed path, if any
	if (D->successor(target) != target)
	  {
	    std::vector<TI> path = D->get_path(target);
	    tools_for_kssp::CandidatePath<TI,TV> *C = new tools_for_kssp::CandidatePath<TI,TV>(path, D->weight(target), 0);
	    heap_sorted_paths->push(C, 0);
	  }
	delete D;
      }
  }


  // Destructor
  template<typename TI, typename TV>
  Yen<TI, TV>::~Yen()
  {
    while (not heap_sorted_paths->empty())
      {
	tools_for_kssp::CandidatePath<TI,TV> *C = heap_sorted_paths->top_item();
	delete C;
	heap_sorted_paths->pop();
      }
    delete heap_sorted_paths;
    while (not yielded_paths.empty())
      {
	tools_for_kssp::CandidatePath<TI,TV> *C = yielded_paths.back();
	delete C;
	yielded_paths.pop_back();
      }
    forbidden_vertices.clear();
    forbidden_edges.clear();
  }



  // Check whether there is a candidate path
  template<typename TI, typename TV>
  inline bool Yen<TI, TV>::empty()
  {
    return heap_sorted_paths->empty();
  }


  // Get next path
  template<typename TI, typename TV>
  std::pair<std::vector<TI>, TV> Yen<TI, TV>::next_path()
  {

    if (heap_sorted_paths->empty())
      {
	// no more paths
	std::vector<TI> path;
	path.clear();
	return std::make_pair(path, 0);
      }

    // Extract the next best path from the heap
    tools_for_kssp::CandidatePath<TI,TV> *prev_path = heap_sorted_paths->top_item();
    heap_sorted_paths->pop();
    
    if ((prev_path->weight == 0) || (prev_path->weight == MAX_WEIGHT))
      {
	return std::make_pair(prev_path->path, 0);
      }

    // Compute deviations
    compute_deviations(prev_path);
    
    // Store and return the extracted path
    yielded_paths.push_back(prev_path);
    cpt_yielded_paths++;
    return std::make_pair(prev_path->path, prev_path->weight);
  }    

  
  // Compute and store deviations of this path
  template<typename TI, typename TV>
  void Yen<TI, TV>::compute_deviations(tools_for_kssp::CandidatePath<TI,TV> *candidate_path)
  {

    /* Compute the sets of forbidden edges for each prefix.
     * We care only of edges incident to prev_path[i], as the edges incident
     * to other vertices in the prefix are forbidden by forbidden vertices
     */
    std::vector<TI> prev_path = candidate_path->path;
    size_t prev_path_size = prev_path.size();
    size_t dev_idx = candidate_path->dev_idx;

    forbidden_edges.clear();
    forbidden_edges.resize(prev_path_size);

    for (size_t i = dev_idx; i < prev_path_size - 1; i++)
	forbidden_edges[i].insert(std::make_pair(prev_path[i], prev_path[i + 1]));
    for (auto const& c_path: yielded_paths)
      {
	size_t jmax = (prev_path_size < c_path->path.size() ? prev_path_size : c_path->path.size()) - 1;
	size_t j = 1;
	while ((j < jmax) && (prev_path[j] == c_path->path[j]))
	  j++;
	if ((dev_idx < j) && (j < c_path->path.size()))
	  {
	    // j is the first index at which the paths differ
	    forbidden_edges[j - 1].insert(std::make_pair(c_path->path[j - 1], c_path->path[j]));
	  }
      }

    // Initialize the set of forbidden vertices
    TV root_weight = 0;
    forbidden_vertices.clear();
    for (size_t i = 1; i < dev_idx; i++)
      {
	root_weight += graph->edges[std::make_pair(prev_path[i - 1], prev_path[i])];
	forbidden_vertices.insert(prev_path[i - 1]);
      }

    /*
     * Deviate from the previous path to find the candidate paths
     */
    for (size_t i = dev_idx; i < prev_path_size - 1; i++)
      {
	TI spur_node = prev_path[i];

	// update root part of the previous path
	if (i > 0)
	  {
	    root_weight += graph->edges[std::make_pair(prev_path[i - 1], prev_path[i])];
	    forbidden_vertices.insert(prev_path[i - 1]);
	  }

	// Compute shortest path from spur node to target
	dijkstra::Dijkstra<TI,TV> *D = new dijkstra::Dijkstra<TI,TV>(graph, forbidden_vertices, forbidden_edges[i],
								     spur_node, false);
	D->run(target);
	cpt_used_trees++;

	if (D->successor(target) == target)
	  { // no path found
	    continue;
	  }

	// Extract the new candidate path
	std::vector<TI> new_path;
	new_path.clear();
	for (size_t j = 0; j < i; j++)
	  new_path.push_back(prev_path[j]);
	for (TI u: D->get_path(target))
	  new_path.push_back(u);

	tools_for_kssp::CandidatePath<TI,TV> *C = new tools_for_kssp::CandidatePath<TI,TV>(new_path, root_weight + D->weight(target), i);
	heap_sorted_paths->push(C, C->weight);
	delete D;
      }
  }

}
#endif
