/*
 * Sidetrack Based algotihm for the k shortest simple paths problem.
 *
 * Implements the k shortest simple paths algorithm as described in [1,2].
 * Our implementation does not assume that k is known in advance. Instead,
 * method next_path() acts as an iterator and returns the next shortest path.
 * When no more path can be found, the method returns an empty path with null
 * weight.
 *
 * [1] Denis Kurz, Petra Mutzel. "A Sidetrack-Based Algorithm for Finding the
 *     k Shortest Simple Paths in a Directed Graph". In Proc. of International
 *     Symposium on Algorithms and Computation (ISAAC), LIPIcs 64:1-13, 2016.
 *     doi:10.4230/LIPIcs.ISAAC.2016.49
 *
 * [2] Denis Kurz. "K-Best Enumeration -- Theory and Application". PhD Thesis,
 *     TU Dortmund, 2018. doi:10.17877/DE290R-19814
 */


#ifndef SIDETRACK_BASED_H
#define SIDETRACK_BASED_H


#include <iostream>
#include <limits>       /* used to get maximum possible value of type TV */
#include "digraph.h"
#include "tools.h"
#include "pairing_heap.h"
#include "lazy_dijkstra.h"

//using namespace std;
using namespace tools_for_kssp;

namespace kssp {
  
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class SidetrackBased
  {
  public:
    SidetrackBased() = default;

    // Constructor
    explicit SidetrackBased(directed_graph::DirectedGraph<TI,TV> *g,
			    const TI ssource,
			    const TI ttarget,
			    const size_t version,
			    bool verbose);

    // Destructor
    virtual ~SidetrackBased();

    // Check whether there is a candidate path
    bool empty();

    // Get path from/to u
    std::pair<std::vector<TI>,TV> next_path();

    // report some informations
    void report();

    size_t cpt_yielded_paths;
    size_t cpt_used_trees;
    bool verbose;

  private:
    directed_graph::DirectedGraph<TI,TV> *graph;
    TI source;
    TI target;
    TI n;
    size_t version;

    // Storage of all trees (pointers to tree)
    std::vector<dijkstra::LazyDijkstra<TI,TV> *> idx_to_tree;

    // Priority queue of candidate paths
    heap::PairingHeap<CandidateSimple<TI,TV> *,TV> *heap_candidate_paths;
    
    // Priority queue of hard candidate paths
    heap::PairingHeap<CandidateHard<TI,TV> *,TV> *heap_hard_candidate_paths;

    // to ease some tests
    size_t *seen;
    size_t this_round;


    // data structure for the last extracted path
    std::vector<TI> extracted_path;
    TV extracted_path_weight;

    void extract_path(CandidateSimple<TI,TV> *candidate_path);
    

    // Compute and store deviations of this path
    void compute_deviations(CandidateSimple<TI,TV> *candidate_path);


    // data structures for colors
    size_t *color;

    void compute_colors(CandidateSimple<TI,TV> *candidate_path,
			std::vector<std::pair<size_t,TI> > *sidetrack_edges);

    // Compute deviations for one hard candidate
    void treat_one_hard_candidate();

    //

    size_t MAX_COLOR;
    TV MAX_WEIGHT;
    size_t MAX_TREE_IDX;
  };


  /*
   * Constructor
   */
  template<typename TI, typename TV>
  SidetrackBased<TI, TV>::SidetrackBased(directed_graph::DirectedGraph<TI,TV> *g,
					 const TI ssource, const TI ttarget,
					 const size_t version,
					 bool verbose):
    verbose(verbose), graph(g), source(ssource), target(ttarget), version(version)
  {
    n = graph->n;
    idx_to_tree.clear();
    heap_candidate_paths = new heap::PairingHeap<CandidateSimple<TI,TV> *,TV>();
    heap_hard_candidate_paths = new heap::PairingHeap<CandidateHard<TI,TV> *,TV>();
    
    MAX_COLOR = std::numeric_limits<size_t>::max();
    MAX_WEIGHT = std::numeric_limits<TV>::max();
    MAX_TREE_IDX = std::numeric_limits<size_t>::max();
    seen = (size_t*)calloc(n, sizeof(size_t));
    this_round = 0;
    color = (size_t*)malloc(n * sizeof(size_t));
    cpt_yielded_paths = 0;
    cpt_used_trees = 0;

    // We create a shortest path tree to target and store it
    dijkstra::LazyDijkstra<TI,TV> *T = new dijkstra::LazyDijkstra<TI,TV>(graph, target);
    idx_to_tree.push_back(T);
    cpt_used_trees++;
    std::vector<std::pair<TI,TI> > sidetrack_sequence;
    sidetrack_sequence.clear();
    std::vector<size_t> tree_index_sequence;
    tree_index_sequence.clear();
    tree_index_sequence.push_back(0);
    CandidateSimple<TI,TV> *candidate = new CandidateSimple<TI,TV>(T->weight(source), 0,
								   sidetrack_sequence,
								   tree_index_sequence);
    heap_candidate_paths->push(candidate, T->weight(source));
  }


  /*
   * Destructor
   */
  template<typename TI, typename TV>
  SidetrackBased<TI, TV>::~SidetrackBased()
  {
    report();
    while (not idx_to_tree.empty())
      {
	if (idx_to_tree.back() != nullptr)
	  delete idx_to_tree.back();
	idx_to_tree.pop_back();
      }
    while (not heap_candidate_paths->empty())
      {
	delete heap_candidate_paths->top_item();
	heap_candidate_paths->pop();
      }
    delete heap_candidate_paths;
    while (not heap_hard_candidate_paths->empty())
      {
	delete heap_hard_candidate_paths->top_item();
	heap_hard_candidate_paths->pop();
      }
    delete heap_hard_candidate_paths;
    free(seen);
    free(color);
  }


  /*
   * Report some informations
   */
  template<typename TI, typename TV>
  void SidetrackBased<TI, TV>::report()
  {
    if (verbose == true)
      {
	std::cout << "# good = " << heap_candidate_paths->size();
	std::cout << "\t# hard = " << heap_hard_candidate_paths->size();
	size_t cpt = 0;
	for (size_t i = 0; i < idx_to_tree.size(); i++)
	  if (idx_to_tree[i] != nullptr)
	    cpt++;
	std::cout << "\t# trees = " << cpt << " / " << idx_to_tree.size();
	std::cout << "\t# yielded = " << cpt_yielded_paths;
	if (not heap_candidate_paths->empty())
	  std::cout << "\ttop simple = " << heap_candidate_paths->top_value();
	if (not heap_hard_candidate_paths->empty())
	  std::cout << "\ttop hard = " << heap_hard_candidate_paths->top_value();
	std::cout << std::endl;
	for (size_t i = 0; i < idx_to_tree.size(); i++)
	  if (idx_to_tree[i] != nullptr)
	    idx_to_tree[i]->show();
      }
  }

  

  /*
   * Check whether there is a candidate path left
   */
  template<typename TI, typename TV>
  inline bool SidetrackBased<TI, TV>::empty()
  {
    return heap_candidate_paths->empty() && heap_hard_candidate_paths->empty();
  }


  /*
   * Check whether there is a candidate path and return it
   */
  template<typename TI, typename TV>
  std::pair<std::vector<TI>,TV> SidetrackBased<TI, TV>::next_path()
  {
    report();

    std::vector<TI> path;
    path.clear();
    TV path_w = 0;

    while (not empty())
      {
	if (heap_candidate_paths->empty() ||
	    ((not heap_hard_candidate_paths->empty()) &&
	     (heap_candidate_paths->top_value() > heap_hard_candidate_paths->top_value())))
	  {
	    treat_one_hard_candidate();
	  }
	else
	  {
	    // Extract the next best path from the heap
	    CandidateSimple<TI,TV> *candidate_path = heap_candidate_paths->top_item();
	    heap_candidate_paths->pop();
	    cpt_yielded_paths++;

	    if (source == target)
	      path.push_back(source);
	    else
	      {
		extract_path(candidate_path);
		compute_deviations(candidate_path);

		for (TI u: extracted_path)
		  path.push_back(u);
		path_w = extracted_path_weight;
	      }

	    // this candidate is no longer needed
	    delete candidate_path;
	    break;
	  }
      }
    return std::make_pair(path, path_w);
  }


  /*
   * Extract the path encoded into candidate_path and store it in extracted_path
   * and extracted_path_weight.
   * We know that this path is simple !
   */
  template<typename TI, typename TV>
  void SidetrackBased<TI,TV>::extract_path(CandidateSimple<TI,TV> *candidate_path)
  {
    extracted_path.clear();
    extracted_path.push_back(source);
    extracted_path_weight = 0;
    size_t idx_tree = 0;  // index of tree in tree_index_sequence
    size_t current_tree = candidate_path->tree_index_sequence[0];

    for (auto const& edge: candidate_path->sidetrack_sequence)
      {
	// Follow the path to the tail of edge in current tree
	TI u = extracted_path.back();
	TI tail = edge.first;
	extracted_path_weight += idx_to_tree[current_tree]->weight(u);
	while (u != tail)
	  {
	    u = idx_to_tree[current_tree]->successor(u);
	    extracted_path.push_back(u);
	  }

	// update weight, follow edge and go to next tree
	extracted_path_weight -= idx_to_tree[current_tree]->weight(tail);
	extracted_path_weight += graph->edges[edge];  // weight of edge
	extracted_path.push_back(edge.second);  // head of edge
	idx_tree++;
	current_tree = candidate_path->tree_index_sequence[idx_tree];
      }
    
    // We follow the last part of the path to target
    TI u = extracted_path.back();
    extracted_path_weight += idx_to_tree[current_tree]->weight(u);
    while (u != target)
      {
	u = idx_to_tree[current_tree]->successor(u);
	extracted_path.push_back(u);
      }
  }


  
  /*
   * Compute deviations for this candidate path
   */
  template<typename TI, typename TV>
  void SidetrackBased<TI, TV>::compute_deviations(CandidateSimple<TI,TV> *candidate_path)
  {
    // Get current shortest path tree
    dijkstra::LazyDijkstra<TI,TV> *T = idx_to_tree[candidate_path->tree_index_sequence.back()];

    // will be used to avoid sidetrack edges pointing to a vertex in the prefix
    this_round++;
    for (size_t i = 0; i <= candidate_path->dev_idx; i++)
      seen[extracted_path[i]] = this_round;
	   
    // 1. Compute the list of sidetrack edges
    // We store the pair (index in extracted path, head)
    std::vector<std::pair<size_t,TI> > sidetracks;
    sidetracks.clear();
    for (size_t i = candidate_path->dev_idx; i < extracted_path.size() - 1; i++)
      {
	seen[extracted_path[i + 1]] = this_round;
	for (auto const& it: graph->out_neighbors[extracted_path[i]])
	  {
	    TI head = it.first;
	    if (seen[head] != this_round)
	      sidetracks.push_back(std::make_pair(i, head));
	  }
      }

    // 2. Compute colors
    compute_colors(candidate_path, &sidetracks);

    // 3. Determine simple and non-simple deviations
    size_t prev_i = n; // used to reuse same tree for each sidetracks from a spur node
    size_t hard_tree_index = MAX_TREE_IDX;
    for (auto const& edge: sidetracks)
      {
	TI head = edge.second;
	if (color[head] == MAX_COLOR)
	  continue;
	size_t i = edge.first;
	if (color[head] <= i)
	  {
	    if (T->weight(head) != MAX_WEIGHT)
	      { // We have a non simple path to target
		// We create a new hard candidate
		TI u = extracted_path[i];
		TV lbw = extracted_path_weight - T->weight(u) + graph->edges[std::make_pair(u, head)] + T->weight(head);
		if (i != prev_i)
		  { // We need a new tree. Otherwise, we can use the same tree
		    hard_tree_index = idx_to_tree.size();
		    idx_to_tree.push_back(nullptr); // The construction of the tree is postponed
		    prev_i = i;
		  }
		CandidateHard<TI,TV> *hard = new CandidateHard<TI,TV>(lbw, i, *candidate_path, edge, hard_tree_index);
		heap_hard_candidate_paths->push(hard, lbw);
	      }
	  }
	else
	  { // We have a simple path to target in T, so a new simple path
	    TI u = extracted_path[i];
	    TV weight_path = extracted_path_weight - T->weight(u) + graph->edges[std::make_pair(u, head)] + T->weight(head);
	    CandidateSimple<TI,TV> *new_candidate = new CandidateSimple<TI,TV>(candidate_path);
	    new_candidate->weight = weight_path;
	    new_candidate->dev_idx = i + 1;
	    new_candidate->sidetrack_sequence.push_back(std::make_pair(u, head));
	    new_candidate->tree_index_sequence.push_back(new_candidate->tree_index_sequence.back());
	    heap_candidate_paths->push(new_candidate, weight_path);
	  }
      }
  }

  
  /* Compute colors
   *
   * Given a candidate path and a list of sidetrack edges for this path, this
   * method performs a node classification à la Feng. That is, it assigns colors
   * red (MAX_COLOR), yellow (0) and green (>= 1) to the heads of sidetrack
   * edges so that:
   * - red: the head is a forbidden vertex
   * - yellow: the path from head to target is not simple in current tree
   * - green: we have a shortest path from head to target in current tree
   *
   * In practice, we initialize the colors to 0, then set color i to the ith
   * vertex of the extracted_path, except the source which is red. Then, we
   * follow the path from each head of a sidetrack to a vertex in extracted_path
   * and assign each vertex along that path index i if we reached the ith vertex.
   */
  template<typename TI, typename TV>
  void SidetrackBased<TI, TV>::compute_colors(CandidateSimple<TI,TV> *candidate_path,
					      std::vector<std::pair<size_t,TI> > *sidetrack_edges)
  {
    // Get current shortest path tree
    dijkstra::LazyDijkstra<TI,TV> *T = idx_to_tree[candidate_path->tree_index_sequence.back()];

    // Initialize colors
    std::memset(color, 0, n * sizeof(size_t));
    for (size_t i = 1; i < extracted_path.size(); i++)
      color[extracted_path[i]] = i;
    color[source] = MAX_COLOR;

    // For each head of a sidetrack edge, we follow the path in T until reaching
    // a colored vertex. We then assign that color to each vertex along that path
    std::vector<TI> path;
    for (auto const& edge: *sidetrack_edges)
      {
	TI u = edge.second;
	path.clear();
	path.push_back(u);
	while ((color[u] == 0) && (u != T->successor(u)))
	  {
	    u = T->successor(u);
	    path.push_back(u);
	  }
	size_t col = color[u];
	if ((col != 0) && (col != MAX_COLOR))
	  for (TI v: path)
	    color[v] = col;
      }
  }


  /*
   * Treat one hard candidate
   *
   * This method searches for the shortest simple path from the head of the
   * sidetrack to target in the graph minus the prefix (path from source to the
   * tail of this sidetrack). If such a path is found, we push a new simple
   * candidate path. Otherwise, we discard this sidetrack.
   * Observe that the tree used here is reused for other sidetracks with same
   * tail and same prefix.
   */
  template<typename TI, typename TV>
  void SidetrackBased<TI, TV>::treat_one_hard_candidate()
  {
    if (heap_hard_candidate_paths->empty())
      return;

    CandidateHard<TI,TV> *hard_candidate = heap_hard_candidate_paths->top_item();
    heap_hard_candidate_paths->pop();
    
    extract_path(&(hard_candidate->candidate));

    std::vector<TI> forbidden;
    for (size_t i = 0; i <= hard_candidate->sidetrack.first; i++)
      forbidden.push_back(extracted_path[i]);

    dijkstra::LazyDijkstra<TI,TV> *T_orig = idx_to_tree[hard_candidate->candidate.tree_index_sequence.back()];
    dijkstra::LazyDijkstra<TI,TV> *T_hard = idx_to_tree[hard_candidate->tree_index];
    if (T_hard == nullptr)
      { // We need to build the tree
	if (version == 1)
	  { // Build a fresh new tree and disable some vertices
	    T_hard = new dijkstra::LazyDijkstra<TI,TV>(graph, &forbidden, target);
	  }
	else
	  { // Better way to copy and disable
	    T_hard = new dijkstra::LazyDijkstra<TI,TV>(T_orig, &forbidden);
	  }
	idx_to_tree[hard_candidate->tree_index] = T_hard;
	cpt_used_trees++;
      }

    TI head = hard_candidate->sidetrack.second;
    if (T_hard->successor(head) != head)
      {
	// We have a simple path. We push it to heap_candidate_paths
	CandidateSimple<TI,TV> *new_candidate = new CandidateSimple<TI,TV>(&(hard_candidate->candidate));
	size_t idx = hard_candidate->sidetrack.first;
	// TI u = extracted_path[idx];
	// new_candidate->weight = extracted_path_weight - T_orig->weight(u) + graph->edges[std::make_pair(u, head)] + T_hard->weight(head);
	new_candidate->weight = hard_candidate->weight - T_orig->weight(head) + T_hard->weight(head);
	new_candidate->dev_idx = idx + 1;
	new_candidate->sidetrack_sequence.push_back(std::make_pair(extracted_path[idx], head));
	new_candidate->tree_index_sequence.push_back(hard_candidate->tree_index);
	heap_candidate_paths->push(new_candidate, new_candidate->weight);
      }

    // We are done with this hard candidate
    delete hard_candidate;
  }
  
} // end namespace
#endif

