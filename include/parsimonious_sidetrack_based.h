/*
 * Parsimonious Sidetrack Based algotihm for the k shortest simple paths problem
 *
 * This is a variant proposed by Al Zoobi, Coudert and Nisse of the Sidetrack
 * Based (SB) algorithm by Kurz and Mutzel [1].
 *
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
 */


#ifndef PARSIMONIOUS_SIDETRACK_BASED_H
#define PARSIMONIOUS_SIDETRACK_BASED_H

#include <iostream>
#include <limits>       /* used to get maximum possible value of type TV */
#include "digraph.h"
#include "tools.h"
#include "pairing_heap.h"
#include "lazy_dijkstra.h"

using namespace tools_for_kssp;

namespace kssp {
  
  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class ParsimoniousSidetrackBased
  {
  public:
    ParsimoniousSidetrackBased() = default;

    // Constructor
    explicit ParsimoniousSidetrackBased(directed_graph::DirectedGraph<TI,TV> *g,
					const TI ssource,
					const TI ttarget,
					const size_t version,
					bool verbose,
					float tthreshold=1);

    // Destructor
    virtual ~ParsimoniousSidetrackBased();

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
    float threshold;

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
  ParsimoniousSidetrackBased<TI, TV>::ParsimoniousSidetrackBased(directed_graph::DirectedGraph<TI,TV> *g,
								 const TI ssource, const TI ttarget,
								 const size_t version,
								 bool verbose,
								 float tthreshold):
    verbose(verbose), graph(g), source(ssource), target(ttarget), version(version), threshold(tthreshold)
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
  ParsimoniousSidetrackBased<TI, TV>::~ParsimoniousSidetrackBased()
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
  void ParsimoniousSidetrackBased<TI, TV>::report()
  {
    if (verbose == true)
      {
	std::cout << "# good = " << heap_candidate_paths->size();
	std::cout << "\t# hard = " << heap_hard_candidate_paths->size();
	std::cout << "\t# trees = " << cpt_used_trees << " / " << idx_to_tree.size();
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
  inline bool ParsimoniousSidetrackBased<TI, TV>::empty()
  {
    return heap_candidate_paths->empty() && heap_hard_candidate_paths->empty();
  }


  /*
   * Check whether there is a candidate path and return it
   */
  template<typename TI, typename TV>
  std::pair<std::vector<TI>,TV> ParsimoniousSidetrackBased<TI, TV>::next_path()
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
	    // Extract the next best simple path from the heap
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
  void ParsimoniousSidetrackBased<TI,TV>::extract_path(CandidateSimple<TI,TV> *candidate_path)
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

	// Update weight, follow edge and go to next tree
	extracted_path_weight -= idx_to_tree[current_tree]->weight(tail);
	extracted_path_weight += graph->edges[edge];  // weight of edge
	extracted_path.push_back(edge.second);  // head of edge
	idx_tree++;
	current_tree = candidate_path->tree_index_sequence[idx_tree];
      }
    
    // We follow the last part of the path to target
    
    if (idx_to_tree[current_tree] == nullptr)
      { // first time we use this tree. We create it
	TI tmp = extracted_path.back();
	extracted_path.pop_back();
	if (version == 1)
	  { // Build a fresh new tree and disable prefix
	    idx_to_tree[current_tree] = new dijkstra::LazyDijkstra<TI,TV>(graph, &extracted_path, target);
	  }
	else
	  {
	    size_t ti = (idx_tree == 0 ? 0 : candidate_path->tree_index_sequence[idx_tree - 1]);
	    idx_to_tree[current_tree] = new dijkstra::LazyDijkstra<TI,TV>(idx_to_tree[ti], &extracted_path);
	  }
	cpt_used_trees++;
	extracted_path.push_back(tmp);
      }
    
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
  void ParsimoniousSidetrackBased<TI, TV>::compute_deviations(CandidateSimple<TI,TV> *candidate_path)
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
    std::vector<std::pair<size_t,TI> > bad_sidetracks;
    bad_sidetracks.clear();
    for (auto const& edge: sidetracks)
      {
	TI head = edge.second;
	if (color[head] == MAX_COLOR)
	  continue;
	size_t i = edge.first;
	if (color[head] <= i)
	  {
	    if (T->weight(head) != MAX_WEIGHT)
	      bad_sidetracks.push_back(edge);
	  }
	else
	  { // We have a simple path to target in T, so a new simple path
	    TI u = extracted_path[i];
	    CandidateSimple<TI,TV> *new_candidate = new CandidateSimple<TI,TV>(candidate_path);
	    new_candidate->weight = extracted_path_weight - T->weight(u) + graph->edges[std::make_pair(u, head)] + T->weight(head);
	    new_candidate->dev_idx = i + 1;
	    new_candidate->sidetrack_sequence.push_back(std::make_pair(u, head));
	    new_candidate->tree_index_sequence.push_back(new_candidate->tree_index_sequence.back());
	    heap_candidate_paths->push(new_candidate, new_candidate->weight);
	  }
      }

    // 4. Deal with non simple paths
    if (not bad_sidetracks.empty())
      {
	// We determine the lowerbound over all bad sidetracks of candidate path weight
	TV lbw = MAX_WEIGHT;
	size_t idx = bad_sidetracks.back().first;
	for (auto const& edge: bad_sidetracks)
	  {
	    TI u = extracted_path[edge.first];
	    TI v = edge.second;
	    TV w = extracted_path_weight - T->weight(u) + graph->edges[std::make_pair(u, v)] + T->weight(v);
	    if (w < lbw)
	      {
		lbw = w;
		idx = edge.first;
	      }
	  }
	// We create a new hard candidate
	CandidateHard<TI,TV> *hard = new CandidateHard<TI,TV>(lbw, idx, *candidate_path, bad_sidetracks, MAX_TREE_IDX);
	heap_hard_candidate_paths->push(hard, lbw);
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
  void ParsimoniousSidetrackBased<TI, TV>::compute_colors(CandidateSimple<TI,TV> *candidate_path,
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
   * Starting from the last sidetrack, this method computes the exact weight
   * of a simple path using that sidetrack. If such a path exists, a new simple
   * candidate is created and stored. Otherwise, that sidetrack is discarded.
   * The method ends either when all sidetracks have been considered, or we
   * reached the sidetrack used to weight this hard candidate. In such case, we
   * compute a new lowerbound using remaining sidetracks and push this modified
   * hard candidate to heap_hard_candidate_paths.
   */
  template<typename TI, typename TV>
  void ParsimoniousSidetrackBased<TI, TV>::treat_one_hard_candidate()
  {
    if (heap_hard_candidate_paths->empty())
      return;  // should not happen

    CandidateHard<TI,TV> *hard_candidate = heap_hard_candidate_paths->top_item();
    heap_hard_candidate_paths->pop();
    
    extract_path(&(hard_candidate->candidate));

    std::vector<TI> prefix;
    prefix.clear();
    prefix.push_back(extracted_path[0]);
    for (size_t i = 0; i <= hard_candidate->sidetracks.back().first; i++)
      prefix.push_back(extracted_path[i]);

    dijkstra::LazyDijkstra<TI,TV> *T_orig = idx_to_tree[hard_candidate->candidate.tree_index_sequence.back()];
    dijkstra::LazyDijkstra<TI,TV> *T_hard;
    if (hard_candidate->tree_index == MAX_TREE_IDX)
      { // We need a new (working) tree

	if (version == 1)
	  { // Build a fresh new tree and disable prefix
	    T_hard = new dijkstra::LazyDijkstra<TI,TV>(graph, &prefix, target);
	  }
	else
	  { // Better way to copy and disable
	    T_hard = new dijkstra::LazyDijkstra<TI,TV>(T_orig, &prefix);
	  }
	hard_candidate->tree_index = idx_to_tree.size();
	idx_to_tree.push_back(T_hard);
	cpt_used_trees++;
      }
    else
      { // We reuse the last working tree. It is ready for use
	T_hard = idx_to_tree[hard_candidate->tree_index];
      }
    bool T_hard_used = false; // Set to true is T_hard is used in a simple path

    // get vertex used to set weight of this hard_candidate
    TI u_min = extracted_path[hard_candidate->idx_min];
    TI last_spur = n;
    size_t current_tree_index = hard_candidate->tree_index; // initialization not needed

    while ((not hard_candidate->sidetracks.empty()) &&
	   (hard_candidate->sidetracks.back().first >= hard_candidate->idx_min))
      {
	size_t idx = hard_candidate->sidetracks.back().first;
	TI u = extracted_path[idx];
	TI head = hard_candidate->sidetracks.back().second;
	hard_candidate->sidetracks.pop_back();
	if (T_hard->successor(head) != head)
	  {
	    // We have a simple path. We push it to heap_candidate_paths
	    CandidateSimple<TI,TV> *new_candidate = new CandidateSimple<TI,TV>(&(hard_candidate->candidate));
	    new_candidate->weight = extracted_path_weight - T_orig->weight(u) + graph->edges[std::make_pair(u, head)] + T_hard->weight(head);
	    new_candidate->dev_idx = idx + 1;
	    new_candidate->sidetrack_sequence.push_back(std::make_pair(u, head));
	    if (u == u_min)
	      { // We use the current working tree
		new_candidate->tree_index_sequence.push_back(hard_candidate->tree_index);
		T_hard_used = true;
	      }
	    else
	      {
		if (u == last_spur)
		  { // we reuse the last tree (same index)
		    new_candidate->tree_index_sequence.push_back(current_tree_index);
		  }
		else
		  { // We "create" a new tree
		    current_tree_index = idx_to_tree.size();
		    new_candidate->tree_index_sequence.push_back(current_tree_index);
		    idx_to_tree.push_back(nullptr);
		    last_spur = u;
		  }
	      }
	    heap_candidate_paths->push(new_candidate, new_candidate->weight);
	  }

	// Should we update the current working tree ?
	if ((not T_hard_used) && (not hard_candidate->sidetracks.empty()))
	  {
	    size_t next_idx = hard_candidate->sidetracks.back().first;
	    if (next_idx != idx)
	      {
		// We enable some vertices
		TI v = extracted_path[next_idx];
		while (prefix.back() != v)
		  {
		    T_hard->enable(prefix.back());
		    prefix.pop_back();
		  }
	      }
	  }
      }

    if (not hard_candidate->sidetracks.empty())
      { // update the hard candidate and push it back to heap
	
	// We determine the lowerbound over all bad sidetracks of path weight
	TV lbw = MAX_WEIGHT;
	size_t idx = hard_candidate->sidetracks.back().first;
	for (auto const& edge: hard_candidate->sidetracks)
	  {
	    TI u = extracted_path[edge.first];
	    TI v = edge.second;
	    TV w = extracted_path_weight - T_orig->weight(u) + graph->edges[std::make_pair(u, v)] + T_orig->weight(v);
	    if (w < lbw)
	      {
		lbw = w;
		idx = edge.first;
	      }
	  }
	hard_candidate->weight = lbw;
	hard_candidate->idx_min = idx;
	if (T_hard_used)
	  { // We will need a new tree
	    hard_candidate->tree_index = MAX_TREE_IDX;
	  }
	heap_hard_candidate_paths->push(hard_candidate, lbw);
      }
    else
      { // We are done with this hard candidate
	if (not T_hard_used)
	  { // This tree has never been used, so we can delete it to save memory
	    delete T_hard;
	    idx_to_tree[hard_candidate->tree_index] = nullptr;
	  }
	delete hard_candidate;
      }
  }
  
} // end namespace
#endif

