/*
 * Dijkstra
 *
 * Implement standard Dijkstra single source shortest path algorithm.
 * By default, it build shortest path tree from specified source vertex encoded
 * with predecessor function. Setting parameter reverse to true, it builds
 * shortest path tree to target with successor function.
 */


#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <cstdint>
#include <set>
#include <assert.h>     /* assert */
#include <limits>       /* used to get maximum possible value of type TV */
#include "digraph.h"
#include "pairing_heap.h"


namespace dijkstra {

  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class Dijkstra
  {
  public:
    Dijkstra() = default;

    // Constructor
    explicit Dijkstra(directed_graph::DirectedGraph<TI,TV> *g,
		      const TI &source,
		      const bool &reverse);

    // Constructor
    explicit Dijkstra(directed_graph::DirectedGraph<TI,TV> *g,
		      std::set<TI> forbidden_vertices,
		      std::set<std::pair<TI,TI> > forbidden_edges,
		      const TI &source,
		      const bool &reverse);

    // Destructor
    virtual ~Dijkstra();

    // Compute shortest paths until target is reached
    bool run(const TI &target);

    // Get distance from/to source
    TV weight(const TI &u);

    // Get predecessor of u in shortest path from source
    TI predecessor(const TI &u);

    // Get successor of u in shortest path to source
    TI successor(const TI &u);

    // Get path from/to u
    std::vector<TI> get_path(const TI &u);

  private:
    directed_graph::DirectedGraph<TI,TV> *graph;
    TI n;
    TI source;
    TI *_predecessor;
    TV *_weight;
    std::set<TI> f_vertices;
    std::set<std::pair<TI,TI> > f_edges;
    bool reverse;
    heap::PairingHeap<TI,TV> *pq;

    // Initialize data structures
    void initialize();

    // Compute shortest paths from source
    void runFrom(const TI &target);

    // Compute shortest paths to source
    void runTo(const TI &target);
  };


  // Constructor
  template<typename TI, typename TV>
  Dijkstra<TI, TV>::Dijkstra(directed_graph::DirectedGraph<TI,TV> *g,
			     const TI &source,
			     const bool &reverse):
    graph(g), n(g->n), source(source), reverse(reverse)
  {
    initialize();
  }

 
  // Constructor
  template<typename TI, typename TV>
  Dijkstra<TI, TV>::Dijkstra(directed_graph::DirectedGraph<TI,TV> *g,
			     std::set<TI> forbidden_vertices,
			     std::set<std::pair<TI,TI> > forbidden_edges,
			     const TI &source,
			     const bool &reverse):
    graph(g), n(g->n), source(source),
    f_vertices(forbidden_vertices), f_edges(forbidden_edges), reverse(reverse)
  {
    initialize();
  }

  
  // Destructor
  template<typename TI, typename TV>
  Dijkstra<TI, TV>::~Dijkstra()
  {
    free(_predecessor);
    free(_weight);
    delete pq;
  }

 
  // Initialize data structures
  template<typename TI, typename TV>
  void Dijkstra<TI, TV>::initialize()
  {
    _predecessor = (TI*)malloc(n * sizeof(TI));
    _weight = (TV*)malloc(n * sizeof(TV));
    TV MAX_WEIGHT = std::numeric_limits<TV>::max();
    for (TI u = 0; u < n; u++)
      {
	_predecessor[u] = u;
	_weight[u] = MAX_WEIGHT;
      }
    pq = new heap::PairingHeap<TI,TV>();
    pq->push(source, 0);
    _weight[source] = 0;
  }


  // Compute shortest paths
  template<typename TI, typename TV>
  bool Dijkstra<TI, TV>::run(const TI &target)
  {
    if (reverse)
      runTo(target);
    else
      runFrom(target);
    return true;
  }


  // Compute shortest paths from source until target is reached
  template<typename TI, typename TV>
  void Dijkstra<TI, TV>::runFrom(const TI &target)
  {
    while (not pq->empty())
      {
	TI u = pq->top_item();
	if (u == target)
	  break;

	pq->pop();

	for(auto const& value: graph->out_neighbors[u])
	  {
	    TI v = value.first;
	    TV w = value.second;
	    if ((f_vertices.count(v) > 0) || (f_edges.count(std::make_pair(u, v)) > 0))
	      continue;
	    if (_weight[u] + w < _weight[v])
	      {
		_weight[v] = _weight[u] + w;
		if (_predecessor[v] == v)
		  pq->push(v, _weight[v]);
		else
		  pq->decrease(v, _weight[v]);		  
		_predecessor[v] = u;
	      }
	  }
      }
  }

  // Compute shortest paths to source until target is reached
  template<typename TI, typename TV>
  void Dijkstra<TI, TV>::runTo(const TI &target)
  {
    while (not pq->empty())
      {
	TI u = pq->top_item();
	if (u == target)
	  break;

	pq->pop();

	for(auto const& value: graph->in_neighbors[u])
	  {
	    TI v = value.first;
	    TV w = value.second;
	    if ((f_vertices.count(v) > 0) || (f_edges.count(std::make_pair(u, v)) > 0))
	      continue;
	    if (_weight[u] + w < _weight[v])
	      {
		_weight[v] = _weight[u] + w;
		if (_predecessor[v] == v)
		  pq->push(v, _weight[v]);
		else
		  pq->decrease(v, _weight[v]);		  
		_predecessor[v] = u;
	      }
	  }
      }
  }

  // Get distance from/to source
  template<typename TI, typename TV>
  inline TV Dijkstra<TI, TV>::weight(const TI &u)
  {
    return _weight[u];
  }

  // Get predecessor of u in shortest path from source
  template<typename TI, typename TV>
  inline TI Dijkstra<TI, TV>::predecessor(const TI &u)
  {
    assert(not reverse);
    return _predecessor[u];
  }

  // Get successor of u in shortest path to source
  template<typename TI, typename TV>
  inline TI Dijkstra<TI, TV>::successor(const TI &u)
  {
    assert(reverse);
    return _predecessor[u];
  }

  // Get path from/to u
  template<typename TI, typename TV>
  inline std::vector<TI> Dijkstra<TI, TV>::get_path(const TI &u)
  {
    std::vector<TI> path;
    path.clear();

    path.push_back(u);
    TI v = u;
    while (_predecessor[v] != v)
      {
        path.push_back(_predecessor[v]);
        v = _predecessor[v];
      }

    if (not reverse)
      {
        // We revert the path
        TI begin = 0;
        TI end = path.size() -1;
        while (begin < end)
	  {
	    TI tmp = path[begin];
	    path[begin] = path[end];
	    path[end] = tmp;
            begin++;
            end--;
	  }
      }

    return path;
  }


}
#endif
