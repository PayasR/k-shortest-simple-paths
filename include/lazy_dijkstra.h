/*
 * Lazy Dijkstra
 *
 * Implement a lazy Dijkstra single source shortest path algorithm **to** target
 * (i.e., we build successor array) with updates. It supports the following
 * operations:
 * - copy and disable vertices: make a copy of the data structure and disable 
 *   vertices on the way. Update the shortest paths tree after the "removal" of
 *   a set of vertices, from the graph
 * - enable a vertex: update shortest paths tree after the "addition" of a
 *   vertex to the graph
 */


#ifndef LAZY_DIJKSTRA_H
#define LAZY_DIJKSTRA_H

#include <cstdint>
#include <cstring>
#include <iostream>
#include <assert.h>     /* assert */
#include <limits>       /* used to get maximum possible value of type TV */
#include "digraph.h"
#include "pairing_heap.h"

namespace dijkstra {

  template<
    typename TI,  // type of node indices (unsigned int, uint32_t, etc.)
    typename TV   // type of edge weights (int, unsigned int, float, etc.)
    >
  class LazyDijkstra
  {
  public:
    LazyDijkstra() = default;

    // Constructor
    explicit LazyDijkstra(directed_graph::DirectedGraph<TI,TV> *g,
			  const TI target);

    // Constructor with forbidden vertices
    explicit LazyDijkstra(directed_graph::DirectedGraph<TI,TV> *g,
			  std::vector<TI> *forbidden_vertices,
			  const TI target);

    // Copy constructor
    LazyDijkstra(LazyDijkstra const *other);

    // Copy constructor that disables vertices on the way
    LazyDijkstra(LazyDijkstra const *other, std::vector<TI> *forbidden_vertices);

    // Destructor
    virtual ~LazyDijkstra();

    // Enable a vertex that was forbidden before
    void enable(const TI &vertex);

    // Get distance to target
    TV weight(const TI &u);

    // Get successor of u in shortest path to target
    TI successor(const TI &u);

    // Get path from u to target
    std::vector<TI> get_path(const TI &u);

    // Get successor of u in shortest path to target
    void show();

  private:
    directed_graph::DirectedGraph<TI,TV> *graph;
    TI n;
    TI target;
    heap::PairingHeap<TI, TV> *pq;
    bool *seen;
    bool *f_vertices;
    TI *_successor;
    TV *_weight;

    // Compute shortest paths to target until vertex is found or the weight bound is reached
    void runTo(const TI vertex, const TV weight_bound);

    TV MAX_WEIGHT;
  };


  /*
   * Constructor
   */
  template<typename TI, typename TV>
  LazyDijkstra<TI, TV>::LazyDijkstra(directed_graph::DirectedGraph<TI,TV> *g,
				     const TI target):
    graph(g), n(g->n), target(target)
  {
    MAX_WEIGHT = std::numeric_limits<TV>::max();
    seen = (bool*)malloc(n * sizeof(bool));
    f_vertices = (bool*)malloc(n * sizeof(bool));
    _successor = (TI*)malloc(n * sizeof(TI));
    _weight = (TV*)malloc(n * sizeof(TV));
    for (TI i = 0; i < n; i++)
      {
	seen[i] = false;
	f_vertices[i] = false;
	_successor[i] = i;
	_weight[i] = MAX_WEIGHT;
      }

    pq = new heap::PairingHeap<TI, TV>();
    pq->push(target, 0);
    _weight[target] = 0;
  }

  /*
   * Constructor with forbidden vertices
   */
  template<typename TI, typename TV>
  LazyDijkstra<TI, TV>::LazyDijkstra(directed_graph::DirectedGraph<TI,TV> *g,
				     std::vector<TI> *forbidden_vertices,
				     const TI target):
    graph(g), n(g->n), target(target)
  {
    MAX_WEIGHT = std::numeric_limits<TV>::max();
    seen = (bool*)malloc(n * sizeof(bool));
    f_vertices = (bool*)malloc(n * sizeof(bool));
    _successor = (TI*)malloc(n * sizeof(TI));
    _weight = (TV*)malloc(n * sizeof(TV));
    std::memset(seen, 0, n * sizeof(bool));
    std::memset(f_vertices, 0, n * sizeof(bool));
    for (size_t i = 0; i < n; i++)
      {
	_successor[i] = i;
	_weight[i] = MAX_WEIGHT;
      }
    for (TI u: *forbidden_vertices)
      f_vertices[u] = true;
    pq = new heap::PairingHeap<TI, TV>();
    pq->push(target, 0);
    _weight[target] = 0;
  }



  /*
   * Copy constructor
   */
  template<typename TI, typename TV>
  LazyDijkstra<TI, TV>::LazyDijkstra(LazyDijkstra const *other):
    graph(other->graph), n(other->n), target(other->target)
  {
    MAX_WEIGHT = std::numeric_limits<TV>::max();
    seen = (bool*)malloc(n * sizeof(bool));
    f_vertices = (bool*)malloc(n * sizeof(bool));
    _successor = (TI*)malloc(n * sizeof(TI));
    _weight = (TV*)malloc(n * sizeof(TV));
    std::memcpy(seen, other->seen, sizeof(bool) * n);
    std::memcpy(f_vertices, other->f_vertices, sizeof(bool) * n);
    std::memcpy(_successor, other->_successor, sizeof(TI) * n);
    std::memcpy(_weight, other->_weight, sizeof(TV) * n);
    pq = new heap::PairingHeap<TI, TV>(other->pq);
  }


  /*
   * Copy constructor that disables vertices on the way
   */
  template<typename TI, typename TV>
  LazyDijkstra<TI, TV>::LazyDijkstra(LazyDijkstra const *other, std::vector<TI> *forbidden_vertices):
    graph(other->graph), n(other->n), target(other->target)
  {
    // 1. create data structures
    MAX_WEIGHT = std::numeric_limits<TV>::max();
    seen = (bool*)malloc(n * sizeof(bool));
    f_vertices = (bool*)malloc(n * sizeof(bool));
    _successor = (TI*)malloc(n * sizeof(TI));
    _weight = (TV*)malloc(n * sizeof(TV));
    std::memcpy(seen, other->seen, sizeof(bool) * n);
    std::memcpy(f_vertices, other->f_vertices, sizeof(bool) * n);
    std::memcpy(_successor, other->_successor, sizeof(TI) * n);
    std::memcpy(_weight, other->_weight, sizeof(TV) * n);

    // 2. Identify impacted vertices
    std::vector<TI> todo;
    todo.clear();
    std::vector<TI> stack;
    stack.clear();
    for (TI const& u: *forbidden_vertices)
      {
	if (f_vertices[u])
	  continue;

	stack.push_back(u);
	while (not stack.empty())
	  {
	    TI v = stack.back();
	    stack.pop_back();
	    for (auto const& it: graph->in_neighbors[v])
	      {
		TI x = it.first;
		if (_successor[x] == v)
		  {
		    todo.push_back(x);
		    seen[x] = false;
		    _successor[x] = x;
		    _weight[x] = MAX_WEIGHT;
		    stack.push_back(x);
		  }
	      }
	  }

	f_vertices[u] = true;
	_successor[u] = u;
	_weight[u] = MAX_WEIGHT;
	seen[u] = false;
      }

    // 3. Update the status of impacted vertices
    for (TI const& u: todo)
      {
	// If u has a seen out neighbor, it must be in the heap
	for (auto const& it: graph->out_neighbors[u])
	  {
	    TI v = it.first;
	    if (seen[v])
	      {
		TV w = it.second;
		if (w + _weight[v] < _weight[u])
		  {
		    _weight[u] = w + _weight[v];
		    _successor[u] = v;
		  }
	      }
	  }
      }

    // 4. Build the heap
    pq = new heap::PairingHeap<TI, TV>();
    for (TI u = 0; u < n; u++)
      if ((not seen[u]) && (_successor[u] != u))
	pq->push(u, _weight[u]);
  }



  /*
   * Destructor
   */
  template<typename TI, typename TV>
  LazyDijkstra<TI, TV>::~LazyDijkstra()
  {
    free(seen);
    free(f_vertices);
    free(_successor);
    free(_weight);
    delete pq;
  }



  /*
   * Enable a vertex
   * After this method, some weights are no longer correct. A call to runTo is needed
   */
  template<typename TI, typename TV>
  void LazyDijkstra<TI, TV>::enable(const TI &vertex)
  {
    if (f_vertices[vertex])
      {
	f_vertices[vertex] = false;
	_weight[vertex] = MAX_WEIGHT;
	_successor[vertex] = vertex;
	seen[vertex] = false;
	for (auto const& it: graph->out_neighbors[vertex])
	  {
	    TI v = it.first;
	    if (seen[v])
	      {
		TV w = it.second;
		if (w + _weight[v] < _weight[vertex])
		  {
		    _weight[vertex] = w + _weight[v];
		    _successor[vertex] = v;
		  }
	      }
	  }
	if (_successor[vertex] != vertex)
	  pq->push(vertex, _weight[vertex]);
      }
  }


  /*
   * Compute shortest paths to target until vertex is found or the weight bound is reached
   */
  template<typename TI, typename TV>
  void LazyDijkstra<TI, TV>::runTo(const TI vertex, const TV weight_bound)
  {
    while (not pq->empty())
      {
	TI u = pq->top_item();
	seen[u] = true;
	pq->pop();

	for(auto const& it: graph->in_neighbors[u])
	  {
	    TI v = it.first;
	    if (f_vertices[v])
	      continue;
	    TV w = it.second;
	    if (w + _weight[u] < _weight[v])
	      {
		_weight[v] = w + _weight[u];
		_successor[v] = u;
		pq->decrease(v, _weight[v]);
		seen[v] = false;
	      }
	  }

	if ((u == vertex) || (_weight[u] > weight_bound))
	  break;
      }
  }


  /*
   * Get weight of a shortest path from u to target
   */
  template<typename TI, typename TV>
  inline TV LazyDijkstra<TI, TV>::weight(const TI &u)
  {
    // If u is not seen, we compute its path.
    // If the top of the heap has a small weight, an update has been done and might change this weight
    if (not seen[u])
      runTo(u, MAX_WEIGHT);
    else if ((not pq->empty()) && (pq->top_value() < _weight[u]))
      runTo(u, _weight[u]);
    return _weight[u];
  }


  /*
   * Get successor of u in shortest path to target
   */
  template<typename TI, typename TV>
  inline TI LazyDijkstra<TI, TV>::successor(const TI &u)
  {
    // If u is not seen, we compute its path.
    // If the top of the heap has a small weight, an update has been done and might change this weight
    if (not seen[u])
      runTo(u, MAX_WEIGHT);
    else if ((not pq->empty()) && (pq->top_value() < _weight[u]))
      runTo(u, _weight[u]);
    return _successor[u];
  }


  /*
   * Get path from u to target
   */
  template<typename TI, typename TV>
  inline std::vector<TI> LazyDijkstra<TI, TV>::get_path(const TI &u)
  {
    // If u is not seen, we compute its path.
    // If the top of the heap has a small weight, an update has been done and might change this weight
    if (not seen[u])
      runTo(u, MAX_WEIGHT);
    else if ((not pq->empty()) && (pq->top_value() < _weight[u]))
      runTo(u, _weight[u]);

    TI v;
    std::vector<TI> path;
    path.clear();

    path.push_back(u);
    v = u;
    while (_successor[v] != v)
      {
        path.push_back(_successor[v]);
        v = _successor[v];
      }

    return path;
  }

  /*
   * Get successor of u in shortest path to target
   */
  template<typename TI, typename TV>
  void LazyDijkstra<TI, TV>::show()
  {
    std::cout << "succ:";
    for (TI u = 0; u < n; u++)
      if (f_vertices[u])
	std::cout << "  f("<<u<<")";
      else
	{
	  if (seen[u])
	    std::cout << "  "<<u<<" --> "<<_successor[u];
	  else
	    std::cout << "  ("<<u<<" --> "<<_successor[u]<<")";
	}
    std::cout << std::endl;
  }


}
#endif
