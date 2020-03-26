# -*- coding: utf-8 -*-
# cython: binding=True
# distutils: language = c++
r"""
MyDiGraph


Methods
-------
"""
# ****************************************************************************
# Copyright (C) 2019 Ali Al Zoobi <ali.al-zoobi@inria.fr>
#                    David Coudert <david.coudert@inria.fr>
#                    Nicolas Nisse <nicolas.nisse@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.set cimport set as cpp_set

from cysignals.signals cimport sig_on, sig_off, sig_check

from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph


from sagemath.types cimport vertex_t, weight_t

from sagemath.digraph cimport DirectedGraph
from sagemath.digraph cimport digraph_from_dimacs_file
from sagemath.digraph cimport digraph_from_edgelist

from sagemath.dijkstra cimport Dijkstra, LazyDijkstra


from sagemath.kssp cimport Yen
from sagemath.kssp cimport NodeClassification  # Feng
from sagemath.kssp cimport SidetrackBased      # Kurz and Mutzel
from sagemath.kssp cimport ParsimoniousSidetrackBased  # Al Zoobi, Coudert and Nisse

# ==============================================================================
# Class MyDiGraph
# ==============================================================================

cdef class MyDiGraph:
    """
    """

    def __cinit__(self, data=None, by_weight=False, weight_function=None):
        """
        Constructor

        INPUT:

        - ``data`` -- one of the following format:

          - a string encoding the full path to a file in DIMACS format

          - a Sage Graph or DiGraph

        - ``by_weight`` -- boolean (default: ``False``); whether to consider the
          graph as weighted or unweighted

        - ``weight_function`` -- function (default: ``None``); a function that
          inputs an edge and outputs a weight (unsigned long). Ignored when
          ``by_weight`` is ``False``.
        """
        self.n = 0
        self.m = 0
        self.int_to_vertex = list()
        self.vertex_to_int = dict()

        if isinstance(data, (str, bytes)):
            self._read_from_file(data)
        elif isinstance(data, DiGraph):
            self._from_sage_graph(data, by_weight=by_weight, weight_function=weight_function, directed=True)
        elif isinstance(data, Graph):
            self._from_sage_graph(data, by_weight=by_weight, weight_function=weight_function, directed=False)
        else:
            raise ValueError("unrecognized input format")
        
        self.n = self.g.n
        self.m = self.g.m

    def __dealloc__(self):
        """
        Destructor
        """
        self.g.reset(0)

    def __repr__(self):
        """
        Return a string representation of ``self``.
        """
        return "MyDiGraph from {}".format(self.name)


    def _read_from_file(self, filename):
        """
        Load the graph from file
        """
        try:
            self.int_to_vertex, self.vertex_to_int = digraph_from_dimacs_file(&(self.g), filename)
            self.name = "file {}".format(filename)
        except:
            raise ValueError("unrecognized file format")

    def _from_sage_graph(self, G, by_weight=False, weight_function=None, directed=True):
        """
        """
        if by_weight:
            if weight_function is None:
                E = G.edge_iterator(labels=True)
            else:
                E = [(e[0], e[1], weight_function(e)) for e in G.edge_iterator(labels=True)]
        else:
            E = [(u, v, 1) for u, v in G.edge_iterator(labels=False)]
        self.int_to_vertex, self.vertex_to_int = digraph_from_edgelist(&(self.g), G, E, directed=False)
        self.name = "Sage {}".format(G.name())


    def __len__(self):
        """
        Return the number of vertices
        """
        return self.n

    def order(self):
        """
        Return the number of vertices
        """
        return self.n

    def size(self):
        """
        Return the number of edges
        """
        return self.m

    def __iter__(self):
        """
        Iterator over the vertices of self
        """
        yield from self.int_to_vertex

    def __bool__(self):
        """
        Check whether self is empty
        """
        return bool(self.n)

    def __contains__(self, u):
        """
        Check wheter u is a vertex of self
        """
        return u in self.vertex_to_int

    def has_vertex(self, u):
        """
        Check wheter u is a vertex of self
        """
        return u in self.vertex_to_int

    def _error_if_vertex_not_in(self, u):
        """
        Raise Error if u not in self
        """
        if u not in self.vertex_to_int:
            raise ValueError("unknown vertex {}".format(u))

    def has_edge(self, u, v=None, label=None):
        """
        Check whether self has edge (u, v, label).
        """
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except Exception:
                    u, v = u
                    label = None
        cdef vertex_t i, j
        try:
            i = self.vertex_to_int[u]
            j = self.vertex_to_int[v]
        except:
            return False
        if not self.g.has_edge(i, j):
            return False
        if label is None:
            return True
        return self.g.weight(i, j) == label

    def vertices(self):
        """
        Iterator over the vertices of self
        """
        yield from self

    def edges(self, labels=True):
        """
        Iterator of the edges of self
        """
        cdef vertex_t u, uu, v
        cdef weight_t w

        if labels:
            for u, uu in enumerate(self.int_to_vertex):
                for v, w in self.g.out_neighbors[u]:
                    yield uu, self.int_to_vertex[v], w
        else:
            for u, uu in enumerate(self.int_to_vertex):
                for v, w in self.g.out_neighbors[u]:
                    yield uu, self.int_to_vertex[v]

    def out_degree(self, u=None):
        """
        Return the out-degree of u.
        """
        if u is None:
            u = self
        elif u in self:
            self._error_if_vertex_not_in(u)
            return self.g.out_neighbors[self.vertex_to_int[u]].size()                
        else:
            u = [v for v in u if v in self]
        return {v: self.g.out_degree(v) for v in u}

    def in_degree(self, u=None):
        """
        Return the in-degree of u.
        """
        if u is None:
            u = self
        elif u in self:
            self._error_if_vertex_not_in(u)
            return self.g.in_neighbors[self.vertex_to_int[u]].size()
        else:
            u = [v for v in u if v in self]
        return {v: self.g.in_degree(v) for v in u}
        
    def out_neighbors(self, u):
        """
        Iterator over the out-neighbors of u.
        """
        self._error_if_vertex_not_in(u)
        cdef vertex_t v, w
        for v, w in self.g.out_neighbors[self.vertex_to_int[u]]:
            yield self.int_to_vertex[v]

    def in_neighbors(self, u):
        """
        Iterator over the in-neighbors of u.
        """
        self._error_if_vertex_not_in(u)
        cdef vertex_t v
        cdef weight_t w
        for v, w in self.g.in_neighbors[self.vertex_to_int[u]]:
            yield self.int_to_vertex[v]

    def to_sage_graph(self):
        """
        """
        cdef vertex_t u, v
        cdef weight_t w
        G = DiGraph()
        for u in range(self.n):
            for v, w in self.g.out_neighbors[u]:
                G.add_edge(self.int_to_vertex[u], self.int_to_vertex[v], w)
        return G


    # ==========================================================================
    # Dijkstra
    # ==========================================================================

    def dijkstra(self, source, target=None, forbidden_vertices=None, forbidden_edges=None, reverse=False):
        """
        """
        self._error_if_vertex_not_in(source)
        cdef vertex_t s = self.vertex_to_int[source]
        cdef vertex_t t = -1
        if self.has_vertex(target):
            t = self.vertex_to_int[target]

        cdef vertex_t u, v
        cdef cpp_set[vertex_t] f_vertices
        cdef cpp_set[pair[vertex_t,vertex_t]] f_edges
        f_vertices.clear()
        f_edges.clear()

        if forbidden_vertices:
            for u in forbidden_vertices:
                f_vertices.insert(self.vertex_to_int[u])
        if forbidden_edges:
            for u, v in forbidden_edges:
                f_edges.insert((self.vertex_to_int[u], self.vertex_to_int[u]))

        cdef Dijkstra[vertex_t,weight_t] *D = new Dijkstra[vertex_t,weight_t](&(self.g), f_vertices, f_edges, s, <bint>reverse)
        sig_on()
        D.run(t)
        sig_off()

        if target is None:
            if reverse:
                L = {self.int_to_vertex[u]: self.int_to_vertex[D.successor(u)] for u in range(self.n)
                         if not f_vertices.count(u)}
            else:
                L = {self.int_to_vertex[u]: self.int_to_vertex[D.predecessor(u)] for u in range(self.n)
                         if not f_vertices.count(u)}
            return (L, {self.int_to_vertex[u]: D.weight(u) for u in range(self.n)
                            if not f_vertices.count(u)})

        path = D.get_path(t)
        res = [self.int_to_vertex[u] for u in path], D.weight(t)
        del D
        return res


    def lazy_dijkstra(self, source, target=None, forbidden_vertices=None):
        """
        """
        self._error_if_vertex_not_in(source)
        cdef vertex_t s = self.vertex_to_int[source]
        cdef vertex_t t = -1
        if self.has_vertex(target):
            t = self.vertex_to_int[target]

        cdef size_t u, v
        cdef vector[vertex_t] f_vertices
        cdef LazyDijkstra[vertex_t,weight_t] *D

        if forbidden_vertices:
            f_vertices.clear()
            for u in forbidden_vertices:
                f_vertices.push_back(self.vertex_to_int[u])
            D = new LazyDijkstra[vertex_t,weight_t](&(self.g), &f_vertices, s)

        else:
            D = new LazyDijkstra[vertex_t,weight_t](&(self.g), s)

        if target is None:
            s_vertices = set(f_vertices)
            return ({self.int_to_vertex[u]: self.int_to_vertex[D.successor(u)] for u in range(self.n)
                         if u not in s_vertices},
                    {self.int_to_vertex[u]: D.weight(u) for u in range(self.n)
                             if u not in s_vertices})

        path = D.get_path(t)
        res = [self.int_to_vertex[u] for u in path], D.weight(t)
        del D
        return res

    # ==========================================================================
    # k shortest simple paths iterators
    # ==========================================================================

    def yen_iterator(self, source, target, version=1):
        """
        """
        self._error_if_vertex_not_in(source)
        self._error_if_vertex_not_in(target)
        cdef vertex_t s = self.vertex_to_int[source]
        cdef vertex_t t = self.vertex_to_int[target]
        cdef Yen[vertex_t,weight_t] *Y = new Yen[vertex_t,weight_t](&(self.g), s, t)
        while not Y.empty():
            sig_on()
            path, cost = Y.next_path()
            sig_off()
            if cost:
                yield [self.int_to_vertex[u] for u in path], cost
        del Y

    def node_classification_iterator(self, source, target):
        """
        """
        self._error_if_vertex_not_in(source)
        self._error_if_vertex_not_in(target)
        cdef vertex_t s = self.vertex_to_int[source]
        cdef vertex_t t = self.vertex_to_int[target]
        cdef NodeClassification[vertex_t,weight_t] *NC = new NodeClassification[vertex_t,weight_t](&(self.g), s, t)
        while not NC.empty():
            sig_on()
            path, cost = NC.next_path()
            sig_off()
            if cost:
                yield [self.int_to_vertex[u] for u in path], cost
        del NC

    def sidetrack_based_iterator(self, source, target, version=1, verbose=False):
        """
        """
        self._error_if_vertex_not_in(source)
        self._error_if_vertex_not_in(target)
        cdef vertex_t s = self.vertex_to_int[source]
        cdef vertex_t t = self.vertex_to_int[target]
        cdef SidetrackBased[vertex_t,weight_t] *SB = new SidetrackBased[vertex_t,weight_t](&(self.g), s, t,
                                                                                               version, verbose)
        while not SB.empty():
            sig_on()
            path, cost = SB.next_path()
            sig_off()
            if cost:
                yield [self.int_to_vertex[u] for u in path], cost
        del SB

    def parsimonious_sidetrack_based_iterator(self, source, target, version=1, verbose=False):
        """
        """
        self._error_if_vertex_not_in(source)
        self._error_if_vertex_not_in(target)
        cdef vertex_t s = self.vertex_to_int[source]
        cdef vertex_t t = self.vertex_to_int[target]
        cdef ParsimoniousSidetrackBased[vertex_t,weight_t] *PSB = new ParsimoniousSidetrackBased[vertex_t,weight_t](&(self.g), s, t, version, verbose)
        while not PSB.empty():
            sig_on()
            path, cost = PSB.next_path()
            sig_off()
            if cost:
                yield [self.int_to_vertex[u] for u in path], cost
        del PSB
