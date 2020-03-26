# -*- coding: utf-8 -*-
# cython: binding=True
# distutils: language = c++
r"""
Static digraph

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

from cysignals.signals cimport sig_on, sig_off, sig_check

from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.set cimport set as cpp_set
from libcpp.queue cimport priority_queue

from sagemath.types cimport vertex_t, weight_t

cdef tuple digraph_from_edgelist(DirectedGraph[vertex_t,weight_t] *g, V, E, directed=False):
    """
    Feed the digraph with input lists.

    INPUT:

    - ``g`` -- a ``DirectedGraph``

    - ``V`` -- list of vertices

    - ``E`` -- list of edges. Edge ``(u, v, w)``, with ``u`` and ``v`` in ``V``,
      indicates a direct edge from ``u`` to ``v`` with weight ``w``.

    - ``directed`` -- boolean (default: ``False``); whether the input list of
      edges is directed or not. If not, a pair of symmetric edges between ``u``
      to ``v`` with weight ``w`` is added to the graph for each edge ``(u, v,
      w)`` in ``E``.

    OUTPUT: a tuple ``(int_to_vertex, vertex_to_int)`` where

    - ``int_to_vertex`` is a mapping from integer in ``0..n-1`` to the vertices
      in ``V``.

    - ``vertex_to_int`` is a mapping from vertices in ``V`` to integers in
      ``0..n-1``.
    """
    cdef size_t i, ui, vi
    cdef list int_to_vertex = sorted(V)
    cdef dict vertex_to_int = {u: i for i, u in enumerate(int_to_vertex)}

    sig_on()
    g.reset(len(int_to_vertex))
    for u, v, w in set(E):  # use set to remove duplicates
        ui = vertex_to_int[u]
        vi = vertex_to_int[v]
        g.add_edge(ui, vi, w)
        if not directed:
            g.add_edge(vi, ui, w)
    sig_off()

    return int_to_vertex, vertex_to_int


cdef tuple digraph_from_dimacs_file(DirectedGraph[vertex_t,weight_t] *g, filename):
    """
    Read the digraph from file.

    En error is raised if something goes wrong, generally indicating that the
    specificied file in not in the DIMACS file format.  See
    http://users.diag.uniroma1.it/challenge9/format.shtml for more details on
    this format.

    INPUT:

    - ``g`` -- a ``DirectedGraph``

    - ``filename`` -- full path to a file encoding a digraph in DIMACS file
      format

    OUTPUT: a tuple ``(int_to_vertex, vertex_to_int)`` where

    - ``int_to_vertex`` is a mapping from integer in ``0..n-1`` to the vertices
      in ``V``.

    - ``vertex_to_int`` is a mapping from vertices in ``V`` to integers in
      ``0..n-1``.
    """
    cdef list int_to_vertex = []
    cdef dict vertex_to_int = dict()

    cdef size_t n, m, i
    cdef size_t u, v, w
    cdef vector[pair[size_t,size_t]] E
    cdef vector[size_t] W
    cdef set V = set()

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            elif line[0] == 'a':
                _, us, vs, ws = line.split()
                u, v, w = int(us), int(vs), int(ws)
                E.push_back((u, v))
                W.push_back(w)
                V.add(u)
                V.add(v)
            elif line[0] == 'p':
                _, _, ns, ms = line.split()
                n = int(ns)
                m = int(ms)
            elif line[0] == 'c':
                continue
            else:
                raise ValueError("unexpected line format: {}".format(line))

    if <size_t>len(V) != n or E.size() != m:
        raise ValueError("something goes wrong")

    int_to_vertex.extend(V)
    vertex_to_int.update({u: i for i, u in enumerate(int_to_vertex)})
    
    sig_on()
    g.reset(n)
    for i in range(m):
        u, v = E[i]
        if not g.has_edge(u, v):
            g.add_edge(vertex_to_int[u], vertex_to_int[v], W[i])
    sig_off()

    return int_to_vertex, vertex_to_int




# ==============================================================================
# Internal test methods
# ==============================================================================

cpdef my_test_read_dimacs(filename):
    """
    """
    cdef DirectedGraph[vertex_t,weight_t] g
    cdef list int_to_vertex
    cdef dict vertex_to_int

    int_to_vertex, vertex_to_int = digraph_from_dimacs_file(&g, filename)

    cdef size_t u, v, w
    for u in range(g.n):
        for v, w in g.out_neighbors[u]:
            print(u, v, w, g.weight(u, v))

