"""
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
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set

from sagemath.types cimport vertex_t, weight_t

cdef extern from "../include/digraph.h":
    pass

# Declare the class with cdef
cdef extern from "../include/digraph.h" namespace "directed_graph":
    cdef cppclass DirectedGraph[TI, TV]:
        DirectedGraph()

        TI n
        TI m
        cpp_map[pair[TI, TI], TV] edges
        vector[vector[pair[TI,TV]]] out_neighbors
        vector[vector[pair[TI,TV]]] in_neighbors

        DirectedGraph(TI) except +
        void reset(TI)
        TI order()
        TI size()
        void add_edge(TI, TI, TV)
        bint has_edge(TI, TI)
        TI out_degree(TI)
        TI in_degree(TI)
        TV weight(TI, TI)


cdef tuple digraph_from_edgelist(DirectedGraph[vertex_t,weight_t] *g, V, E, directed=?)

cdef tuple digraph_from_dimacs_file(DirectedGraph[vertex_t,weight_t] *g, filename)


cpdef my_test_read_dimacs(filename)


