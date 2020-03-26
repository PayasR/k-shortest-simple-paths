"""
Interface for Dijkstra
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

from sagemath.types cimport vertex_t, weight_t

from sagemath.digraph cimport DirectedGraph

cdef extern from "../include/dijkstra.h":
    pass

# Declare the class with cdef
cdef extern from "../include/dijkstra.h" namespace "dijkstra":
    cdef cppclass Dijkstra[TI, TV]:
        Dijkstra()
        Dijkstra(DirectedGraph[TI,TV] *, TI, bint) except +
        Dijkstra(DirectedGraph[TI,TV] *, cpp_set[TI], cpp_set[pair[TI,TI]], TI, bint) except +
        bint run(TI)
        TV weight(TI)
        TI predecessor(TI)
        TI successor(TI)
        vector[TI] get_path(TI)

ctypedef Dijkstra[vertex_t,weight_t] StandardDijkstra



cdef extern from "../include/lazy_dijkstra.h":
    pass

# Declare the class with cdef
cdef extern from "../include/lazy_dijkstra.h" namespace "dijkstra":
    cdef cppclass LazyDijkstra[TI, TV]:
        LazyDijkstra()
        LazyDijkstra(DirectedGraph[TI,TV] *, TI) except +
        LazyDijkstra(DirectedGraph[TI,TV] *, vector[TI] *, TI) except +

        enable(TI)

        TV weight(TI)
        TI successor(TI)
        vector[TI] get_path(TI)
