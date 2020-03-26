# -*- coding: utf-8 -*-
# cython: binding=True
# distutils: language = c++
r"""
Interface for Dijkstra

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

from sagemath.digraph cimport DirectedGraph

cdef dummy_dijkstra():
    cdef DirectedGraph[vertex_t,weight_t] g
    cdef cpp_set[vertex_t] FV
    cdef cpp_set[pair[vertex_t,vertex_t]] FE
    cdef Dijkstra[vertex_t,weight_t] D = Dijkstra[vertex_t,weight_t](&g, FV, FE, 0, True)

cdef dummy_lazy_dijkstra():
    cdef DirectedGraph[vertex_t,weight_t] g
    cdef cpp_set[vertex_t] FV
    cdef LazyDijkstra[vertex_t,weight_t] D = LazyDijkstra[vertex_t,weight_t](&g, 0)

