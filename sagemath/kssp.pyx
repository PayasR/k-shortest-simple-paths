# -*- coding: utf-8 -*-
# cython: binding=True
# distutils: language = c++
r"""
Interface for k shortest simple paths algorithms

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

from sagemath.types cimport vertex_t, weight_t

from sagemath.digraph cimport DirectedGraph

cdef dummy_yen():
    cdef DirectedGraph[vertex_t,weight_t] g
    cdef Yen[vertex_t,weight_t] D = Yen[vertex_t,weight_t](&g, 0, 0)

cdef dummy_node_classification():  # Feng
    cdef DirectedGraph[vertex_t,weight_t] g
    cdef NodeClassification[vertex_t,weight_t] D = NodeClassification[vertex_t,weight_t](&g, 0, 0)

cdef dummy_sidetrack_based():      # Kurz and Mutzel
    cdef DirectedGraph[vertex_t,weight_t] g
    cdef SidetrackBased[vertex_t,weight_t] D = SidetrackBased[vertex_t,weight_t](&g, 0, 0, 1, False)

cdef dummy_parsimonious_sidetrack_based():  # Al Zoobi, Coudert and Nisse
    cdef DirectedGraph[vertex_t,weight_t] g
    cdef ParsimoniousSidetrackBased[vertex_t,weight_t] D = ParsimoniousSidetrackBased[vertex_t,weight_t](&g, 0, 0, 1, False)

