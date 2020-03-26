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

from sagemath.digraph cimport DirectedGraph
from sagemath.types cimport vertex_t, weight_t

cdef class MyDiGraph:
    """
    """
    cdef size_t n
    cdef size_t m
    cdef DirectedGraph[vertex_t,weight_t] g
    cdef list int_to_vertex
    cdef dict vertex_to_int
    cdef str name
