"""
Interface for k shortest simple paths algorithms
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

cdef extern from "../include/yen.h":
    pass

# Declare the class with cdef
cdef extern from "../include/yen.h" namespace "kssp":
    cdef cppclass Yen[TI, TV]:
        Yen()
        Yen(DirectedGraph[TI,TV] *, TI, TI) except +
        bint empty()
        pair[vector[TI],TV] next_path()


# Feng

cdef extern from "../include/node_classification.h":
    pass

# Declare the class with cdef
cdef extern from "../include/node_classification.h" namespace "kssp":
    cdef cppclass NodeClassification[TI, TV]:
        NodeClassification()
        NodeClassification(DirectedGraph[TI,TV] *, TI, TI) except +
        bint empty()
        pair[vector[TI],TV] next_path()


# Kurz and Mutzel

cdef extern from "../include/sidetrack_based.h":
    pass

# Declare the class with cdef
cdef extern from "../include/sidetrack_based.h" namespace "kssp":
    cdef cppclass SidetrackBased[TI, TV]:
        SidetrackBased()
        SidetrackBased(DirectedGraph[TI,TV] *, TI, TI, TI, bint) except +
        bint empty()
        pair[vector[TI],TV] next_path()


# Al Zoobi, Coudert and Nisse

cdef extern from "../include/parsimonious_sidetrack_based.h":
    pass

# Declare the class with cdef
cdef extern from "../include/parsimonious_sidetrack_based.h" namespace "kssp":
    cdef cppclass ParsimoniousSidetrackBased[TI, TV]:
        ParsimoniousSidetrackBased()
        ParsimoniousSidetrackBased(DirectedGraph[TI,TV] *, TI, TI, TI, bint) except +
        bint empty()
        pair[vector[TI],TV] next_path()

