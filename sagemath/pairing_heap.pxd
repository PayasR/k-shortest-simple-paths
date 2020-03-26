"""
Interface for Pairing Heap

[1] M. L. Fredman, R. Sedgewick, D. D. Sleator, and R. E. Tarjan.
    "The pairing heap: a new form of self-adjusting heap".
    Algorithmica. 1 (1): 111-129, 1986. doi:10.1007/BF01840439.

[2] See https://en.wikipedia.org/wiki/Pairing_heap for more details
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

cdef extern from "../include/pairing_heap.h":
    pass


# Declare the class with cdef
cdef extern from "../include/pairing_heap.h" namespace "heap":
    cdef cppclass PairingHeap[TI, TV]:
        PairingHeap() except +
        PairingHeap(PairingHeap[TI,TV]) except +
        bint empty()
        void reset()
        void push(TI, TV)
        pair[TI, TV] top()
        TI top_item()
        TV top_value()
        void pop()
        void decrease(TI, TV)
        bint contains(TI)
        TV value(TI)

