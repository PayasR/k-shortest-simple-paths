# distutils: language = c++
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

from cysignals.signals cimport sig_on, sig_off, sig_check

from libcpp.pair cimport pair

from sage.misc.prandom import shuffle

def test_simple():
    cdef PairingHeap[size_t, size_t] *PH = new PairingHeap[size_t, size_t]()
    cdef size_t i

    for i in range(10):
        PH.push(i, 10 + i)
        print(PH.top(), PH.top_item(), PH.top_value())

    print("-------")

    while(not PH.empty()):
        print(PH.top())
        PH.pop()

    print("-------")

    for i in range(10):
        PH.push(i, 10 + i)

    print(PH.top())
    print("-------")

    for i in range(10):
        PH.decrease(i, PH.value(i) - 1)

    while(not PH.empty()):
        print(PH.top())
        PH.pop()

    print("-------")
    
    for i in range(10):
        PH.push(i, 10 + i)

    print(PH.top())
    print("-------")


    sig_on()
    while(not PH.empty()):
        print(PH.top())
        PH.pop()
    sig_off()

    print("-------")

    
def test_random(n, k=10):
    cdef PairingHeap[size_t, size_t] *PH = new PairingHeap[size_t, size_t]()
    cdef size_t i
    cdef list L = list(range(n))
    
    shuffle(L)
    for i in L:
        PH.push(i, i + k)

    for _ in range(k):
        shuffle(L)
        for i in L:
            PH.decrease(i, PH.value(i) - 1)

    L = []
    while not PH.empty():
        L.append(PH.top_item())
        PH.pop()
    print(L)
    
