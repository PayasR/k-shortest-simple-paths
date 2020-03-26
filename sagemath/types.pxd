"""
Define types to use

See e.g. https://en.cppreference.com/w/cpp/types/integer
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

from libc.stdint cimport uint32_t, uint_fast32_t, UINT_FAST32_MAX

ctypedef uint_fast32_t vertex_t
ctypedef uint_fast32_t weight_t

cdef vertex_t MAX_VERTEX_ID = UINT_FAST32_MAX
cdef weight_t MAX_WEIGHT = UINT_FAST32_MAX

"""

from libc.stdint cimport uint64_t, uint_fast64_t, UINT_FAST64_MAX

ctypedef uint_fast64_t vertex_t
ctypedef uint_fast64_t weight_t

cdef vertex_t MAX_VERTEX_ID = UINT_FAST64_MAX
cdef weight_t MAX_WEIGHT = UINT_FAST64_MAX
"""
