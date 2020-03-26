"""
Setup
"""

from distutils.core import setup
from Cython.Build import cythonize
import sys

directives = {
    'optimize.inline_defnode_calls': True,
    'language_level': sys.version_info[0] # "2" or "3"
}

extensions = [
    "sagemath/types.pyx",
    "sagemath/pairing_heap.pyx",
    "sagemath/digraph.pyx",
    "sagemath/dijkstra.pyx",
    "sagemath/kssp.pyx",
    "my_digraph.pyx"
              ]

setup(ext_modules=cythonize(extensions,
                            language_level=sys.version_info[0],  # 2 or 3
                            annotate=True,
                            compiler_directives=directives)
          )
