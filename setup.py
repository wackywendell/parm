#! /usr/bin/env python

# System imports
from distutils.core import Extension, setup

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# ezrange extension module
_d2 = Extension("_sim2d",
    ["src/sim.i", "src/array.i", "src/collection.hpp", "src/constraints.hpp",
    "src/interaction.hpp", "src/trackers.hpp", "src/box.hpp",
    "src/vecrand.hpp", "src/collection.cpp", "src/constraints.cpp",
    "src/interaction.cpp", "src/trackers.cpp", "src/box.cpp",
    "src/vecrand.cpp"],
    include_dirs=[numpy_include],
    swig_opts=["-py3", "-shadow", "-c++", "-DVEC2D"]
)

# ezrange setup
setup(
    name="pyparm",
    description="None",
    author="Wendell Smith",
    version="0.2",
    ext_modules=[_d2]
)
