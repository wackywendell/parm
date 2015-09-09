#! /usr/bin/env python

# System imports
from setuptools import Extension, setup

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

module_opts = [
    ("_sim2d", "pyparm/sim_wrap2d.cxx", ["-DVEC2D"]),
    ("_sim3d", "pyparm/sim_wrap3d.cxx", ["-DVEC3D"]),
    ("_sim2dlong", "pyparm/sim_wrap2dlong.cxx", ["-DVEC2D", "-DLONGFLOAT"]),
    ("_sim3dlong", "pyparm/sim_wrap3dlong.cxx", ["-DVEC3D", "-DLONGFLOAT"])
]
hpp_files = ["src/collection.hpp", "src/constraints.hpp",
    "src/interaction.hpp", "src/trackers.hpp", "src/box.hpp",
    "src/vecrand.hpp"]
cpp_files = ["src/collection.cpp", "src/constraints.cpp",
    "src/interaction.cpp", "src/trackers.cpp", "src/box.cpp",
    "src/vecrand.cpp"]

swigged_modules = [
    Extension(
        name,
        [swig_file],
        include_dirs=[numpy_include, "src"],
        extra_compile_args=compile_opts,
    ) for name, swig_file, compile_opts in module_opts
]

setup(
    name="pyparm",
    description="None",
    author="Wendell Smith",
    version="0.2",
    ext_modules=swigged_modules
)
