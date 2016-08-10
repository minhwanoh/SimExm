from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

extensions = []

#database
modules = ["connectomics", "conversion"]
for module in modules:
	extensions.append(Extension("database/methods/" + module, ["database/methods/" + module + ".pyx"], include_dirs = [numpy.get_include(), "database/methods/", "./"], extra_compile_args=["-O3", "-w"]))

#simulation
extensions.append(Extension("simulation/methods", ["simulation/methods.pyx"], include_dirs = [numpy.get_include(), "simulation/", "./"], extra_compile_args=["-O3", "-w"]))

sim_modules = ["labeling", "optics"]
for sim_module in sim_modules:
	extensions.append(Extension("simulation/models/" + sim_module, ["simulation/models/" + sim_module + ".pyx"], include_dirs = [numpy.get_include(), "simulation/models/", "./"], extra_compile_args=["-O3", "-w"]))#"-fopenmp"]))

# Visualization
visualization_modules = ["methods"]
for visualization_module in visualization_modules:
	extensions.append(Extension("visualization/" + visualization_module, ["visualization/" + visualization_module + ".pyx"], include_dirs = [numpy.get_include(), "visualization/", "./"], extra_compile_args=["-O3", "-w"]))

setup(
    ext_modules = cythonize(extensions, compiler_directives={'extra_compile_args': ["-03"] }, annotate=False),
)
