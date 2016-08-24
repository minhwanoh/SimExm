from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

extensions = []

sim_modules = ["labeling", "optics"]
for sim_module in sim_modules:
	extensions.append(Extension("src/simulation/models/" + sim_module, ["src/simulation/models/" + sim_module + ".pyx"],\
	 include_dirs = [numpy.get_include(), "."], extra_compile_args=["-O3", "-w"]))#"-fopenmp"]))

setup(
    ext_modules = cythonize(extensions, compiler_directives={'extra_compile_args': ["-03"] }, annotate=False),
)