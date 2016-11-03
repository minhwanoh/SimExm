from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy

extensions = []

#Compile sim extensions
sim_modules = ["expansion","labeling", "optics"]
for sim_module in sim_modules:
	extensions.append(Extension("src/simulation/models/" + sim_module, ["src/simulation/models/" + sim_module + ".pyx"],\
	 include_dirs = [numpy.get_include(), "."], extra_compile_args=["-O3", "-w"]))#"-fopenmp"]))

setup(ext_modules = cythonize(extensions, compiler_directives={'extra_compile_args': ["-03"] }, annotate=False),)

#Compile psf extension, currently unsused
setup(name = "_psf", ext_modules = [Extension("src/simulation/_psf",["src/simulation/psf.c"], include_dirs = [numpy.get_include()])])

#Compile tiffile extension
setup(name = "_tifffile", ext_modules = [Extension("src/database/models/_tifffile",["src/database/models/tifffile.c"], include_dirs = [numpy.get_include()])])
