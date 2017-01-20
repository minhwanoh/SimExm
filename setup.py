from distutils.core import setup, Extension
import numpy

setup(name='_psf', ext_modules=[Extension('src/_psf', ['src/psf.c'], include_dirs=[numpy.get_include()])])