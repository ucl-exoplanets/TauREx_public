#! /usr/bin/python

#compilation script turning all classes and libraries into shared (.so) C# libraries. 
#
# To compile type: ./c_setup.py build_ext --inplace

import glob,os,shutil,numpy
from distutils.core import setup
from Cython.Build import cythonize


# classlist = glob.glob('classes/*.py')
librarylist = glob.glob('library/*.pyx')

print 'Compiling code...'

setup(
    ext_modules = cythonize(librarylist),
    include_dirs=[numpy.get_include()]
)

