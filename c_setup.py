#! /usr/bin/python

#compilation script turning all classes and libraries into shared (.so) C# libraries. 

import glob
from distutils.core import setup
from Cython.Build import cythonize




classlist = glob.glob('classes/*.py')
librarylist = glob.glob('library/*.py')

print 'Compiling code...'

setup(
    ext_modules = cythonize(classlist+librarylist)
)
