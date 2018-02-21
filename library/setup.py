
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np


ext_modules = [Extension("library_cythonised_functions", ["library_cythonised_functions.pyx"])]

setup(
    name = 'library_cythonised_functions',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
    include_dirs=[np.get_include()]
)
