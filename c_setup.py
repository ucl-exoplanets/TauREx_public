#! /usr/bin/python

#compilation script turning all classes and libraries into shared (.so) C# libraries. 
#
# To compile type: ./c_setup.py build_ext --inplace

import glob,os,shutil,numpy
from distutils.core import setup
from Cython.Build import cythonize


classlist = glob.glob('classes/*.py')
librarylist = glob.glob('library/*.py')

print 'Compiling code...'

setup(
    ext_modules = cythonize(classlist+librarylist),
    include_dirs=[numpy.get_include()]
)


print ''
delc = raw_input('Delete .c files (y/n)? [y]: ')
if delc == '' or delc == 'y':
    delc = True
else:
    delc = False
    
delpy = raw_input('Delete .py files (y/n)? [n]: ')
if delpy == '' or delpy == 'n':
    delpy = False
elif delpy == 'y':
    delpy = True
    
if delc:
    print 'deleting .c files'
    c_classlist = glob.glob('classes/*.c')
    c_librarylist = glob.glob('library/*.c')
    for file in c_classlist:
        try:
            os.remove(file)
        except OSError:
            pass
    for file in c_librarylist:
        try:
            os.remove(file)
        except OSError:
            pass
    
if delpy:
    print 'deleting .py files'
    for file in classlist:
        try:
            os.remove(file)
        except OSError:
            pass
        
    for file in librarylist:
        try:
            os.remove(file)
        except OSError:
            pass
    
    pyc_classlist = glob.glob('classes/*.pyc')
    pyc_librarylist = glob.glob('library/*.pyc')

    print 'deleting .pyc files'
    for file in pyc_classlist:
        try:
            os.remove(file)
        except OSError:
            pass
    for file in pyc_librarylist:
        try:
            os.remove(file)
        except OSError:
            pass   
        
print 'removing build directory'
try:
    shutil.rmtree('./build')
except OSError:
    pass
    
print 'done'
    
     