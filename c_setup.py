#! /usr/bin/python

#compilation script turning all classes and libraries into shared (.so) C# libraries. 

import glob,os
from distutils.core import setup
from Cython.Build import cythonize




classlist = glob.glob('classes/*.py')
librarylist = glob.glob('library/*.py')

print classlist
exit()
print 'Compiling code...'

setup(
    ext_modules = cythonize(classlist+librarylist)
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
    c_librarylist = glob.glob('classes/*.c')
    os.remove(c_classlist+c_librarylist)
    
if delpy:
    print 'deleting .py files'
    os.remove(classlist_librarylist)
    
print 'done'
    
     