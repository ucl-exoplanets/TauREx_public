#! /usr/bin/python
#
# script compiling all necessary c/cpp files for main code
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, May 2013  
########################################################

import os
# import glob

HOMEPATH = os.getcwd()
LIBPATH = '/library/'
MULTNESTPATH = '/PyMultiNest/'

#compiling pathintegral.cpp
print 'compiling pathintegral.cpp'
os.chdir(HOMEPATH+LIBPATH)
if os.path.isfile('pathintegral.so'):
    os.system('rm pathintegral.so')

os.system('g++ -fPIC -shared -o pathintegral.so pathintegral.cpp')


print '--------'
print 'all done'