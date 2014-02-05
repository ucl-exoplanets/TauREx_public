#! /usr/bin/python
#
# script compiling all necessary c/cpp files for main code
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, May 2013
#   - v1.1 : included script creation and first run full install option, IPW, Jan 2014
########################################################

import os
from sys import platform as _platform
import glob

HOMEPATH = os.getcwd()
LIBPATH = '/library/'
MULTNESTPATH = '/PyMultiNest/multinest/'

#compiling pathintegral.cpp
print 'compiling pathintegral.cpp'
os.chdir(HOMEPATH+LIBPATH)
if os.path.isfile('pathintegral.so'):
    os.system('rm pathintegral.so')

os.system('g++ -fPIC -shared -o pathintegral.so pathintegral.cpp')


########################################
#Multinest installation part

#generating source script for python path
print 'writing setup script. Before running code please source script:'
print '>> source setup.sh'

SCRIPT = open('setup.sh','w')
SCRIPT.write('echo \"setting MultiNest path\" \n')
SCRIPT.write('export LD_LIBRARY_PATH='+HOMEPATH+MULTNESTPATH+'lib:$LD_LIBRARY_PATH \n')
SCRIPT.close()

print ''
print '--------'
NESTINSTALL = raw_input('Automatically install/recompile/fix MultiNest? Y/[N]: ')
if NESTINSTALL == '' or NESTINSTALL == 'N' or NESTINSTALL == 'n':
    print 'no MultiNest installation'
    print '--------'
    print 'all done'
    exit()

print 'trying to install pymultinest library'
pipstatus = os.system('sudo pip install pymultinest')
if pipstatus != 0:
    print 'ERROR: please make sure you have PIP installed'
    print 'on MAC OSX (using macports) type: sudo port install py27-pip'
    print 'on LINUX (using apt-get) type: sudo apt-get install python-pip'


print 'trying to compile the multinest libraries...'
print 'please make sure you have the correct C++, Fortran compilers '
print 'and CMAKE installed'
os.chdir(HOMEPATH+MULTNESTPATH+'build')

if _platform == 'darwin':
    multistatus1 = os.system('cmake -DCMAKE_{C,CXX}_FLAGS="-arch x86_64" -DCMAKE_Fortran_FLAGS="-m64" ..')
else:
    multistatus1 = os.system('cmake ..')
if multistatus1 == 0:
    multistatus2 = os.system('make')
if multistatus2 == 0:
    multistatus3 = os.system('sudo make install')
if multistatus3 == 0 and _platform == 'darwin':
    os.chdir(HOMEPATH+MULTNESTPATH+'lib')
    libraries = glob.glob('*dylib')
    for file in libraries:
        os.system('sudo ln -s '+file+' '+file[:-5]+'so')


print '--------'
print 'all done'