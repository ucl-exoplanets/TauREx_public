#! /usr/bin/python

################################################a
# Small wrapper script executing line_by_line profiling.
#
# Requirements: line_profiling package
# URL: https://pythonhosted.org/line_profiler/
# pip: sudo pip install --pre line_profiler
#
# Usage: set @profile decorator in front of functions that need to be profiled
#        then just run this wrapper and it will create a .lprofile file in /Profile dir.
#
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013
#
################################################


import os, sys

PROFNAME = 'taurex_emission'
PROFDIR  = 'Profiling/'

print 'Running line-by-line runtime profiling:'
print '-----------------------------------------'

os.chdir('../')
# os.system('./'+PROFDIR+'kernprof.py -l -v taurex.py -p Parfiles/taurex_superearth.par')
os.system('./'+PROFDIR+'kernprof.py -l -v create_spectrum.py -p Output/newgrid/newgrid.par')

print '-----------------------------------------'
print 'Finished profiling'
print''
print 'coallating results into '+PROFNAME+'.lprofile'

# sys.stdout = open(PROFDIR+PROFNAME+'.lprofile','wb')
os.system('python -m line_profiler create_spectrum.py.lprof > '+PROFDIR+PROFNAME+'.lprofile')
# os.system('rm exonest.py.lprof')



################################################
# Additional script running memory line-by-line profiling
# this requires the memory profiler
# URL: https://pypi.python.org/pypi/memory_profiler
# pip: sudo pip install memory_profiler
#
# USAGE: same as runtime line-by-line profiler
#
# OUTPUTS: .mprofile file
#
################################################

# print ''
# print 'Running line-by-line memory profiling:'
# print '-----------------------------------------'

# os.system('python -m memory_profiler exonest.py > '+PROFDIR+PROFNAME+'.mprofile')