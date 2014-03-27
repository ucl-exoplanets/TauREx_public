#! /usr/bin/python -W ignore

###########################################################
# Script generating transmission spectra from absorption coeff
# files for spectral libarary needed at preselector stage.
#
# Modification History:
#   - v1.0 : first definition, Ingo Waldmann, Apr 2013
#
###########################################################

#loading libraries
import numpy, pylab, optparse, glob,string,os
from numpy import * #nummerical array library
from pylab import * #science and plotting library for python
from ConfigParser import SafeConfigParser



#loading classes
from classes.parameters import *

from classes.transmission import *

from classes.profile import *
from classes.data import *


#loading libraries
from library.library_transmission import *
from library.library_general import *



parser = optparse.OptionParser()
parser.add_option('-p', '--parfile',
                  dest="param_filename",
                  default="exonest.par",
)
parser.add_option('-v', '--verbose',
                  dest="verbose",
                  default=False,
                  action="store_true",
)
options, remainder = parser.parse_args()

#Initialise parameters object
params = parameters(options.param_filename)
if params.verbose == True:
    print 'ARGV      :', sys.argv[1:]
    print 'VERBOSE   :', params.verbose
    print 'PARFILE    :', options.param_filename
    print 'REMAINING :', remainder

#####################################################################


#initialising data object
dataob = data(params)

#adding some molecules to the atmosphere
dataob.add_molecule('H2', 2.0, 2.0e-9, 1.0001384, 0.85)

#initialising TP profile object
profileob = profile(params, dataob)

#initialising transmission radiative transfer code object
transob = transmission(params, dataob)


########
###example of how to manually reading in ABS file and computing transmission spectrum
# dataob.set_ABSfile(path='/Users/ingowaldmann/UCLlocal/REPOS/exonest/exonestpy/test-code/crosssections/',
#                    filelist=['1H2-16O_0-29999_600K_0.010000.sigma.abs'])


#reading available cross section lists in PATH
PATH = '/Users/ingowaldmann/UCLlocal/REPOS/exonest/exonestpy/test-code/crosssections/'
OUTPATH = '/Users/ingowaldmann/UCLlocal/REPOS/exonest/exonestpy/test-code/speclib/'
globlist = glob.glob(PATH+'*.abs')

if os.path.isdir(OUTPATH) == False:
    os.mkdir(OUTPATH)


for FILE in globlist:
    fname = string.rsplit(FILE,'/',1)[1] #splitting the name
    temp  = float(string.rsplit(fname,'_',2)[1][:-1]) #getting temperature from file name

    print fname

    for exponent in range(3,7):
        print 1.0/10**exponent

        X_in   = zeros((profileob.nlayers,profileob.ngas)) #setting up mixing ratio array
        X_in  += 1.0/10**exponent #setting mixing ratio
        rho_in = profileob.get_rho(T=temp) #calculating T-P profile

        dataob.set_ABSfile(path=PATH,filelist=[fname]) #reading in cross section file
        transob.reset(dataob) #resets transob to reflect changes in dataob

        #manually setting mixing ratio and T-P profile
        MODEL = transob.cpath_integral(rho=rho_in,X=X_in) #computing transmission
        # print OUTPATH+fname[:-4]+'_'+str(1.0/10**exponent)+'d.spec'
        savetxt(OUTPATH+fname[:-4]+'_'+str(1.0/10**exponent)+'d.spec',column_stack((dataob.wavegrid,MODEL)))

        # print shape(dataob.wavegrid), shape(MODEL), shape(column_stack([dataob.wavegrid,MODEL]))

        # ion()
        # clf()
        # figure(1)
        # plot(dataob.wavegrid,MODEL)
        # draw()
        # exit()