# TauREx version 2.5

TauREx (Tau Retrieval for Exoplanets) is a fully bayesian inverse atmospheric retrieval framework. A full documentation of the code is still to come but here are installation instructions and worked examples for now. 

For any questions on how to run the code, features and bugs, please email Ingo Waldmann (ingo@star.ucl.ac.uk).


### References:
Waldmann et al. (2015a), “Tau-REx I: A Next Generation Retrieval Code for Exoplanetary Atmospheres”, ApJ, 802, 107
Waldmann et al. (2015b), “Tau-REx II: Retrieval of Emission Spectra, ApJ, 813, 13

### License: 
This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

### Installation:

These are some preliminary sets of instructions to install TauREx on your system. They are very Mac centric but most of it should apply to a linux installation as well. To install things easily and smoothly for Mac, I strongly recommend to use macports (https://www.macports.org). For ubuntu, most of the work will be the same using apt-get but the individual commands will obviously be different. 

The main issue here is the installation of multinest. You only need to do this once for your system but it can be a bit of a hassle first time you try. There is quite a lot of documentation on the web though if you google it. 

First, download multinest (you will need to register to access the download page): 

https://ccpforge.cse.rl.ac.uk/gf/project/multinest/

Make sure you download the cmake installation version. This sorts most of the makefile and unpack it somewhere. The default installation directory (unless you change it) is
/usr/local/multinest
which is the default installation path for the further description here. You can also choose a more local directory. Before installing, read the README file in the multinest folder. This should explain quite well how to do that.  

Before compiling multinest, I found it helpful to make sure the latest version of gcc is installed on your system (do not install 4.8 as that has a bug w.r.t. some lapack libraries). Anything above version 4.9 is fine.  gcc 5 or above is recommended and everything will be linked to this from here on 

To install gcc on a Mac using macports, type

```
sudo port install gcc5
```

now make sure this gcc is the default option

```
sudo port select --set gcc mp-gcc5
```

install cmake if you haven’t already got it 

```
sudo port install cmake
```

install the OpenBLAS library (this can take an hour, not 100% sure you really need that. It depends whether you have a set of lapack libraries already installed on your system. You may want to skip this step on first try and then return to it later if needed. )

```
sudo port install openblas +gcc5
```

Install openmpi with gcc5 compiler linking

```
sudo port install openmpi +gcc5
sudo port select --set mpi openmpi-mp-fortran
```

now install needed python libraries using macports: The minimal set of standard libraries you need are: numpy, matplotlib, pymc, configparser, optparse 

```
sudo port install py27-pymc, py27-sklearn, py27-cython, py27-configparser, py27-optparse, py27-threading, py27-subprocess, pip
```

now install the mpi linking to python with the correct links . This is important: make always sure that you are installing everything with openmpi and not mpich linking

```
sudo port install py27-mpi4py +openmpi
```

now install pymultinest using pip 
```
sudo pip install pymultinest 
```

Now we get to compiling multinest itself. Go to the folder where you unpacked the multinest files into. There is a very helpful little readme file but just in case here is what you have to do: 

go to the multinest directory and make a folder ‘build’
```
$ cd $MULTNEST_DIR
$ mkdir build 
```
now go into the build folder and make and install multinest 

```
$ cd build
$ sudo cmake ..
$ sudo make
$ sudo make install 
```

Make sure that it correctly finds the gcc, gfortran (or icc, ifort) compilers, BLAS and LAPACK libraries and OpenMPI libraries. Cmake will complain if it doesn’t. 

If you are using a Mac, you need to soft-link the .dylib libraries to .so version. You do not need to do that on a linux system. Go to the library folder and /lib and set up soft links from .dylib to .so libraries 

these files need to have softlinks to .so endings
```
libmultinest.dylib
libmultinest.3.9.dylib
libmultinest_mpi.dylib
libmultinest_mpi.3.9.dylib
```
```
sudo ln -s libmultinest.dylib libmultinest.so
sudo ln -s libmultinest_mpi.dylib libmultines_mpi.so
```

Finally, add the multinest libraries to your LD_LIBRARY_PATH variable. Best is if you do that in your .bashrc or .bash_profile or .profile (whichever your system uses)

```
export LD_LIBRARY_PATH=/usr/local/multinest/lib:$LD_LIBRARY_PATH
```

That’s it. Now you should be able to run TauREx. 


