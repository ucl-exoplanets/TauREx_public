# TauREx version 2.5

TauREx (Tau Retrieval for Exoplanets) is a fully bayesian inverse atmospheric retrieval framework. 
TauREx is a very extensive retrieval framework with a wide range of functionalities. Here are installation instructions and worked examples whilst we write a more exhaustive manual.  
Already existing documentation can be found in the *Documentation* folder. 

For any questions on how to run the code, features and bugs, please email Ingo Waldmann (ingo@star.ucl.ac.uk).


## References:
Waldmann et al. (2015a), “Tau-REx I: A Next Generation Retrieval Code for Exoplanetary Atmospheres”, ApJ, 802, 107
Waldmann et al. (2015b), “Tau-REx II: Retrieval of Emission Spectra, ApJ, 813, 13

## Active developers:
- Ingo Waldmann 
- Marco Rocchetto 
- Tiziano Zingales

---

## License: 
This work is licensed under the Creative Commons Attribution 4.0 International License. 

We would like to draw your attention to Section 3 of the license and to:
- retain identification of the creators by including the above listed references in future work and publications.
- indicate if You modified the Licensed Material and retain an indication of any previous modifications

To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

---

## Installation:

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

Finally, the TauREx C++ and fortran libraries need to be compiled. In the /library folder 

```
sh comipile.sh
```

Optional: if the chemical equilibrium model is required, it must be compiled separately first. In the /library/ACE folder

```
gfortran -shared -fPIC  -o ACE.so Md_ACE.f90 Md_Constantes.f90 Md_Types_Numeriques.f90 Md_Utilitaires.f90 Md_numerical_recipes.f90
```

That’s it. Now you should be able to run TauREx. 

## Adding Input folder data

TauREx requires input data such as ktables/cross-sections/cia/mie etc. These cannot be provided on GitHub due to size. The required Input folder can be downloaded here: 
http://bit.ly/2y7XkKq

## Downloading from Dropbox

Alternatively, you can download the full code with its Input files in place from here: http://bit.ly/2wmPjMV

---

## Running TauREx

TauREx has two modes in which it can be run: Forward model and Retrieval. 

### Forward model
In the forward model mode, we can create individual spectra through either the provision of a 
parameter file and listed abundances or C/O and metallicity, or by providing external temperature-pressure and chemical profile files. 

For a standard example, we provide example parameter files in the /tests folder. 

#### For a transmisison/emission case:

```python
python create_spectrum.py -p tests/test_0_transmission/test_fm.par --plot
```
The spectrum is saved as ascii in SPECTRUM_out.dat in the Output folder specified in the parameter file. 
Note these parameter files are minimal files and modify the default parameter in Parfiles/default.par. You can add as many paramters in your minimal par file
as required.

Flags:
-p : path to parameter file 
--plot : creates plot of spectrum 
--save_instance : saves the spectrum, together with other data (e.g. temperature-pressure profiles, transmissivity, etc) in a numpy pickle file. 

To get a list of parameters that can be set on command line:

```python
python create_spectrum.py --help
```

For the emission example, use 

```python
python create_spectrum.py -p tests/test_0_emission/test_fm.par --plot
```

#### For C/O equilibrium chemistry models 

TauREx uses the Atmospheric Chemical Equilibrium (ACE) model (Agundez et al. 2012, A&A, 548, A73) using the thermochemical data by Venot et al (2012, A&A,546,A43).
To run the C/O forward model:

```python
python create_spectrum.py -p tests/test_0_transmission_ace/test_fm.par --plot
```

#### External input files 

External input files for temperature-pressure profiles and chemical abundance profiles can be provided using the following routine

```python
create_spectrum_chemprofile.py
```

We will add documentation for this mode soon. If you want to use it, please just email ingo@star.ucl.ac.uk for now.


### Retrieval
Retrievals can be run in transmission and emission modes (currently this version only supports transiting planets but directly imaged planets can also be modeled)

To run a retrieval on a spectrum, specify the spectrum path in the parameter file 

```
[Input]
spectrum_file = PATH_TO_SPECTRUM
```

The spectrum file should have 3 to 4 columns: wavelength (microns), (Rp/Rs)^2 or Fp/Fs, Error, wavelength bin width (microns, optional)
The parameter file should contain all necessary attributes to model the spectrum. We have provided examples for 30 planets in transmission in the Tsiaras et al. (2017, arXiv: ) paper. 
The data, parameter files and corresponding retrievals can be found here: http://bit.ly/HSTDATA

The /tests folder provides transmission and emission examples. To run the tranmission retrieval case: 

```python
python taurex.py -p tests/test_0_transmission/test_retrieval.par --plot
```

The --plot flag will plot a standard set of output results (posterior distirbutions, fitted spectrum, opacity contribution breakdown). To plot all available data (e.g. TP-profiles), plots can be 
generated after the retrieval has terminated using the created nest_out.pickle files. All parameters can be specified in the parameter file or given as flags in the command line. To print available paramters:

```python
python taurex.py --help
```

To generate all available plots: 

```python
python tools/taurex_plots.py --db_dir tests/tests_0_transmission/ --plot_all
```

where --db_dir is the directory containing the nest_out.pickle file. For the data structure of the nest_out.pickle file and the data contained within, we refer you to the documentation in the *Documentation* folder.

There are obviously many many topics that are not covered here, like the various temperature-pressure profile parameterisations, the various Mie and grey cloud models, Ktables and cross-sections, etc. 
More info will be available in the upcoming documentation. Until then, please just ask me (ingo@star.ucl.ac.uk). Your questions will actually prompt me to get on with writing the manual. 


