===================================================================
TAU v1.8 README
===================================================================

Directories
-------------------------------------------
	./doc/		directory containing program documentation
	./out/		directory containing output files
	./run/		directory containing required and optional input files
	./src/		directory containing source code


./doc/:
	User_Guide_v1.8c.pdf			User Guide for version 1.8c
	readme.txt				Readme file

./out/:
	sample_tau_output.dat		absorption as a function of wavelength for run mode 9

./run/:
	h2o_800K.abs			absorption cross-section as a function of wavelength for water (at T=800K from the BT2 line list)
	h2o_1000K.abs			absorption cross-section as a function of wavelength for water (at T=1000K from the BT2 line list)
	h2o_1500K.abs			absorption cross-section as a function of wavelength for water (at T=1500K from the BT2 line list)
	profile.atm			sample atmospheric pressure-temperature profile for an isothermal atmosphere (T=1500K) of HD189733b

./src/:
	functions.h			header file containing external functions and definitions
	tau.cpp				program source code file


Installation
-------------------------------------------
Step 1: Download the tar ball 'tau.tar.gz' from "http://www.ucl.ac.uk/exoplanets/".

Step 2: In the directory where the file is downloaded, type into the terminal: 
	tar -xvf tau.tar.gz

	This will extract all of the source files into a directory named 'tau/'.

Step 3: To compile and make the executable, whilst in the 'tau' directory type: 
	g++ -lgomp -fopenmp -o tau ./src/tau.cpp

	This creates the executable file 'tau' in the current directory. The header file containing the required functions and scientific constants is located in the directory './tau/src/', output files are placed in the directory './tau/out/', and the required input files should be in the directory './tau/run/'. Documentation for the code is in './tau/doc/'.

Note that for the compilation, it is assumed that g++ from GCC, the Gnu Compiler Collection is present on the system. If not, it may be downloaded from "http://gcc.gnu.org". Additionally, the library 'omp.h' from the OpenMP API specification for parallel programming must be installed, obtainable from "http://openmp.org/".


Execution
-------------------------------------------
The current version of the code has the option to run in 3 modes, as detailed below. Run Mode 0 is used to calculate the spectrum for a single trace molecular absorber, Run Mode 1 is used to model multiple molecular absorbers, and Run Mode 9 is used for a quick run with custom settings. Run Modes 0 and 1 allow the '.atm' and/or '.abs' file inputs to be specified as command line arguments, and Run Mode 9 provides a quick testing mode, allowing all interactive elements to be bypassed and all input files and mixing ratio values hard-coded prior to compilation. The usage instructions and mode formats can be displayed by typing './tau' at the command line, and the option sorting and default file definitions can be found in the function 'optionSort_H' in the header file.

Run Mode 0: The run command for this mode is './tau 0 [atmfile] [[absfile]]'
This mode allows the user to specify the atmospheric profile file and a single absorption cross-section file (i.e. to model a single absorber molecule in the atmosphere). If run in this mode, both files are optional command line arguments, and if nothing is entered for either the code reverts to the default files, './run/profile.atm' and './run/h2o_1500K.abs'. It is possible to only enter a filename for the atmospheric profile, and leave the code to take the default absorption cross-section file, but if the '.abs' file is to be entered, so the '.atm' file must also be. Mixing ratios are read in from the atmospheric profile file.

Run Mode 1: The run command for this mode is './tau 1 [atmfile]'
This mode allows the user to specify on the command line the atmospheric profile file and (interactively) multiple absorption cross-section files (i.e. to model an arbitrary number of absorber molecules in the atmosphere). If run in this mode, the file is an optional command line argument, and if nothing is entered the code reverts to the default file, './run/profile.atm'. Mixing ratios are read in from the atmospheric profile file (one column per absorption cross-section file).

Run Mode 9: The run command for this mode is './tau 9'
This mode allows the user to bypass all interactive input, and to specify all input files and mixing ratios in the program code pre-compilation. This mode is useful for example when comparing models where the atmospheric constituents and most of the mixing ratios remain constant between runs but one or a few are altered. The default atmospheric profile and absorption coefficient files are set in the file 'functions.h', and mixing ratios are set to 10^{-5} in the file 'tau.cpp'. 
