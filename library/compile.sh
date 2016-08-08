#!/bin/bash

# Compile supporting libraries


COMPILER=g++
OPENMP_FLAG=-fopenmp

# COMPILER=icc
# OPENMP_FLAG=-openmp


# transmission, single core version (used in retrieval or create_spectrum)
$COMPILER -fPIC -shared -o ctypes_pathintegral_transmission_xsec.so ctypes_pathintegral_transmission_xsec.cpp

# transmission, openmp version (only used in create_spectrum)
$COMPILER -fPIC -shared $OPENMP_FLAG -o ctypes_pathintegral_transmission_parallel_xsec.so ctypes_pathintegral_transmission_xsec.cpp

# emission, single core version (used in retrieval or create_spectrum)
$COMPILER -fPIC -shared -o ctypes_pathintegral_emission.so ctypes_pathintegral_emission.cpp

# emission, openmp version (only used in create_spectrum)
$COMPILER -fPIC -shared $OPENMP_FLAG -o ctypes_pathintegral_emission_parallel.so ctypes_pathintegral_emission.cpp