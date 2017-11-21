#!/bin/bash

# Compile supporting libraries
# compilers: g++, pgc++

# openACC flag: -fopenacc
# openMP flag: -fopenmp

COMPILER=pgc++
FLAG=-fopenacc

# COMPILER=icc
# OPENMP_FLAG=-openmp


# transmission, ktab
# $COMPILER -fPIC -shared $OPENMP_FLAG -o ctypes_pathintegral_transmission_ktab.so ctypes_pathintegral_transmission_ktab.cpp

# transmission, xsec
$COMPILER -fPIC -shared -o ctypes_pathintegral_transmission_xsec.so ctypes_pathintegral_transmission_xsec.cpp

# emission, xsec
# $COMPILER -fPIC -shared $OPENMP_FLAG -o ctypes_pathintegral_emission.so ctypes_pathintegral_emission.cpp

# emission, ktab
# $COMPILER -fPIC -shared $OPENMP_FLAG -o ctypes_pathintegral_emission_ktab.so ctypes_pathintegral_emission_ktab.cpp
