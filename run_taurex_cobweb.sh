#! /bin/bash
#PBS -S /bin/bash

#################################
## Script to run Tau-REx on Cobweb
## 
## DATADIR    : data directory (e.g. /share/data/taurex) 
## SCRATCHDIR : temporary scratch dir. gets deleted at the end
## RUNDIR     : same as SCRATCHDIR but can include subfolders
## OUTDIR     : final results folder (DATADIR=OUTDIR usually)
## 
## NP         : number of total cpus, i.e. nodes x ppn
##################################

## Name of job
#PBS -N taurex_run

## Queue to submit to: compute/test
#PBS -q compute

##run parameters
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=24
#PBS -l mem=100gb

##error output and path variables
#PBS -j oe
#PBS -V

##cd ${PBS_O_WORKDIR}
##cat $PBS_NODEFILE > nodes

##Number of total cores for mpirun command
NP=24 

##setting up scratch
DATADIR=/share/data/ingo/repos/taurex_cluster
SCRATCHDIR=/scratch/ingo/run/single 
RUNDIR=$SCRATCHDIR
OUTDIR=/share/data/ingo/taurex/output

##coping data
mkdir -p $SCRATCHDIR
cp -rf $DATADIR/* $SCRATCHDIR
cd $RUNDIR  

#creating hostfile
cat $PBS_NODEFILE > nodes

##run job
mpirun -np $NP -hostfile nodes /share/apps/anaconda/2.2.0/bin/python taurex.py -p Parfiles/taurex_emission_wasp76.par

##copy stuff from scratch to /share/data
cp -rf $SCRATCHDIR/* $OUTDIR
rm -rf $SCRATCHDIR

