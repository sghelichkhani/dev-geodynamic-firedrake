#!/bin/bash
#PBS -N SIMUNAME 
#PBS -P xd2
#PBS -q normal 
#PBS -l walltime=48:00:00
#PBS -l mem=192GB
#PBS -l ncpus=16 
#PBS -l storage=scratch/xd2+gdata/xd2
#PBS -l wd
#PBS -m abe

#  Loading firedrake
module use /g/data/xd2/modulefiles/
module load firedrake

# Adding python path
export PYTHONPATH="/g/data/xd2/sg8812/local/lib/python3.8/site-packages"

# Where we are
echo ALL_PRMS 

# Run the code
mpiexec -np 16 python3 optimisation.py 

