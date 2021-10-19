#!/bin/bash
#PBS -N testio
#PBS -j oe
#PBS -q single
#PBS -l nodes=1:ppn=10
#PBS -m n

#################################################
# batch script for more.cfca.nao.ac.jp
#################################################

module load ompi4

# Go to this job's working directory
cd  $PBS_O_WORKDIR
export OMP_NUM_THREADS=1

date  >& log.$PBS_JOBID
time mpiexec -n 8 ./iotest  >> log.$PBS_JOBID
date  >> log.$PBS_JOBID
