#!/bin/bash
#SBATCH --job-name=testio
#SBATCH --partition=M-large-cfca
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100G
#SBATCH -o slmlog%J.out

# Go to this job's working directory
cd ${SLURM_SUBMIT_DIR}

export OMP_NUM_THREADS=1

date >& log.$SLURM_JOB_ID
time srun ./iotest >> log.$SLURM_JOB_ID
date >> log.$SLURM_JOB_ID




