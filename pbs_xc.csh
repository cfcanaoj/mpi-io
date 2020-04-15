#PBS -N testio
#PBS -j oe
#PBS -q large-cfca
#PBS -l nodes=1
#PBS -m b

# Go to this job's working directory
cd  $PBS_O_WORKDIR
setenv OMP_NUM_THREADS 1

date  >& log.$PBS_JOBID
time aprun -cc none -n 8 ./iotest  >> log.$PBS_JOBID
date  >> log.$PBS_JOBID
