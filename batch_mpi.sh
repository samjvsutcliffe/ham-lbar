#!/bin/bash

# Request resources:
#SBATCH --time=12:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH -p shared
#SBATCH -n 4                    # number of MPI ranks
#SBATCH --ntasks-per-socket=1   # number of MPI ranks per CPU socket
#SBATCH --cpus-per-task=16   # number of MPI ranks per CPU socket
#SBATCH --mem-per-cpu=1G
#SBATCH -N 2-4                    # number of compute nodes. 

module load gcc
module load intelmpi/2021.10
module load aocl

echo "Running code"
rm output/*

sbcl --dynamic-space-size 16000  --disable-debugger --load "build_step.lisp" --quit

#cp ~/quicklisp/local-projects/cl-mpm-worker/mpi-worker ./

#echo $OMP_NUM_THREADS
mpirun ./mpi-worker --dynamic-space-size 16000

