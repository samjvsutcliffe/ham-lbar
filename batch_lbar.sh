#!/bin/bash

# Request resources:
#SBATCH --time=8:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH -p shared
#SBATCH -n 4                    # number of MPI ranks
#SBATCH --ntasks-per-socket=1   # number of MPI ranks per CPU socket
#SBATCH --cpus-per-task=16   # number of MPI ranks per CPU socket
#SBATCH --mem-per-cpu=1G
#SBATCH -N 1-4                    # number of compute nodes. 

module load gcc
#module load intelmpi
module load mvapich2
module load aocl

echo "Running code"
#rm -r output-*


#sbcl --dynamic-space-size 16000  --disable-debugger --load "build_step.lisp" --quit

#cp ~/quicklisp/local-projects/cl-mpm-worker/mpi-worker ./

export MV2_ENABLE_AFFINITY=0
echo $OMP_NUM_THREADS
#export REFINE=6.0
#export KAPPA=1.0
#export lc=1.0
#2.65
#rm -r output-*
#export lc=0.0
#export REFINE=0.5
mpirun ./mpi-worker --dynamic-space-size 16000
