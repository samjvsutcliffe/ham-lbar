#!/bin/bash

# Request resources:
#SBATCH -c 32     # 1 entire node
#SBATCH --time=10:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=32G      # 1 GB RAM
#SBATCH -p shared

module load gcc
module load intelmpi
module load aocl

echo "Running code"
rm output/*

sbcl --dynamic-space-size 32000  --disable-debugger --load "lbar.lisp" --quit
