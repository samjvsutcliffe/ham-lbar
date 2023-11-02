#!/bin/bash

# Request resources:
#SBATCH -c 32     # 1 entire node
#SBATCH --time=12:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=64G      # 1 GB RAM
#SBATCH -p shared

module load gcc
module load aocl

echo "Running code"
rm output/*

sbcl --dynamic-space-size 64000  --disable-debugger --load "template.lisp" --quit
