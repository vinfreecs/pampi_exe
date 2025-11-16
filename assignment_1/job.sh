#!/bin/bash -l                 
#
#SBATCH --nodes=1
#SBATCH --time=00:05:00            
#SBATCH --partition=singlenode
#SBATCH --output=output.txt
#SBATCH --job-name=job_assign_1 
#SBATCH --export=NONE

module load intel
make
./exe-ICX "$@"
make clean


