#!/bin/bash -l
#SBATCH --job-name=assignment_2_1
#SBATCH --output=output_2_1.txt
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=00:10:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV


module load intel intelmpi
make
srun -n 4 ./exe-ICX 
make clean
