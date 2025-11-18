#!/bin/bash -l
#SBATCH --job-name=bench_internode
#SBATCH --output=Fritz_ICX_DMVM_internode_non_blocking
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=72
#SBATCH --time=00:30:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0
export I_MPI_PIN_PROCESSOR_LIST=0-71

FILENAME="result_one_node_sor_domain.csv"


make

rm $FILENAME
touch $FILENAME

echo "scaling of SOR of different domain sizes in 1 node across 1-72 cores part 1 a)">>$FILENAME
for i in {1..72}
do
    echo "Running SOR on $i cores" >> $FILENAME
    mpirun -n $i ./exe-ICC poisson.par | tail -2 >> $FILENAME
done