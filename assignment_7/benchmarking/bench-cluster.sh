#!/bin/bash -l
#SBATCH --job-name=bench_internode
#SBATCH --output=Fritz_ICX_DMVM_internode
#SBATCH --partition=multinode
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=72
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0
export I_MPI_PIN_PROCESSOR_LIST=0-71

FILENAME="result_bench_internode.csv"

cd /home/hpc/pavl/pavl166v/pampi_exe/assignment_7/benchmarking

rm $FILENAME
touch $FILENAME
echo "Ranks,Perfoemance" >>$FILENAME

_iterate() {
    for np in $(seq 1 4); do
        npn=$(($np * $NPM))

        result="$(mpirun -n $npn ./exe-ICX canal_cluster.par | awk '/Performance total:/ {sum+=$14; count++} END {if (count > 0) print sum/count}')"

        echo $result
        # result="$(echo $result | sed 's/MPI startup(): Warning: I_MPI_PMI_LIBRARY will be ignored since the hydra process manager was found//g')"

        echo $npn $result >>$FILENAME
    done
}

NPM=72
_iterate


sed -i 's/ /,/g' $FILENAME
