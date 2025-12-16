#!/bin/bash -l
#SBATCH --job-name=bench_intranode
#SBATCH --output=Fritz_ICX_DMVM_intranode
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0

FILENAME="result_bench_intranode.csv"

cd /home/hpc/pavl/pavl166v/pampi_exe/assignment_7/benchmarking

rm $FILENAME
touch $FILENAME
echo "Ranks,Performance MLUPs/s" >>$FILENAME

NPM=18
_iterate() {
    for np in $(seq 1 4); do
        npn=$(($np * $NPM))
        np_1=$(($npn - 1))
        export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

        result="$(mpirun -n $npn ./exe-ICX canal_node.par | awk '/Performance total:/ {sum+=$14; count++} END {if (count > 0) print sum/count}')"
        
        echo $npn $result >>$FILENAME
    done
}

_iterate

