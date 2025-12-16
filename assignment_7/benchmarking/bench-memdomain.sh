#!/bin/bash -l
#SBATCH --job-name=bench_memdomain
#SBATCH --output=Fritz_ICX_DMVM_memdomain
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=05:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0

FILENAME="result_bench_memdomain.csv"

cd /home/hpc/pavl/pavl166v/pampi_exe/assignment_7/benchmarking
make distclean
make

rm $FILENAME
touch $FILENAME
echo "Ranks,Performance MLUPs/s" >>$FILENAME

_iterate() {
    for np in $(seq 1 18); do
        np_1=$(($np - 1))
        export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

        result="$(mpirun -n $np ./exe-ICX canal_memdomain.par | awk '/Performance total:/ {sum+=$14; count++} END {if (count > 0) print sum/count}')"

        echo $np $result >>$FILENAME
    done
}

NPM=18

_iterate

sed -i 's/ /,/g' $FILENAME
