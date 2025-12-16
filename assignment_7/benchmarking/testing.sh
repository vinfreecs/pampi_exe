#!/bin/bash -l
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

        echo $npn
        result="$(mpirun -n $npn ./exe-ICX canal_node.par | awk '/Performance total:/ {sum+=$14; count++} END {if (count > 0) print sum/count}')"

        echo $result
        
        echo $npn $result >>$FILENAME
    done
}

_iterate
