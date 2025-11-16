#!/bin/bash -l
#SBATCH --job-name=bench_internode
#SBATCH --output=Fritz_ICX_DMVM_internode_non_blocking
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

FILENAME="result_bench_internode_non_blocking.csv"


make distclean
make

rm $FILENAME
touch $FILENAME
echo "Ranks,NITER,N,MFlops,Time" >>$FILENAME

_iterate() {
    for np in $(seq 1 4); do
        npn=$(($np * $NPM))

        result="$(mpirun -n $npn ./exe-ICX $N $NITER)"
        result="$(echo $result | sed 's/MPI startup(): Warning: I_MPI_PMI_LIBRARY will be ignored since the hydra process manager was found//g')"

        echo $npn $result >>$FILENAME
    done
}

NPM=72

# For domain of 1000x1000
NITER=10000
N=1000
_iterate

# For domain of 4000x4000
NITER=10000
N=4000
_iterate

# For domain of 4000x4000
NITER=10000
N=10000
_iterate

# For domain of 4000x4000
NITER=5000
N=20000
_iterate

sed -i 's/ /,/g' $FILENAME
