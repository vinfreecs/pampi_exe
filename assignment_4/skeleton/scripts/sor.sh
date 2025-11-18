
module load intel intelmpi likwid

export I_MPI_PIN=1 I_MPI_DEBUG=1 I_MPI_PIN_PROCESSOR_LIST=0-71

make clean

make

echo "scaling of SOR in 1 node across 1-72 cores part 1 a)">>sor_scaling_one_node_w1.2
for i in {1..72}
do
    echo "Running SOR on $i cores">>sor_scaling_one_node
    mpirun -n $i ./exe-ICC poisson.par | tail -2 >> sor_scaling_one_node_2_w1.2
done