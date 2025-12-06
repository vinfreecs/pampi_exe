#!/bin/bash

t=72

make clean
make

# Correct arithmetic expansion:
export I_MPI_PIN=1
export I_MPI_DEBUG=1
export I_MPI_PIN_PROCESSOR_LIST=0-$((t-1))

mpirun -n $t ./exe-ICX canal.par

gnuplot vector.plot
gnuplot surface.plot
