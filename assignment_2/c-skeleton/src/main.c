/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

#include "affinity.h"
#include "allocate.h"
#include "timing.h"

int main(int argc, char** argv) {
    int rank,size; 
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Hello World from %d of %d \n", rank, size);
    MPI_Finalize();
    return EXIT_SUCCESS; 
}
