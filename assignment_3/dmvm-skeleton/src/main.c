/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include<mpi.h>
#include <string.h>

#include "allocate.h"
#include "timing.h"

extern double dmvm(double* restrict y,
    const double* restrict a,
    double* restrict x,
    int N,
    int iter,int rank,int size);

double dmvm_non_blocking(double* restrict y,
    const double* restrict a,
    double* restrict x,
    int N,
    int iter,
    int rank,
    int size);

int main(int argc, char** argv)
{
    MPI_Init(&argc,&argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    size_t bytesPerWord = sizeof(double);
    size_t N            = 0;
    size_t iter         = 1;
    double *a, *x, *y;
    double t0, t1;
    double walltime;

    if (argc > 2) {
        N    = atoi(argv[1]);
        iter = atoi(argv[2]);
    } else {
        printf("Usage: %s <N> <iter>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
int num = N / size;
    int rest = N % size;
    int Nlocal = num + ((rank < rest) ? 1 : 0);
    int startRow = rank * num + ((rank < rest) ? rank : rest);

    // Allocate ONLY local rows for matrix a
    a = (double*)allocate(ARRAY_ALIGNMENT, Nlocal * N * bytesPerWord);
    // Full x vector on all ranks (needed for ring shift)
    x = (double*)allocate(ARRAY_ALIGNMENT, N * bytesPerWord);
    // Local y vector
    y = (double*)allocate(ARRAY_ALIGNMENT, Nlocal * bytesPerWord);

    for (int i = 0; i < N; i++) {
        x[i] = (double)i;
    }

    for (int i = 0; i < Nlocal; i++) {
        y[i] = 0.0;
    }

    for (int i = 0; i < Nlocal; i++) {
        for (int j = 0; j < N; j++) {
            // Global row index is (startRow + i)
            a[i * N + j] = (double)j + (startRow + i);
        }
    }

    walltime = dmvm(y, a, x, N, iter,rank,size);

    double flops = (double)2.0 * N * N * iter;
    // # iterations, problem size, flop rate, walltime
    if (rank==0) printf("%zu %zu %.2f %.2f\n", iter, N, 1.0E-06 * flops / walltime, walltime);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
