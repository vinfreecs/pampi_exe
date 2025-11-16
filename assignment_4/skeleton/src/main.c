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

#include "likwid-marker.h"
#include "parameter.h"
#include "solver.h"
#include "timing.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc,&argv);
    double startTime, endTime;
    Parameter params;
    Solver solver;
    initParameter(&params);
    int rank;
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (argc < 2) {
        printf("Usage: %s <configFile>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    readParameter(&params, argv[1]);
    printParameter(&params);

    initSolver(&solver, &params, 2);
    startTime = getTimeStamp();
    solve(&solver);
    getResult(&solver);
    endTime = getTimeStamp();
    // writeResult(&solver, "p.dat");

    MPI_Reduce(MPI_IN_PLACE, &startTime,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &endTime,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
    if(rank == 0) printf("Walltime %.2fs\n", endTime - startTime);
    MPI_Finalize();
    return EXIT_SUCCESS;
}
