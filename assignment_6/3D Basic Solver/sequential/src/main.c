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

#include "allocate.h"
#include "parameter.h"
#include "progress.h"
#include "solver.h"
#include "timing.h"
#include "vtkWriter.h"

#define G(v, i, j, k) v[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]

static void createBulkArrays(Solver* s, double* pg, double* ug, double* vg, double* wg)
{
    int imax = s->grid.imax;
    int jmax = s->grid.jmax;
    int kmax = s->grid.kmax;
    int idx  = 0;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                pg[idx++] = G(s->p, i, j, k);
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                ug[idx++] = (G(s->u, i, j, k) + G(s->u, i - 1, j, k)) / 2.0;
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                vg[idx++] = (G(s->v, i, j, k) + G(s->v, i, j - 1, k)) / 2.0;
            }
        }
    }

    idx = 0;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                wg[idx++] = (G(s->w, i, j, k) + G(s->w, i, j, k - 1)) / 2.0;
            }
        }
    }
}

int main(int argc, char** argv)
{
    double timeStart, timeStop;
    Parameter p;
    Solver s;
    initParameter(&p);

    if (argc != 2) {
        printf("Usage: %s <configFile>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }

    readParameter(&p, argv[1]);
    printParameter(&p);
    initSolver(&s, &p);
#ifndef VERBOSE
    initProgress(s.te);
#endif

    double tau = s.tau;
    double te  = s.te;
    double t   = 0.0;

    timeStart = getTimeStamp();
    while (t <= te) {
        if (tau > 0.0) computeTimestep(&s);
        setBoundaryConditions(&s);
        setSpecialBoundaryCondition(&s);
        computeFG(&s);
        computeRHS(&s);
        solveRB(&s);
        adaptUV(&s);
        t += s.dt;

#ifdef VERBOSE
        printf("TIME %f , TIMESTEP %f\n", t, s.dt);
#else
        printProgress(t);
#endif
    }
    timeStop = getTimeStamp();
#ifndef VERBOSE
    stopProgress();
#endif
    printf("Solution took %.2fs\n", timeStop - timeStart);

    double *pg, *ug, *vg, *wg;

    size_t bytesize = (size_t)(s.grid.imax * s.grid.jmax * s.grid.kmax) * sizeof(double);

    pg = allocate(64, bytesize);
    ug = allocate(64, bytesize);
    vg = allocate(64, bytesize);
    wg = allocate(64, bytesize);

    createBulkArrays(&s, pg, ug, vg, wg);
    VtkOptions opts = { .grid = s.grid };
    vtkOpen(&opts, s.problem);
    vtkScalar(&opts, "pressure", pg);
    vtkVector(&opts, "velocity", (VtkVector) { ug, vg, wg });
    vtkClose(&opts);
    return EXIT_SUCCESS;
}
