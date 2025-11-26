/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "parameter.h"
#include "solver.h"
#include "util.h"

#define P(i, j)   p[(j) * (imaxLocal + 2) + (i)]
#define F(i, j)   f[(j) * (imaxLocal + 2) + (i)]
#define G(i, j)   g[(j) * (imaxLocal + 2) + (i)]
#define U(i, j)   u[(j) * (imaxLocal + 2) + (i)]
#define V(i, j)   v[(j) * (imaxLocal + 2) + (i)]
#define RHS(i, j) rhs[(j) * (imaxLocal + 2) + (i)]
#define Pall(i, j)   Pall[(j) * (imax + 2) + (i)]
#define Uall(i, j)   Uall[(j) * (imax + 2) + (i)]
#define Vall(i, j)   Vall[(j) * (imax + 2) + (i)]

static int sizeOfRank(int rank, int size, int N)
{
    return N / size + ((N % size > rank) ? 1 : 0);
}

/* Use this function to check if your exchange function works or not. */
/* Assigns every ranks local pressure array with rank number.         */
/* After exchange, you should see your neighbours rank in ghost cells */
/* of the local pressure array.                                       */
static void printExchange(Solver* solver, double* grid)
{
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double* p = grid;
    for (int j = 0; j < jmaxLocal + 2; j++) {
        for (int i = 0; i < imaxLocal + 2; i++) {
            P(i, j) = solver->rank;
        }
    }

    exchange(solver, p);

    for (int i = 0; i < solver->size; i++) {
        if (i == solver->rank) {
            printf("### RANK %d "
                   "#######################################################\n",
                solver->rank);
            for (int j = 0; j < jmaxLocal + 2; j++) {
                printf("%02d: ", j);
                for (int i = 0; i < imaxLocal + 2; i++) {
                    printf("%12.2f  ", grid[j * (imaxLocal + 2) + i]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(solver->comm);
    }
}

/* Use this function to check if your shift function works or not. */
/* Assigns every ranks local pressure array with rank number.      */
/* After shift, you should see your neighbours rank in ghost cells */
/* of the local F and G array.                                     */
static void printShift(Solver* solver)
{
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double* f = solver->f;
    double* g = solver->g;


    for (int j = 0; j < jmaxLocal + 2; j++) {
        for (int i = 0; i < imaxLocal + 2; i++) {
            F(i, j) = solver->rank;
            G(i, j) = solver->rank;
        }
    }

    shift(solver);

    printf("Vector F\n");
    for (int i = 0; i < solver->size; i++) {
        if (i == solver->rank) {
            printf("### RANK %d "
                   "#######################################################\n",
                solver->rank);
            for (int j = 0; j < jmaxLocal + 2; j++) {
                printf("%02d: ", j);
                for (int i = 0; i < imaxLocal + 2; i++) {
                    printf("%12.2f  ", F(i, j));
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    printf("Vector G\n");
    for (int i = 0; i < solver->size; i++) {
        if (i == solver->rank) {
            printf("### RANK %d "
                   "#######################################################\n",
                solver->rank);
            for (int j = 0; j < jmaxLocal + 2; j++) {
                printf("%02d: ", j);
                for (int i = 0; i < imaxLocal + 2; i++) {
                    printf("%12.2f  ", G(i, j));
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

/* subroutines local to this module */
static int sum(int* sizes, int init, int offset, int coord)
{
    int sum = 0;

    for (int i = init - offset; coord > 0; i -= offset, --coord) {
        sum += sizes[i];
    }

    return sum;
}

void exchange(Solver* solver, double* grid)
{
    /* TODO ghost cell exchange with respective neighbors   */
    /* Remember this is a 2D domain decomposition and you   */
    /* have your displacements in solver->sdispls and       */
    /* solver->rdispls. Use solver->bufferTypes to indicate */
    /* datatype to send. You have initialized them during   */
    /* initSolver() call and you can use them directly here.*/
    /* Refer to lecture slides for more information.        */
    int counts[NDIRS] = { 1, 1, 1, 1 };
    MPI_Neighbor_alltoallw(...);
}

void shift(Solver* solver)
{
    /* TODO shift (send) f ghost layer from left to right */
    /* TODO shift (send) g ghost layer from bottom to top */
    double* f = solver->f;
    double* g = solver->g;

    MPI_Request requests[4] = { MPI_REQUEST_NULL,
        MPI_REQUEST_NULL,
        MPI_REQUEST_NULL,
        MPI_REQUEST_NULL };

    /* shift G                                              */
    /* receive ghost cells from bottom neighbor             */
    /* use solver->bufferTypes and solver->neighbours array */
    /* to receive from BOTTOM rank                          */
    double* buf = g + 1;
    MPI_Irecv(...);

    /* send ghost cells to top neighbor                     */
    /* use solver->bufferTypes and solver->neighbours array */
    /* to send to TOP rank                                  */
    buf = g + (solver->jmaxLocal) * (solver->imaxLocal + 2) + 1;
    MPI_Isend(...);

    /* shift F                                              */
    /* receive ghost cells from left neighbor               */
    /* use solver->bufferTypes and solver->neighbours array */
    /* to receive from LEFT rank                            */
    buf = f + (solver->imaxLocal + 2);
    MPI_Irecv(...);

    /* send ghost cells to right neighbor                   */
    /* use solver->bufferTypes and solver->neighbours array */
    /* to send to RIGHT rank                                */
    buf = f + (solver->imaxLocal + 2) + (solver->imaxLocal);
    MPI_Isend(...);

    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
}

static int isBoundary(Solver* solver, int direction)
{
    /* Use solver->coords to check if the rank has boundaries */
    /* Fill in the rest of the boundaries. One example given  */
    switch (direction) {
    case LEFT:
        return solver->coords[IDIM] == 0;
        break;
    case RIGHT:
        return solver->coords[IDIM] == ...;
        break;
    case BOTTOM:
        return solver->coords[JDIM] == ...;
        break;
    case TOP:
        return solver->coords[JDIM] == ...;
        break;
    }
    return 1;
}

static void assembleResult(Solver* solver, double* src, double* dst)
{
    /* do not adapt these to local imaxLocal and jmaxLocal */
    int imax = solver->imax;
    int jmax = solver->jmax;

    MPI_Request* requests;
    int numRequests = 1;

    if (solver->rank == 0) {
        numRequests = solver->size + 1;
    } else {
        numRequests = 1;
    }

    requests = (MPI_Request*)malloc(numRequests * sizeof(MPI_Request));

    /* all ranks send their bulk array, including the external boundary layer */
    MPI_Datatype bulkType;
    int oldSizes[NDIMS] = { solver->jmaxLocal + 2, solver->imaxLocal + 2 };
    int newSizes[NDIMS] = { solver->jmaxLocal, solver->imaxLocal };
    int starts[NDIMS]   = { 1, 1 };

    if (isBoundary(solver, LEFT)) {
        newSizes[CIDIM] += 1;
        starts[CIDIM] = 0;
    }
    if (isBoundary(solver, RIGHT)) {
        newSizes[CIDIM] += 1;
    }
    if (isBoundary(solver, BOTTOM)) {
        newSizes[CJDIM] += 1;
        starts[CJDIM] = 0;
    }
    if (isBoundary(solver, TOP)) {
        newSizes[CJDIM] += 1;
    }

    /* Refer to lecture slides for more information                       */
    /* https://rookiehpc.org/mpi/docs/mpi_type_create_subarray/index.html */
    /* Create a subarray in bulkType object declared above                */
    /* Use oldSizes, newSizes and starts defined above                    */
    MPI_Type_create_subarray(...);
    MPI_Type_commit(&bulkType);

    /* Send the data using bulkType subarray to rank 0    */
    MPI_Isend(...);

    int newSizesI[solver->size];
    int newSizesJ[solver->size];
    MPI_Gather(&newSizes[CIDIM], 1, MPI_INT, newSizesI, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&newSizes[CJDIM], 1, MPI_INT, newSizesJ, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* rank 0 assembles the subdomains */
    if (solver->rank == 0) {
        for (int i = 0; i < solver->size; i++) {
            MPI_Datatype domainType;
            int oldSizes[NDIMS] = { jmax + 2, imax + 2 };
            int newSizes[NDIMS] = { newSizesJ[i], newSizesI[i] };
            int coords[NDIMS];
            MPI_Cart_coords(solver->comm, i, NDIMS, coords);
            int starts[NDIMS] = { sum(newSizesJ, i, 1, coords[JDIM]),
                sum(newSizesI, i, solver->dims[JDIM], coords[IDIM]) };
            printf(
                "Rank: %d, Coords(i,j): %d %d, Size(i,j): %d %d, Target Size(i,j): %d %d "
                "Starts(i,j): %d %d\n",
                i,
                coords[IDIM],
                coords[JDIM],
                oldSizes[CIDIM],
                oldSizes[CJDIM],
                newSizes[CIDIM],
                newSizes[CJDIM],
                starts[CIDIM],
                starts[CJDIM]);

            /* Refer to lecture slides for more information                       */
            /* Create a subarray in domainType object declared above              */
            /* Use oldSizes, newSizes and starts defined above                    */
            MPI_Type_create_subarray(...);

            MPI_Type_commit(&domainType);

            /* Recv the data using domainType subarray from rest of the ranks */
            MPI_Irecv(...);
            MPI_Type_free(&domainType);
        }
    }

    MPI_Waitall(numRequests, requests, MPI_STATUSES_IGNORE);
}

void collectResult(Solver* solver)
{
    size_t bytesize = (solver->imax + 2) * (solver->jmax + 2) * sizeof(double);

    double* Pall = allocate(64, bytesize);
    double* Uall = allocate(64, bytesize);
    double* Vall = allocate(64, bytesize);
    /* collect P */
    assembleResult(solver, solver->p, Pall);

    /* collect U */
    assembleResult(solver, solver->u, Uall);

    /* collect V */
    assembleResult(solver, solver->v, Vall);

    if (solver->rank == 0) {
        writeResult(solver, Pall, Uall, Vall);
    }

    free(Pall);
    free(Uall);
    free(Vall);
}

static void printConfig(Solver* solver)
{
    if (solver->rank == 0) {
        printf("Parameters for #%s#\n", solver->problem);
        printf("Boundary conditions Left:%d Right:%d Bottom:%d Top:%d\n",
            solver->bcLeft,
            solver->bcRight,
            solver->bcBottom,
            solver->bcTop);
        printf("\tReynolds number: %.2f\n", solver->re);
        printf("\tGx Gy: %.2f %.2f\n", solver->gx, solver->gy);
        printf("Geometry data:\n");
        printf("\tDomain box size (x, y): %.2f, %.2f\n",
            solver->xlength,
            solver->ylength);
        printf("\tCells (x, y): %d, %d\n", solver->imax, solver->jmax);
        printf("Timestep parameters:\n");
        printf("\tDefault stepsize: %.2f, Final time %.2f\n", solver->dt, solver->te);
        printf("\tdt bound: %.6f\n", solver->dtBound);
        printf("\tTau factor: %.2f\n", solver->tau);
        printf("Iterative solver parameters:\n");
        printf("\tMax iterations: %d\n", solver->itermax);
        printf("\tepsilon (stopping tolerance) : %f\n", solver->eps);
        printf("\tgamma factor: %f\n", solver->gamma);
        printf("\tomega (SOR relaxation): %f\n", solver->omega);
        printf("Communication parameters:\n");
    }
    for (int i = 0; i < solver->size; i++) {
        if (i == solver->rank) {
            printf("\tRank %d of %d\n", solver->rank, solver->size);
            printf("\tNeighbours (bottom, top, left, right): %d %d, %d, %d\n",
                solver->neighbours[BOTTOM],
                solver->neighbours[TOP],
                solver->neighbours[LEFT],
                solver->neighbours[RIGHT]);
            printf("\tIs boundary:\n");
            printf("\t\tLEFT: %d\n", isBoundary(solver, LEFT));
            printf("\t\tRIGHT: %d\n", isBoundary(solver, RIGHT));
            printf("\t\tBOTTOM: %d\n", isBoundary(solver, BOTTOM));
            printf("\t\tTOP: %d\n", isBoundary(solver, TOP));
            printf("\tCoordinates (i,j) %d %d\n", solver->coords[IDIM], solver->coords[JDIM]);
            printf("\tDims (i,j) %d %d\n", solver->dims[IDIM], solver->dims[JDIM]);
            printf("\tLocal domain size (i,j) %dx%d\n", solver->imaxLocal, solver->jmaxLocal);
            fflush(stdout);
        }
        MPI_Barrier(solver->comm);
    }
}

void initSolver(Solver* solver, Parameter* params)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &(solver->rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(solver->size));
    solver->problem   = params->name;
    solver->bcLeft    = params->bcLeft;
    solver->bcRight   = params->bcRight;
    solver->bcBottom  = params->bcBottom;
    solver->bcTop     = params->bcTop;
    solver->imax      = params->imax;
    solver->jmax      = params->jmax;
    solver->xlength = params->xlength;
    solver->ylength = params->ylength;
    solver->dx      = params->xlength / params->imax;
    solver->dy      = params->ylength / params->jmax;
    solver->eps     = params->eps;
    solver->omega   = params->omg;
    solver->itermax = params->itermax;
    solver->re      = params->re;
    solver->gx      = params->gx;
    solver->gy      = params->gy;
    solver->dt      = params->dt;
    solver->te      = params->te;
    solver->tau     = params->tau;
    solver->gamma   = params->gamma;

    int imax = solver->imax;
    int jmax = solver->jmax;

    /* Use virtual topologies to get dimension and coordinates of the rank */
    int dims[NDIMS]    = { 0, 0 };
    int periods[NDIMS] = { 0, 0 };

    /* Retrieves best possible #ranks in X and Y dimensions from the given communicator size. */
    /* Say for 6 ranks, best possible #ranks in in X and Y dimension would be 3 and 2         */
    /* Store in dims array                                                                    */
    /* https://rookiehpc.org/mpi/docs/mpi_dims_create/index.html                              */
    MPI_Dims_create(...);
    memcpy(&solver->dims, &dims, NDIMS);

    /* Create a new communicator based on the dimensions             */
    /* Store the new communicator in &solver->comm                   */
    /* Use NDIMS to specify dimension size                           */
    /* https://rookiehpc.org/mpi/docs/mpi_cart_create/index.html     */
    MPI_Cart_create(...);

    /* Use newly created communicator solver->comm                   */
    /* For left and right neighbors                                  */
    /* Store in solver->neighbors[LEFT] and solver->neighbors[RIGHT] */
    /* https://rookiehpc.org/mpi/docs/mpi_cart_shift/index.html      */
    MPI_Cart_shift(..);

    /* Use newly created communicator solver->comm                   */
    /* For bottom and top neighbors                                  */
    /* Store in solver->neighbors[BOTTOM] and solver->neighbors[TOP] */
    MPI_Cart_shift(...);

    /* Retrieve the cartesion topology of the rank in solver->coords */
    /* https://rookiehpc.org/mpi/docs/mpi_cart_get/index.html        */
    MPI_Cart_get(...);

    solver->imaxLocal = sizeOfRank(solver->coords[IDIM], dims[IDIM], imax);
    solver->jmaxLocal = sizeOfRank(solver->coords[JDIM], dims[JDIM], jmax);

    int imaxLocal   = solver->imaxLocal;
    int jmaxLocal   = solver->jmaxLocal;
    size_t bytesize = (imaxLocal + 2) * (jmaxLocal + 2) * sizeof(double);

    /* Refer to lecture slides for more information                  */
    /* https://rookiehpc.org/mpi/docs/mpi_type_contiguous/index.html */
    MPI_Datatype jBufferType;
    MPI_Type_contiguous(...);
    MPI_Type_commit(...);

    /* Refer to lecture slides for more information                  */
    /* https://rookiehpc.org/mpi/docs/mpi_type_vector/index.html     */
    MPI_Datatype iBufferType;
    MPI_Type_vector(...);
    MPI_Type_commit(...);

    /* Refer to lecture slides for more information */
    solver->bufferTypes[LEFT]   = iBufferType;
    solver->bufferTypes[RIGHT]  = ...;
    solver->bufferTypes[BOTTOM] = ...;
    solver->bufferTypes[TOP]    = ...;

    size_t dblsize     = sizeof(double);
    solver->sdispls[LEFT]   = ...;
    solver->sdispls[RIGHT]  = ...;
    solver->sdispls[BOTTOM] = ...;
    solver->sdispls[TOP]    = ...;

    solver->rdispls[LEFT]   = ...;
    solver->rdispls[RIGHT]  = ...;
    solver->rdispls[BOTTOM] = ...;
    solver->rdispls[TOP]    = ...;

    solver->u       = allocate(64, bytesize);
    solver->v       = allocate(64, bytesize);
    solver->p       = allocate(64, bytesize);
    solver->rhs     = allocate(64, bytesize);
    solver->f       = allocate(64, bytesize);
    solver->g       = allocate(64, bytesize);

    /* TODO adapt to local size imaxLocal & jmaxLocal */
    for (int i = 0; i < (imax + 2) * (jmax + 2); i++) { 
        solver->u[i]   = params->u_init;
        solver->v[i]   = params->v_init;
        solver->p[i]   = params->p_init;
        solver->rhs[i] = 0.0;
        solver->f[i]   = 0.0;
        solver->g[i]   = 0.0;
    }

    double dx        = solver->dx;
    double dy        = solver->dy;
    double invSqrSum = 1.0 / (dx * dx) + 1.0 / (dy * dy);
    solver->dtBound  = 0.5 * solver->re * 1.0 / invSqrSum;
#ifdef VERBOSE
    printConfig(solver);
#endif
}

void computeRHS(Solver* solver)
{
    int imax    = solver->imax;
    int jmax    = solver->jmax;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double idx  = 1.0 / solver->dx;
    double idy  = 1.0 / solver->dy;
    double idt  = 1.0 / solver->dt;
    double* rhs = solver->rhs;
    double* f   = solver->f;
    double* g   = solver->g;

    /* Shift function is similar to exchange. Just that you need to exchange values for F and G array */
    /* A simple Send Recv function for F and G array can be implemented. */
    /* You can look at the formulation below to guess which values will be required. */
    /* For F, you just have to send the values from left neighbor to right neighbor's ghost cells */
    /* For G, you just have to send the values from bottom neighbor to top neighbor's ghost cells */

    shift(solver);

    /* TODO adapt to local size imaxLocal & jmaxLocal */
    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            RHS(i, j) = ((F(i, j) - F(i - 1, j)) * idx + (G(i, j) - G(i, j - 1)) * idy) *
                        idt;
        }
    }
}

void solve(Solver* solver)
{
    int imax      = solver->imax;
    int jmax      = solver->jmax;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double eps    = solver->eps;
    int itermax   = solver->itermax;
    double dx2    = solver->dx * solver->dx;
    double dy2    = solver->dy * solver->dy;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double factor = solver->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double* p     = solver->p;
    double* rhs   = solver->rhs;
    double epssq  = eps * eps;
    int it        = 0;
    double res    = 1.0;

    while ((res >= epssq) && (it < itermax)) {
        res = 0.0;

        /* Exchange the ghost cells with respective neighbors */
        /* Remember this is a 2D domain decomposition.        */
        exchange(solver, p);

        /* TODO adapt to local size imaxLocal & jmaxLocal     */
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {

                double r = RHS(i, j) -
                           ((P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) * idx2 +
                               (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) * idy2);

                P(i, j) -= (factor * r);
                res += (r * r);
            }
        }

        /* Extend boundary condition handling from Poisson assignment   */
        /* use isBoundary(solver, LEFT/RIGHT/TOP/BOTTOM) to check       */
        /* whether a rank has boundary and apply BC only for those rank */
        /* Use indices for more intuition.                              */
        for (int i = 1; i < imax + 1; i++) {
            P(i, 0)        = P(i, 1);
            P(i, jmax + 1) = P(i, jmax);
        }

        for (int j = 1; j < jmax + 1; j++) {
            P(0, j)        = P(1, j);
            P(imax + 1, j) = P(imax, j);
        }

        /* TODO add global MPI_SUM reduction for res */
        MPI_Allreduce(...);
        
        res = res / (double)(imax * jmax);
#ifdef DEBUG
        printf("%d Residuum: %e\n", it, res);
#endif
        it++;
    }

#ifdef VERBOSE
    if(solver->rank == 0)
    printf("Solver took %d iterations to reach %f\n", it, sqrt(res));
#endif
}

static double maxElement(Solver* solver, double* m)
{
    /* TODO adapt to local size imaxLocal & jmaxLocal */
    int size      = (solver->imax + 2) * (solver->jmax + 2);
    double maxval = DBL_MIN;

    for (int i = 0; i < size; i++) {
        maxval = MAX(maxval, fabs(m[i]));
    }

    /* TODO add global MPI_MAX reduction for maxval */
    MPI_Allreduce(...);

    return maxval;
}

void normalizePressure(Solver* solver)
{
    /* TODO adapt to local size imaxLocal & jmaxLocal */
    int size    = (solver->imax + 2) * (solver->jmax + 2);
    double* p   = solver->p;
    double avgP = 0.0;

    for (int i = 0; i < size; i++) {
        avgP += p[i];
    }
    avgP /= size;
    
    /* TODO add global MPI_SUM reduction for avgP */
    MPI_Allreduce(...);

    /* do not adapt to local size imaxLocal & jmaxLocal */
    avgP /= (solver->imax + 2) * (solver->jmax + 2);

    for (int i = 0; i < size; i++) {
        p[i] = p[i] - avgP;
    }
}

void computeTimestep(Solver* solver)
{
    double dt   = solver->dtBound;
    double dx   = solver->dx;
    double dy   = solver->dy;
    double umax = maxElement(solver, solver->u);
    double vmax = maxElement(solver, solver->v);

    if (umax > 0) {
        dt = (dt > dx / umax) ? dx / umax : dt;
    }
    if (vmax > 0) {
        dt = (dt > dy / vmax) ? dy / vmax : dt;
    }

    solver->dt = dt * solver->tau;
}

/* TODO adapt to local size imaxLocal & jmaxLocal                                           */
/* TODO Use isBoundary(solver, TOP/BOTTOM/LEFT/RIGHT) to check                              */
/* whether the rank has any boundaries and apply NOSLIP, SLIP or OUTFLOW boundary condition */
void setBoundaryConditions(Solver* solver)
{
    int imax = solver->imax;
    int jmax = solver->jmax;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;
    double* u = solver->u;
    double* v = solver->v;

    // Left boundary
    switch (solver->bcLeft) {
    case NOSLIP:
        for (int j = 1; j < jmax + 1; j++) {
            U(0, j) = 0.0;
            V(0, j) = -V(1, j);
        }
        break;
    case SLIP:
        for (int j = 1; j < jmax + 1; j++) {
            U(0, j) = 0.0;
            V(0, j) = V(1, j);
        }
        break;
    case OUTFLOW:
        for (int j = 1; j < jmax + 1; j++) {
            U(0, j) = U(1, j);
            V(0, j) = V(1, j);
        }
        break;
    case PERIODIC:
        break;
    }

    // Right boundary
    switch (solver->bcRight) {
    case NOSLIP:
        for (int j = 1; j < jmax + 1; j++) {
            U(imax, j)     = 0.0;
            V(imax + 1, j) = -V(imax, j);
        }
        break;
    case SLIP:
        for (int j = 1; j < jmax + 1; j++) {
            U(imax, j)     = 0.0;
            V(imax + 1, j) = V(imax, j);
        }
        break;
    case OUTFLOW:
        for (int j = 1; j < jmax + 1; j++) {
            U(imax, j)     = U(imax - 1, j);
            V(imax + 1, j) = V(imax, j);
        }
        break;
    case PERIODIC:
        break;
    }

    // Bottom boundary
    switch (solver->bcBottom) {
    case NOSLIP:
        for (int i = 1; i < imax + 1; i++) {
            V(i, 0) = 0.0;
            U(i, 0) = -U(i, 1);
        }
        break;
    case SLIP:
        for (int i = 1; i < imax + 1; i++) {
            V(i, 0) = 0.0;
            U(i, 0) = U(i, 1);
        }
        break;
    case OUTFLOW:
        for (int i = 1; i < imax + 1; i++) {
            U(i, 0) = U(i, 1);
            V(i, 0) = V(i, 1);
        }
        break;
    case PERIODIC:
        break;
    }

    // Top boundary
    switch (solver->bcTop) {
    case NOSLIP:
        for (int i = 1; i < imax + 1; i++) {
            V(i, jmax)     = 0.0;
            U(i, jmax + 1) = -U(i, jmax);
        }
        break;
    case SLIP:
        for (int i = 1; i < imax + 1; i++) {
            V(i, jmax)     = 0.0;
            U(i, jmax + 1) = U(i, jmax);
        }
        break;
    case OUTFLOW:
        for (int i = 1; i < imax + 1; i++) {
            U(i, jmax + 1) = U(i, jmax);
            V(i, jmax)     = V(i, jmax - 1);
        }
        break;
    case PERIODIC:
        break;
    }
}

/* TODO adapt to local size imaxLocal & jmaxLocal */
void setSpecialBoundaryCondition(Solver* solver)
{
    int imax   = solver->imax;
    int jmax   = solver->jmax;
    double mDy = solver->dy;
    double* u  = solver->u;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;


    /* In case of dcavity, apply special boundary condition to those */
    /* ranks who has top boundary. Check the indices for intuition.  */
    /* Use isBoundary(solver, TOP/BOTTOM/LEFT/RIGHT) to check if a   */
    /* rank has TOP boundary and the apply to it.                    */
    if (strcmp(solver->problem, "dcavity") == 0) {
        for (int i = 1; i < imax; i++) {
            U(i, jmax + 1) = 2.0 - U(i, jmax);
        }
    } 
    
    /* In case of canal, apply special boundary condition to those        */
    /* ranks who has left boundary. Check the indices for intuition.      */
    /* Use isBoundary(solver, TOP/BOTTOM/LEFT/RIGHT) to check if a        */
    /* rank has LEFT boundary and the apply to it.                        */
    /* Also you have to adapt the special boundary condition for each     */
    /* rank with LEFT boundary because in sequential case, it a           */
    /* parabola. Use solver->dims to adjust the formula according to rank */
    else if (strcmp(solver->problem, "canal") == 0) {
        double ylength = solver->ylength;
        double y;

        for (int j = 1; j < jmax + 1; j++) {
            y       = mDy * (j - 0.5);
            U(0, j) = y * (ylength - y) * 4.0 / (ylength * ylength);
        }
    }
}

/* TODO adapt to local size imaxLocal & jmaxLocal */
void computeFG(Solver* solver)
{
    double* u        = solver->u;
    double* v        = solver->v;
    double* f        = solver->f;
    double* g        = solver->g;
    int imax         = solver->imax;
    int jmax         = solver->jmax;
    double gx        = solver->gx;
    double gy        = solver->gy;
    double gamma     = solver->gamma;
    double dt        = solver->dt;
    double inverseRe = 1.0 / solver->re;
    double inverseDx = 1.0 / solver->dx;
    double inverseDy = 1.0 / solver->dy;
    double du2dx, dv2dy, duvdx, duvdy;
    double du2dx2, du2dy2, dv2dx2, dv2dy2;
    int imaxLocal = solver->imaxLocal;
    int jmaxLocal = solver->jmaxLocal;

    exchange(solver, u);
    exchange(solver, v);

    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            du2dx = inverseDx * 0.25 *
                        ((U(i, j) + U(i + 1, j)) * (U(i, j) + U(i + 1, j)) -
                            (U(i, j) + U(i - 1, j)) * (U(i, j) + U(i - 1, j))) +
                    gamma * inverseDx * 0.25 *
                        (fabs(U(i, j) + U(i + 1, j)) * (U(i, j) - U(i + 1, j)) +
                            fabs(U(i, j) + U(i - 1, j)) * (U(i, j) - U(i - 1, j)));

            duvdy = inverseDy * 0.25 *
                        ((V(i, j) + V(i + 1, j)) * (U(i, j) + U(i, j + 1)) -
                            (V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j) + U(i, j - 1))) +
                    gamma * inverseDy * 0.25 *
                        (fabs(V(i, j) + V(i + 1, j)) * (U(i, j) - U(i, j + 1)) +
                            fabs(V(i, j - 1) + V(i + 1, j - 1)) *
                                (U(i, j) - U(i, j - 1)));

            du2dx2  = inverseDx * inverseDx * (U(i + 1, j) - 2.0 * U(i, j) + U(i - 1, j));
            du2dy2  = inverseDy * inverseDy * (U(i, j + 1) - 2.0 * U(i, j) + U(i, j - 1));
            F(i, j) = U(i, j) + dt * (inverseRe * (du2dx2 + du2dy2) - du2dx - duvdy + gx);

            duvdx = inverseDx * 0.25 *
                        ((U(i, j) + U(i, j + 1)) * (V(i, j) + V(i + 1, j)) -
                            (U(i - 1, j) + U(i - 1, j + 1)) * (V(i, j) + V(i - 1, j))) +
                    gamma * inverseDx * 0.25 *
                        (fabs(U(i, j) + U(i, j + 1)) * (V(i, j) - V(i + 1, j)) +
                            fabs(U(i - 1, j) + U(i - 1, j + 1)) *
                                (V(i, j) - V(i - 1, j)));

            dv2dy = inverseDy * 0.25 *
                        ((V(i, j) + V(i, j + 1)) * (V(i, j) + V(i, j + 1)) -
                            (V(i, j) + V(i, j - 1)) * (V(i, j) + V(i, j - 1))) +
                    gamma * inverseDy * 0.25 *
                        (fabs(V(i, j) + V(i, j + 1)) * (V(i, j) - V(i, j + 1)) +
                            fabs(V(i, j) + V(i, j - 1)) * (V(i, j) - V(i, j - 1)));

            dv2dx2  = inverseDx * inverseDx * (V(i + 1, j) - 2.0 * V(i, j) + V(i - 1, j));
            dv2dy2  = inverseDy * inverseDy * (V(i, j + 1) - 2.0 * V(i, j) + V(i, j - 1));
            G(i, j) = V(i, j) + dt * (inverseRe * (dv2dx2 + dv2dy2) - duvdx - dv2dy + gy);
        }
    }

    /* Adapt boundary conditons for F only for ranks with LEFT and RIGHT boundary */
    /* Check the indices for more intuition. Use isBoundary(solver, LEFT/RIGHT)   */
    /* to check if a rank has LEFT or RIGHT boundary.                             */
    /* ----------------------------- boundary of F --------------------------- */
    for (int j = 1; j < jmax + 1; j++) {
        F(0, j)    = U(0, j);
        F(imax, j) = U(imax, j);
    }

    /* Adapt boundary conditons for G only for ranks with TOP and BOTTOM boundary */
    /* Check the indices for more intuition. Use isBoundary(solver, TOP/BOTTOM)   */
    /* to check if a rank has TOP or BOTTOM boundary.                             */
    /* ----------------------------- boundary of G --------------------------- */
    for (int i = 1; i < imax + 1; i++) {
        G(i, 0)    = V(i, 0);
        G(i, jmax) = V(i, jmax);
    }
}

/* TODO adapt to local size imaxLocal & jmaxLocal */
void adaptUV(Solver* solver)
{
    int imax       = solver->imax;
    int jmax       = solver->jmax;
    int imaxLocal  = solver->imaxLocal;
    int jmaxLocal  = solver->jmaxLocal;
    double* p      = solver->p;
    double* u      = solver->u;
    double* v      = solver->v;
    double* f      = solver->f;
    double* g      = solver->g;
    double factorX = solver->dt / solver->dx;
    double factorY = solver->dt / solver->dy;

    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            U(i, j) = F(i, j) - (P(i + 1, j) - P(i, j)) * factorX;
            V(i, j) = G(i, j) - (P(i, j + 1) - P(i, j)) * factorY;
        }
    }
}

void writeResult(Solver* solver, double* Pall, double* Uall, double* Vall)
{
    
    int imax  = solver->imax;
    int jmax  = solver->jmax;
    double dx = solver->dx;
    double dy = solver->dy;
    double x = 0.0, y = 0.0;

    FILE* fp;
    fp = fopen("pressure.dat", "w");

    if (fp == NULL) {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j < jmax + 1; j++) {
        y = (double)(j - 0.5) * dy;
        for (int i = 1; i < imax + 1; i++) {
            x = (double)(i - 0.5) * dx;
            fprintf(fp, "%.2f %.2f %f\n", x, y, Pall(i, j));
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

    fp = fopen("velocity.dat", "w");

    if (fp == NULL) {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 1; j < jmax + 1; j++) {
        y = dy * (j - 0.5);
        for (int i = 1; i < imax + 1; i++) {
            x           = dx * (i - 0.5);
            double velU = (Uall(i, j) + Uall(i - 1, j)) / 2.0;
            double velV = (Vall(i, j) + Vall(i, j - 1)) / 2.0;
            double len  = sqrt((velU * velU) + (velV * velV));
            fprintf(fp, "%.2f %.2f %f %f %f\n", x, y, velU, velV, len);
        }
    }

    fclose(fp);
}
