/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocate.h"
#include "parameter.h"
#include "solver.h"
#include "util.h"

#define P(i, j, k)   p[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define F(i, j, k)   f[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define G(i, j, k)   g[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define H(i, j, k)   h[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define U(i, j, k)   u[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define V(i, j, k)   v[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define W(i, j, k)   w[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]
#define RHS(i, j, k) rhs[(k) * (imax + 2) * (jmax + 2) + (j) * (imax + 2) + (i)]

static void printConfig(Solver* s)
{
    printf("Parameters for #%s#\n", s->problem);
    printf("BC Left:%d Right:%d Bottom:%d Top:%d Front:%d Back:%d\n",
        s->bcLeft,
        s->bcRight,
        s->bcBottom,
        s->bcTop,
        s->bcFront,
        s->bcBack);
    printf("\tReynolds number: %.2f\n", s->re);
    printf("\tGx Gy: %.2f %.2f %.2f\n", s->gx, s->gy, s->gz);
    printf("Geometry data:\n");
    printf("\tDomain box size (x, y, z): %.2f, %.2f, %.2f\n",
        s->grid.xlength,
        s->grid.ylength,
        s->grid.zlength);
        printf("\tCells (x, y, z): %d, %d, %d\n",
            s->grid.imax,
            s->grid.jmax,
            s->grid.kmax);
        printf("\tCell size (dx, dy, dz): %f, %f, %f\n",
            s->grid.dx,
            s->grid.dy,
            s->grid.dz);
    printf("Timestep parameters:\n");
    printf("\tDefault stepsize: %.2f, Final time %.2f\n", s->dt, s->te);
    printf("\tdt bound: %.6f\n", s->dtBound);
    printf("\tTau factor: %.2f\n", s->tau);
        printf("Iterative parameters:\n");
    printf("\tMax iterations: %d\n", s->itermax);
    printf("\tepsilon (stopping tolerance) : %f\n", s->eps);
    printf("\tgamma factor: %f\n", s->gamma);
    printf("\tomega (SOR relaxation): %f\n", s->omega);
}

void initSolver(Solver* s, Parameter* params)
{
    s->problem      = params->name;
    s->bcLeft       = params->bcLeft;
    s->bcRight      = params->bcRight;
    s->bcBottom     = params->bcBottom;
    s->bcTop        = params->bcTop;
    s->bcFront      = params->bcFront;
    s->bcBack       = params->bcBack;

    s->grid.imax    = params->imax;
    s->grid.jmax    = params->jmax;
    s->grid.kmax    = params->kmax;
    s->grid.xlength = params->xlength;
    s->grid.ylength = params->ylength;
    s->grid.zlength = params->zlength;
    s->grid.dx      = params->xlength / params->imax;
    s->grid.dy      = params->ylength / params->jmax;
    s->grid.dz      = params->zlength / params->kmax;

    s->eps          = params->eps;
    s->omega        = params->omg;
    s->itermax      = params->itermax;
    s->re           = params->re;
    s->gx           = params->gx;
    s->gy           = params->gy;
    s->gz           = params->gz;
    s->dt           = params->dt;
    s->te           = params->te;
    s->tau          = params->tau;
    s->gamma        = params->gamma;

    int imax        = s->grid.imax;
    int jmax        = s->grid.jmax;
    int kmax        = s->grid.kmax;
    size_t bytesize = (imax + 2) * (jmax + 2) * (kmax + 2) * sizeof(double);
    s->u            = allocate(64, bytesize);
    s->v            = allocate(64, bytesize);
    s->w            = allocate(64, bytesize);
    s->p            = allocate(64, bytesize);
    s->rhs          = allocate(64, bytesize);
    s->f            = allocate(64, bytesize);
    s->g            = allocate(64, bytesize);
    s->h            = allocate(64, bytesize);

    for (int i = 0; i < (imax + 2) * (jmax + 2) * (kmax + 2); i++) {
        s->u[i]   = params->u_init;
        s->v[i]   = params->v_init;
        s->w[i]   = params->w_init;
        s->p[i]   = params->p_init;
        s->rhs[i] = 0.0;
        s->f[i]   = 0.0;
        s->g[i]   = 0.0;
        s->h[i]   = 0.0;
    }

    double dx = s->grid.dx;
    double dy = s->grid.dy;
    double dz = s->grid.dz;

    double invSqrSum = 1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);
    s->dtBound       = 0.5 * s->re * 1.0 / invSqrSum;

#ifdef VERBOSE
    printConfig(s);
#endif /* VERBOSE */
}

void computeRHS(Solver* s)
{
    int imax    = s->grid.imax;
    int jmax    = s->grid.jmax;
    int kmax    = s->grid.kmax;
    double idx  = 1.0 / s->grid.dx;
    double idy  = 1.0 / s->grid.dy;
    double idz  = 1.0 / s->grid.dz;
    double idt  = 1.0 / s->dt;

    double* rhs = s->rhs;
    double* f   = s->f;
    double* g   = s->g;
    double* h   = s->h;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                RHS(i, j, k) = ((F(i, j, k) - F(i - 1, j, k)) * idx +
                                   (G(i, j, k) - G(i, j - 1, k)) * idy +
                                   (H(i, j, k) - H(i, j, k - 1)) * idz) *
                               idt;
            }
        }
    }
}

void solve(Solver* s)
{
    int imax      = s->grid.imax;
    int jmax      = s->grid.jmax;
    int kmax      = s->grid.kmax;

    double eps    = s->eps;
    int itermax   = s->itermax;
    double dx2    = s->grid.dx * s->grid.dx;
    double dy2    = s->grid.dy * s->grid.dy;
    double dz2    = s->grid.dz * s->grid.dz;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double idz2   = 1.0 / dz2;

    double factor = s->omega * 0.5 * (dx2 * dy2 * dz2) /
                    (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);
    double* p    = s->p;
    double* rhs  = s->rhs;
    double epssq = eps * eps;
    int it       = 0;
    double res   = 1.0;

    while ((res >= epssq) && (it < itermax)) {
        res = 0.0;

        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                for (int i = 1; i < imax + 1; i++) {

                    double r = RHS(i, j, k) -
                               ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) *
                                       idx2 +
                                   (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) *
                                       idy2 +
                                   (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) *
                                       idz2);

                    P(i, j, k) -= (factor * r);
                    res += (r * r);
                }
            }
        }

        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, j, 0)        = P(i, j, 1);
                P(i, j, kmax + 1) = P(i, j, kmax);
            }
        }

        for (int k = 1; k < kmax + 1; k++) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, 0, k)        = P(i, 1, k);
                P(i, jmax + 1, k) = P(i, jmax, k);
            }
        }

        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                P(0, j, k)        = P(1, j, k);
                P(imax + 1, j, k) = P(imax, j, k);
            }
        }

        res = res / (double)(imax * jmax * kmax);
#ifdef DEBUG
        printf("%d Residuum: %e\n", it, res);
#endif
        it++;
    }

#ifdef VERBOSE
    printf("Solver took %d iterations to reach %f\n", it, sqrt(res));
#endif
}

void solveRB(Solver* s)
{
    int imax      = s->grid.imax;
    int jmax      = s->grid.jmax;
    int kmax      = s->grid.kmax;
    double eps    = s->eps;
    int itermax   = s->itermax;
    double dx2    = s->grid.dx * s->grid.dx;
    double dy2    = s->grid.dy * s->grid.dy;
    double dz2    = s->grid.dz * s->grid.dz;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double idz2   = 1.0 / dz2;
    double factor = s->omega * 0.5 * (dx2 * dy2 * dz2) /
                    (dy2 * dz2 + dx2 * dz2 + dx2 * dy2);
    double* p    = s->p;
    double* rhs  = s->rhs;
    double epssq = eps * eps;
    int it       = 0;
    double res   = 1.0;
    int pass, ksw, jsw, isw;

    while ((res >= epssq) && (it < itermax)) {
        res = 0.0;
        ksw = 1;

        for (pass = 0; pass < 2; pass++) {
            jsw = ksw;

            for (int k = 1; k < kmax + 1; k++) {
                isw = jsw;
                for (int j = 1; j < jmax + 1; j++) {
                    for (int i = isw; i < imax + 1; i += 2) {

                        double r =
                            RHS(i, j, k) -
                            ((P(i + 1, j, k) - 2.0 * P(i, j, k) + P(i - 1, j, k)) * idx2 +
                                (P(i, j + 1, k) - 2.0 * P(i, j, k) + P(i, j - 1, k)) *
                                    idy2 +
                                (P(i, j, k + 1) - 2.0 * P(i, j, k) + P(i, j, k - 1)) *
                                    idz2);

                        P(i, j, k) -= (factor * r);
                        res += (r * r);
                    }
                    isw = 3 - isw;
                }
                jsw = 3 - jsw;
            }
            ksw = 3 - ksw;
        }

        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, j, 0)        = P(i, j, 1);
                P(i, j, kmax + 1) = P(i, j, kmax);
            }
        }

        for (int k = 1; k < kmax + 1; k++) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, 0, k)        = P(i, 1, k);
                P(i, jmax + 1, k) = P(i, jmax, k);
            }
        }

        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                P(0, j, k)        = P(1, j, k);
                P(imax + 1, j, k) = P(imax, j, k);
            }
        }

        res = res / (double)(imax * jmax * kmax);
#ifdef DEBUG
        printf("%d Residuum: %e\n", it, res);
#endif
        it++;
    }

#ifdef VERBOSE
    printf("Solver took %d iterations to reach %f\n", it, sqrt(res));
#endif
}

static double maxElement(Solver* s, double* m)
{
    int size      = (s->grid.imax + 2) * (s->grid.jmax + 2) * (s->grid.kmax + 2);
    double maxval = DBL_MIN;

    for (int i = 0; i < size; i++) {
        maxval = MAX(maxval, fabs(m[i]));
    }

    return maxval;
}

void normalizePressure(Solver* s)
{
    int size    = (s->grid.imax + 2) * (s->grid.jmax + 2) * (s->grid.kmax + 2);
    double* p   = s->p;
    double avgP = 0.0;

    for (int i = 0; i < size; i++) {
        avgP += p[i];
    }
    avgP /= size;

    for (int i = 0; i < size; i++) {
        p[i] = p[i] - avgP;
    }
}

void computeTimestep(Solver* s)
{
    double dt   = s->dtBound;
    double dx   = s->grid.dx;
    double dy   = s->grid.dy;
    double dz   = s->grid.dz;

    double umax = maxElement(s, s->u);
    double vmax = maxElement(s, s->v);
    double wmax = maxElement(s, s->w);

    if (umax > 0) {
        dt = (dt > dx / umax) ? dx / umax : dt;
    }
    if (vmax > 0) {
        dt = (dt > dy / vmax) ? dy / vmax : dt;
    }
    if (wmax > 0) {
        dt = (dt > dz / wmax) ? dz / wmax : dt;
    }

    s->dt = dt * s->tau;
}

void setBoundaryConditions(Solver* s)
{
    int imax = s->grid.imax;
    int jmax = s->grid.jmax;
    int kmax = s->grid.kmax;

    double* u = s->u;
    double* v = s->v;
    double* w = s->w;

    switch (s->bcTop) {
    case NOSLIP:
        for (int k = 1; k < kmax + 1; k++) {
            for (int i = 1; i < imax + 1; i++) {
                V(i, jmax, k)     = 0.0;
                U(i, jmax + 1, k) = -U(i, jmax, k);
                W(i, jmax + 1, k) = -W(i, jmax, k);
            }
        }
        break;
    case SLIP:
        for (int k = 1; k < kmax + 1; k++) {
            for (int i = 1; i < imax + 1; i++) {
                V(i, jmax, k)     = 0.0;
                U(i, jmax + 1, k) = U(i, jmax, k);
                W(i, jmax + 1, k) = W(i, jmax, k);
            }
        }
        break;
    case OUTFLOW:
        for (int k = 1; k < kmax + 1; k++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, jmax + 1, k) = U(i, jmax, k);
                V(i, jmax, k)     = V(i, jmax - 1, k);
                W(i, jmax + 1, k) = W(i, jmax, k);
            }
        }
        break;
    case PERIODIC:
        break;
    }

    switch (s->bcBottom) {
    case NOSLIP:
        for (int k = 1; k < kmax + 1; k++) {
            for (int i = 1; i < imax + 1; i++) {
                V(i, 0, k) = 0.0;
                U(i, 0, k) = -U(i, 1, k);
                W(i, 0, k) = -W(i, 1, k);
            }
        }
        break;
    case SLIP:
        for (int k = 1; k < kmax + 1; k++) {
            for (int i = 1; i < imax + 1; i++) {
                V(i, 0, k) = 0.0;
                U(i, 0, k) = U(i, 1, k);
                W(i, 0, k) = W(i, 1, k);
            }
        }
        break;
    case OUTFLOW:
        for (int k = 1; k < kmax + 1; k++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, 0, k) = U(i, 1, k);
                V(i, 0, k) = V(i, 1, k);
                W(i, 0, k) = W(i, 1, k);
            }
        }
        break;
    case PERIODIC:
        break;
    }

    switch (s->bcLeft) {
    case NOSLIP:
        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                U(0, j, k) = 0.0;
                V(0, j, k) = -V(1, j, k);
                W(0, j, k) = -W(1, j, k);
            }
        }
        break;
    case SLIP:
        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                U(0, j, k) = 0.0;
                V(0, j, k) = V(1, j, k);
                W(0, j, k) = W(1, j, k);
            }
        }
        break;
    case OUTFLOW:
        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                U(0, j, k) = U(1, j, k);
                V(0, j, k) = V(1, j, k);
                W(0, j, k) = W(1, j, k);
            }
        }
        break;
    case PERIODIC:
        break;
    }

    switch (s->bcRight) {
    case NOSLIP:
        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                U(imax, j, k)     = 0.0;
                V(imax + 1, j, k) = -V(imax, j, k);
                W(imax + 1, j, k) = -W(imax, j, k);
            }
        }
        break;
    case SLIP:
        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                U(imax, j, k)     = 0.0;
                V(imax + 1, j, k) = V(imax, j, k);
                W(imax + 1, j, k) = W(imax, j, k);
            }
        }
        break;
    case OUTFLOW:
        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                U(imax, j, k)     = U(imax - 1, j, k);
                V(imax + 1, j, k) = V(imax, j, k);
                W(imax + 1, j, k) = W(imax, j, k);
            }
        }
        break;
    case PERIODIC:
        break;
    }

    switch (s->bcFront) {
    case NOSLIP:
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, j, 0) = -U(i, j, 1);
                V(i, j, 0) = -V(i, j, 1);
                W(i, j, 0) = 0.0;
            }
        }
        break;
    case SLIP:
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, j, 0) = U(i, j, 1);
                V(i, j, 0) = V(i, j, 1);
                W(i, j, 0) = 0.0;
            }
        }
        break;
    case OUTFLOW:
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, j, 0) = U(i, j, 1);
                V(i, j, 0) = V(i, j, 1);
                W(i, j, 0) = W(i, j, 1);
            }
        }
        break;
    case PERIODIC:
        break;
    }

    switch (s->bcBack) {
    case NOSLIP:
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, j, kmax + 1) = -U(i, j, kmax);
                V(i, j, kmax + 1) = -V(i, j, kmax);
                W(i, j, kmax + 1) = 0.0;
            }
        }
        break;
    case SLIP:
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, j, kmax + 1) = U(i, j, kmax);
                V(i, j, kmax + 1) = V(i, j, kmax);
                W(i, j, kmax + 1) = 0.0;
            }
        }
        break;
    case OUTFLOW:
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, j, kmax + 1) = U(i, j, kmax);
                V(i, j, kmax + 1) = V(i, j, kmax);
                W(i, j, kmax)     = W(i, j, kmax - 1);
            }
        }
        break;
    case PERIODIC:
        break;
    }
}

void setSpecialBoundaryCondition(Solver* s)
{
    int imax = s->grid.imax;
    int jmax = s->grid.jmax;
    int kmax = s->grid.kmax;

    double mDy = s->grid.dy;
    double* u  = s->u;

    if (strcmp(s->problem, "dcavity") == 0) {
        for (int k = 1; k < kmax; k++) {
            for (int i = 1; i < imax; i++) {
                U(i, jmax + 1, k) = 2.0 - U(i, jmax, k);
            }
        }
    } else if (strcmp(s->problem, "canal") == 0) {
        double ylength = s->grid.ylength;
        double y;

        for (int k = 1; k < kmax + 1; k++) {
            for (int j = 1; j < jmax + 1; j++) {
                y          = mDy * (j - 0.5);
                U(0, j, k) = y * (ylength - y) * 4.0 / (ylength * ylength);
            }
        }
    }
}

void computeFG(Solver* s)
{
    int imax = s->grid.imax;
    int jmax = s->grid.jmax;
    int kmax = s->grid.kmax;

    double* u = s->u;
    double* v = s->v;
    double* w = s->w;
    double* f = s->f;
    double* g = s->g;
    double* h = s->h;

    double gx    = s->gx;
    double gy    = s->gy;
    double gz    = s->gz;
    double dt    = s->dt;

    double gamma     = s->gamma;
    double inverseRe = 1.0 / s->re;
    double inverseDx = 1.0 / s->grid.dx;
    double inverseDy = 1.0 / s->grid.dy;
    double inverseDz = 1.0 / s->grid.dz;
    double du2dx, dv2dy, dw2dz;
    double duvdx, duwdx, duvdy, dvwdy, duwdz, dvwdz;
    double du2dx2, du2dy2, du2dz2;
    double dv2dx2, dv2dy2, dv2dz2;
    double dw2dx2, dw2dy2, dw2dz2;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                du2dx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i + 1, j, k)) *
                                    (U(i, j, k) + U(i + 1, j, k)) -
                                (U(i, j, k) + U(i - 1, j, k)) *
                                    (U(i, j, k) + U(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i + 1, j, k)) *
                                    (U(i, j, k) - U(i + 1, j, k)) +
                                fabs(U(i, j, k) + U(i - 1, j, k)) *
                                    (U(i, j, k) - U(i - 1, j, k)));

                duvdy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i + 1, j, k)) *
                                    (U(i, j, k) + U(i, j + 1, k)) -
                                (V(i, j - 1, k) + V(i + 1, j - 1, k)) *
                                    (U(i, j, k) + U(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i + 1, j, k)) *
                                    (U(i, j, k) - U(i, j + 1, k)) +
                                fabs(V(i, j - 1, k) + V(i + 1, j - 1, k)) *
                                    (U(i, j, k) - U(i, j - 1, k)));

                duwdz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i + 1, j, k)) *
                                    (U(i, j, k) + U(i, j, k + 1)) -
                                (W(i, j, k - 1) + W(i + 1, j, k - 1)) *
                                    (U(i, j, k) + U(i, j, k - 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i + 1, j, k)) *
                                    (U(i, j, k) - U(i, j, k + 1)) +
                                fabs(W(i, j, k - 1) + W(i + 1, j, k - 1)) *
                                    (U(i, j, k) - U(i, j, k - 1)));

                du2dx2 = inverseDx * inverseDx *
                         (U(i + 1, j, k) - 2.0 * U(i, j, k) + U(i - 1, j, k));
                du2dy2 = inverseDy * inverseDy *
                         (U(i, j + 1, k) - 2.0 * U(i, j, k) + U(i, j - 1, k));
                du2dz2 = inverseDz * inverseDz *
                         (U(i, j, k + 1) - 2.0 * U(i, j, k) + U(i, j, k - 1));
                F(i, j, k) = U(i, j, k) + dt * (inverseRe * (du2dx2 + du2dy2 + du2dz2) -
                                                   du2dx - duvdy - duwdz + gx);

                duvdx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i, j + 1, k)) *
                                    (V(i, j, k) + V(i + 1, j, k)) -
                                (U(i - 1, j, k) + U(i - 1, j + 1, k)) *
                                    (V(i, j, k) + V(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i, j + 1, k)) *
                                    (V(i, j, k) - V(i + 1, j, k)) +
                                fabs(U(i - 1, j, k) + U(i - 1, j + 1, k)) *
                                    (V(i, j, k) - V(i - 1, j, k)));

                dv2dy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i, j + 1, k)) *
                                    (V(i, j, k) + V(i, j + 1, k)) -
                                (V(i, j, k) + V(i, j - 1, k)) *
                                    (V(i, j, k) + V(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i, j + 1, k)) *
                                    (V(i, j, k) - V(i, j + 1, k)) +
                                fabs(V(i, j, k) + V(i, j - 1, k)) *
                                    (V(i, j, k) - V(i, j - 1, k)));

                dvwdz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i, j + 1, k)) *
                                    (V(i, j, k) + V(i, j, k + 1)) -
                                (W(i, j, k - 1) + W(i, j + 1, k - 1)) *
                                    (V(i, j, k) + V(i, j, k + 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i, j + 1, k)) *
                                    (V(i, j, k) - V(i, j, k + 1)) +
                                fabs(W(i, j, k - 1) + W(i, j + 1, k - 1)) *
                                    (V(i, j, k) - V(i, j, k + 1)));

                dv2dx2 = inverseDx * inverseDx *
                         (V(i + 1, j, k) - 2.0 * V(i, j, k) + V(i - 1, j, k));
                dv2dy2 = inverseDy * inverseDy *
                         (V(i, j + 1, k) - 2.0 * V(i, j, k) + V(i, j - 1, k));
                dv2dz2 = inverseDz * inverseDz *
                         (V(i, j, k + 1) - 2.0 * V(i, j, k) + V(i, j, k - 1));
                G(i, j, k) = V(i, j, k) + dt * (inverseRe * (dv2dx2 + dv2dy2 + dv2dz2) -
                                                   duvdx - dv2dy - dvwdz + gy);

                duwdx = inverseDx * 0.25 *
                            ((U(i, j, k) + U(i, j, k + 1)) *
                                    (W(i, j, k) + W(i + 1, j, k)) -
                                (U(i - 1, j, k) + U(i - 1, j, k + 1)) *
                                    (W(i, j, k) + W(i - 1, j, k))) +
                        gamma * inverseDx * 0.25 *
                            (fabs(U(i, j, k) + U(i, j, k + 1)) *
                                    (W(i, j, k) - W(i + 1, j, k)) +
                                fabs(U(i - 1, j, k) + U(i - 1, j, k + 1)) *
                                    (W(i, j, k) - W(i - 1, j, k)));

                dvwdy = inverseDy * 0.25 *
                            ((V(i, j, k) + V(i, j, k + 1)) *
                                    (W(i, j, k) + W(i, j + 1, k)) -
                                (V(i, j - 1, k + 1) + V(i, j - 1, k)) *
                                    (W(i, j, k) + W(i, j - 1, k))) +
                        gamma * inverseDy * 0.25 *
                            (fabs(V(i, j, k) + V(i, j, k + 1)) *
                                    (W(i, j, k) - W(i, j + 1, k)) +
                                fabs(V(i, j - 1, k + 1) + V(i, j - 1, k)) *
                                    (W(i, j, k) - W(i, j - 1, k)));

                dw2dz = inverseDz * 0.25 *
                            ((W(i, j, k) + W(i, j, k + 1)) *
                                    (W(i, j, k) + W(i, j, k + 1)) -
                                (W(i, j, k) + W(i, j, k - 1)) *
                                    (W(i, j, k) + W(i, j, k - 1))) +
                        gamma * inverseDz * 0.25 *
                            (fabs(W(i, j, k) + W(i, j, k + 1)) *
                                    (W(i, j, k) - W(i, j, k + 1)) +
                                fabs(W(i, j, k) + W(i, j, k - 1)) *
                                    (W(i, j, k) - W(i, j, k - 1)));

                dw2dx2 = inverseDx * inverseDx *
                         (W(i + 1, j, k) - 2.0 * W(i, j, k) + W(i - 1, j, k));
                dw2dy2 = inverseDy * inverseDy *
                         (W(i, j + 1, k) - 2.0 * W(i, j, k) + W(i, j - 1, k));
                dw2dz2 = inverseDz * inverseDz *
                         (W(i, j, k + 1) - 2.0 * W(i, j, k) + W(i, j, k - 1));
                H(i, j, k) = W(i, j, k) + dt * (inverseRe * (dw2dx2 + dw2dy2 + dw2dz2) -
                                                   duwdx - dvwdy - dw2dz + gz);
            }
        }
    }

    /* ----------------------------- boundary of F ---------------------------
     */
    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            F(0, j, k)    = U(0, j, k);
            F(imax, j, k) = U(imax, j, k);
        }
    }

    /* ----------------------------- boundary of G ---------------------------
     */
    for (int k = 1; k < kmax + 1; k++) {
        for (int i = 1; i < imax + 1; i++) {
            G(i, 0, k)    = V(i, 0, k);
            G(i, jmax, k) = V(i, jmax, k);
        }
    }

    /* ----------------------------- boundary of G ---------------------------
     */
    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            H(i, j, 0)    = W(i, j, 0);
            H(i, j, kmax) = W(i, j, kmax);
        }
    }
}

void adaptUV(Solver* s)
{
    int imax = s->grid.imax;
    int jmax = s->grid.jmax;
    int kmax = s->grid.kmax;

    double* p = s->p;
    double* u = s->u;
    double* v = s->v;
    double* w = s->w;
    double* f = s->f;
    double* g = s->g;
    double* h = s->h;

    double factorX = s->dt / s->grid.dx;
    double factorY = s->dt / s->grid.dy;
    double factorZ = s->dt / s->grid.dz;

    for (int k = 1; k < kmax + 1; k++) {
        for (int j = 1; j < jmax + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                U(i, j, k) = F(i, j, k) - (P(i + 1, j, k) - P(i, j, k)) * factorX;
                V(i, j, k) = G(i, j, k) - (P(i, j + 1, k) - P(i, j, k)) * factorY;
                W(i, j, k) = H(i, j, k) - (P(i, j, k + 1) - P(i, j, k)) * factorZ;
            }
        }
    }
}
