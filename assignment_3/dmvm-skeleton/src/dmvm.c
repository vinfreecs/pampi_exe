/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include "allocate.h"
#include "timing.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int sizeOfRank(int rank, int size, int N) {

  int num = N / size;
  int rest = N % size;
  return (rank < rest) ? (num + 1) : num;
}

double dmvm(double *restrict y, const double *restrict a, double *restrict x,
            int N, int iter, int rank, int size) {
  double ts, te;
  size_t bytesPerWord = sizeof(double);

  MPI_Status status;
  double *x_recv = (double *)allocate(ARRAY_ALIGNMENT, N * bytesPerWord);

  int upperNeighbor = (rank - 1) % size;
  if (upperNeighbor < 0)
    upperNeighbor = size - 1;
  int lowerNeighbor = (rank + 1) % size;

  int num = N / size;
  int rest = N % size;
  int cs = rank * num + ((rank < rest) ? rank : rest);
  int Nlocal = (N / size) + ((N % size > rank) ? 1 : 0);
  int Ncurrent = Nlocal;

  ts = getTimeStamp();

  for (int j = 0; j < iter; j++) {
    // Reset for each iteration
    cs = rank * num + ((rank < rest) ? rank : rest);
    Ncurrent = Nlocal;

    for (int rot = 0; rot < size; rot++) {
      for (int r = 0; r < Nlocal; r++) {
        for (int c = cs; c < cs + Ncurrent; c++) {
          y[r] += a[r * N + c] * x[c - cs];
        }
      }

      cs += Ncurrent;
      if (cs >= N)
        cs = 0; // wrap around

      // Get size of next chunk from left neighbor
      int nextSourceRank = (rank + rot + 1) % size;
      Ncurrent = sizeOfRank(nextSourceRank, size, N);

      // Ncurrent = sizeOfRank(lowerNeighbor, size, N);

      if (rot != size - 1) {
        // sendrecvreplace is bugghy
        MPI_Sendrecv(x, N, MPI_DOUBLE, upperNeighbor, 0, x_recv, N, MPI_DOUBLE,
                     lowerNeighbor, 0, MPI_COMM_WORLD, &status);
      }
      memcpy(x, x_recv, N * sizeof(double));
    }

#ifdef CHECK
    {
      double sum = 0.0;
      for (int i = 0; i < Nlocal; i++) {
        sum += y[i];
        y[i] = 0.0;
      }
      double global_sum;
      MPI_Reduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (rank == 0) {
        fprintf(stderr, "Sum: %f\n", global_sum);
      }
    }
#endif
  }

  te = getTimeStamp();

  // Return maximum time across all ranks
  double walltime = te - ts;
  double max_walltime;
  MPI_Reduce(&walltime, &max_walltime, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    return max_walltime;
  } else {
    return walltime;
  }
}
