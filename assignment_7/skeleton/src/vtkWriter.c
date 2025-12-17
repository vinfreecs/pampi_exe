/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "allocate.h"
#include "comm.h"
#include "mpi.h"
#include "mpio.h"
#include "vtkWriter.h"

#define G(v, i, j, k) v[(k) * imax * jmax + (j) * imax + (i)]

// reset fileview for output of string headers
static void resetFileview(VtkOptions *o) {
  MPI_Offset disp;
  MPI_File_sync(o->handle);
  MPI_Barrier(o->comm.comm);
  MPI_File_get_size(o->handle, &disp);
  MPI_File_set_view(o->handle, disp, MPI_CHAR, MPI_CHAR, "native",
                    MPI_INFO_NULL);
}

static double floatSwap(double f) {
  union {
    double f;
    char b[8];
  } dat1, dat2;

  dat1.f = f;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.f;
}

static void writeHeader(VtkOptions *o) {

  // Only rank 0 should write this header to the file.
  // Use MPI IO along with rank check to write this header to the vtk file.

  if (o->comm.rank != 0) {
    return;
  }

  fprintf(o->fh, "# vtk DataFile Version 3.0\n");
  fprintf(o->fh, "PAMPI cfd solver output\n");
  if (o->fmt == ASCII) {
    fprintf(o->fh, "ASCII\n");
  } else if (o->fmt == BINARY) {
    fprintf(o->fh, "BINARY\n");
  }

  fprintf(o->fh, "DATASET STRUCTURED_POINTS\n");
  fprintf(o->fh, "DIMENSIONS %d %d %d\n", o->grid.imax, o->grid.jmax,
          o->grid.kmax);
  fprintf(o->fh, "ORIGIN %f %f %f\n", o->grid.dx * 0.5, o->grid.dy * 0.5,
          o->grid.dz * 0.5);
  fprintf(o->fh, "SPACING %f %f %f\n", o->grid.dx, o->grid.dy, o->grid.dz);
  fprintf(o->fh, "POINT_DATA %d\n", o->grid.imax * o->grid.jmax * o->grid.kmax);
}

void vtkOpen(VtkOptions *o, char *problem) {
  char filename[50];
  snprintf(filename, 50, "/lustre/pavl/pavl166v/%s.vtk", problem);

  MPI_Info info;
  MPI_Info_create(&info);
  // Hint: stripe over 10 I/O devices
  MPI_Info_set(info, "striping_factor", "10");
  MPI_File_open(o->comm.comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, info,
                &o->handle);
  MPI_Info_free(&info);
  // Rank 0 writes header using MPI_File_write
  if (o->comm.rank == 0) {
    char header[512];
    int len =
        snprintf(header, 512,
                 "# vtk DataFile Version 3.0\n"
                 "PAMPI cfd solver output\n"
                 "%s\n"
                 "DATASET STRUCTURED_POINTS\n"
                 "DIMENSIONS %d %d %d\n"
                 "ORIGIN %f %f %f\n"
                 "SPACING %f %f %f\n"
                 "POINT_DATA %d\n",
                 (o->fmt == BINARY) ? "BINARY" : "ASCII", o->grid.imax,
                 o->grid.jmax, o->grid.kmax, o->grid.dx * 0.5, o->grid.dy * 0.5,
                 o->grid.dz * 0.5, o->grid.dx, o->grid.dy, o->grid.dz,
                 o->grid.imax * o->grid.jmax * o->grid.kmax);

    MPI_File_write(o->handle, header, len, MPI_CHAR, MPI_STATUS_IGNORE);
  }

  MPI_Barrier(o->comm.comm);

  if (commIsMaster(&o->comm)) {
    printf("Writing VTK output for %s\n", problem);
  }
}

// TODO verfiy if this also should be changed
static void writeScalar(VtkOptions *o, double *s) {

  int imax = o->grid.imax;
  int jmax = o->grid.jmax;
  int kmax = o->grid.kmax;

  for (int k = 0; k < kmax; k++) {
    for (int j = 0; j < jmax; j++) {
      for (int i = 0; i < imax; i++) {
        if (o->fmt == ASCII) {
          fprintf(o->fh, "%f\n", G(s, i, j, k));
        } else if (o->fmt == BINARY) {
          fwrite((double[1]){floatSwap(G(s, i, j, k))}, sizeof(double), 1,
                 o->fh);
        }
      }
    }
  }
  if (o->fmt == BINARY)
    fprintf(o->fh, "\n");
}

static bool isInitialized(FILE *ptr) {
  if (ptr == NULL) {
    printf("vtkWriter not initialize! Call vtkOpen first!\n");
    return false;
  }
  return true;
}

void vtkScalar(VtkOptions *o, char *name, double *s) {
  // Write MPI IO code to parallelise writing to a single file.
  // Steps to perform MPI IO
  // 1. Use resetFileview(o); here before starting MPI IO
  resetFileview(o);
  if (o->comm.rank == 0) {
    printf("Register scalar %s\n", name);
  }
  // if (!isInitialized(o->fh))
  //   return;

  // Make sure this is only written by rank 0.
  if (o->comm.rank == 0) {
    char header[256];
    int len = snprintf(header, 256,
                       "SCALARS %s double 1\nLOOKUP_TABLE default\n", name);
    MPI_File_write(o->handle, header, len, MPI_CHAR, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(o->comm.comm);

  // 2. Get offsets from all ranks using commGetOffsets from comm.c. This will
  // be

  int offsets[NDIMS];
  commGetOffsets(&o->comm, offsets, o->grid.kmax, o->grid.jmax, o->grid.imax);

  // 3. Use the MPI_File_sync and MPI_File_get_size (along with MPI_Offset)
  // to get file size.
  MPI_File_sync(o->handle);
  // MPI_Barrier(o->comm.comm);
  MPI_Offset file_size;
  MPI_File_get_size(o->handle, &file_size);
  // 4. Use offsets gathered in 3 to create subarray. Create subarray from
  // global problem (Imax, Jmax, Kmax)
  MPI_Datatype fileViewType;
  MPI_Type_create_subarray(
      NDIMS, (int[NDIMS]){o->grid.kmax, o->grid.jmax, o->grid.imax},
      (int[NDIMS]){o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal},
      offsets, MPI_ORDER_C, MPI_DOUBLE, &fileViewType);
  MPI_Type_commit(&fileViewType);
  // 5. Finally use the MPI_File_set_view to set the view in file for each
  // rank.
  MPI_File_set_view(o->handle, file_size, MPI_DOUBLE, fileViewType,
                    "external32", MPI_INFO_NULL);
  // 6. Instead of using writeScalar(), create another subarray excluding
  // the ghost cells (imaxlocal+2 etc.)
  MPI_Datatype bulkType;
  MPI_Type_create_subarray(
      NDIMS,
      (int[NDIMS]){o->comm.kmaxLocal + 2, o->comm.jmaxLocal + 2,
                   o->comm.imaxLocal + 2},
      (int[NDIMS]){o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal},
      (int[NDIMS]){1, 1, 1}, MPI_ORDER_C, MPI_DOUBLE, &bulkType);
  MPI_Type_commit(&bulkType);
// 7. Now use the MPI_File_write using the subarray created in 6 so that
// each rank can write to the same file at the given displacements.
#define GG(v, i, j, k)                                                         \
  v[(k) * (o->comm.kmaxLocal + 2) * (o->comm.jmaxLocal + 2) +                  \
    (j) * (o->comm.imaxLocal + 2) + (i)]
  size_t cnt = o->comm.imaxLocal * o->comm.jmaxLocal * o->comm.kmaxLocal;
  double *tmp = allocate(64, cnt * NDIMS * sizeof(double));
  int idx = 0;
  for (int k = 1; k < o->comm.kmaxLocal + 1; k++) {
    for (int j = 1; j < o->comm.jmaxLocal + 1; j++) {
      for (int i = 1; i < o->comm.imaxLocal + 1; i++) {
        tmp[idx++] = GG(s, i, j, k);
      }
    }
  }
  MPI_File_write(o->handle, tmp, 1, bulkType, MPI_STATUS_IGNORE);
  // writeScalar(o, s);

  MPI_Type_free(&fileViewType);
  MPI_Type_free(&bulkType);

  // Binary segment must be terminated with newline character
  resetFileview(o);
  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->handle, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
  }
}

static void writeVector(VtkOptions *o, VtkVector vec) {
  int imax = o->grid.imax;
  int jmax = o->grid.jmax;
  int kmax = o->grid.kmax;

  for (int k = 0; k < kmax; k++) {
    for (int j = 0; j < jmax; j++) {
      for (int i = 0; i < imax; i++) {
        if (o->fmt == ASCII) {
          fprintf(o->fh, "%f %f %f\n", G(vec.u, i, j, k), G(vec.v, i, j, k),
                  G(vec.w, i, j, k));
        } else if (o->fmt == BINARY) {
          fwrite((double[3]){floatSwap(G(vec.u, i, j, k)),
                             floatSwap(G(vec.v, i, j, k)),
                             floatSwap(G(vec.w, i, j, k))},
                 sizeof(double), 3, o->fh);
        }
      }
    }
  }
  if (o->fmt == BINARY)
    fprintf(o->fh, "\n");
}

void vtkVector(VtkOptions *o, char *name, VtkVector vec) {
  // Write MPI IO code to parallelise writing to a single file.
  // Steps to perform MPI IO
  // 1. Use resetFileview(o); here before starting MPI IO
  resetFileview(o);

  if (o->comm.rank == 0) {
    printf("Register vector %s\n", name);
  }
  // if (!isInitialized(o->fh))
  //   return;

  // Make sure this is only written by rank 0.
  if (o->comm.rank == 0) {
    char header[256];
    int len = snprintf(header, 256, "VECTORS %s double\n", name);
    MPI_File_write(o->handle, header, len, MPI_CHAR, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(o->comm.comm);

  // 2. Get offsets from all ranks using commGetOffsets from comm.c. This will
  // be
  int offsets[NDIMS];
  commGetOffsets(&o->comm, offsets, o->grid.kmax, o->grid.jmax, o->grid.imax);

  // 3. Use the MPI_File_sync and MPI_File_get_size (along with MPI_Offset) to
  // get file size.
  MPI_Offset disp;
  MPI_Datatype fileViewType, vectorType;
  MPI_File_get_size(o->handle, &disp);
  // 4. Use offsets gathered in 3 to create subarray. Create subarray from
  // global problem (Imax, Jmax, Kmax) Note : Feel free to use a temporary array
  // to gather the x, y and z velocities in contiguous fashion to write vectors.
  MPI_Type_contiguous(NDIMS, MPI_DOUBLE, &vectorType);
  MPI_Type_commit(&vectorType);
  MPI_Type_create_subarray(
      NDIMS, (int[NDIMS]){o->grid.kmax, o->grid.jmax, o->grid.imax},
      (int[NDIMS]){o->comm.kmaxLocal, o->comm.jmaxLocal, o->comm.imaxLocal},
      offsets, MPI_ORDER_C, vectorType, &fileViewType);
  MPI_Type_commit(&fileViewType);
  // 5. Finally use the MPI_File_set_view to set the view in file for each rank.
  MPI_File_set_view(o->handle, disp, MPI_DOUBLE, fileViewType, "external32",
                    MPI_INFO_NULL);
// 6. Instead of using writeVector(), create another subarray excluding the
// ghost cells (imaxlocal+2 etc.)
#define GGG(v, i, j, k)                                                        \
  v[(k) * (o->comm.kmaxLocal + 2) * (o->comm.jmaxLocal + 2) +                  \
    (j) * (o->comm.imaxLocal + 2) + (i)]
  size_t cnt = o->comm.imaxLocal * o->comm.jmaxLocal * o->comm.kmaxLocal;
  double *tmp = allocate(64, cnt * NDIMS * sizeof(double));
  int idx = 0;
  for (int k = 1; k < o->comm.kmaxLocal + 1; k++) {
    for (int j = 1; j < o->comm.jmaxLocal + 1; j++) {
      for (int i = 1; i < o->comm.imaxLocal + 1; i++) {
        tmp[idx++] = (GGG(vec.u, i, j, k) + GGG(vec.u, i - 1, j, k)) / 2.0;
        tmp[idx++] = (GGG(vec.v, i, j, k) + GGG(vec.v, i, j - 1, k)) / 2.0;
        tmp[idx++] = (GGG(vec.w, i, j, k) + GGG(vec.w, i, j, k - 1)) / 2.0;
      }
    }
  }
  // 7. Now use the MPI_File_write using the subarray created in 6 so that each
  // rank can write to the same file at the given displacements.
  MPI_File_write(o->handle, tmp, cnt, vectorType, MPI_STATUS_IGNORE);
  MPI_Type_free(&fileViewType);
  MPI_Type_free(&vectorType);

  // writeVector(o, vec);

  // Binary segment must be terminated with newline character
  resetFileview(o);
  if (commIsMaster(&o->comm)) {
    MPI_File_write(o->handle, "\n", 1, MPI_CHAR, MPI_STATUS_IGNORE);
  }
}

void vtkClose(VtkOptions *o) {
  // fclose(o->fh);
  // o->fh = NULL;
  MPI_File_close(&o->handle);
}
