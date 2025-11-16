#include "timing.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "allocate.h"
#include <string.h>

// Helper function to calculate size of rank
int size_of_rank(int rank, int size, int N)
{
    int num = N / size;
    int rest = N % size;
    return (rank < rest) ? (num + 1) : num;
}

double dmvm_non_blocking(double* restrict y,
    const double* restrict a,
    double* restrict x,
    int N,
    int iter,
    int rank,
    int size)
{
    // new buffers to sen and recv and while that happens do the dmvm and wait and repeat
    MPI_Request send_req, recv_req;
    size_t bytesPerWord = sizeof(double);
    double *x_send = (double*)allocate(ARRAY_ALIGNMENT, N * bytesPerWord);
    x_send = x;
    double *x_recv = (double*)allocate(ARRAY_ALIGNMENT, N * bytesPerWord);
    double ts, te;
    MPI_Status status;
    int Nnext;
    
    // Calculate neighbors
    int upperNeighbor = (rank - 1) % size;
    if (upperNeighbor < 0) upperNeighbor = size - 1;
    int lowerNeighbor = (rank + 1) % size;
    
    // Calculate starting column for this rank
    int num = N / size;
    int rest = N % size;
    int cs = rank * num + ((rank < rest) ? rank : rest);
    int Nlocal = (N/size) + ((N%size>rank)?1:0);
    int Ncurrent = Nlocal;
    
    ts = getTimeStamp();
    
    for (int j = 0; j < iter; j++) {
        // Reset for each iteration
        cs = rank * num + ((rank < rest) ? rank : rest);
        Ncurrent = Nlocal;
        
        // Loop over RHS ring shifts
        for (int rot = 0; rot < size; rot++) {
            // Compute local contribution
            //MPI_Isend(x_send, Ncurrent, MPI_DOUBLE, lowerNeighbor, 0, MPI_COMM_WORLD,)
            //MPI_Irecv(x_recv)

            if (rot < size - 1) {
                // Calculate next segment size
                int nextSourceRank = (rank + rot + 1) % size;
                Nnext = size_of_rank(nextSourceRank, size, N);
                // int Nnext = sizeOfRank_block_non(lowerNeighbor, size, N);
                // Non-blocking send current segment up
                MPI_Isend(x_send, Ncurrent, MPI_DOUBLE, 
                         upperNeighbor, rot, MPI_COMM_WORLD, &send_req);
                
                // Non-blocking receive next segment from below
                MPI_Irecv(x_recv, Nnext, MPI_DOUBLE,
                         lowerNeighbor, rot, MPI_COMM_WORLD, &recv_req);
            }
            
            for (int r = 0; r < Nlocal; r++) {
                for (int c = cs; c < cs + Ncurrent; c++) {
                    y[r] += a[r * N + c] * x[c - cs];
                }
            }

            if (rot < size - 1) {
                MPI_Wait(&send_req, &status);
                MPI_Wait(&recv_req, &status);
                
                // Swap buffers (received data becomes send data for next round)
                double *temp = x_send;
                x_send = x_recv;
                x_recv = temp;
                
                // Update for next rotation
                cs += Ncurrent;
                if (cs >= N) cs = 0;
                Ncurrent = Nnext;
                memcpy(x, x_recv, N * sizeof(double));
            }

            //a barrier here to sync the sends and recvieves
            //TODO: replace the x with recieved x_recv
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
    // this is not necessary but makes its more consistent
    double walltime = te - ts;
    double max_walltime;
    MPI_Reduce(&walltime, &max_walltime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        return max_walltime;
    } else {
        return walltime;
    }
}
