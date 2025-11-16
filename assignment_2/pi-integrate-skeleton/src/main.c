/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

#include "timing.h"

double integrate(double, double, int);

int main (int argc, char** argv) {

    //for serial while executing the code just run it on 1 core ==> 1 rank so serial

    double wcs, wce;
    double  Pi;

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double partition = 1.0/size; //TODO change to command line arg or not verify

    double a = (double)(rank)/size;
    double b = (double)(rank+1)/size;
    //printf("This is of the rank %d and size is %d and the a %f and b %f and the partitions are %f are \n ", rank,size,a,b,partition);
    wcs = getTimeStamp();
    Pi = integrate(a, b, size);
    wce = getTimeStamp();
    //here we are recieving n-1 ranks pi data and add them to the rank 0(Master) and the at the multiplying with 4(for four quadrants)
    if(rank == 0 ){
        double recieve = Pi;
        for(int i=1;i<size;i++){
            
            MPI_Recv(&recieve,1,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //printf("I have recived %f \n ", recieve);
            Pi+=recieve;
        }
    }else{
        //a simple mpi_send , where all the rank send data to master(rank 0)
        MPI_Send(&Pi,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        //printf("i am % d I have sent %f Pi \n",rank, Pi);
    }
    //as the values are accumulated in the rank 0 we will print the data only if the rank is of the master
    Pi *= 4;
    if(rank == 0) printf("Pi=%.15lf in %.3lf ms \n", Pi,(wce-wcs)*1e3);
    MPI_Finalize();

    return EXIT_SUCCESS;
}

double integrate(double a, double b, int size) {
	
	/*
	
		Your logic to integrate between given interval a to b.
		Declare SLICES here and calculate delta x using a, b and SLICES.
		Iterate over number of slices, calculate the area and sum them.
		Return sum * delta x to get the value of PI.
		
	*/
    //the mid point rule and 10^6 works for the 8 point accuracy
    int iter = 1e6/size;
    double delta_x = (b-a)/iter;
    //printf("the number of iter : %d \n", iter);
    double sum=0;
    // a Midpoint method is better but normal also works
    for(int i=0;i<iter;i++){
        double x = (a+((i+0.5)*delta_x));
        double f_x = sqrt(1-x*x);
        double area = f_x * delta_x;
       sum+=area;
    }
    // sum *= 4;
	return sum;
}
