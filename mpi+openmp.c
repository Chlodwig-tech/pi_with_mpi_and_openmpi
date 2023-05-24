#include "utility.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>

#define THREADNUM 8
#define RESULT 1
double precision = 1000000000;
double step;

void calculate (int inputArgument, int threadsNumber)
{
  int value;
  #pragma omp parallel reduction(+:value) 
  {
    int r = 10000;
    for (int i = 0; i < inputArgument; i++){
      long x = rand() % r + 1;
      long y = rand() % r + 1;
      if(x * x + y * y <= r * r){
        value++;
      }
    }
  }
  double PI = 0.0;
  PI = (float) value / (inputArgument * threadsNumber) * 4;
  MPI_Send (&PI, 1, MPI_DOUBLE, 0, RESULT, MPI_COMM_WORLD);
}

int main(int argc,char **argv) {
  double pi_final = 0;
  int threadsupport;
  MPI_Status status;
  Args ins__args;
  parseArgs(&ins__args, &argc, argv);
  int threadsNumber = ins__args.n_thr;

  //set number of threads
  omp_set_num_threads(threadsNumber);
  
  //program input argument
  long inputArgument = ins__args.arg; 

  struct timeval ins__tstart, ins__tstop;

  int myrank,nproc;

  // Initialize MPI with desired support for multithreading -- state your desired support level

  MPI_Init_thread(&argc, &argv,MPI_THREAD_FUNNELED,&threadsupport); 

  if (threadsupport<MPI_THREAD_FUNNELED) {
    printf("\nThe implementation does not support MPI_THREAD_FUNNELED, it supports level %d\n",threadsupport);
    MPI_Finalize();
    return -1;
  }
  
  // obtain my rank
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  // and the number of processes
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  if (!myrank) 
    gettimeofday(&ins__tstart, NULL);
  // run your computations here (including MPI communication and OpenMP stuff)

  calculate(inputArgument / nproc, threadsNumber);
  if(myrank == 0){
    double resulttemp;
	   for (int i = 0; i < nproc; i++)
      {
        MPI_Recv (&resulttemp, 1, MPI_DOUBLE, i, RESULT, MPI_COMM_WORLD, &status);
        printf ("\nReceived result %f for process %d\n",resulttemp, i);
        fflush (stdout);
        pi_final += resulttemp;
	   }

    pi_final /= nproc;
    printf ("\npi=%f\n", pi_final);
  }

  // synchronize/finalize your computations

  if (!myrank) {
    gettimeofday(&ins__tstop, NULL);
    ins__printtime(&ins__tstart, &ins__tstop, ins__args.marker);
  }
    
  MPI_Finalize();
  
}
