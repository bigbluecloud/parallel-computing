#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

#define NUM_BODY 100
#define ITERATIONS 100
#define BODY_DATA_COLS 7

#define MASS 0
#define XPOS 1
#define YPOS 2
#define ZPOS 3
#define XVEL 4
#define YVEL 5
#define ZVEL 6

#define TIMESTEP 0.005
#define MAX_MASS 1000
#define SPACE_SIZE 1000
#define BODY_VEL_START 200

const double GRAV_CONST = 1;

void print_data(double bodyData[][BODY_DATA_COLS]) {
  int i, j;
  for(i=0;i<NUM_BODY;i++) {
    printf("[%d]\t", i+1);
    for(j=0;j<BODY_DATA_COLS;j++) {
      if(j == 0) printf("%-6.0f", bodyData[i][j]);
      else printf("%-16.4f", bodyData[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void init_bodies(double bodyData[][BODY_DATA_COLS], int rank, int worldSize) {
  int i, low, high, remaining, work, dataSize, nodes = worldSize - 1;
  srand((time(NULL) >> rank));

  work = NUM_BODY / (worldSize - 1);
  dataSize = BODY_DATA_COLS * work;
  if(rank == 0) {
    MPI_Status stat;
    remaining = NUM_BODY % worldSize;

    if(remaining != 0) {
      low = NUM_BODY - remaining;
      for(i=low;i<NUM_BODY;i++) {
        bodyData[i][MASS] = rand()%MAX_MASS + 100;
        bodyData[i][XPOS] = rand()%SPACE_SIZE;
        bodyData[i][YPOS] = rand()%SPACE_SIZE;
        bodyData[i][ZPOS] = rand()%SPACE_SIZE;
        bodyData[i][XVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2);
        bodyData[i][YVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2);
        bodyData[i][ZVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2);
      }
    }
    for(i=0;i<nodes;i++) {
      MPI_Recv(&low, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
      if(stat.MPI_SOURCE!=1) MPI_Recv(&bodyData[low], dataSize, MPI_DOUBLE, stat.MPI_SOURCE, 0, MPI_COMM_WORLD, &stat);
      else MPI_Recv(&bodyData[low], dataSize << 1, MPI_DOUBLE, stat.MPI_SOURCE, 0, MPI_COMM_WORLD, &stat);
    }
  }
  else {
    if(rank == 1) {
      low = 0;
      high = work << 1;
      dataSize = dataSize << 1;
    }
    else {
      low = work * rank;
      high = low + work;
    }
    for(i=low;i<high;i++) {
      bodyData[i][MASS] = rand()%MAX_MASS + 100;
      bodyData[i][XPOS] = rand()%SPACE_SIZE;
      bodyData[i][YPOS] = rand()%SPACE_SIZE;
      bodyData[i][ZPOS] = rand()%SPACE_SIZE;
      bodyData[i][XVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2);
      bodyData[i][YVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2);
      bodyData[i][ZVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2);
    }
    MPI_Send(&low, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&bodyData[low], dataSize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}

void run_simulation(double body_data[][BODY_DATA_COLS], int rank, int worldSize) {
  int i, step, newBody, nodes = worldSize - 1;

  for(step=0;step<ITERATIONS;step++) {
    MPI_Bcast(body_data, NUM_BODY * BODY_DATA_COLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank==0) {
      MPI_Status stat;
      for(i=0;i<NUM_BODY;i++) {
        MPI_Recv(&newBody, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
        MPI_Recv(&(body_data[newBody]), BODY_DATA_COLS, MPI_DOUBLE, stat.MPI_SOURCE, 0, MPI_COMM_WORLD, &stat);
      }
      printf("Iteration %d\n\n", step+1);
      print_data(body_data);
    } //End master node operations

    else {
      int workLoad = NUM_BODY % nodes;
      double fx, fy, fz, r;
      double data_copy[BODY_DATA_COLS];
      newBody = rank - 1;

      while(1) {
        if(newBody > NUM_BODY) //Exit condition
          break;
        data_copy[MASS] = body_data[newBody][MASS];
        data_copy[XPOS] = body_data[newBody][XPOS];
        data_copy[YPOS] = body_data[newBody][YPOS];
        data_copy[ZPOS] = body_data[newBody][ZPOS];
        data_copy[XVEL] = body_data[newBody][XVEL];
        data_copy[YVEL] = body_data[newBody][YVEL];
        data_copy[ZVEL] = body_data[newBody][ZVEL];

        for(i=0;i<NUM_BODY;i++) {
          if(i==newBody)
            continue;
          r = sqrt((pow(body_data[i][XPOS] - data_copy[XPOS], 2) + pow(body_data[i][YPOS] - data_copy[YPOS], 2) + pow(body_data[i][ZPOS] - data_copy[ZPOS], 2)));
          fx = ((GRAV_CONST * body_data[i][MASS] * data_copy[MASS]) / (pow(r, 2))) * ((data_copy[XPOS] - body_data[i][XPOS]) / r);
          fy = ((GRAV_CONST * body_data[i][MASS] * data_copy[MASS]) / (pow(r, 2))) * ((data_copy[YPOS] - body_data[i][YPOS]) / r);
          fz = ((GRAV_CONST * body_data[i][MASS] * data_copy[MASS]) / (pow(r, 2))) * ((data_copy[ZPOS] - body_data[i][ZPOS]) / r);
          data_copy[XVEL] += (fx * TIMESTEP / data_copy[MASS]);
          data_copy[YVEL] += (fy * TIMESTEP / data_copy[MASS]);
          data_copy[ZVEL] += (fz * TIMESTEP / data_copy[MASS]);
        } //End force routine for
        data_copy[XPOS] += data_copy[XVEL] * TIMESTEP;
        data_copy[YPOS] += data_copy[YVEL] * TIMESTEP;
        data_copy[ZPOS] += data_copy[ZVEL] * TIMESTEP;
        MPI_Send(&newBody, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&data_copy, BODY_DATA_COLS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        newBody += workLoad;
      }
    } //End slave node operations
  } //End ITERATION for
}

int main(int argc, char* argv[]) {
  int rank, size;
  double time;
  double body_data[NUM_BODY][BODY_DATA_COLS];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  
  if(rank == 0)
    time = MPI_Wtime();
  init_bodies(body_data, rank, size);
  if(rank == 0) { printf("Intial state\n"); print_data(body_data); }
  run_simulation(body_data, rank, size); //Master node rank is changed after running this for some reason
  if(rank <= 0) {
    rank = 0; //Hack method to restore the master node rank so it doesn't segfault when calling Finalize()
    time = MPI_Wtime() - time;
    printf("Master: Simulation finished\nExecuted in %f seconds\n", time);
  }
  MPI_Finalize();
  return 0;
}
