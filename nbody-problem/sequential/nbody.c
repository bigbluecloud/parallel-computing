#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h> //For high resolution timer (UNIX only)
#include <time.h> //For getting a seed for srand()
#include <math.h> //For sqrt() and pow()

#define NUM_BODY 100
#define ITERATIONS 100
#define BODY_DATA_COLS 7	//Number of columns in the body_data array

#define MASS 0	                //Macros for readability when addressing the array
#define XPOS 1
#define YPOS 2
#define ZPOS 3
#define XVEL 4
#define YVEL 5
#define ZVEL 6

#define TIMESTEP 0.005          //Determines accuracy of position updates
#define MAX_MASS 1000           //Maximum mass of body when randomly generating
#define SPACE_SIZE 1000         //Size of the space to randomly place the bodies in
#define BODY_VEL_START 200	//Maximum for body starting velocity

const double GRAV_CONST = 1;    //Simulation gravitational constant

void init_bodies(double bodyData[][BODY_DATA_COLS]) { //Sets the data for each individual body in the initial state
  int i;
  srand(time(NULL)); //Seed random
  for(i=0;i<NUM_BODY;i++) {
    bodyData[i][MASS] = rand()%MAX_MASS + 100; //Bodies may have mass in the range of 100 to 1100
    bodyData[i][XPOS] = rand()%SPACE_SIZE;
    bodyData[i][YPOS] = rand()%SPACE_SIZE;
    bodyData[i][ZPOS] = rand()%SPACE_SIZE;
    bodyData[i][XVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2); //Range of -100 to 100
    bodyData[i][YVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2);
    bodyData[i][ZVEL] = (rand()%BODY_VEL_START + 1) - (BODY_VEL_START / 2);
  }
}

void print_data(double bodyData[][BODY_DATA_COLS]) { //Prints out the data for each body
  int i, j;
  for(i=0;i<NUM_BODY;i++) {
    printf("[%d]\t", i+1);
    for(j=0;j<BODY_DATA_COLS;j++) {
	  if(j == 0)
	    printf("%-6.0f", bodyData[i][j]);
	  else
        printf("%-16.4f", bodyData[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void run_simulation(double body_data[][BODY_DATA_COLS]) {
  double data_copy[NUM_BODY][3] = {0}; //[0] is the x velocity, [1] is the y, [2] is z. Other data is not copied out of the body_data array
  double fx, fy, fz, r; //Force xyz, distance
  int i, j, step; //i is the current body, j is other body being calculated against

  for(step=0;step<ITERATIONS;step++) { //For each iteration
	for(i=0;i<NUM_BODY;i++) { //For every body
	  data_copy[i][0] = body_data[i][XVEL];
	  data_copy[i][1] = body_data[i][YVEL];
	  data_copy[i][2] = body_data[i][ZVEL];
      for(j=0;j<NUM_BODY;j++) { //Calculate the force on it from every other body
        if(j == i)	//If the body is being calculated against is itself
          continue;
        r = sqrt((pow(body_data[i][XPOS] - body_data[j][XPOS], 2) + pow(body_data[i][YPOS] - body_data[j][YPOS], 2) + pow(body_data[i][ZPOS] - body_data[j][ZPOS], 2))); //Calculate distance between the bodies i and j
        fx = ((GRAV_CONST * body_data[i][MASS] * body_data[j][MASS]) / (pow(r, 2))) * ((body_data[j][XPOS] - body_data[i][XPOS]) / r); //Calculate forces
        fy = ((GRAV_CONST * body_data[i][MASS] * body_data[j][MASS]) / (pow(r, 2))) * ((body_data[j][YPOS] - body_data[i][YPOS]) / r);
	    fz = ((GRAV_CONST * body_data[i][MASS] * body_data[j][MASS]) / (pow(r, 2))) * ((body_data[j][ZPOS] - body_data[i][ZPOS]) / r);
	    data_copy[i][0] += (fx * TIMESTEP / body_data[i][MASS]); //New velocities
	    data_copy[i][1] += (fy * TIMESTEP / body_data[i][MASS]);
	    data_copy[i][2] += (fz * TIMESTEP / body_data[i][MASS]);
      }
      body_data[i][XPOS] += data_copy[i][0] * TIMESTEP; //Update the original array
      body_data[i][YPOS] += data_copy[i][1] * TIMESTEP;
      body_data[i][ZPOS] += data_copy[i][2] * TIMESTEP;
      body_data[i][XVEL] = data_copy[i][0];
      body_data[i][YVEL] = data_copy[i][1];
      body_data[i][ZVEL] = data_copy[i][2];
    }
    printf("Iteration %d\n\n", step+1);
    print_data(body_data); //Print data array
  }
}

int main(int argc, char* argv[]) {
  double body_data[NUM_BODY][BODY_DATA_COLS];
  struct timespec start, end, elapsed;
  time_t seconds;
  long milliseconds;

  clock_gettime(CLOCK_MONOTONIC, &start);
  init_bodies(body_data);
  printf("Initial state\n");
  print_data(body_data);
  run_simulation(body_data);
  clock_gettime(CLOCK_MONOTONIC, &end);

  if((end.tv_nsec - start.tv_nsec) < 0) {
    seconds = end.tv_sec - start.tv_sec - 1;
    milliseconds = 1000000000 + end.tv_nsec - start.tv_nsec;
  }
  else {
    seconds = end.tv_sec - start.tv_sec;
    milliseconds = end.tv_nsec - start.tv_nsec;
  }
  printf("\nSimulation finished\nExecuted in %ld.%ld seconds\n", (long)seconds, (milliseconds / 1.0e6));
  return 0;
}
