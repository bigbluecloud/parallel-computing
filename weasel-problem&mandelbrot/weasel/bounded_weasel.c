#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

#define MAX_NEGATIVE -999999
#define POOLSIZE 1000
#define VALIDLEN strlen(validChars)
#define TARGETLEN strlen(targetStr)

#define MUTATION_RATE 5

#define FOUND_TAG 0
#define MASTER_SEND_TAG 1
#define MASTER_RECV_TAG 2
#define MASTER 0

const char* validChars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ ";
const char* targetStr = "WEASEL";

int min(int a, int b) {
  return a < b ? a : b;
}


int findChar(char key, const char* characterSet) { //Finds the index of the character in the valid chars list
  int i, len = strlen(characterSet);
  for(i=0;i<len;++i)
    if(characterSet[i] == key)
      break;
  return i;
}


int getFitness(const char* candidate) { //Returns the fitness of the passed string with regard to the target string
  int i, targetPos, candidatePos, currentFitness = 0;
  for(i=0;i<TARGETLEN;++i) {
    targetPos = findChar(targetStr[i], validChars);
    candidatePos = findChar(candidate[i], validChars);
    int difference = abs(targetPos - candidatePos);
    currentFitness -= min(difference, (int)(VALIDLEN) - difference);
  }
  return currentFitness;
}


void generateString(char* string) {
  int j;
  memset(string, 0, TARGETLEN + 1);
  for(j=0;j<TARGETLEN;++j)
    string[j] = validChars[rand()%VALIDLEN];
}


int selection(int* fitnesses) { //Roulette selection implementation
  int i, currentIndex = 0, selectionFitness = MAX_NEGATIVE;
  for(i=0;i<POOLSIZE;++i) {
    int random = 0.9 * fitnesses[i] - rand()%100;
    if(random >= selectionFitness) { //Greater than or equal to to increase the chances an unfit string will be sent
      currentIndex = i;
      selectionFitness = random;
    }
  }
  return currentIndex;
}


int findWorstFitIndex(int* poolFitness) { //Finds the index of the least fit string in the pool
  int i, index, fitness = 0;
  for(i=0;i<POOLSIZE;++i) {
    if(poolFitness[i] < fitness) {
      fitness = poolFitness[i];
      index = i;
    }
  }
  return index;
}

int findBestFitIndex(int* poolFitness) { //Finds the index of the fittest string in the pool
  int i, index, fitness = MAX_NEGATIVE;
  for(i=0;i<POOLSIZE;++i) {
    if(poolFitness[i] > fitness) {
      fitness = poolFitness[i];
      index = i;
    }
  }
  return index;
}


int similarity(char* str1, char* str2) {
  int i, similarity = TARGETLEN;
  for(i=0;i<TARGETLEN+1;++i)
    if(str1[i] != str2[i])
      --similarity;
  return similarity;
}


void masterLogic(int size) {
  int i, nodes = size - 1, iteration = 1, worstFitIndex = 0, worstFitness = 0, bestFitIndex = 0, bestFitness = MAX_NEGATIVE, found = 0, prevBest = MAX_NEGATIVE;
  char pool[POOLSIZE][TARGETLEN + 1];
  char recvPool[POOLSIZE][TARGETLEN+1];
  int poolFitness[POOLSIZE] = {MAX_NEGATIVE};
  MPI_Status stat;
  FILE* bestStringsFile = fopen("strings.txt", "w");

  for(i=0;i<POOLSIZE;++i)
    memset(recvPool[i], 0, TARGETLEN + 1);

  for(i=0;i<POOLSIZE;++i) { //Initial generation of the pool
    generateString(pool[i]);
    poolFitness[i] = getFitness(pool[i]);

    if(poolFitness[i] < worstFitness) {
      worstFitness = poolFitness[i];
      worstFitIndex = i;
    }
    else if(poolFitness[i] > bestFitness) {
      bestFitness = poolFitness[i];
      bestFitIndex = i;
    }
  } //End initial pool generation for loop

  if(bestFitness==0) //Optimistically the solution is randomly generated
    found = 1;

  while(1) {
    MPI_Bcast(&found, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(found==1) break;

    for(i=0;i<nodes*2;++i)
      MPI_Recv(&recvPool[i], TARGETLEN+1, MPI_CHAR, MPI_ANY_SOURCE, MASTER_RECV_TAG, MPI_COMM_WORLD, &stat); //Recv every string pair from the nodes

    for(i=1;i<size;++i) {
      MPI_Send(&pool[selection(poolFitness)], TARGETLEN+1, MPI_CHAR, i, MASTER_SEND_TAG, MPI_COMM_WORLD); //Send the node a new pair of strings to operate on with bias selection
      MPI_Send(&pool[selection(poolFitness)], TARGETLEN+1, MPI_CHAR, i, MASTER_SEND_TAG, MPI_COMM_WORLD);
    }

    for(i=0;i<POOLSIZE-1;++i) {
      int strFit;
      if(recvPool[i][0] == '\0') //The end of the received strings is reached, exit the loop
        break;
      strFit = getFitness(recvPool[i]);

      //If the current string being looked at is fitter or passes a diversity check
      if(strFit > worstFitness || (strFit + 5 > worstFitness && similarity(pool[worstFitIndex], recvPool[i]) < TARGETLEN/2+1)) { 
        strncpy(pool[worstFitIndex], recvPool[i], TARGETLEN); //Add the string into the pool
        poolFitness[worstFitIndex] = strFit; //Change the string fitness to match

        worstFitIndex = findWorstFitIndex(poolFitness); //Find the next least fit string
        worstFitness = poolFitness[worstFitIndex];
      }
    }
    bestFitIndex = findBestFitIndex(poolFitness); //Re-find the best fit string
    bestFitness = poolFitness[bestFitIndex];

    if(bestFitness == 0) {
      printf("Target string found in %d iterations!\n", iteration);
      printf("Target string: %s\n", pool[bestFitIndex]);
      fprintf(bestStringsFile, "[%d]\t%s\n", iteration, pool[bestFitIndex]);
      found = 1;
    }
    else if(prevBest != bestFitness) {
      prevBest = bestFitness;
      printf("[%d]\t%s %d\n", iteration, pool[bestFitIndex], bestFitness);
      fprintf(bestStringsFile, "[%d]\t%s\n", iteration, pool[bestFitIndex]);
    }
    ++iteration;
  }
  fclose(bestStringsFile);
}


void crossOverStrings(char* str1, char* str2) {
  char temp;
  int i, split = (TARGETLEN/2) + (rand()%(TARGETLEN/4)+1); //Split in the middle of the string plus or minus a quarter of the entire string length

  for(i=split;i<TARGETLEN+1;++i) {
    temp = str1[i];
    str1[i] = str2[i];
    str2[i] = temp;
  }
}


void mutateString(char* str) {
  int i;
  for(i=0;i<TARGETLEN+1;++i)
    if(rand()%100 < MUTATION_RATE) //Check if the character will mutate
      str[i] = validChars[rand()%VALIDLEN]; //Mutate the character to another random character in validChars
}


void nodeLogic(int rank) {
  int i, mode = (rank % 2 == 0 ? 0 : 1), found = 0;
  char strings [2][TARGETLEN + 1];

  memset(strings[0], 0, TARGETLEN + 1);
  memset(strings[1], 0, TARGETLEN + 1);

  generateString(strings[0]); //Generate two random strings
  generateString(strings[1]);
  
  while(1) {
    MPI_Bcast(&found, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(found==1) break;

    for(i=0;i<2;++i)
      MPI_Send(&strings[i], TARGETLEN+1, MPI_CHAR, MASTER, MASTER_RECV_TAG, MPI_COMM_WORLD); //Send string pair to master node
    for(i=0;i<2;++i)
      MPI_Recv(&strings[i], TARGETLEN+1, MPI_CHAR, MASTER, MASTER_SEND_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Recv new string pair
    
    switch(mode) {
      case 0:
        crossOverStrings(strings[0], strings[1]); //Single split crossover
        mode = 1;
        break;
      case 1:
        mutateString(strings[0]); //Random mutate
        mutateString(strings[1]);
        mode = 0;
        break;
    }
  }
}


int main(int argc, char** argv) {
  int rank, worldSz;
  char randStr[32] = {0};

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &worldSz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  srand(time(NULL) >> rank);
  if(rank==MASTER)
    masterLogic(worldSz);
  else
    nodeLogic(rank);
  
  MPI_Finalize();
  return 0;
}
