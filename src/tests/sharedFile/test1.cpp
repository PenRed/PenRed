
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
//
//    This file is part of PenRed: Parallel Engine for Radiation Energy Deposition.
//
//    PenRed is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    PenRed is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.
//
//    You should have received a copy of the GNU Affero General Public License
//    along with PenRed.  If not, see <https://www.gnu.org/licenses/>. 
//
//    contact emails:
//
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//


//Simple multiple read test 

#include <time.h>
#include <stdlib.h>
#include <thread>
#include "sharedFile.hh"

#define BUFFERSIZE 5000

char* readAndCheck(pen_sharedFile* sharedF, const unsigned ID, const char* data, const size_t toRead, size_t* pos, int* err){

  char* buffer = nullptr;
  size_t nread;
  buffer = (char*) malloc(sizeof(char)*toRead);
  *err = sharedF->read(ID,buffer,
		       sharedF->isBinary() ? toRead : toRead+1,
		       nread);

  if((*err) != SHARED_FILE_SUCCESS){
    free(buffer);
    return nullptr;
  }
  
  if(strncmp(buffer,&(data[(*pos)]),nread) != 0){
    free(buffer);
    (*err) = -1;
    return nullptr;
  }

  (*pos) += nread;
  
  (*err) = SHARED_FILE_SUCCESS;
  return buffer;
}

void readAndCheckAll(pen_sharedFile* sharedF, const unsigned ID, const char* data, const size_t toRead, size_t* pos, int* err){

  size_t remaining = toRead;
  (*err) = SHARED_FILE_SUCCESS;
  
  while(remaining > 0){
    size_t prevPos = (*pos);
    char* pbuffer = readAndCheck(sharedF,ID,data,remaining,pos,err);
    if(pbuffer != nullptr)
      free(pbuffer);
    
    if((*err) != SHARED_FILE_SUCCESS)
      return;
    if(prevPos == (*pos)){
      printf("remaining: %lu/%lu\n",remaining,toRead);
      printf("prevPos: %lu\n",prevPos);
      printf("pos: %lu\n",(*pos));
      (*err) = -2;
      return;
    }
    remaining -= (*pos) - prevPos;
  }
  
}

void readTwice(pen_sharedFile* sharedF, const unsigned ID, const char* data, const size_t toRead, const size_t toRead2, int* err){

  char* buffer1 = nullptr;
  char* buffer2 = nullptr;
  
  //Set file to origin
  (*err) = sharedF->seek(ID, 0, SEEK_SET);
  if((*err) != SHARED_FILE_SUCCESS)
    return;
  
  //Read first data chunk
  size_t pos = 0;
  readAndCheckAll(sharedF,ID,data,toRead,&pos,err);
  if((*err) != SHARED_FILE_SUCCESS)
    return;

  //Get current position
  fpos_t fpos1;
  (*err) = sharedF->getPos(ID,fpos1);
  if((*err) != SHARED_FILE_SUCCESS){
    free(buffer1);
    return;
  }

  //Read second data chunk
  size_t pos2 = pos;

  buffer1 = readAndCheck(sharedF,ID,data,toRead2,&pos2,err);
  if((*err) != SHARED_FILE_SUCCESS){
    if(buffer1 != nullptr)
      free(buffer1);
    return;
  }

  if(ID == 0 && !(sharedF->isBinary())){
    printf("Initial pos reader 0: %lu\n",toRead);
    printf("%s\n",buffer1);
  }
  
  //Return to previous position
  (*err) = sharedF->setPos(ID,fpos1);

  //Read and check using buffer1
  size_t pos3 = 0;
  buffer2 = readAndCheck(sharedF,ID,buffer1,toRead2,&pos3,err);
  if((*err) != SHARED_FILE_SUCCESS){
    if(buffer2 != nullptr)
      free(buffer2);
    free(buffer1);
    return;
  }

  if(ID == 0 && !(sharedF->isBinary())){
    printf("%s\n",buffer2);
  }

  //Free memory
  free(buffer1);
  free(buffer2);

}

unsigned test1(pen_sharedFile& sharedF, const unsigned nreaders, const char* data, const size_t dataSize){

  int nerr = 0;

  
  //Threads Vector
  std::vector<std::thread> threads;

  //Error vector
  std::vector<int> errors;
  errors.resize(nreaders,-1);

  //Positions vector
  std::vector<size_t> positions;
  positions.resize(nreaders,0);
  
  
  for(unsigned i = 0; i < nreaders; i++){
    threads.push_back(std::thread(readAndCheckAll,&sharedF, i, data, dataSize, &(positions[i]), &(errors[i])));
  }

  //Join threads
  for(unsigned i = 0; i < nreaders; i++){
    threads[i].join();
  }

  //Check errors
  for(unsigned i = 0; i < nreaders; i++){
    if(errors[i] != 0){
      nerr++;
      printf("Error on test 1 reader %u: %d\n",i,errors[i]);
    }
  }

  //Check number of read data
  for(unsigned i = 0; i < nreaders; i++){
    if(positions[i] != dataSize){
      nerr++;
      printf("Error on test 1 reader %u. Read data is not the expected.\n",i);
      printf("   Expected: %lu\n",dataSize);
      printf("       Read: %lu\n",positions[i]);
    }
  }
  
  //Clear threads
  threads.clear();
    
  
  return nerr;
}

unsigned test2(pen_sharedFile& sharedF, const unsigned nreaders, const char* data, const size_t dataSize){

  int nerr = 0;
  
  //Threads Vector
  std::vector<std::thread> threads;

  //Error vector
  std::vector<int> errors;
  errors.resize(nreaders,-1);

  //Position vectors
  std::vector<size_t> positions1;
  positions1.resize(nreaders,0);
  std::vector<size_t> positions2;
  positions2.resize(nreaders,0);

  for(unsigned i = 0; i < nreaders; i++){
    double toReadD = (double)dataSize*((double) rand() / (RAND_MAX));
    positions1[i] = (size_t)toReadD;    
  }
  for(unsigned i = 0; i < nreaders; i++){
    double toReadD = (double)(dataSize-positions1[i])*((double) rand() / (RAND_MAX));
    positions2[i] = (size_t)toReadD;    
  }
  
  for(unsigned i = 0; i < nreaders; i++){
    threads.push_back(std::thread(readTwice,&sharedF,i,data, positions1[i], positions2[i], &(errors[i])));
  }

  //Join threads
  for(unsigned i = 0; i < nreaders; i++){
    threads[i].join();
  }

  //Check errors
  for(unsigned i = 0; i < nreaders; i++){
    if(errors[i] != 0){
      nerr++;
      printf("Error on test 2 reader %u: %d\n",i,errors[i]);
    }
  }
  
  //Clear threads
  threads.clear();
    
  
  return nerr;
}

int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s n-readers\n",argv[0]);
    return 1;
  }

  int nthreads = atoi(argv[1]);

  if(nthreads > 10){
    printf("Warning: 10 is the maximum number of readers for this test.\n");
    nthreads = 10;
  }
  
  if(nthreads < 1){
    printf("Error: number of readers must be greater than 0.\n");
    return -1;
  }

  unsigned nreaders = (unsigned) nthreads;

  //Create dummy data
  size_t dataSize = BUFFERSIZE*nreaders;
  size_t bufferSize = dataSize + 1;
  size_t minLines = dataSize*0.01;
  char* data = (char*) malloc(bufferSize*sizeof(char));

  srand (time(NULL));
  for(size_t i = 0; i < dataSize; i++){
    data[i] = rand() % 94 + 32;
  }

  for(size_t i = 0; i < minLines; i++){
    size_t pos = dataSize*((double) rand() / (RAND_MAX));
    data[pos] = '\n';
  }

  //Set end of string character
  data[dataSize] = '\0';

  //Save dummy data in binary and ASCII format
  //******************************************
  
  // *** Binary
  FILE* fdummy = nullptr;
  fdummy = fopen("binaryT1.bin","wb");
  if(fdummy == nullptr){
    printf("Unable to write binary dummy data file.\n");
    return -1;
  }

  //Write binary data
  size_t nwrote = fwrite(data,sizeof(char),dataSize,fdummy);
  if((unsigned)nwrote != dataSize){
    printf("Unable to write all data.\n");
    printf(" Expected: %lu\n",dataSize);
    printf("   Wroten: %lu\n",nwrote);
    return -2;
  }

  fclose(fdummy);
  fdummy = nullptr;
  
  // *** ASCII
  fdummy = fopen("ASCIIT1.dat","w");
  if(fdummy == nullptr){
    printf("Unable to write binary dummy data file.\n");
    return -1;
  }
  int nwrote2 = 0;
  fprintf(fdummy,"%s%n",data,&nwrote2);
  if((unsigned) nwrote2 != dataSize){
    printf("Unable to write all data.\n");
    printf(" Expected: %lu\n",dataSize);
    printf("   Wroten: %d\n",nwrote2);
    return -3;
  }
  
  fclose(fdummy);

  // Read test
  //********************

  //Create shared files
  pen_sharedFile sharedASCII;
  pen_sharedFile sharedBinary;

  //Open files
  sharedASCII.open("ASCIIT1.dat",false);
  sharedBinary.open("binaryT1.bin",true);

  //Create readers
  for(unsigned i = 0; i < nreaders; i++){
    int errCR = -1;
    errCR = sharedASCII.createReader(i);
    if(errCR != SHARED_FILE_SUCCESS){
      printf("Error: Unable to create reader %u for ASCII shared file: %d",i,errCR);
      return -4;
    }
    errCR = sharedBinary.createReader(i);
    if(errCR != SHARED_FILE_SUCCESS){
      printf("Error: Unable to create reader %u for Binary shared file: %d",i,errCR);
      return -4;
    }    
  }
  
  int nerr = 0;
  int prevErr = 0;
  // TEST 1
  //*********

  printf("\n *** Complete read test ***\n\n");

  // *** ASCII
  printf("\n *** ASCII\n");
  
  nerr += test1(sharedASCII, nreaders, data, dataSize);
  
  if(nerr != prevErr)
    printf("\nErrors on ASCII complete read test.\n");
  else
    printf(" Done!\n");
  prevErr = nerr;

  // *** Binary
  printf("\n *** Binary\n");
  nerr += test1(sharedBinary, nreaders, data, dataSize);

  if(nerr != prevErr)
    printf("\nErrors on Binary complete read test.\n");
  else
    printf(" Done!\n");
  

  // TEST 2
  //*********

  printf("\n *** Read twice test ***\n\n");

  // *** ASCII
  printf("\n *** ASCII\n");

  prevErr = nerr;
  nerr += test2(sharedASCII, nreaders, data, dataSize);
  
  if(nerr != prevErr)
    printf("\nErrors on ASCII read twice test.\n");
  else
    printf(" Done!\n");

  // *** Binary
  printf("\n *** Binary\n");
  prevErr = nerr;
  nerr += test2(sharedBinary, nreaders, data, dataSize);

  if(nerr != prevErr)
    printf("\nErrors on Binary read twice test.\n");
  else
    printf(" Done!\n");
  
  prevErr = nerr;
  
  if(nerr != 0)
    printf("Errors on tests!\n");
  else
    printf("Tests passed!\n");

  return 0;
}
