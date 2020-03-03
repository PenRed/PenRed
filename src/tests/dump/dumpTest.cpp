
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

 
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>

#include "pen_dump.hh"


unsigned compare(const double* p1, const double* p2, const size_t n){
  unsigned diffs = 0;
  for(size_t i = 0; i < n; i++){
    if(p1[i] == 0.0 || p2[i] == 0.0){
      if(p1[i] != p2[i]){
	printf("%.22E =? %.22E\n", p1[i], p2[i]);
	diffs++;
      }
    }
    else{
      if(fabs(1.0-p1[i]/p2[i]) > 1.0e-14){
	printf("%.22E =? %.22E\n", p1[i], p2[i]);
	diffs++;
      }
    }
  }
  return diffs;
}

unsigned compare(const int* p1, const int* p2, const size_t n){
  unsigned diffs = 0;
  for(size_t i = 0; i < n; i++){
    if(p1[i] != p2[i])
      diffs++;
  }
  return diffs;
}

unsigned compare(const long int* p1, const long int* p2, const size_t n){
  unsigned diffs = 0;
  for(size_t i = 0; i < n; i++){
    if(p1[i] != p2[i])
      diffs++;
  }
  return diffs;
}

unsigned compare(const unsigned* p1, const unsigned* p2, const size_t n){
  unsigned diffs = 0;
  for(size_t i = 0; i < n; i++){
    if(p1[i] != p2[i])
      diffs++;
  }
  return diffs;
}

unsigned compare(const long unsigned* p1, const long unsigned* p2, const size_t n){
  unsigned diffs = 0;
  for(size_t i = 0; i < n; i++){
    if(p1[i] != p2[i])
      diffs++;
  }
  return diffs;
}

unsigned compare(const unsigned char* p1, const unsigned char* p2, const size_t n){
  unsigned diffs = 0;
  for(size_t i = 0; i < n; i++){
    if(p1[i] != p2[i])
      diffs++;
  }
  return diffs;
}

void fill(double* p, size_t n, double max){
  
  for(size_t i = 0; i < n; i++){
    p[i] = ((double)rand() / RAND_MAX) - 0.5e0;
    p[i] *= max;
  }  
}

void fill(int* p, size_t n, double max){
  
  for(size_t i = 0; i < n; i++){
    double random = (double)rand() / RAND_MAX;
    p[i] = int((random-0.5)*max);
  }  
}

void fill(long int* p, size_t n, double max){
  
  for(size_t i = 0; i < n; i++){
    double random = (double)rand() / RAND_MAX;
    p[i] = (long int)((random-0.5)*max);
  }  
}

void fill(unsigned* p, size_t n, double max){
  
  for(size_t i = 0; i < n; i++){
    double random = (double)rand() / RAND_MAX;
    p[i] = unsigned(random*max);
  }  
}

void fill(long unsigned* p, size_t n, double max){
  
  for(size_t i = 0; i < n; i++){
    double random = (double)rand() / RAND_MAX;
    p[i] = (long unsigned)(random*max);
  }  
}

void fill(unsigned char* p, size_t n, double /*max*/){
  
  for(size_t i = 0; i < n; i++){
    double random = (double)rand() / RAND_MAX;
    p[i] = random*(126-32) + 32;
  }  
}

template <class T>
int test_dumpRead(const char* TypeName){

  pen_dump dump1;
  pen_dump dump2;

  const size_t dim1 = 10000;
  const size_t dim2 = 246733;
  const size_t dim3 = 200;
  const size_t dim4 = 1;

  //Error variable
  int err;
  
  //Verbose mode
  unsigned verbose = 3;
  
  //Double pointers
  T  pdin1[dim1];
  T* pdin2 = (T*) malloc(sizeof(T)*dim2);
  T  pdin3[dim3];
  T  pdin4;

  T  pdout1[dim1];
  T* pdout2 = (T*) malloc(sizeof(T)*dim2);
  T  pdout3[dim3];
  T  pdout4;

  unsigned char* pdump = nullptr;
  size_t dimDump;
  size_t nread = 0;

  //Print header
  printf("\n\n **** %s read after dump\n\n",TypeName);
  
  /* initialize random seed: */
  srand (time(NULL));

  //Fill arrays
  fill(pdin1 ,dim1,1.0e6);
  fill(pdin2 ,dim2,1.0e7);
  fill(pdin3 ,dim3,1.0e8);
  fill(&pdin4,dim4,1.0e9);

  //Register input arrays at dump1
  err = dump1.toDump(pdin1 ,dim1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 1 to dump1: %d\n",err);
    return -1;
  }
  err = dump1.toDump(pdin2 ,dim2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 2 to dump1: %d\n",err);
    return -1;
  }
  err = dump1.toDump(pdin3 ,dim3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 3 to dump1: %d\n",err);
    return -1;
  }
  err = dump1.toDump(&pdin4,dim4);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 4 to dump1: %d\n",err);
    return -1;
  }

  //Dump dump1
  err = dump1.dump(pdump,dimDump,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping dump1: %d\n",err);
    return -2;
  }

  //Register input arrays at dump2
  err = dump2.toDump(pdout1 ,dim1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 1 to dump2: %d\n",err);
    return -3;
  }
  err = dump2.toDump(pdout2 ,dim2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 2 to dump2: %d\n",err);
    return -3;
  }
  err = dump2.toDump(pdout3 ,dim3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 3 to dump2: %d\n",err);
    return -3;
  }
  err = dump2.toDump(&pdout4,dim4);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 4 to dump2: %d\n",err);
    return -3;
  }

  //Read dumped array in dump2
  err = dump2.read(pdump,nread,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading dumped data: %d\n",err);
    return -5;
  }

  printf("bytes printed: %lu\n",dimDump);
  printf("bytes read : %lu\n",nread);
  
  if(dimDump != nread){
    printf("bytes dumpend doesn't match read");
    return -6;
  }
  
  //for(unsigned i = 0; i < dim3; i++)
  //  printf("%.22E =? %.22E\n", pdin3[i], pdout3[i]);
  
  printf("Differenses at array 1: %u\n",compare(pdin1 ,pdout1 ,dim1));
  printf("Differenses at array 2: %u\n",compare(pdin2 ,pdout2 ,dim2));
  printf("Differenses at array 3: %u\n",compare(pdin3 ,pdout3 ,dim3));
  printf("Differenses at array 4: %u\n",compare(&pdin4,&pdout4,dim4));

  //Free memory
  free(pdump);
  free(pdin2);
  free(pdout2);
  
  return 0;
}

template <class T>
int test_remove(const char* TypeName){

  pen_dump dump;
  pen_dump read;

  unsigned verbose = 3;
  
  T p1[100];
  T p2[100];
  T p3[100];
  T p4[100];

  T p5[100];
  T p6[100];
  T p7[100];
  
  unsigned char* pdump1 = nullptr;
  size_t dimDump1;
  unsigned char* pdump2 = nullptr;
  size_t dimDump2;
  size_t dimRead;
  
  int err;

  unsigned diffs;
  
  //Print header
  printf("\n\n **** %s remove\n\n",TypeName);
  
  //fill arrays
  fill(p1,100,50.0);
  fill(p2,100,50.0);
  fill(p3,100,50.0);
  fill(p4,100,50.0);
    

  //Register input arrays at dump
  err = dump.toDump(p1 ,100);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 1 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(p2 ,100);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 2 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(p3 ,100);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 3 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(p4,100);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 4 to dump: %d\n",err);
    return -1;
  }

  printf("Arrays registered to dump %lu\n",dump.nRegistered());

  //Register input arrays at dump
  err = read.toDump(p5 ,100);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 5 to read: %d\n",err);
    return -2;
  }
  err = read.toDump(p6 ,100);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 6 to read: %d\n",err);
    return -2;
  }
  err = read.toDump(p7 ,100);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering array 7 to read: %d\n",err);
    return -2;
  }

  printf("Arrays registered to read %lu\n",read.nRegistered());
  
  //Dump entire arrays
  dump.dump(pdump1,dimDump1,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping with 4 arrays: %d\n",err);
    return -3;
  }  
  printf("Dumped %lu bytes\n",dimDump1);
  
  printf("Remove array 3 from dump...\n");

  //Remove one array
  err = dump.remove(p3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error removing array 3 in dump: %d\n",err);
    return -4;
  }

  //Confirm that we can't remove the same array
  err = dump.remove(p3);
  if(err != PEN_DUMP_ELEMENT_NOT_FOUND){
    printf("Error, array 3 still registered at dump: %d\n",err);
    return -5;
  }  

  printf("Arrays registered to dump %lu\n",dump.nRegistered());
  
  //Dump remainding arrays
  err = dump.dump(pdump2,dimDump2,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping with 3 arrays: %d\n",err);
    return -6;
  }  
  printf("Dumped %lu bytes\n",dimDump2);

  //Read arrays
  printf("Read dump with 3 arrays.\n");
  dimRead = 0;
  err = read.read(pdump2,dimRead,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading dumped data: %d\n",err);
    return -7;
  }
  
  //Check if dump1 array 4 is equal to dump2 array 3
  diffs = compare(p4,p7,100);

  if(diffs > 0){
    printf("Last arrays are not equal, %u diferences found.\n", diffs);
    return -8;
  }
  else{
    printf("Last arrays of both dumps match!\n");
  }

  //Try to read dump with 4 arrays
  printf("Read dump with 4 arrays.\n");
  dimRead = 0;
  err = read.read(pdump1,dimRead,verbose);
  if(err == PEN_DUMP_SUCCESS){
    printf("Error, doesn't fail reading 4 arrays: %d\n",err);
    return -9;
  }

  printf("Test passed!\n");
  
  return 0;
}

int test_generalDump(){

  pen_dump dump;
  pen_dump read;

  int err;

  unsigned verbose = 3;
  
  //Create double arrays
  const size_t dimpd1 = 432;
  double pd1_1[dimpd1];
  double pd1_2[dimpd1];
  const size_t dimpd2 = 21;
  double pd2_1[dimpd2];
  double pd2_2[dimpd2];

  //Create integer arrays
  const size_t dimpi1 = 10432;
  int pi1_1[dimpi1];
  int pi1_2[dimpi1];
  const size_t dimpi2 = 56;
  int pi2_1[dimpi2];
  int pi2_2[dimpi2];
  const size_t dimpi3 = 765;
  int pi3_1[dimpi3];
  int pi3_2[dimpi3];

  //Create unsigned arrays
  const size_t dimpu1 = 123;
  unsigned pu1_1[dimpu1];
  unsigned pu1_2[dimpu1];
  const size_t dimpu2 = 98;
  unsigned pu2_1[dimpu2];
  unsigned pu2_2[dimpu2];

  //Create long unsigned arrays
  const size_t dimplu1 = 1232;
  long unsigned plu1_1[dimplu1];
  long unsigned plu1_2[dimplu1];
  const size_t dimplu2 = 346;
  long unsigned plu2_1[dimplu2];
  long unsigned plu2_2[dimplu2];
  
  //Create char arrays
  const size_t dimpuc1 = 645;
  unsigned char puc1_1[dimpuc1];
  unsigned char puc1_2[dimpuc1];
  const size_t dimpuc2 = 12;
  unsigned char puc2_1[dimpuc2];
  unsigned char puc2_2[dimpuc2];

  //Dump pointer
  unsigned char* pdump = nullptr;
  size_t dimDump = 0;
  size_t dimRead = 0;

  //Print header
  printf("\n\n **** Generic dump\n\n");
  
  //Fill arrays
  fill(pd1_1,dimpd1,12.32e30);
  fill(pd2_1,dimpd2,10.0);

  fill(pi1_1,dimpi1,12500000);
  fill(pi2_1,dimpi2,35);
  fill(pi3_1,dimpi3,35);

  fill(pu1_1,dimpu1,80);
  fill(pu2_1,dimpu2,120);

  fill(puc1_1,dimpuc1,0);
  fill(puc2_1,dimpuc2,0);

  //   Register dump arrays
  // ** Doubles
  err = dump.toDump(pd1_1,dimpd1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array 1 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(pd2_1,dimpd2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array 2 to dump: %d\n",err);
    return -1;
  }

  // ** Int
  err = dump.toDump(pi1_1,dimpi1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering integer array 1 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(pi2_1,dimpi2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering integer array 2 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(pi3_1,dimpi3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering integer array 3 to dump: %d\n",err);
    return -1;
  }

  // ** Unsigned
  err = dump.toDump(pu1_1,dimpu1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array 1 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(pu2_1,dimpu2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array 2 to dump: %d\n",err);
    return -1;
  }

  // ** Long Unsigned
  err = dump.toDump(plu1_1,dimplu1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering long unsigned integer array 1 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(plu2_1,dimplu2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering long unsigned integer array 2 to dump: %d\n",err);
    return -1;
  }

  // ** Chars
  err = dump.toDump(puc1_1,dimpuc1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned char array 1 to dump: %d\n",err);
    return -1;
  }
  err = dump.toDump(puc2_1,dimpuc2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned char array 2 to dump: %d\n",err);
    return -1;
  }

  //   Register read arrays
  // ** Doubles
  err = read.toDump(pd1_2,dimpd1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array 1 to read: %d\n",err);
    return -2;
  }
  err = read.toDump(pd2_2,dimpd2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array 2 to read: %d\n",err);
    return -2;
  }

  // ** Int
  err = read.toDump(pi1_2,dimpi1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering integer array 1 to read: %d\n",err);
    return -2;
  }
  err = read.toDump(pi2_2,dimpi2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering integer array 2 to read: %d\n",err);
    return -2;
  }
  err = read.toDump(pi3_2,dimpi3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering integer array 3 to read: %d\n",err);
    return -2;
  }

  // ** Unsigned
  err = read.toDump(pu1_2,dimpu1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array 1 to read: %d\n",err);
    return -2;
  }
  err = read.toDump(pu2_2,dimpu2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array 2 to read: %d\n",err);
    return -2;
  }

  // ** Long Unsigned
  err = read.toDump(plu1_2,dimplu1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering long unsigned integer array 1 to read: %d\n",err);
    return -2;
  }
  err = read.toDump(plu2_2,dimplu2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering long unsigned integer array 2 to read: %d\n",err);
    return -2;
  }
  
  // ** Chars
  err = read.toDump(puc1_2,dimpuc1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned char array 1 to read: %d\n",err);
    return -2;
  }
  err = read.toDump(puc2_2,dimpuc2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned char array 2 to read: %d\n",err);
    return -2;
  }

  //   Dump arrays
  err = dump.dump(pdump,dimDump,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping data: %d\n",err);
    return -3;
  }

  //   Read dumped arrays
  err = read.read(pdump,dimRead,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading data from dump array: %d\n",err);
    return -4;
  }

  //Check read arrays
  printf("Double arrays:\n");
  printf("Differenses at array 1: %u\n",compare(pd1_1, pd1_2,dimpd1));
  printf("Differenses at array 2: %u\n",compare(pd2_1, pd2_2,dimpd2));

  printf("Integer arrays:\n");
  printf("Differenses at array 1: %u\n",compare(pi1_1, pi1_2,dimpi1));
  printf("Differenses at array 2: %u\n",compare(pi2_1, pi2_2,dimpi2));
  printf("Differenses at array 3: %u\n",compare(pi3_1, pi3_2,dimpi3));
  
  printf("Unsigned integer arrays:\n");
  printf("Differenses at array 1: %u\n",compare(pu1_1, pu1_2,dimpu1));
  printf("Differenses at array 2: %u\n",compare(pu2_1, pu2_2,dimpu2));

  printf("Long Unsigned integer arrays:\n");
  printf("Differenses at array 1: %u\n",compare(plu1_1, plu1_2,dimplu1));
  printf("Differenses at array 2: %u\n",compare(plu2_1, plu2_2,dimplu2));
  
  printf("Unsigned char arrays:\n");
  printf("Differenses at array 1: %u\n",compare(puc1_1, puc1_2,dimpuc1));
  printf("Differenses at array 2: %u\n",compare(puc2_1, puc2_2,dimpuc2));

  return 0;
}

int test_splitDump(){

  pen_dump dump1,dump2,dump3;
  pen_dump read1,read2,read3;

  // double arrays
  const size_t dimpd1_1 = 571;
  const size_t dimpd1_2 = 27;
  const size_t dimpd1_3 = 83;
  double pd1_1[dimpd1_1];
  double pd1_2[dimpd1_2];
  double pd1_3[dimpd1_3];

  double pd1_1_read[dimpd1_1];
  double pd1_2_read[dimpd1_2];
  double pd1_3_read[dimpd1_3];
  
  // Unsigned arrays
  const size_t dimpu1_1 = 234;
  const size_t dimpu1_2 = 1000;
  const size_t dimpu1_3 = 5;
  unsigned pu1_1[dimpu1_1];
  unsigned pu1_2[dimpu1_2];
  unsigned pu1_3[dimpu1_3];  

  unsigned pu1_1_read[dimpu1_1];
  unsigned pu1_2_read[dimpu1_2];
  unsigned pu1_3_read[dimpu1_3];  
  
  unsigned char* pdump;
  size_t dumpSize;
  size_t readSize;
  unsigned char* pdump1;
  size_t dump1Size;
  unsigned char* pdump2;
  size_t dump2Size;
  unsigned char* pdump3;
  size_t dump3Size;

  unsigned verbose = 3;
  int err;
  
  //Print header
  printf("\n\n **** Splitted dump\n\n");
  
  // Fill arrays
  fill(pd1_1,dimpd1_1,1.5e1);
  fill(pd1_2,dimpd1_2,1.5e5);
  fill(pd1_3,dimpd1_3,1.5e15);

  fill(pu1_1,dimpu1_1,2.5e1);
  fill(pu1_2,dimpu1_2,2.5e5);
  fill(pu1_3,dimpu1_3,2.5e7);

  printf("Arrays filled\n");
  
  // Register arrays
  // ** dump1
  err = dump1.toDump(pd1_1,dimpd1_1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array to dump 1: %d\n",err);
    return -1;
  }
  err = dump1.toDump(pu1_1,dimpu1_1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array to dump 1: %d\n",err);
    return -1;
  }  

  // ** dump2
  err = dump2.toDump(pd1_2,dimpd1_2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array to dump 2: %d\n",err);
    return -1;
  }
  err = dump2.toDump(pu1_2,dimpu1_2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array to dump 2: %d\n",err);
    return -1;
  }  

  // ** dump3
  err = dump3.toDump(pd1_3,dimpd1_3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array to dump 3: %d\n",err);
    return -1;
  }
  err = dump3.toDump(pu1_3,dimpu1_3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array to dump 3: %d\n",err);
    return -1;
  }  

  // ** read1
  err = read1.toDump(pd1_1_read,dimpd1_1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array to read 1: %d\n",err);
    return -1;
  }
  err = read1.toDump(pu1_1_read,dimpu1_1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array to read 1: %d\n",err);
    return -1;
  }  

  // ** read2
  err = read2.toDump(pd1_2_read,dimpd1_2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array to read 2: %d\n",err);
    return -1;
  }
  err = read2.toDump(pu1_2_read,dimpu1_2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned integer array to read 2: %d\n",err);
    return -1;
  }

  // ** read3
  err = read3.toDump(pd1_3_read,dimpd1_3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering double array to read 3: %d\n",err);
    return -1;
  }
  err = read3.toDump(pu1_3_read,dimpu1_3);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering unsigned array to read 3: %d\n",err);
    return -1;
  }  

  // Dump partial dumps
  printf("Dump partial dumps.\n");
  err = dump1.dump(pdump1,dump1Size,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping from dump 1: %d.\n",err);
    return -2;
  }
  err = dump2.dump(pdump2,dump2Size,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping from dump 2: %d.\n",err);
    return -2;
  }
  err = dump3.dump(pdump3,dump3Size,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping from dump 3: %d.\n",err);
    return -2;
  }

  printf("Dumped. Join partial dumps\n");
  
  //Store to global dump
  dumpSize = dump1Size + dump2Size + dump3Size;
  pdump = (unsigned char*) malloc(dumpSize*sizeof(unsigned char));

  memcpy(pdump,pdump1,dump1Size);
  free(pdump1);
  memcpy(&pdump[dump1Size],pdump2,dump2Size);
  free(pdump2);
  memcpy(&pdump[dump1Size+dump2Size],pdump3,dump3Size);
  free(pdump3);

  pdump1 = pdump2 = pdump3 = nullptr;

  printf("Joined.\n");
  
  //Read from global dump
  readSize = 0;
  err = read1.read(pdump,readSize,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading dumped data to read1: %d\n",err);
    return -3;
  }
  err = read2.read(pdump,readSize,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading dumped data to read1: %d\n",err);
    return -3;
  }
  err = read3.read(pdump,readSize,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading dumped data to read1: %d\n",err);
    return -3;
  }

  if(readSize != dumpSize){
    printf("Error: Read and dumped bytes doesn't match.\n");
    return -4;
  }

  printf("Dumped and read data size match: %lu bytes\n",readSize);

  //Check if arrays match

  printf("\nDump 1:\n");
  printf("Differences between dumped and read double array: %u\n",compare(pd1_1,pd1_1_read,dimpd1_1));
  printf("Differences between dumped and read unsigned integer array: %u\n",compare(pu1_1,pu1_1_read,dimpu1_1));  

  printf("\nDump 2:\n");
  printf("Differences between dumped and read double array: %u\n",compare(pd1_2,pd1_2_read,dimpd1_2));
  printf("Differences between dumped and read unsigned integer array: %u\n",compare(pu1_2,pu1_2_read,dimpu1_2));  

  printf("\nDump 3:\n");
  printf("Differences between dumped and read double array: %u\n",compare(pd1_3,pd1_3_read,dimpd1_3));
  printf("Differences between dumped and read unsigned integer array: %u\n",compare(pu1_3,pu1_3_read,dimpu1_3));  
  
  
  return 0;
}

template <class T>
int test_redimension(const char* TypeName){

  pen_dump dump;
  pen_dump read1, read2, read3;  

  //Array to dump
  const size_t dim1 = 100;
  const size_t dim2 = 90;  
  T pin[dim1];

  T pout1[dim1];
  T pout2[dim1];
  T pout3[dim2];

  unsigned char* pdump1;
  size_t dimDump1;
  unsigned char* pdump2;
  size_t dimDump2;

  size_t dimRead1 = 0;
  size_t dimRead2 = 0;
  size_t dimRead3 = 0;
  
  int err;
  unsigned verbose = 3;
  
  //Print header
  printf("\n\n **** %s resize array to dump\n\n",TypeName);

  //Fill the array
  fill(pin,dim1,2.3e4);

  printf("Input array filled\n");
  
  //Register array to dump
  err = dump.toDump(pin,dim1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering input array: %d\n",err);
    return -2;
  }

  //Register read array
  err = read1.toDump(pout1,dim1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering output array to read 1: %d\n",err);
    return -3;
  }
  err = read2.toDump(pout2,dim1);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering output array to read 2: %d\n",err);
    return -3;
  }
  err = read3.toDump(pout3,dim2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registering output array to read 3: %d\n",err);
    return -3;
  }

  printf("Input and output arrays registered\n");
  
  //Dump entire array
  err = dump.dump(pdump1,dimDump1,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping entire array: %d.\n",err);
    return -4;
  }
  
  printf("Dumped entire array\n");
  
  //Resize dumping array
  err = dump.toDump(pin,dim2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registrering dump array with less dimension: %d\n",err);
    return -5;
  }
  //Resize read 2 array
  err = read2.toDump(pout2,dim2);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error registrering read2 array with less dimension: %d\n",err);
    return -5;
  }

  printf("Dump and read2 array resized\n");

  //Dump cutted array
  err = dump.dump(pdump2,dimDump2,0,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error dumping cutted array: %d\n",err);
    return -6;
  }

  printf("Dumped cutted array\n");
  
  //Read entire array
  err = read1.read(pdump1,dimRead1,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading read 1: %d\n",err);
    return -7;
  }
  printf("Read1 done\n");
  
  //Read cutted array
  printf("Try to read the complete dumped array on read2. This call should fail.\n");
  err = read2.read(pdump1,dimRead2,verbose);
  if(err == PEN_DUMP_SUCCESS){
    printf("Error reading entire array in read2 doesn't fail\n");
    return -7;
  }

  printf("Failed! Now, read the cutted one.\n");
  
  dimRead2 = 0;
  err = read2.read(pdump2,dimRead2,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading cutted array to read2: %d\n",err);
    return -7;
  }
  printf("Read2 done\n");
  
  err = read3.read(pdump2,dimRead3,verbose);
  if(err != PEN_DUMP_SUCCESS){
    printf("Error reading small array to read3: %d\n",err);
    return -7;
  }
  printf("Read3 done\n");

  //Check dimensions
  if(dimRead1 != dimDump1){
    printf("Dumped and read dimensions of entire array doesn't match:\n");
    printf("  Dumped bytes: %lu\n",dimDump1);
    printf("    Read bytes: %lu\n",dimRead1);
    return -8;
  }

  printf("Dumped and read bytes for entire array match!\n");
  
  if(dimRead2 != dimDump2){
    printf("Dumped and read dimensions of cutted array doesn't match:\n");
    printf("  Dumped bytes: %lu\n",dimDump2);
    printf("    Read bytes: %lu\n",dimRead2);
    return -8;
  }

  printf("Dumped and read bytes for cutted array match!\n");

  if(dimRead3 != dimDump2){
    printf("Dumped and read dimensions of cutted array doesn't match:\n");
    printf("  Dumped bytes: %lu\n",dimDump2);
    printf("    Read bytes: %lu\n",dimRead3);
    return -8;
  }

  printf("Dumped and read bytes for small array match!\n");

  //Compare arrays
  printf("Differences for entire array:\n");
  printf("Differences after read: %u\n",compare(pin,pout1,dim1));

  printf("Differences for cutted array:\n");
  printf("Differences with entire array (%lu expected): %u\n",dim1-dim2,compare(pin,pout2,dim1));
  printf("Differences with cutted array: %u\n",compare(pin,pout2,dim2));

  printf("Differences for small array:\n");
  printf("Differences after read: %u\n",compare(pin,pout3,dim2));
  
  return 0;
}

int main(){

  //Error variable
  int err;
  unsigned testsFailed = 0;

  // Check double read after dump
  //*******************************

  err = test_dumpRead<double>("Double");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Double: read after dump\n");
    testsFailed++;
  }

  err = test_dumpRead<int>("Int");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Int: read after dump\n");
    testsFailed++;
  }

  err = test_dumpRead<long int>("Long Int");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Long Int: read after dump\n");
    testsFailed++;
  }
  
  err = test_dumpRead<unsigned>("Unsigned");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Unsigned: read after dump\n");
    testsFailed++;
  }

  err = test_dumpRead<long unsigned>("Long Unsigned");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Long Unsigned: read after dump\n");
    testsFailed++;
  }
  
  err = test_dumpRead<unsigned char>("Unsigned Char");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Unsigned Char: read after dump\n");
    testsFailed++;
  }

  // Test remove
  //*****************
  err = test_remove<double>("Double");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Double: 'remove'\n");
    testsFailed++;
  }

  err = test_remove<int>("Int");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Int: 'remove'\n");
    testsFailed++;
  }

  err = test_remove<long int>("Long Int");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Long Int: 'remove'\n");
    testsFailed++;
  }
  
  err = test_remove<unsigned>("Unsigned");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Unsigned: 'remove'\n");
    testsFailed++;
  }

  err = test_remove<long unsigned>("Long Unsigned");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Long Unsigned: 'remove'\n");
    testsFailed++;
  }
  
  err = test_remove<unsigned char>("Unsigned Char");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: Unsigned Char: 'remove'\n");
    testsFailed++;
  }

  // Test generic dump
  //*********************
  err = test_generalDump();
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: general dump\n");
    testsFailed++;
  }

  // Test splitted dump
  //*********************
  err = test_splitDump();
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: splitted dump\n");
    testsFailed++;
  }

  // Test redimension dump
  //**************************
  err = test_redimension<double>("Double");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: redimension registered array\n");
    testsFailed++;
  }

  err = test_redimension<int>("Integer");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: redimension registered array\n");
    testsFailed++;
  }
  
  err = test_redimension<long int>("Long Integer");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: redimension registered array\n");
    testsFailed++;
  }
  
  err = test_redimension<unsigned>("Unsigned");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: redimension registered array\n");
    testsFailed++;
  }

  err = test_redimension<long unsigned>("Long Unsigned");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: redimension registered array\n");
    testsFailed++;
  }
  
  err = test_redimension<unsigned char>("Unsigned Char");
  if(err != PEN_DUMP_SUCCESS){
    printf("Error on test: redimension registered array\n");
    testsFailed++;
  }


  printf("\n\n");
  printf("*********************\n");
  printf("  Tests failed: %u\n",testsFailed);
  printf("*********************\n");
  
  return 0;
}
