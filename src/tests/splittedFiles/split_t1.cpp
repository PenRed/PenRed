
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


#include <stdlib.h>
#include <thread>
#include "splittedFile.hh"

#define BUFFERSIZE 500

void generateASCII(pen_splittedFile* splitted, unsigned ID, unsigned nIter){

  splitted->createPartition(ID);
  
  char c = ID + '0';

  for(unsigned iter = 0; iter < nIter; iter++){
    for(unsigned i = 0; i < ID+1; i++){
      char buffer[2*BUFFERSIZE+1];

      //Fill buffer
      for(unsigned j = 0; j < BUFFERSIZE; j++){
	buffer[2*j  ] = c;
	buffer[2*j+1] = '\n';
      }
      buffer[2*BUFFERSIZE] = '\0';
      
      //Write to partition
      splitted->write(ID,buffer,2*BUFFERSIZE);
    }
  }
  splitted->flush(ID);  
}

void generateBinary(pen_splittedFile* splitted, unsigned ID, unsigned nIter){

  splitted->createPartition(ID);

  char c = ID + '0';
  char buffer[BUFFERSIZE];
  //Fill buffer
  memset(buffer,c,BUFFERSIZE);

  for(unsigned iter = 0; iter < nIter; iter++){
    for(unsigned i = 0; i < ID+1; i++){
      splitted->write(ID,buffer,BUFFERSIZE);
    }
  }
  splitted->flush(ID);
}

unsigned checkPartitions(const char* prefix, const bool binary, const unsigned npartitions, unsigned nIter){

  //Test splitted file partitions
  
  unsigned testsFailed = 0;

  //getc(stdin);
  
  char format[3] = {'r','\0','\0'};
  if(binary) format[1] = 'b';
  
  for(unsigned i = 0; i < npartitions; i++){

    printf(" Testing partition %u...",i);
    
    char filename[300];
    sprintf(filename,"%s-%u",prefix,i);
    
    FILE* fin = nullptr;
    fin = fopen(filename,format);
    if(fin == nullptr){
      printf("Error: Unable to open partition file '%s'\n",filename);
      testsFailed++;
      continue;
    }

    //Read and count lines/elements
    unsigned elements = 0;
    if(binary){
      char c;
      char cobj = i + '0';
      while(fread(&c,sizeof(char),1,fin) == 1){
	if(c != cobj){
	  printf("\nError on partition %u:\n",i);
	  printf("   Number expected: %c\n",cobj);
	  printf("       Number read: %c\n\n",c);
	  testsFailed++;
	  break;	  
	}
	elements++;
      }
    }
    else{
      int num;
      while(fscanf(fin," %d ",&num) == 1){
	if(num != (int)i){
	  printf("\n\nError on partition %u:\n",i);
	  printf("   Number expected: %u\n",i);
	  printf("       Number read: %u\n\n",num);
	  testsFailed++;
	  break;
	}
	elements++;
      }
    }

    if(elements != (unsigned)BUFFERSIZE*(i+1)*nIter){
      printf("\n\nError on partition %u:\n",i);
      printf("  Expected elements: %u\n",BUFFERSIZE*(i+1)*nIter);
      printf("      Read elements: %u\n\n",elements);
      testsFailed++;
    }
    else{
      printf(" %u/%u\n",elements,BUFFERSIZE*(i+1)*nIter);
    }

    //Close file
    fclose(fin);    
  }
  
  
  return testsFailed;
}

unsigned checkMerged(const char* prefix, const bool binary, const unsigned npartitions, unsigned nIter){

  //Test partitions with specified prefix
  
  unsigned testsFailed = 0;

  char format[3] = {'r','\0','\0'};
  if(binary) format[1] = 'b';

  char filename[300];
  sprintf(filename,"%s-merged.dat",prefix);
  
  FILE* fmerged = nullptr;
  fmerged = fopen(filename,format);
  if(fmerged == nullptr){
    printf("Error: Unable to open ASCII merged file ('%s').\n",filename);
    testsFailed++;
  }
  else{

    //Get expected number of elements
    unsigned expectedElements = 0;
    for(unsigned ipart = 0; ipart < npartitions; ipart++){
      expectedElements += BUFFERSIZE*(ipart+1);
    }
    expectedElements *= nIter;
    
    unsigned elements = 0;

    //Iterate over partitions
    for(unsigned ipart = 0; ipart < npartitions; ipart++){

      //Create buffers for binary files
      char cobj = ipart + '0';
      char bufferComp[BUFFERSIZE];
      char bufferRead[BUFFERSIZE];
      memset(bufferComp,cobj,BUFFERSIZE);
      
      //Iterate over chunks
      for(unsigned iter = 0; iter < nIter; iter++){
	for(unsigned i = 0; i < ipart+1; i++){
	
	  if(binary){
	    // Binary files	  

	    //Read data chunck
	    size_t nread = fread(bufferRead,sizeof(char),BUFFERSIZE,fmerged);
	    elements += nread;
	    if(nread != BUFFERSIZE){
	      printf("Error: Unable to read chunk %u.\n",i);
	      printf("         elements expected: %u\n",expectedElements);
	      printf("         elements read  : %u\n",elements);
	      testsFailed++;
	      return testsFailed;	    
	    }

	    //Compare with expected data
	    if(memcmp(bufferComp,bufferRead,BUFFERSIZE) != 0){
	      printf("Error: Data in partition %u, chunk %u is not the expected one.\n",ipart,i);
	      printf("     Expected: %s\n",bufferComp);
	      printf("     Read  : %s\n",bufferRead);
	      testsFailed++;
	      return testsFailed;	    
	    }
	  
	  }
	  else{
	    // ASCII files

	    //Iterate over chunk data
	    for(unsigned j = 0; j < BUFFERSIZE; j++){

	      int num;
	      if(fscanf(fmerged," %d ",&num) != 1){
		printf("Error: Unable to read line %u.\n",elements+1);
		printf("         lines expected: %u\n",expectedElements);
		testsFailed++;
		return testsFailed;
	      }else{
		if(num != (int)ipart){
		  printf("Error: Unexpected value on line %u.\n",elements+1);
		  printf("         Expected: %d\n",ipart);
		  printf("         Read  : %d\n",num);
		  testsFailed++;
		  return testsFailed;
		}
		elements++;
	      }	  
	    }
	  }
	}
      }
    }

    if(elements != expectedElements){
      printf("Error: Expected and read elements mismatch.\n");
      printf("                  Expected: %u\n",expectedElements);
      printf("                      Read: %u\n",elements);
      testsFailed++;
    }
    fclose(fmerged);
  }
  
  
  return testsFailed;
}

void createData(pen_splittedFile& splittedASCII, pen_splittedFile& splittedBinary, unsigned npartitions, unsigned nIter){
  
  //Threads Vector
  std::vector<std::thread> threads;
  
  printf(" **** Data creation \n\n");

  // ASCII
  printf(" Create ASCII data\n");
  //Creat partitions with independent threads
  for(unsigned i = 0; i < npartitions; i++){
    threads.push_back(std::thread(generateASCII,&splittedASCII,i,nIter));
  }

  //Join threads
  for(unsigned i = 0; i < npartitions; i++){
    threads[i].join();
  }

  //Clear threads
  threads.clear();

  // Binary
  printf(" Create binary data\n");
  //Creat partitions with independent threads
  for(unsigned i = 0; i < npartitions; i++){
    threads.push_back(std::thread(generateBinary,&splittedBinary,i,nIter));
  }

  //Join threads
  for(unsigned i = 0; i < npartitions; i++){
    threads[i].join();
  }

  //Clear threads
  threads.clear();  
}

unsigned checkWrotePartitions(const char* prefixASCII, const char* prefixBinary, unsigned npartitions, unsigned nIter){

  unsigned testsFailed = 0;
  unsigned prevErr = 0;
  
  printf("\n **** Partitions check\n");

  // ASCII
  printf("\n - Splitted ASCII:\n");

  testsFailed += checkPartitions(prefixASCII,false,npartitions,nIter);

  if(prevErr == testsFailed){
    printf("\n\n Partitions check success!\n");
  }
  else{
    printf("\n\n Partitions check failed\n");
  }

  // Binary
  printf("\n - Splitted binary:\n");
  prevErr = testsFailed;

  testsFailed += checkPartitions(prefixBinary,true,npartitions,nIter);

  if(prevErr == testsFailed){
    printf("\n\n Partitions check success!\n");
  }
  else{
    printf("\n\n Partitions check failed\n");
  }

  return testsFailed;
}

unsigned checkWroteMerges(const char* prefixASCII, const char* prefixBinary, unsigned npartitions, unsigned nIter){

  unsigned testsFailed = 0;
  unsigned prevErr = 0;
  
  printf("\n **** Merge check\n");

  // ASCII
  printf("\n - Splitted ASCII:\n");
  prevErr = testsFailed;
  
  testsFailed += checkMerged(prefixASCII, false, npartitions,nIter);
  
  if(prevErr == testsFailed){
    printf("\n\n Merge check success!\n");
  }
  else{
    printf("\n\n Merge check failed\n");
  }

  // Binary
  printf("\n - Splitted binary:\n");
  prevErr = testsFailed;
  
  testsFailed += checkMerged(prefixBinary, true, npartitions,nIter);

  
  if(prevErr == testsFailed){
    printf("\n\n Merge check success!\n");
  }
  else{
    printf("\n\n Merge check failed\n");
  }  

  return testsFailed;
}

unsigned test1(pen_splittedFile& splittedASCII, const char* prefixASCII, pen_splittedFile& splittedBinary, const char* prefixBinary, unsigned npartitions){

  unsigned testsFailed = 0;
  
  //*********************
  //* Test merge at end *
  //*********************

  printf("\n----------------------------\n");
  printf("      Test merge at end");
  printf("\n----------------------------\n");
  
  // Create data files
  //********************

  createData(splittedASCII, splittedBinary, npartitions, 1);
  
  // Check partitions
  //********************

  testsFailed += checkWrotePartitions(prefixASCII, prefixBinary, npartitions, 1);
  
  // Check merged file
  //********************

  int errMergeASCII = splittedASCII.merge();
  int errMergeBinary = splittedBinary.merge();

  if(errMergeASCII != SPLITTED_FILE_SUCCESS){
    printf("Unable to merge ASCII file. Error code: %d\n",errMergeASCII);
    testsFailed++;
  }

  if(errMergeBinary != SPLITTED_FILE_SUCCESS){
    printf("Unable to merge Binary file. Error code: %d\n",errMergeBinary);
    testsFailed++;
  }
  
  testsFailed += checkWroteMerges(prefixASCII, prefixBinary, npartitions, 1);

  return testsFailed;
}

unsigned test2(pen_splittedFile& splittedASCII, const char* prefixASCII, pen_splittedFile& splittedBinary, const char* prefixBinary, unsigned npartitions){

  unsigned testsFailed = 0;

  //********************************
  //*   Test checkpoint merging    *
  //********************************

  printf("\n----------------------------\n");
  printf("  Test checkpoint merging");
  printf("\n----------------------------\n");

  // Create data (first iteration)
  //************************************
  
  createData(splittedASCII, splittedBinary, npartitions, 2);

  // Check partitions (first iteration)
  //************************************

  testsFailed += checkWrotePartitions(prefixASCII, prefixBinary, npartitions, 2);

  // Check merged file (first iteration)
  //**************************************

  int errMergeASCII = splittedASCII.merge();
  int errMergeBinary = splittedBinary.merge();

  if(errMergeASCII != SPLITTED_FILE_SUCCESS){
    printf("Unable to merge ASCII file. Error code: %d\n",errMergeASCII);
    testsFailed++;
  }

  if(errMergeBinary != SPLITTED_FILE_SUCCESS){
    printf("Unable to merge Binary file. Error code: %d\n",errMergeBinary);
    testsFailed++;
  }
  
  testsFailed += checkWroteMerges(prefixASCII, prefixBinary, npartitions, 2);
  
  // Create data (second iteration)
  //************************************
  
  createData(splittedASCII, splittedBinary, npartitions, 3);

  // Check partitions (second iteration)
  //************************************

  testsFailed += checkWrotePartitions(prefixASCII, prefixBinary, npartitions, 5);

  // Check merged file (second iteration)
  //**************************************

  errMergeASCII = splittedASCII.merge();
  errMergeBinary = splittedBinary.merge();

  if(errMergeASCII != SPLITTED_FILE_SUCCESS){
    printf("Unable to merge ASCII file. Error code: %d\n",errMergeASCII);
    testsFailed++;
  }

  if(errMergeBinary != SPLITTED_FILE_SUCCESS){
    printf("Unable to merge Binary file. Error code: %d\n",errMergeBinary);
    testsFailed++;
  }
  
  testsFailed += checkWroteMerges(prefixASCII, prefixBinary, npartitions, 5);

  
  
  return testsFailed;
}

unsigned test3(pen_splittedFile* splittedASCII, const char* prefixASCII, pen_splittedFile* splittedBinary, const char* prefixBinary, unsigned npartitions){

  unsigned testsFailed = 0;
  
  //************************
  //* Test merge at detroy *
  //************************

  printf("\n----------------------------\n");
  printf("     Test merge at detroy");
  printf("\n----------------------------\n");
  
  // Create data files
  //********************

  createData(*splittedASCII, *splittedBinary, npartitions, 1);
  
  // Check partitions
  //********************

  testsFailed += checkWrotePartitions(prefixASCII, prefixBinary, npartitions, 1);

  // Destroy split handlers
  //************************

  delete splittedASCII;
  delete splittedBinary;
  
  // Check merged file
  //********************
  
  testsFailed += checkWroteMerges(prefixASCII, prefixBinary, npartitions, 1);

  return testsFailed;
}

int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s n-partitions\n",argv[0]);
    return 1;
  }

  int npartitions = atoi(argv[1]);

  if(npartitions > 10){
    printf("Warning: 10 is the maximum number of partitions for this test.\n");
    npartitions = 10;
  }
  
  if(npartitions < 1){
    printf("Error: number of partitions must be greater than 0.\n");
    return -1;
  }

  //Splitted files handlers
  pen_splittedFile splittedASCII("ASCII",false);
  pen_splittedFile splittedBinary("Binary",true);

  pen_splittedFile splittedASCII2("ASCII2",false);
  pen_splittedFile splittedBinary2("Binary2",true);

  pen_splittedFile* splittedASCII3 = new pen_splittedFile("ASCII3",false);
  pen_splittedFile* splittedBinary3 = new pen_splittedFile("Binary3",true);
  
  unsigned testsFailed = 0;
  
  //*********************
  //* Test merge at end *
  //*********************

  testsFailed += test1(splittedASCII,"ASCII",splittedBinary,"Binary",(unsigned)npartitions);

  //********************************
  //*   Test checkpoint merging    *
  //********************************
  
  testsFailed += test2(splittedASCII2,"ASCII2",splittedBinary2,"Binary2",(unsigned)npartitions);

  //********************************
  //*   Test merging at delete     *
  //********************************
  
  testsFailed += test3(splittedASCII3,"ASCII3",splittedBinary3,"Binary3",(unsigned)npartitions);
  
  //*****************
  //   End report
  //*****************

  printf("\n----------------------------\n");
  printf("         Tests ended");
  printf("\n----------------------------\n");
  
  if(testsFailed == 0){
    printf("Test success!\n");
  }
  else{
    printf("Test with errors\n");    
  }
  
  return 0;
}
