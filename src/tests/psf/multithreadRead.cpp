
//
//
//    Copyright (C) 2019-2022 Universitat de València - UV
//    Copyright (C) 2019-2022 Universitat Politècnica de València - UPV
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


#include <cmath>
#include "pen_samplers.hh"
#include "splittedFile.hh"

void processPSF(pen_splittedFile* splittedFile, const unsigned ID, psf_specificSampler* sampler){

  pen_particleState state;
  pen_KPAR genKpar = PEN_ELECTRON;
  unsigned long long dhist;
  pen_rand random;

  char buffer[2000];
  
  sampler->sample(state,genKpar,dhist,random);
  while(genKpar != ALWAYS_AT_END){
    snprintf(buffer,2000,
	     "  %5llu  %4u %s\n",
	     dhist,genKpar,state.stringifyBaseNoGeo().c_str());
    splittedFile->write(ID,buffer,strlen(buffer));
    sampler->sample(state,genKpar,dhist,random);
  }
  
}

int main(int argc, char** argv){

 if(argc < 4){
    printf("usage: %s filename nreaders npartitions\n",argv[0]);
    return 1;
  }
  
  //Get nreaders and npartitions
  int nreaders    = atoi(argv[2]);
  int npartitions = atoi(argv[3]);

  if(nreaders < 1 || npartitions < 1){
    printf("Error: The number of readers and partitions must be greater than zero.\n");
    printf("     readers: %d\n",nreaders);
    printf("  partitions: %d\n",npartitions);
    return -1;
  }

  printf(" Number of partitions: %d\n",npartitions);
  printf(" Number of readers: %d\n",nreaders);
  
  //Create a phase space file
  pen_psfreader psf;

  //Create a vector of phase space file sources
  std::vector<psf_specificSampler> samplers;
  samplers.resize(nreaders);
  //Create a vector to store errors on source configurations
  std::vector<int> configErr;
  configErr.resize(nreaders,0);

  // Single thread read
  //**********************
  printf("Convert psf to ASCII format using one single thread.\n");

  //Read phase space file and write it to ASCII
  FILE* fin = nullptr;
  fin = fopen(argv[1],"rb");
  if(fin == nullptr){
    printf("Error: unable to open file '%s'.\n",argv[1]);
    return -2;
  }

  //Open output file
  FILE* fout = nullptr;
  fout = fopen("ascii.psf","w");
  if(fout == nullptr){
    printf("Error: unable to open file '%s'.\n","ascii.psf");
    return -3;
  }

  //Print header
  fprintf(fout,"# DHIST  KPAR %s\n",baseStateHeaderNoGeo());

  //Read input file until the end
  unsigned nchunks = 0;
  while(psf.read(fin,1) == PEN_PSF_SUCCESS){
    nchunks++;
    //Iterate over read states
    pen_particleState state;
    unsigned long dhist;
    unsigned kpar;
    while(psf.get(dhist,kpar,state) > 0){
      //Print base state without body and material information
      fprintf(fout,"  %5lu  %4u %s\n",dhist,kpar,state.stringifyBaseNoGeo().c_str());
    }
  }

  //Close files
  fclose(fin);
  fclose(fout);  

  // Single thread read
  //**********************
  printf("Convert psf to ASCII format using specified number of threads.\n");

  //Create splitted file
  pen_splittedFile splittedFile("sourcePSF-ascii",false);
  //Create one partition for each reader
  for(size_t i = 0; i < (unsigned)nreaders; i++){
    splittedFile.createPartition(i);
  }
  
  //Create configuration for sampler
  pen_parserSection config;

  config.set("filename",std::string(argv[1]));
  config.set("Emax",1.0e9);
  config.set("npartitions",npartitions);
  pen_parserArray wghtWindow;
  wghtWindow.append(0.0);
  wghtWindow.append(2.0);
  config.set("wght-window",wghtWindow);

  printf("Generated configuration:\n%s\n",config.stringify().c_str());
  
  //Set thread to each source
  for(size_t i = 0; i < (unsigned)nreaders; i++){
    samplers[i].setThread(i);
  }

  //Config all sources
  double Emax;
  configErr[0] = samplers[0].configure(Emax,config,3);
  for(size_t i = 1; i < (unsigned)nreaders; i++){
    configErr[i] = samplers[i].configure(Emax,config,1);
  }

  //Check configuration errors
  bool errors = false;
  for(unsigned i = 0; i < (unsigned)nreaders; i++){
    if(configErr[i] != 0){
      printf("Error at source configuration of thread %d. Error code: %d.\n",
	     i,configErr[i]);
      errors = true;
    }
  }

  if(errors)
    return -4;

  //Print header
  char buffer[2000];
  snprintf(buffer,2000,"# DHIST  KPAR %s\n",baseStateHeaderNoGeo());
  splittedFile.write(0,buffer,strlen(buffer));
  
  //Create a array of threads
  std::vector<std::thread> threads;
  //Run partial psf processing on each thread
  for(unsigned i = 0; i < (unsigned)nreaders; i++){
    threads.push_back(std::thread(processPSF,&splittedFile, i, &samplers[i]));
  }

  //Join threads
  for(unsigned i = 0; i < (unsigned)nreaders; i++){
    threads[i].join();
  }  

  //Flush partitions
  for(unsigned i = 0; i < (unsigned)nreaders; i++){
    splittedFile.flush(i);
  }    
  
  //Clear threads
  threads.clear();

  //Check each partition
  FILE* fascii = nullptr;
  fascii = fopen("ascii.psf","r");
  if(fascii == nullptr){
    printf("Error: unable to open file '%s'.\n","ascii.psf");
    return -5;
  }

  printf("\n *** Compare partitions with previous processed ASCII psf.\n");
  fflush(stdout);
  unsigned diffs = 0;
  unsigned nline1 = 0;
  for(unsigned i = 0; i < (unsigned)npartitions; i++){
    
    char buffer1[2000];
    char buffer2[2000];
    unsigned nline2 = 0;

    FILE* fsource = nullptr;
    char filename[200];
    snprintf(filename,200,"sourcePSF-ascii-%u",i);
    fsource = fopen(filename,"r");
    if(fsource == nullptr){
      printf("Error: unable to open file '%s'.\n",filename);
      return -6;
    }

    //Read lines until end of source file    
    while(fgets(buffer2,2000,fsource) != nullptr){

      if(fgets(buffer1,2000,fascii) == nullptr){
	printf("Error: End of psf ascii file reached\n");
	fclose(fsource);
	fclose(fascii);
	return -7;
      }
      nline1++;
      nline2++;
      
      if(strcmp(buffer1,buffer2) != 0){
	printf("\nlines differ!\n");
	printf("%10u:%s",nline1,buffer1);
	printf("%10u:%s",nline2,buffer2);
	printf("\n");
	diffs++;
      }
    }

    fclose(fsource);
  }

  fclose(fascii);
  
  if(diffs > 0){
    printf("Processed ASCII files are not equal!\n");
  }
  else{
    printf("Processed ASCII files match!\n");
  }
  
  printf("Read psf chunks: %u (%lu B)\n",nchunks,size_t(nchunks)*psf.memory());

  
  return 0;
}
