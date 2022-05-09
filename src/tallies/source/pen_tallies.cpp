
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


#include "pen_tallies.hh"


instantiator<pen_genericTally<pen_particleState>>& pen_commonTallyCluster::genericTallies(){
  static instantiator<pen_genericTally<pen_particleState>>* ans =
    new instantiator<pen_genericTally<pen_particleState>>;
  return *ans;  
}

#ifdef _PEN_USE_THREADS_
std::thread pen_commonTallyCluster::configure_async(const wrapper_geometry* geometry,
						    const abc_material* const materials[constants::MAXMAT],
						    const unsigned threadNum,
						    const pen_parserSection& config,
						    const unsigned verbose){
  return std::thread(&pen_commonTallyCluster::configure,this,
		     geometry,materials,threadNum,config,verbose);
}
#endif

void pen_commonTallyCluster::configure(const wrapper_geometry* geometry,
				       const abc_material* const materials[constants::MAXMAT],
				       const unsigned threadNum,
				       const pen_parserSection config,
				       const unsigned verbose){

  int err = 0;
    
  //Clear previous configuration
  clear();

  if(verbose > 1){
    printf("\n------------------------------------\n");
    printf("\n **** Tally group '%s'\n",name.c_str());
  }

  //Set thread
  nthread = threadNum;
  
  //Extract tally names
  std::vector<std::string> tallyNames;
  config.ls(tallyNames);
  
  //Iterate for each tally name
  for(unsigned i = 0; i < tallyNames.size(); i++){

    if(verbose > 1){
      printf("\n------------------------------------\n\n");
      printf("\nTally '%s':\n\n",tallyNames[i].c_str());
    }
    
    //Get subsection for tally 'i'
    pen_parserSection tallySec;
    if(config.readSubsection(tallyNames[i],tallySec) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("commonTallyCluster: configure: Error: unable to read section '%s' to configure tally.\n",tallyNames[i].c_str());
      }
      err++;
      continue;
    }

    //Try to read tally type ID
    std::string tallyID;
    if(tallySec.read("type",tallyID) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("commonTallyCluster: configure: Error: unable to read field %s/type. String expected\n",tallyNames[i].c_str());
	err++;
	continue;
      }
    }

    if(verbose > 1){
      printf("Tally type: '%s'\n\n",tallyID.c_str());
    } 

    //Try to read tally output dir
    std::string tallyOutDir;
    if(tallySec.read("outputdir",tallyOutDir) != INTDATA_SUCCESS){
      if(verbose > 2){
	printf("commonTallyCluster: configure: Warning: Unable to read field %s/outputdir. Assumed default output dir path ./\n",tallyNames[i].c_str());
      }
      tallyOutDir.clear();
    }

    if(verbose > 1){
      printf("Tally outputdir: '%s'\n\n",tallyOutDir.c_str());
    }
   
    
    //Try to create and configure tally 'i'
    if(createTally(tallyOutDir.c_str(),tallyID.c_str(),tallyNames[i].c_str(),*geometry,materials,tallySec,verbose) != 0){
      if(verbose > 0){
	printf("commonTallyCluster: configure: Error: Unable to create and configure tally '%s' of type '%s'.\n",tallyNames[i].c_str(),tallyID.c_str());
	err++;
	continue;
      }
    }    
  }
  if(verbose > 1){printf("\n------------------------------------\n\n");}

  if(err > 0){
    if(verbose > 0){
      printf("commonTallyCluster: configure: Error: %d tallies creation failed.\n",err);
    }
  }

  //Print created tallies
  if(verbose > 1){
    printf("\nCreated common tallies (name and type):\n\n");
    for(unsigned i = 0; i < tallies.size(); i++){
      printf("  %20s -> %20s  %-60s\n",
	     tallies[i]->readName().c_str(),
	     tallies[i]->readID(),
       tallies[i]->readOutputDirPath().c_str());
    }
  }

  configStatus = err;

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  //If this tally cluster runs on thread 0, we need a buffer cluster
  //to reduce the results of all MPI processes
  if(nthread == 0){
    mpiBuffer = new pen_commonTallyCluster;
    mpiBuffer->name = std::string("MPI_buffer_") + name;
    //Configure tally cluster with no verbose and a non zero thread
    //to avoid recursive creation of mpi buffers
    mpiBuffer->configure(geometry,
			 materials,
			 1,
			 config,
			 std::min(verbose,unsigned(1)));
    if(mpiBuffer->configureStatus() > 0){
      if(verbose > 0){
	printf("commonTallyCluster: configure: Error creating %d "
	       "tallies on MPI buffer cluster.\n",mpiBuffer->configureStatus());
      }
      configStatus = mpiBuffer->configureStatus();
    }
    else if(verbose > 1){
      printf("\nMPI cluster tally buffer created and configured.\n");
    }
  }
#endif
  // ***************************** MPI END ********************************** //
}

int pen_commonTallyCluster::createTally(const char* tallyOutDir, const char* ID,
					const char* tallyname,
					const wrapper_geometry& geometry,
					const abc_material* const materials[constants::MAXMAT],
					const pen_parserSection& config,
					const unsigned verbose){

  
  //Check if tally name already exists
  std::string strName;
  if(tallyname == nullptr || strName.assign(tallyname).length() == 0){
    if(verbose > 0){
      printf("commonTallyCluster: createTally: Error: empty tally name.\n");
    }
    return -3;
  }
  for(unsigned i = 0; i < tallies.size(); i++){
    if(strName.compare(tallies[i]->readName()) == 0){
      if(verbose > 0){
	printf("commonTallyCluster: createTally: Error: Tally name '%s' already used.\n", tallyname);
      }
      return -4;
    }
  }
  
  //Create specified tally by ID
  pen_genericTally<pen_particleState>* ptally = nullptr;  
  ptally = genericTallies().createInstance(ID);
  if(ptally == nullptr){
    if(verbose > 0){
      printf("commonTallyCluster: createTally: Error: unable to create a tally of type '%s'.\n",ID);
    }
    return -1;
  }

  //Set tally name
  ptally->setName(tallyname);

  //Set thread
  ptally->setThread(nthread);

  //Set OutDir
  ptally->setOutputDirPath(tallyOutDir);
  
  //Configure tally
  int errConfig = ptally->configure(geometry,materials,config,verbose);
  if(errConfig != 0){
    delete ptally;
    if(verbose > 0){
      printf("commonTallyCluster: createTally: Error: tally '%s' of type '%s' failed on configuration step.\n",tallyname,ID);
    }
    return -2;
  }
  
  //Append tally pointer to necessary vectors
  tallies.push_back(ptally);


  if(ptally->has_beginSim()){ tallies_beginSim.push_back(ptally);}
  if(ptally->has_endSim()){ tallies_endSim.push_back(ptally);}
  if(ptally->has_sampledPart()){ tallies_sampledPart.push_back(ptally);}
  if(ptally->has_endHist()){ tallies_endHist.push_back(ptally);}
  if(ptally->has_move2geo()){ tallies_move2geo.push_back(ptally);}
  if(ptally->has_beginPart()){ tallies_beginPart.push_back(ptally);}
  if(ptally->has_endPart()){ tallies_endPart.push_back(ptally);}
  if(ptally->has_localEdep()){ tallies_localEdep.push_back(ptally);}
  if(ptally->has_step()){ tallies_step.push_back(ptally);}
  if(ptally->has_interfCross()){ tallies_interfCross.push_back(ptally);}
  if(ptally->has_matChange()){ tallies_matChange.push_back(ptally);}
  if(ptally->has_jump()){ tallies_jump.push_back(ptally);}
  if(ptally->has_knock()){ tallies_knock.push_back(ptally);}
  if(ptally->has_lastHist()){ tallies_lastHist.push_back(ptally);}

  //Sort vector tallies by name
  std::sort(tallies.begin(), tallies.end(), __tallySort);
  
  return 0;
  
}

int pen_commonTallyCluster::writeDump(unsigned char*& pdump,
				      size_t& dim,
				      const unsigned long long nhist,
				      const int seed1, const int seed2,
				      const int lastSource,
				      const unsigned long long sourceHists,
				      const unsigned verbose) {

  if(tallies.size() < 1){
    if(verbose > 0){
      printf("commontallyCluster: writeDump: Error: Any tally has been registered.\n");
    }
    return -1;
  }
    
  //Create pointers to store dumped data of each tally
  std::vector<unsigned char*> vdumps;
  std::vector<size_t> dumpsDim;
  vdumps.resize(tallies.size());
  dumpsDim.resize(tallies.size());

  size_t dumpSize = 0;
  //To store last random seeds
  dumpSize += 2*sizeof(int32_t);
  //To store the total number of simulated histories
  dumpSize += sizeof(uint64_t); 
  //To store the actual source number identifier
  dumpSize += sizeof(int32_t);
  //To store the histories processed at that source
  dumpSize += sizeof(uint64_t); 

  //Get seeds
  const int32_t seeds[2] = {seed1,seed2}; 
  //Get nhist
  uint64_t nhistUint = static_cast<uint64_t>(nhist);
  //Get source ID
  const int32_t source32 = lastSource; 
  //Get completed source histories
  uint64_t nhistSourceUint = static_cast<uint64_t>(sourceHists);
  
  //Write all dumps and get its size, name and ID size
  for(unsigned i = 0; i < tallies.size(); i++){
    int err;
    err = tallies[i]->writeDump(vdumps[i],dumpsDim[i],verbose);
    if(err != PEN_DUMP_SUCCESS){
      //Free allocated data
      for(unsigned j = 0; j < i; j++){
	free(vdumps[i]);
	vdumps[i] = nullptr;
      }
      if(verbose > 0){
	printf("Error dumping data of tally %s with name %s (position %d)\n",
	       tallies[i]->readID(),tallies[i]->readName().c_str(), i);
      }
      return -2;
    }
    dumpSize += dumpsDim[i];
    dumpSize += 2*sizeof(uint32_t); //To store name and ID length
    dumpSize += tallies[i]->readName().length(); //To store name
    dumpSize += strlen(tallies[i]->readID()); //To store ID      
  }
    
  //Allocate memory
  if(dumpSize < 1){
    if(verbose > 0){
      printf("commontallyCluster: writeDump: Error: No data to dump.\n");
    }      
    return -3;
  }

  pdump = nullptr;
  pdump = (unsigned char*) malloc(dumpSize);
  if(pdump == nullptr){
    if(verbose > 0){
      printf("commontallyCluster: writeDump: Error: Allocation fail.\n");
    }      
    return -4;
  }

  //Write all data to dump buffer
  size_t pos = 0;

  // Write last random seeds
  memcpy(&pdump[pos],seeds,2*sizeof(int32_t));
  pos += 2*sizeof(int32_t);
  // Write number of simulated histories
  memcpy(&pdump[pos],&nhistUint,sizeof(uint64_t));
  pos += sizeof(uint64_t);
  // Write last completed source identifier
  memcpy(&pdump[pos],&source32,sizeof(int32_t));
  pos += sizeof(int32_t);
  // Write number of simulated histories at last source
  memcpy(&pdump[pos],&nhistSourceUint,sizeof(uint64_t));
  pos += sizeof(uint64_t);
      
  for(unsigned i = 0; i < tallies.size(); i++){
    uint32_t nameSize = (uint32_t)tallies[i]->readName().length();
    uint32_t IDsize = (uint32_t)strlen(tallies[i]->readID());

    //Save name size
    memcpy(&pdump[pos],&nameSize,sizeof(uint32_t));
    pos += sizeof(uint32_t);
    //Save name
    if(nameSize > 0){
      memcpy(&pdump[pos],tallies[i]->readName().c_str(),nameSize);
      pos += nameSize;
    }

    //Save ID size
    memcpy(&pdump[pos],&IDsize,sizeof(uint32_t));
    pos += sizeof(uint32_t);
    //Save ID
    if(IDsize > 0){
      memcpy(&pdump[pos],tallies[i]->readID(),IDsize);
      pos += IDsize;
    }

    //Store dump
    memcpy(&pdump[pos],vdumps[i],dumpsDim[i]);
    pos += dumpsDim[i];
  }

  //Check dumped data
  if(pos != dumpSize){
    if(verbose > 0){
      printf("commontallyCluster: writeDump: Error: Dumped data size doesn't match with allocated.\n");
      printf("                                dumped: %lu\n",pos);
      printf("                                dumped: %lu\n",dumpSize);
    }
    free(pdump);
    pdump = nullptr;
    return -5;      
  }

  //Free memory
  for(unsigned i = 0; i < tallies.size(); i++)
    free(vdumps[i]);

  //Store final dimension
  dim = dumpSize;
    
  return 0;
}

int pen_commonTallyCluster::readDump(const unsigned char* const pdump,
				     size_t& pos,
				     unsigned long long& nhist,
				     int& seed1, int& seed2,
				     int& lastSource,
				     unsigned long long& sourceHists,
				     const unsigned verbose){

  if(tallies.size() < 1){
    if(verbose > 0){
      printf("commontallyCluster: readDump: Error: Any tally has been registered.\n");
    }
    return -1;
  }

  //Read seeds
  int32_t seeds[2];
  memcpy(seeds,&pdump[pos],2*sizeof(int32_t));
  pos += 2*sizeof(int32_t);

  seed1 = seeds[0];
  seed2 = seeds[1];

  // Read simulated nhistories
  uint64_t nhistUint;
  memcpy(&nhistUint,&pdump[pos],sizeof(uint64_t));
  pos += sizeof(uint64_t);
  nhist = static_cast<unsigned long long>(nhistUint);

  // Read last used source identifier
  int32_t source32;
  memcpy(&source32,&pdump[pos],sizeof(int32_t));
  pos += sizeof(int32_t);
  lastSource = source32;

  // Read number of completed nhistories at that source
  memcpy(&nhistUint,&pdump[pos],sizeof(uint64_t));
  pos += sizeof(uint64_t);
  sourceHists = static_cast<unsigned long long>(nhistUint);
  
  
  //Read tally dumped data
  for(unsigned i = 0; i < tallies.size(); i++){
    uint32_t nameSize, IDsize;

    //Read name size
    memcpy(&nameSize,&pdump[pos],sizeof(uint32_t));
    pos += sizeof(uint32_t);      
    if(nameSize < 1){
      if(verbose > 0){
	printf("commontallyCluster: readDump: Error: Dumped tally %d has null name size.\n",i);
      }
      return -2;
    }
      
    //Compare names
    if(memcmp(&pdump[pos], tallies[i]->readName().c_str(),nameSize) != 0){
      if(verbose > 0){
	printf("commontallyCluster: readDump: Error: Tally names doesn't match at read tally number %d.\n",i);
	printf("                             Dumped name: %*s\n",nameSize,&pdump[pos]);
	printf("                           Expected name: %s\n",tallies[i]->readName().c_str());
      }
      return -2;	
    }
    pos += nameSize;
      

    //Read ID size
    memcpy(&IDsize,&pdump[pos],sizeof(uint32_t));
    pos += sizeof(uint32_t);      
    if(IDsize < 1){
      if(verbose > 0){
	printf("commontallyCluster: readDump: Error: Dumped tally %d has null ID size.\n",i);
      }
      return -3;
    }

    //Compare IDs
    if(memcmp(&pdump[pos], tallies[i]->readID(),IDsize) != 0){
      if(verbose > 0){
	printf("commontallyCluster: readDump: Error: Tally IDs doesn't match at read tally number %d.\n",i);
	printf("                             Dumped ID: %*s\n",IDsize,&pdump[pos]);
	printf("                           Expected ID: %s\n",tallies[i]->readID());
      }
      return -3;	
    }
    pos += IDsize;      

    int err = tallies[i]->readDump(pdump,pos,verbose);
    if(err != PEN_DUMP_SUCCESS){
      printf("commontallyCluster: readDump: Error reading dumped data of tally %d.\n",i);
      return -4;
    }
  }

  return 0;
}

int pen_commonTallyCluster::dump2file(const char* filename,
				      const unsigned long long nhist,
				      const int seed1, const int seed2,
				      const int lastSource,
				      const unsigned long long sourceHists,
				      const unsigned verbose){

  unsigned char* pdump = nullptr;
  size_t dumpSize;

  int err;
    
  //Create dump
  err = writeDump(pdump,dumpSize,nhist,
		  seed1,seed2,lastSource,
		  sourceHists,verbose);
  if(err != 0){
    if(verbose > 0){
      printf("commontallyCluster: dump2file: Error creating dump.\n");
    }
    return -1;
  }

    std::string auxstr(filename); 
    std::size_t found = auxstr.find_last_of("/\\");

    std::string dumpFilename;
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  int rank;
  //Add the MPI rank number to filename
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
      if(found != std::string::npos){
	dumpFilename = auxstr.substr(0,found+1) + std::string("MPI") + std::to_string(rank) +
	  std::string("-th") + std::to_string(nthread) + auxstr.substr(found+1);
      }else{
	dumpFilename = std::string("MPI") + std::to_string(rank) +
	  std::string("-th") + std::to_string(nthread) + filename;
      }
#else
      // ***************************** MPI END ********************************** //
      if(found != std::string::npos){
	dumpFilename = auxstr.substr(0,found+1) + std::string("th") + std::to_string(nthread) + auxstr.substr(found+1);
      }else{
	dumpFilename = std::string("th") + std::to_string(nthread) + filename;
      }
#endif
  
  FILE* fdump = nullptr;
  fdump = fopen(dumpFilename.c_str(), "wb");
  if(fdump == nullptr){
    if(verbose > 0){
      printf("commontallyCluster: dump2file: Error creating dump file %s.\n",dumpFilename.c_str());
    }
    free(pdump);
    return -2;
  }

  //Write data
  size_t nwrite = fwrite(pdump,1,dumpSize,fdump);
  if(nwrite != dumpSize){
    if(verbose > 0){
      printf("commontallyCluster: dump2file: Failed to write all data.\n");
    }
    free(pdump);
    return -3;
  }

  //Close file
  fclose(fdump);
    
  //Free memory
  free(pdump);
    
  return 0;
}

int pen_commonTallyCluster::readDumpfile(const char* filename,
					 unsigned long long& nhist,
					 int& seed1, int& seed2,
					 int& lastSource,
					 unsigned long long& sourceHists,
					 const unsigned verbose){

  unsigned char* pdump = nullptr;
  size_t dumpSize;

  int err;

  //Read data
  FILE* fdump = nullptr;
  fdump = fopen(filename, "rb");
  if(fdump == nullptr){
    if(verbose > 0){
      printf("commontallyCluster: readDumpfile: Error opening dump file %s.\n",filename);
    }
    return -2;
  }

  // Get file size
  //***************

  unsigned char buffer[5000];
  dumpSize = 0;
  size_t nread = 0;
  while((nread = fread(buffer,sizeof(unsigned char),5000,fdump)) > 0){
    dumpSize += nread;
  }

  //Return to file beginning
  rewind(fdump);

  // Allocate memory to store the entire file
  //*******************************************

  pdump = (unsigned char*) malloc(dumpSize);
  if(pdump == nullptr){
    if(verbose > 0){
      printf("commontallyCluster: readDumpfile: Bad alloc.\n");
    }
    return -4;
  }
    
  // Read data
  //*************
  nread = fread(pdump,1,dumpSize,fdump);
  if(nread != dumpSize){
    if(verbose > 0){
      printf("commontallyCluster: readDumpfile: Failed to read all data file.\n");
    }
    free(pdump);
    return -5;
  }

  //Close file
  fclose(fdump);


  // Extract data
  //***************

  nread = 0;
  err = readDump(pdump,nread,nhist,seed1,seed2,
		 lastSource,sourceHists,verbose);
  if(err != 0){
    if(verbose > 0){
      printf("commontallyCluster: readDumpfile: Failed to extract data.\n");
    }
    free(pdump);
    return -6;
  }

  //Free memory
  free(pdump);
    
  return 0;
}

int pen_commonTallyCluster::sum(pen_commonTallyCluster& cluster,
				const unsigned verbose){
    //Check if both clusters have the same ammount of tallies
    unsigned size1 = tallies.size();
    unsigned size2 = cluster.tallies.size();
    if(size1 != size2){
      if(verbose > 0){
	printf("pen_commonTallyCluster: sum: Number of tallies "
	       "doesn't match on clusters %s and %s of threads %u and %u.\n",
	       name.c_str(),cluster.name.c_str(),
	       getThread() ,cluster.getThread());
	printf("tallies at local cluster (thread %u): %u\n",
	       getThread(),size1);
	printf("tallies at input cluster (thread %u): %u\n",
	       cluster.getThread(),size2);

	printf("Expected tally names:\n");
	for(size_t i = 0; i < size1; ++i){
	  printf("%s\n",tallies[i]->readName().c_str());
	}
	printf("Input tally names:\n");
	for(size_t i = 0; i < size2; ++i){
	  printf("%s\n",cluster.tallies[i]->readName().c_str());
	}
      }
      return -1;
    }

    //Check if both clusters have the same tally types and names
    for(unsigned i = 0; i < size1; i++){
      if(tallies[i]->readName().compare(cluster.tallies[i]->readName()) != 0){
	if(verbose > 0){
	  printf("pen_commonTallyCluster: sum: Tallies names doesn't match:\n");
	  printf("                          %s\n",tallies[i]->readName().c_str());
	  printf("                          %s\n",cluster.tallies[i]->readName().c_str());	  
	}      	
	return i+1;
      }
      if(strcmp(tallies[i]->readID(), cluster.tallies[i]->readID()) != 0){
	if(verbose > 0){
	  printf("pen_commonTallyCluster: sum: Tallies IDs doesn't match:\n");
	  printf("                          %s\n",tallies[i]->readID());
	  printf("                          %s\n",cluster.tallies[i]->readID());
	}
	return i+1;
      }
    }

    //Sum all tallies
    int err = 0;
    for(unsigned i = 0; i < size1; ++i){

      //Flush both tallies
      tallies[i]->flush();
      cluster.tallies[i]->flush();
      
      int err2 = tallies[i]->sum(*(cluster.tallies[i]));
      if(err2 != 0){
	if(verbose > 0){
	  printf("pen_commonTallyCluster: sum: Error adding up tallies %s (position %d)\n",tallies[i]->readName().c_str(),i);
	  printf("                             Error code: %d\n", err2);
	}
	++err;
      }
    }

    return err;
  }

// ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
int pen_commonTallyCluster::reduceMPI(unsigned long long& simulatedHists,
				      const MPI_Comm comm,
				      const unsigned verbose){

  //Get current process "rank" and the total number of processes
  int rank;
  int mpiSize;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpiSize);

  //Get global rank
  int globalRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);

  if(nthread != 0){
    if(verbose > 0){
      printf("MPI process %d: Error: MPI reduce step must be done from thread with ID 0. This thread has the ID %u",globalRank,nthread);
    }
    return -1;
  }
  
  if(mpiSize <= 1){
    // Only one process should be reduced in the current chunk,
    // i.e. nothing to do
    if(verbose > 2){
      printf("Rank %d: Local rank %d/%d: No more processes in this group, nothing to reduce.\n",globalRank,rank,mpiSize-1);
    }    
    return 0;
  }
  
  //Calculate the number of MPI processes to reduce in this step.
  //Each step can reduce a power of two number of processes
  unsigned toReduce = mpiSize;
  while(!isPowerOf2(toReduce)){
    --toReduce;
  }
  unsigned remaining = mpiSize-toReduce;

  if(verbose > 2){
    printf("Rank %d: Local rank %d/%d: Current step will reduce %u processes, %u remaining\n",globalRank,rank,mpiSize-1,toReduce,remaining);
  }

  //Check if this rank must perform the reduction in this step
  if(rank >= 0 && rank < (int)toReduce){
    //This rank will be reduced in this processes chunk
    
    if(remaining > 0){
      //Create a sub communicator for remaining processes
      MPI_Comm subComm;
      MPI_Comm_split(comm, MPI_UNDEFINED, rank, &subComm);
    }

    //Sum up tallies of all MPI processes
    unsigned reduceIter = logb2(toReduce); //Calculate reduce iterations
    unsigned twoPi = 1; //2^i
    for(unsigned ired = 1; ired <= reduceIter; ++ired){
      unsigned twoPim1 = twoPi; //2^(ired-1)
      twoPi *= 2;
      if(rank % twoPi == 0){
	//This process will sumup the results	
	int mpiTag = ired;
	MPI_Status status;
	//Recive dump data from process rank+2^(ired-1)
	int source = rank + twoPim1;

	if(verbose > 2){
	  printf("Rank %d: Local rank %d/%d: Receive results of local rank %d\n",globalRank,rank,mpiSize-1,source);
	}
	
	//First, receive dump size
	unsigned long long dumpSize;
	MPI_Recv(&dumpSize,1,MPI_UNSIGNED_LONG_LONG,source,mpiTag,comm,&status);

	//If no dump information will be send, skip to next iteration
	if(dumpSize == 0){
	  if(verbose > 2){
	    printf("Rank %d: Local rank %d/%d: No dump information will be send, skip to next iteration.\n",globalRank,rank,mpiSize-1);
	  }
	  continue;
	} else if(verbose > 2){
	  printf("Rank %d: Local rank %d/%d: %llu bytes expected.\n",globalRank,rank,mpiSize-1,dumpSize);	  
	}
	
	//Allocate memory to store the entire dump
	size_t dim = dumpSize;
	unsigned char* pdump = nullptr;
	pdump = (unsigned char*) malloc(dim);
	if(pdump == nullptr){
	  if(verbose > 0){
	    printf("MPI process %d: Error: Unable to allocate memory to receive dump (%llu bytes) ",globalRank,(unsigned long long)dim);
	  }
	}
      
	//Next, calculate how many receives are required to receive all data
	const int intmax = std::numeric_limits<int>::max();
	size_t nSends = dim/intmax;
	const int residual = int(dim % size_t(intmax));
      
	size_t pos = 0;
	for(size_t isend = 0; isend < nSends; ++isend){
	  //Receive dump partition
	  if(verbose > 2){
	    printf("             Receive partition %llu with %d bytes.\n",(unsigned long long)isend,intmax);	    
	  }
	  MPI_Recv(&pdump[pos],intmax,MPI_BYTE,source,mpiTag,comm,&status);
	  pos += intmax;
	}

	if(residual > 0){
	  //Receive residual dump bytes
	  if(verbose > 2){
	    printf("             Receive partition %llu with %d bytes.\n",(unsigned long long)nSends,residual);	    
	  }	  
	  MPI_Recv(&pdump[pos],residual,MPI_BYTE,source,mpiTag,comm,&status);
	  pos += residual;
	}
	
	//Check ammount of received data
	if(pos != dim){
	  if(verbose > 0){
	    printf("MPI process %d: Error: Expected and received data mismatch'\n       Expected: %llu\n       Received: %llu\n",globalRank,(unsigned long long)dim,(unsigned long long)pos);
	  }
	}

	//To sumup the results, load dump information into MPI buffer cluster
	unsigned long long dumpHists;
	int dummySeed1,dummySeed2;
	int dummyUint;
	unsigned long long dummyull;
	size_t posDump = 0;
	int err = mpiBuffer->readDump(pdump,posDump,dumpHists,
				      dummySeed1,dummySeed2,
				      dummyUint,dummyull,
				      verbose);
	
	//Free memory
	free(pdump);
	pdump = nullptr;
	
	if(err == 0){
	  //Join results
	  sum(*mpiBuffer,verbose);

	  //Add simulated histories
	  simulatedHists += dumpHists;
	}
	else if(verbose > 0){
	  printf("MPI process %d: Error: Unable to read dump data\n",globalRank);
	}
      }
      else if((rank-twoPim1) % twoPi == 0){
	//The results of this process will be added to the results
	//stored at process with rank-2^(i-1)
	int mpiTag = ired;
	int destination = rank-twoPim1;

	if(verbose > 2){
	  printf("Rank %d: Local rank %d/%d: Send results to local rank %d\n",globalRank,rank,mpiSize-1,destination);
	}
	
	//First, dump the results. Seed information will be filled
	//with dummy data (rank identifier)
	unsigned char* pdump = nullptr;
	size_t dim;
	int err = writeDump(pdump,dim,simulatedHists,
			    rank,rank,0,0,verbose);
	if(err != 0){
	  if(verbose > 0){
	    printf("reduceMPI: Rank %d dump2file: Error creating dump.\n",globalRank);
	  }
	  dim = 0;
	}
	//Convert dump size to unsigned long long to send it
	unsigned long long int dimll = dim;
	if(size_t(dimll) != dim){ //I think this is not possible, but...
	  if(verbose > 0){
	    printf("reduceMPI: Rank %d: Error: Unable to convert 'size_t' to 'unsigned long long'\n       %llu != %llu",rank,dimll,(unsigned long long)dim);
	  }
	  dimll = 0;
	}
	//Send dump size
	MPI_Send(&dimll,1,MPI_UNSIGNED_LONG_LONG,destination,mpiTag,comm);

	if(dimll == 0){
	  if(verbose > 2){
	    printf("Rank %d: Local rank %d/%d: No data to send. Send 0 as data size and skip to next iteration.\n",globalRank,rank,mpiSize-1);
	  }
	  if(pdump != nullptr)
	    free(pdump);
	  pdump = nullptr;
	  continue;
	} else if(verbose > 2){
	  printf("Rank %d: Local rank %d/%d: Send %llu bytes of data.\n",globalRank,rank,mpiSize-1,dimll);
	}

	//Calculate how many sends are required to send all data
	const int intmax = std::numeric_limits<int>::max();
	const size_t nSends = dim/intmax;
	const int residual = int(dim % size_t(intmax));
      
	size_t pos = 0;
	for(size_t isend = 0; isend < nSends; ++isend){
	  //Send dump partition
	  if(verbose > 2){
	    printf("             Send partition %llu with %d bytes.\n",(unsigned long long)isend,intmax);	    
	  }	  
	  MPI_Send(&pdump[pos],intmax,MPI_BYTE,destination,mpiTag,comm);	
	  pos += intmax;
	}

	if(residual > 0){
	  if(verbose > 2){
	    printf("             Send partition %llu with %d bytes.\n",(unsigned long long)nSends,residual);	    
	  }	  
	  //Send residual dump bytes
	  MPI_Send(&pdump[pos],residual,MPI_BYTE,destination,mpiTag,comm);
	  pos += residual;
	}
	
	//Free memory
	free(pdump);
	pdump = nullptr;
	
	//Check ammount of send data
	if(pos != dim){
	  if(verbose > 0){
	    printf("MPI process %d: Error: Expected and send data mismatch'\n       Expected: %llu\n           Send: %llu\n",rank,(unsigned long long)dim,(unsigned long long)pos);
	  }
	}
      
      }
    }
    
  }
  else if(remaining > 0){
    //This rank will be reduced in another processes chunk.
    
    //Create a sub communicator for remaining processes
    MPI_Comm subComm;
    MPI_Comm_split(comm, 1, rank, &subComm);

    if(verbose > 2){
      printf("Rank %d: Local rank %d/%d: Out of range in this reduce step, going to subgroup reduce.\n",globalRank,rank,mpiSize-1);
    }
    
    reduceMPI(simulatedHists,subComm,verbose);
    
    //Free communicator
    MPI_Comm_free(&subComm);
  }
  
  //Wait until all comm processes have been finished their reduces
  MPI_Barrier(comm);

  //If there are some remaining processes, we need to reduce the results
  //of rank 0 and "toReduce"
  if(remaining > 0){

    if(verbose > 2){
      printf("Rank %d: Local rank %d/%d: Reducing first ranks of main and sub reduce groups.\n",globalRank,rank,mpiSize-1);
    }
    
    //Create a sub communicator for process 0 and "toReduce"
    if(rank == 0 || rank == (int)toReduce){
      MPI_Comm subComm;
      MPI_Comm_split(comm, 1, rank, &subComm);

      //Reduce results
      reduceMPI(simulatedHists,subComm,verbose);
    
      //Free communicator
      MPI_Comm_free(&subComm);      
    }
    else{
      //Nothing to do
      MPI_Comm subComm;
      MPI_Comm_split(comm, MPI_UNDEFINED, rank, &subComm);
    }
  }

  //Final results are stored in tally cluster of rank 0 and thread 0
  return 0;
}
#endif
// ***************************** MPI END ********************************** //

void pen_commonTallyCluster::clear(){
  for(std::size_t i = 0; i < tallies.size(); i++){
    delete tallies[i];
  }
  tallies.clear();
  tallies_beginSim.clear();
  tallies_endSim.clear();
  tallies_sampledPart.clear();
  tallies_endHist.clear();
  tallies_move2geo.clear();
  tallies_beginPart.clear();
  tallies_endPart.clear();
  tallies_localEdep.clear();
  tallies_step.clear();
  tallies_interfCross.clear();
  tallies_matChange.clear();
  tallies_jump.clear();
  tallies_knock.clear();
  tallies_lastHist.clear();

  if(mpiBuffer != nullptr)
    delete mpiBuffer;
  mpiBuffer = nullptr;

  nthread = 0;
  configStatus = 0;
}

pen_commonTallyCluster::~pen_commonTallyCluster(){
  clear();
}

bool __tallySort (pen_genericTally<pen_particleState>* i,
		  pen_genericTally<pen_particleState>* j) {
  return (i->readName().compare(j->readName()) < 0);
}


//Include defined tallies
#include "genericTallies.cpp"
#include "specificTallies.cpp"
