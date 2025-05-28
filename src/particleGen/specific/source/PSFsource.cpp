
//
//
//    Copyright (C) 2019-2023 Universitat de València - UV
//    Copyright (C) 2019-2023 Universitat Politècnica de València - UPV
//    Copyright (C) 2024-2025 Vicent Giménez Alventosa
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

#include "PSFsource.hh" 

void psf_specificSampler::skip(const unsigned long long  dhists){

  unsigned long long remaining = dhists;

  // Iterate until all histories has been skipped
  while(remaining > 0){

    unsigned long long skipped = psf.skip(remaining);    

    // Check if end of particles chunk has been reached 
    if(skipped == 0){
      //Check for remaining chunks
      if(remainingChunks <= 0){
	return;
      }

      //Read next chunk
      size_t nread = 0;
      int err = pSF->read(getThread(),buffer,bufferSize,nread);
      if(err != SHARED_FILE_SUCCESS || nread != bufferSize){
	printf("psf_specificSampler:skip: Error on thread %d: Unable to read psf chunk.\n",getThread());
	printf("                Read data: %lu/%lu\n",nread,bufferSize);
	printf("           Returned value: %d\n",err);
	return;      
      }

      //Substract remaining chunks
      remainingChunks--;

      //Process binary data chunk
      size_t pos = 0;
      err = psf.read(buffer,pos,3);
      if(err != PEN_DUMP_SUCCESS){
	printf("psf_specificSampler:skip: Error processing binary psf data. Corrupted file.\n");
	printf("                   Error code: %d\n",err);
	return;
      }
      
    } else{
      remaining -= skipped;
    }
    
  }
  
}

void psf_specificSampler::sample(pen_particleState& state,
				 pen_KPAR& genKpar,
				 unsigned long long& dhist,
				 pen_rand& random){

  //Check for remaining splits
  if(splitted < requiredSplits){
    state.copyBase(splitState);
    genKpar = lastKpar;
    dhist = 0; //Splitted particles counts as secondaries
    ++splitted;
    return;
  }  

  //Read particles until someone survive the Russian roulette
  unsigned sumDHist = 0;
  for(;;){
    //Read next particle from phase space file
    unsigned kpar;
    unsigned long localDHist;
    if(psf.get(localDHist,kpar,state) == 0){
      //Phase space file is empty, refill the buffer
        

      //Check for remaining chunks
      if(remainingChunks <= 0){
	state.reset();
	genKpar = ALWAYS_AT_END;
	dhist = sumDHist;
	splitted = 0;
	requiredSplits = 0;      
	return;
      }
      
      //Read next chunk
      size_t nread = 0;
      int err = pSF->read(getThread(),buffer,bufferSize,nread);
      if(err != SHARED_FILE_SUCCESS || nread != bufferSize){
	printf("psf_specificSampler:sample: Error on thread %d: Unable to read psf chunk.\n",getThread());
	printf("                Read data: %lu/%lu\n",nread,bufferSize);
	printf("           Returned value: %d\n",err);
	state.reset();
	genKpar = ALWAYS_AT_END;
	dhist = sumDHist;
	remainingChunks = 0;
	splitted = 0;
	requiredSplits = 0;      
    return;      
      }

      //Substract remaining chunks
      remainingChunks--;

      //Process binary data chunk
      size_t pos = 0;
      err = psf.read(buffer,pos,3);
      if(err != PEN_DUMP_SUCCESS){
	printf("psf_specificSampler:sample: Error processing binary psf data. Corrupted file.\n");
	printf("                   Error code: %d\n",err);
	state.reset();
	genKpar = ALWAYS_AT_END;
	dhist = sumDHist;
	remainingChunks = 0;
	splitted = 0;
	requiredSplits = 0;  
	return;            
      }

      //Get state from processed chunk
      if(psf.get(localDHist,kpar,state) <= 0){
	//Misterious error
	printf("psf_specificSampler:sample:Error: This should not be happening, revise 'psf_specificSampler::sample' code!\n");
	state.reset();
	genKpar = ALWAYS_AT_END;
	dhist = sumDHist;
	remainingChunks = 0;
	splitted = 0;
	requiredSplits = 0; 
	return;
      }

      //printf("Thread %u: %s\n", getThread(), state.stringify().c_str());      
    }
    //Check if is a valid kpar
    if(isKpar(kpar)){
      genKpar = static_cast<pen_KPAR>(kpar);
    }
    else{
      //Is not a valid KPAR!!
      printf("psf_specificSampler:sample: Error on thread %d: Unknown read kpar (%u).\n",getThread(),kpar);
      state.reset();
      genKpar = ALWAYS_AT_END;
      dhist = sumDHist;
      remainingChunks = 0;
      splitted = 0;
      requiredSplits = 0;
      return;
    }

    //Check energy
    if(state.E > expectedMaxEnergy){
      //Energy is out of range!!
      printf("psf_specificSampler:sample: Error on thread %d: Particle out of energy range (%.5E).\n",getThread(),state.E);
      state.reset();
      genKpar = ALWAYS_AT_END;
      dhist = sumDHist;
      remainingChunks = 0;
      splitted = 0;
      requiredSplits = 0;      
      return;
    }

    //Apply rotation to each particle state
    if (rotation){
        double pos[3] = {state.X,state.Y,state.Z};
        double dir[3] = {state.U,state.V,state.W};
        matmul3D(particleRot,pos);
        matmul3D(particleRot,dir);
        state.X=pos[0];
        state.Y=pos[1];
        state.Z=pos[2];
        state.U=dir[0];
        state.V=dir[1];
        state.W=dir[2];
    }
        
    state.X += dx;
    state.Y += dy;
    state.Z += dz;
        
    //Add history increment
    sumDHist += localDHist;

    //Reset splitted particles
    splitted = 1;
    requiredSplits = 1;
    
    //Check if russian roulette or splitting technics must be applied
    if(state.WGHT < WGHTL){
      // Russian roulette
      if(VRR(state.WGHT,random)){
	//Particle survive, return this state
	dhist = sumDHist;
	return;
      }
      //Particle has been killed, continue to next iteration and get
      //the next particle in psf 
    }
    else if(state.WGHT > WGHTU){
      // Splitting
      //Update dhist
      dhist = sumDHist;

      //Calculate number of required splits
      double auxSplits = state.WGHT*RWGHTL;
      if(auxSplits >= static_cast<double>(MAXsplit)){
          requiredSplits = MAXsplit;
      }
      else{
          requiredSplits = static_cast<unsigned>(auxSplits);
      }

      requiredSplits = std::min(NSPLIT,requiredSplits);

      
      //If required, store current state for future copies 
      if(requiredSplits > 1){
	//Store kpar
	lastKpar = genKpar;
	
	//Modify weight
	state.WGHT /= double(requiredSplits);
	splitState.copyBase(state);
	splitState.ILB[0] += 1; //Split particles are treated as secondaries
	splitState.ILB[3] =  0;
      }
      return;
    }
    else{
      //Particle weight is in the window. Return this state with no
      //Russian roulette nor splitting applied
      dhist = sumDHist;
      return;
    }
  }
}

int psf_specificSampler::configure(double& Emax,
				   const pen_parserSection& config,
				   const unsigned nthreads,
				   const unsigned verbose){

  int err;

  //Clear phase space file
  psf.clear();

  // Read Emax
  //************
  err = config.read("Emax", Emax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("psfSource:configure: Error: Unable to read 'Emax' field. Number expected.\n");
    }
    Emax = 0.0;
    return -4;
  }

  if(Emax <= 0.0){
    if(verbose > 0){
      printf("psfSource:configure: Error: Invalid maximum energy (%14.5EeV). "
	     "Must be greater than zero.\n",Emax);
    }
    Emax = 0.0;
    return -5;
  }

  //Save expected maximum energy
  expectedMaxEnergy = Emax;

  // Read weight window
  //********************
  pen_parserArray window;
  err = config.read("wght-window", window);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("psfSource:configure: Error: Unable to read 'wght-window' in configuration. Array expected.\n");
    }
    return -6;
  }

  err = window[0].read(WGHTL);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("psfSource:configure: Error: Unable to read 'WGHTL' in array at position 0. Double expected\n");
    }
    return -7;
  }
  
  err = window[1].read(WGHTU);
  if(err != INTDATA_SUCCESS){    
    if(verbose > 0){
      printf("psfSource:configure: Error: Unable to read 'WGHTU' in array at position 1. Double expected\n");
    }
    return -8;
  }

  WGHTL = fabs(WGHTL);
  
  if(WGHTL > WGHTU){
    if(verbose > 0){
      printf("psfSource:configure: Error: Inconsistent weight window.\n");
      printf("          WGHTL = %15.4E\n",WGHTL);
      printf("          WGHTU = %15.4E\n",WGHTU);
    }
    return -9;
  }

  if(WGHTL > 0.0){
    RWGHTL = 1.0/WGHTL;
  }
  else{
    RWGHTL = 1.0e35;
  }

  // Read NSPLIT
  //************
  int auxNSplit;
  err = config.read("nsplit", auxNSplit);
  if(err != INTDATA_SUCCESS){
    auxNSplit = 1;
  }

  if(auxNSplit <= 0){
    printf("psfSource:configure: Error: 'nsplit' must be, as least, 1.\n");
    return -10;
  }

  NSPLIT = unsigned(auxNSplit);

  // Save Number of partitions
  //***************************
  npartitions = nthreads;
  if(npartitions <= 0){
    if(verbose > 0){
      printf("psfSource:configure: Error: Invalid number of partitions (%d). "
	     "No threads are used? Please, report this error.\n",npartitions);
    }
    return -12;    
  }

  // Rotation
  //***************************
  if(config.isSection(std::string("rotation"))){
      
    //Read rotation Euler angles 
    double omega, theta, phi;
    err = config.read("rotation/omega", omega);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfSource:configure: Error: Unable to read 'omega' Euler angle\n");
      }
      return -14;
    }
  
    err = config.read("rotation/theta", theta);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfSource:configure: Error: Unable to read 'theta' Euler angle\n");
      }
      return -15;
    }
  
    err = config.read("rotation/phi", phi);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfSource:configure: Error: Unable to read 'phi' Euler angle\n");
      }
      return -16;
    }
  
    //Change to rad
    omega*=M_PI/180.0;
    theta*=M_PI/180.0;
    phi*=M_PI/180.0;  
  
    if(omega != 0.0 || theta != 0.0 || phi != 0.0){
      rotation = true;
      if(verbose > 1){
	printf("Psf rotation has been applied with Euler angles values:\n"
	       "omega = %15.5E\n"
	       "theta = %15.5E\n"
	       "phi   = %15.5E\n", omega*180.0/M_PI, theta*180.0/M_PI, phi*180.0/M_PI);
      }
      //Create particle rotation matrix
      createRotationZYZ(omega,theta,phi,particleRot);
    }
  }else{
    rotation=false;
    if(verbose > 1){
      printf("No psf rotation has been applied.\n");
    }
  }


  // Translation
  //***************************
  if(config.isSection(std::string("translation"))){
    //Read translation 
    err = config.read("translation/dx", dx);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfSource:configure: Error: Unable to read 'dx' translation in x axis\n");
      }
      return -17;
    }
  
    err = config.read("translation/dy", dy);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfSource:configure: Error: Unable to read 'dy' translation in y axis\n");
      }
      return -18;
    }
  
    err = config.read("translation/dz", dz);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfSource:configure: Error: Unable to read 'dz' translation in z axis\n");
      }
      return -19;
    }
  
    if(verbose > 1){
      printf("Psf translation has been applied:\n"
	     "dx = %15.5E\n"
	     "dy = %15.5E\n"
	     "dz = %15.5E\n", dx, dy, dz);
    }
    
  }else
    {
      dx=0.0;
      dy=0.0;
      dz=0.0;
    
      if(verbose > 1){
        printf("No psf translation has been applied.\n");
      }
    }  
  
  //Try to create a shared file with this filename. Only on thread 0, it will be shared later
  if(getThread() == 0){

    // Read psf filename
    //*******************
    std::string psfFilename;
    err = config.read("filename", psfFilename);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfSource:configure: Error: No psf filename specified. String expected\n");
      }
      return -1;
    }

    // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_

    //Get process rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //Add MPI rank to filename 
    psfFilename = std::string("MPI") + std::to_string(rank) +
      std::string("-") + psfFilename;
#endif
    // ***************************** MPI END ********************************** //  
    
  
    //Create a shared sharedFile
    pSF = std::make_shared<pen_sharedFile>();

    //Open the specified file
    err = pSF->open(psfFilename.c_str(),true);
    if(err != SHARED_FILE_SUCCESS){
      if(verbose > 0){
	printf("psfSource:configure: Error: Unable to open file '%s' as shared file\n",
	       psfFilename.c_str());
      }
      return -2;
    }

    //Try to create a reader for each thread
    for(int i = 0; i < npartitions; ++i){
      err = pSF->createReader(i);
      if(err != SHARED_FILE_SUCCESS){
	if(verbose > 0){
	  printf("psfSource:configure: Error: Unable to create reader with ID %d. "
		 "Please, report this error.\n", i);
	}
	return -3;
      }
    }

    //Check file size and number of data chunks
    size_t fsize = pSF->size();

    //Check if specified file is empty
    if(fsize == 0){
      if(verbose > 0){
	printf("psfSource:configure: Error: Specified PSF ('%s') file is empty.\n",psfFilename.c_str());
      }
      return -13;
    }

    //Calculate the number of particle chunks
    nChunks = fsize/psf.memory();

    if(nChunks*psf.memory() != fsize){
      if(verbose > 0){
	printf("psfSource:configure: Error: PSF corrupted or incompatible with current "
	       "version of 'phase space file' library.\n");
      }
      return -14;
    }

    if(nChunks < (unsigned)npartitions){
      if(verbose > 1){
	printf("psfSource:configure: Warning: Number of psf chunks smaller than "
	       "specified partitions. The same number of partitions as chunks will be used\n");
      }
      npartitions = nChunks;
    }

    //Calculate chunks per partition and remaining chunks
    chunksPerPart = nChunks/npartitions;
    offsetChunks = nChunks - chunksPerPart*npartitions;

    if(verbose > 1){
      printf("Maximum psf energy (eV):\n");
      printf(" %14.5E\n",Emax);
      printf("PSF size:\n");
      printf(" %lu B\n",fsize);
      printf("Expected chunk size:\n");
      printf(" %lu B\n",psf.memory());
      printf("PSF partitions:\n");
      printf(" %d\n",npartitions);
      printf("PSF particle chunks:\n");
      printf(" %lu\n",nChunks);
      printf("Offset chunks:\n");
      printf(" %lu\n",offsetChunks);
      printf("NSPLIT:\n");
      printf(" %u\n",NSPLIT);
      printf("Weight window :\n");
      printf(" %15.4E - %15.4E\n\n",WGHTL,WGHTU);
    }

    //Calculate reader start chunks and move them
    long int psfChunkSize = (long int)psf.memory();
    size_t startChunk = 0;
    for(unsigned i = 0; i < (unsigned)npartitions; i++){
      size_t chunks2read = chunksPerPart;
      if(npartitions-i <= offsetChunks)
	chunks2read++;
      if(verbose > 1)
	printf("Partition %u start on chunk %lu and finish at %lu (%lu chunks)\n",
	       i,startChunk,startChunk+chunks2read,chunks2read);

      //Move reader to partition beginning
      for(size_t j = 0; j < startChunk; j++){
	//Skip one data chunk
	int fseekErr;
	err = pSF->seek(i,psfChunkSize,SEEK_CUR,&fseekErr);
	if(err != SHARED_FILE_SUCCESS || fseekErr != 0){
	  if(verbose > 0){
	    printf("psfSource:configure: Error: Error moving file cursor for partition %u.\n"
		   "          Shared file seek error code: %d\n"
		   "                     fseek error code: %d\n",
		   i, err,fseekErr);
	  }
	  return -20;
	}
      }

      startChunk += chunks2read;
    }

    //Set remaining chunks for first thread
    remainingChunks = chunksPerPart;
    
  }
  
  return 0;
}

int psf_specificSampler::sharedConfig(const psf_specificSampler& o){
  //Get missing configured vaues
  pSF = o.pSF;
    
  npartitions = o.npartitions;
  nChunks = o.nChunks;
  chunksPerPart = o.chunksPerPart;
  offsetChunks = o.offsetChunks;
  remainingChunks = chunksPerPart;
  if(npartitions-getThread() <= offsetChunks){
    ++remainingChunks;
  }
    
  return 0;
}

REGISTER_SPECIFIC_SAMPLER(psf_specificSampler,pen_particleState, PSF)
