
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
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

#include "memoryPSFsource.hh" 

void psfMemory_specificSampler::skip(const unsigned long long  dhists){

  unsigned long long remaining = dhists;

  // Iterate until all histories has been skipped
  while(remaining > 0){

    //Get actual chunk
    const pen_psfMemort::chunk& chunk = (*pchunks)[actualChunk];

    //Iterate over chunk particles to skip them
    for(actualPart = 0; actualPart < chunk.nPart; ++actualPart){
      
      //Get actual particle
      const pen_psfMemort::particle& particle = chunk.particles[actualPart];      
      if(particle.dhist < remaining)
	remaining -= particle.dhist;
      else{
	remaining = 0;
	break;
      }
    }

    if(remaining != 0){
      //End of chunk reached, skip it
      actualChunk += npartitions;
    }
  }
}

void psfMemory_specificSampler::sample(pen_particleState& state,
				       pen_KPAR& genKpar,
				       unsigned long long& dhist,
				       pen_rand& random){
  
  //Check for remaining splits
  if(splitted < requiredSplits){
    state.copyBase(splitState);
    genKpar =
      static_cast<pen_KPAR>((*pchunks)[actualChunk].particles[actualPart].kpar);
    dhist = 0; //Splitted particles counts as secondaries
    ++splitted;
    return;
  }
  
  //Read particles until someone survive the Russian roulette
  unsigned sumDHist = 0;
  for(;;){
    
    //Check if the chunk contains more particles
    if(actualPart >= (*pchunks)[actualChunk].nPart){

      //Move to the next chunk if exists
      actualPart = 0;
      actualChunk += npartitions;

      
      if(actualChunk > pchunks->size()){
	//No remaining chunks, finish sampling
	state.reset();
	genKpar = pen_KPAR::ALWAYS_AT_END;
	dhist = sumDHist;
	splitted = 0;
	requiredSplits = 0;      
	return;	
      }

      //Get particle from next chunk
      if((*pchunks)[actualChunk].nPart == 0){
	//Empty chunk, should not happen
	printf("psfMemory_specificSampler:sample:Error: This should not be happening, "
	       "revise 'psfMemory_specificSampler::sample' code!\n");
	printf("    Actual chunk   : %lu\n",
	       static_cast<unsigned long>(actualChunk));
	printf("    Chunk particles: %lu\n",
	       static_cast<unsigned long>((*pchunks)[actualChunk].nPart));
	state.reset();
	genKpar = pen_KPAR::ALWAYS_AT_END;
	dhist = sumDHist;
	splitted = 0;
	requiredSplits = 0; 
	return;	
      }
    }

    //Get final chunk
    const pen_psfMemort::chunk& chunk = (*pchunks)[actualChunk];
    
    //Get actual particle
    const pen_psfMemort::particle& particle = chunk.particles[actualPart];

    //Increase actual particle
    actualPart += 1;
    
    unsigned kpar = particle.kpar;
    unsigned long localDHist = particle.dhist;
    
    //Check if is a valid kpar
    if(isKpar(kpar)){
      genKpar = static_cast<pen_KPAR>(kpar);
    }
    else{
      //Is not a valid KPAR!!
      printf("psf_specificSampler:sample: Error on thread %d: "
	     "Unknown read kpar (%u).\n",getThread(),kpar);
      state.reset();
      genKpar = pen_KPAR::ALWAYS_AT_END;
      dhist = sumDHist;
      splitted = 0;
      requiredSplits = 0;
      return;
    }

    //Save particle state
    state = particle.state;

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

int psfMemory_specificSampler::configure(double& Emax,
					 const pen_parserSection& config,
					 const unsigned nthreads,
					 const unsigned verbose){

  
  int err;

  Emax = 0.0;

  // Read weight window
  //********************
  pen_parserArray window;
  err = config.read("wght-window", window);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("psfMemorySource:configure: Error: Unable to read 'wght-window' in configuration. Array expected.\n");
    }
    return -6;
  }

  err = window[0].read(WGHTL);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("psfMemorySource:configure: Error: Unable to read 'WGHTL' in array at position 0. Double expected\n");
    }
    return -7;
  }
  
  err = window[1].read(WGHTU);
  if(err != INTDATA_SUCCESS){    
    if(verbose > 0){
      printf("psfMemorySource:configure: Error: Unable to read 'WGHTU' in array at position 1. Double expected\n");
    }
    return -8;
  }

  WGHTL = fabs(WGHTL);
  
  if(WGHTL > WGHTU){
    if(verbose > 0){
      printf("psfMemorySource:configure: Error: Inconsistent weight window.\n");
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
    printf("psfMemorySource:configure: Error: 'nsplit' must be, as least, 1.\n");
    return -10;
  }

  NSPLIT = unsigned(auxNSplit);

  // ** Rotation  
  if(config.isSection(std::string("rotation"))){
      
    //Read rotation Euler angles 
    double omega, theta, phi;
    err = config.read("rotation/omega", omega);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfMemorySource:configure: Error: Unable to read 'omega' Euler angle\n");
      }
      return -14;
    }
  
    err = config.read("rotation/theta", theta);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfMemorySource:configure: Error: Unable to read 'theta' Euler angle\n");
      }
      return -15;
    }
  
    err = config.read("rotation/phi", phi);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfMemorySource:configure: Error: Unable to read 'phi' Euler angle\n");
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
      
  // ** Translation  
  if(config.isSection(std::string("translation"))){
    //Read translation 
    err = config.read("translation/dx", dx);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfMemorySource:configure: Error: Unable to read 'dx' translation in x axis\n");
      }
      return -17;
    }
  
    err = config.read("translation/dy", dy);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfMemorySource:configure: Error: Unable to read 'dy' translation in y axis\n");
      }
      return -18;
    }
  
    err = config.read("translation/dz", dz);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfMemorySource:configure: Error: Unable to read 'dz' translation in z axis\n");
      }
      return -19;
    }
  
    if(verbose > 1){
      printf("Psf translation has been applied:\n"
	     "dx = %15.5E\n"
	     "dy = %15.5E\n"
	     "dz = %15.5E\n", dx, dy, dz);
    }
    
  }else{
    dx=0.0;
    dy=0.0;
    dz=0.0;
    
    if(verbose > 1){
      printf("No psf translation has been applied.\n");
    }
  }

  // Save Number of partitions
  //***************************
  npartitions = nthreads;
  if(npartitions <= 0){
    if(verbose > 0){
      printf("psfMemorySource:configure: Error: Invalid number of partitions: %d.\n",npartitions);
    }
    return -12;    
  }
  
  //Process psf file in thread 0
  if(getThread() == 0){

    pen_psfreader psf;

    //Clear phase space file
    psf.clear();

    // Read psf filename
    //*******************
    std::string psfFilename;
    err = config.read("filename", psfFilename);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("psfMemorySource:configure: Error: No psf 'filename' specified. String expected\n");
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


    //Open specified file
    FILE* psfFile = nullptr;
    psfFile = fopen(psfFilename.c_str(), "rb");
    if(psfFile == nullptr){
      if(verbose > 0){
	printf("psfSource:configure: Error: Unable to open file '%s'\n",
	       psfFilename.c_str());
      }
      return -2;
    }

    //Read the psf until the end of file
    for(;;){

      //Read next chunk
      int errPSF = psf.read(psfFile, verbose);
      if(errPSF == PEN_PSF_BAD_READ){
	break; //End of file
      }else if(errPSF == PEN_PSF_BAD_DUMP_READ){
	//The file is corrupted
	if(verbose > 0){
	  printf("psfMemorySource:configure: Error: The PSF file ('%s') is corrupted or "
		 "incompatible with this version.\n",psfFilename.c_str());
	}
	return -3;
      }

      //Create a new chunk
      chunks.emplace_back();
      pen_psfMemort::chunk& chunk = chunks.back();
      
      //Extract all particle states for this chunk
      //and obtain the maximum energy in the psf (Emax)
      unsigned long nRead = 0;      
      do{
	pen_psfMemort::particle& p = chunk.particles[nRead];
	nRead = psf.get(p.dhist, p.kpar, p.state);
	if(p.state.E > Emax)
	  Emax = p.state.E;
      }while(nRead != 0);
      
      chunk.nPart = psf.nread();
    }

    //Save final psf maximum energy
    psfEmax = Emax;

    //Get total number of chunks
    size_t nChunks = chunks.size();

    //Check if we have more partitions than chunks
    if(nChunks < (unsigned)npartitions){
      if(verbose > 0){
	printf("psfMemorySource:configure: Error: Number of psf chunks smaller than specified "
	       "partitions. Use less threads or a bigger PSF\n");
      }
      return -13;
    }

    //Calculate chunks per partition and remaining chunks
    unsigned long chunksPerPart = nChunks/npartitions;
    unsigned long offsetChunks = nChunks - chunksPerPart*npartitions;
    

    if(verbose > 1){
      printf("Maximum psf energy (eV):\n");
      printf(" %14.5E\n",Emax);
      printf("PSF partitions:\n");
      printf(" %d\n",npartitions);
      printf("PSF particle chunks:\n");
      printf(" %lu\n",nChunks);
      printf("Offset chunks:\n");
      printf(" %lu\n",offsetChunks);
      printf("NSPLIT:\n");
      printf(" %u\n",NSPLIT);
      printf("Weight window :\n");
      printf(" %15.4E - %15.4E\n",WGHTL,WGHTU);
    }

    //Save pointer to chunks
    pchunks = &chunks;
  }else{
    //Avoid extra memory usage
    chunks.clear();
    chunks.shrink_to_fit();
  }
  
  //Set actual chunk to thread number
  actualChunk = getThread();

  //Init other sampling variables
  actualPart = 0;
  splitted = 0;
  requiredSplits = 0;
  
  return 0;
}

int psfMemory_specificSampler::sharedConfig(const psfMemory_specificSampler& s0){
  
  //Get chunks from thread 0
  pchunks = s0.pchunks;
  psfEmax = s0.psfEmax;
  return 0;
  
}

REGISTER_SPECIFIC_SAMPLER(psfMemory_specificSampler,pen_particleState, MEMORY_PSF)
