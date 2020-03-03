
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


#include "pen_phaseSpaceFile.hh"

long int pen_psfwriter::store(const double nhist,
			      const unsigned kpar,
			      const pen_particleState& state){
  
  if(nStored >= BUFFERSIZE){
    return -1;
  }
  if(nhist < lastHist){
    printf("pen_psf:store:Error: Trying to store a particle of concluded history\n");
    return -2;
  }

  unsigned long indexD = baseStateDoubles*nStored;
  unsigned long indexI = baseStateInts*nStored;

  //Store state data
  bufferDoubles[indexD  ] = state.E;
  bufferDoubles[indexD+1] = state.X;
  bufferDoubles[indexD+2] = state.Y;
  bufferDoubles[indexD+3] = state.Z;
  bufferDoubles[indexD+4] = state.U;
  bufferDoubles[indexD+5] = state.V;
  bufferDoubles[indexD+6] = state.W;
  bufferDoubles[indexD+7] = state.WGHT;
  bufferDoubles[indexD+8] = state.PAGE;

  long int increment = (unsigned long)(nhist-lastHist+0.5);
  bufferInts[indexI  ] = increment;
  bufferInts[indexI+1] = kpar;
  bufferInts[indexI+2] = state.LAGE;
  //memcpy(&bufferInts[indexI+3],state.ILB,pen_particleState::ILBsize);
  
  bufferInts[indexI+3] = state.ILB[0];
  bufferInts[indexI+4] = state.ILB[1];
  bufferInts[indexI+5] = state.ILB[2];
  bufferInts[indexI+6] = state.ILB[3];
  bufferInts[indexI+7] = state.ILB[4];  
  
  if(increment > 0){
    lastHist = nhist;
    nConcludedPart += particlesLastHist;
    particlesLastHist = 0;
  }
  
  particlesLastHist++;  
  nStored++;
  
  return nStored;
}

int pen_psfwriter::clearBuffer(){

  //Move particles from incomplete history to the buffer beginning
  if(BUFFERSIZE == particlesLastHist){
    return PEN_PSF_HISTORY_OVERFLOW_CHUNK;
  }

  if(particlesLastHist > 0){
  
    size_t first2store = BUFFERSIZE-particlesLastHist;
    size_t doubleInit = first2store*baseStateDoubles;
    size_t intInit = first2store*baseStateInts;

    size_t B2copyDouble = sizeof(double)*baseStateDoubles*particlesLastHist;
    size_t B2copyInt = sizeof(long int)*baseStateInts*particlesLastHist;

    //Copy remaining particles at the beginning of the buffer
    memmove(&bufferDoubles[0],&bufferDoubles[doubleInit],B2copyDouble);
    memmove(&bufferInts[0],&bufferInts[intInit],B2copyInt);

  }
  
  nStored = particlesLastHist;
  nConcludedPart = 0;

  return PEN_PSF_SUCCESS;
}

int pen_psfwriter::write(FILE* fout,
			 const unsigned verbose,
			 const bool dumpAll){

  if(fout == nullptr)
    return PEN_PSF_INVALID_FILE;

  if(dump(verbose,dumpAll) != PEN_DUMP_SUCCESS){
    return PEN_PSF_BAD_DUMP_WRITE;
  }

  if(fwrite(bufferBinary,sizeof(unsigned char), dumpSize, fout) != dumpSize){
    return PEN_PSF_BAD_WRITE;
  }
  
  return PEN_PSF_SUCCESS;
}

double pen_psfreader::skip(const double nHists2Skip){

  double toSkip = floor(nHists2Skip);
  if(toSkip < 0.5)
    return 0.0; //Nothing to do

  //Iterate until specified number of histories or
  //the entire buffer have been skipped
  double skipped = -0.5; // -0.5 to avoid rounding errors
  while(nRead < nStored){

    unsigned long indexI = baseStateInts*nRead;
    unsigned long dhist = (unsigned long)bufferInts[indexI];

    skipped += (double)dhist;
    if(skipped > toSkip){
      //Objective reached, doesn't skip current particle
      return toSkip;
    }
    //Skip current particle
    nRead++;    
  }

  //End of buffer
  return ceil(skipped);  
  
}

unsigned long pen_psfreader::get(unsigned long& dhist,
				 unsigned& kpar,
				 pen_particleState& state){
  
  if(nRead >= nStored){
    return 0;
  }

  unsigned long indexD = baseStateDoubles*nRead;
  unsigned long indexI = baseStateInts*nRead;

  //Store state data
  state.E = bufferDoubles[indexD  ];
  state.X = bufferDoubles[indexD+1];
  state.Y = bufferDoubles[indexD+2];
  state.Z = bufferDoubles[indexD+3];
  state.U = bufferDoubles[indexD+4];
  state.V = bufferDoubles[indexD+5];
  state.W = bufferDoubles[indexD+6];
  state.WGHT = bufferDoubles[indexD+7];
  state.PAGE = bufferDoubles[indexD+8];

  dhist = (unsigned long)bufferInts[indexI];
  kpar  = (unsigned)bufferInts[indexI+1];
  state.LAGE = bufferInts[indexI+2];
  //memcpy(state.ILB,&bufferInts[indexI+3],pen_particleState::ILBsize);
  
  state.ILB[0] = bufferInts[indexI+3];
  state.ILB[1] = bufferInts[indexI+4];
  state.ILB[2] = bufferInts[indexI+5];
  state.ILB[3] = bufferInts[indexI+6];
  state.ILB[4] = bufferInts[indexI+7];
  
  nRead++;
  return nRead;
}

int pen_psfreader::read(FILE* fin, const unsigned verbose){

  clear();
  
  if(fin == nullptr)
    return PEN_PSF_INVALID_FILE;

  if(fread(bufferBinary,sizeof(unsigned char), dumpSize, fin) != dumpSize){
    return PEN_PSF_BAD_READ;
  }

  size_t pos = 0;
  if(ioDump.read(bufferBinary, pos, verbose) != PEN_DUMP_SUCCESS){
    return PEN_PSF_BAD_DUMP_READ;
  }

  return PEN_PSF_SUCCESS;
}
