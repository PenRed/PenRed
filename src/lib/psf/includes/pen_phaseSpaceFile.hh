
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

 
#ifndef __PEN_PHASE_SPACE_FILE__
#define __PEN_PHASE_SPACE_FILE__

#include <cmath>

#include "pen_dump.hh"
#include "pen_states.hh"
#include "pen_constants.hh"

enum psfState{
	      PEN_PSF_SUCCESS = 0,
	      PEN_PSF_INVALID_FILE,
	      PEN_PSF_BAD_READ,
	      PEN_PSF_BAD_DUMP_READ,
	      PEN_PSF_BAD_WRITE,
	      PEN_PSF_BAD_DUMP_WRITE,	      
	      PEN_PSF_HISTORY_OVERFLOW_CHUNK,	      
	      PEN_PSF_UNKNOWN_ERROR
};

class pen_psfwriter{

public:
  static const unsigned baseStateDoubles = 9;
  static const unsigned baseStateInts = 8;
  static const unsigned long BUFFERSIZE = 75000;  
private:

  size_t dumpSize;
  
  double* bufferDoubles;
  long int*  bufferInts;
  unsigned char* bufferBinary;
  unsigned long nStored;

  //Dump handler
  pen_dump ioDump;

  //Save last stored history
  double lastHist;
  //Save stored particles of last history
  unsigned long particlesLastHist;
  //Save number of stored particles from concluded histories
  unsigned long nConcludedPart;
  
public:

  pen_psfwriter() : bufferDoubles(nullptr),
		    bufferInts(nullptr),
		    bufferBinary(nullptr),
		    nStored(0),
		    lastHist(0),
		    particlesLastHist(0),
		    nConcludedPart(0)
  {

    bufferDoubles = (double*) malloc(sizeof(double)*baseStateDoubles*BUFFERSIZE);
    bufferInts    = (long int*) malloc(sizeof(long int)*baseStateInts*BUFFERSIZE);
    
    ioDump.toDump(bufferDoubles  ,baseStateDoubles*BUFFERSIZE);
    ioDump.toDump(bufferInts     ,baseStateInts*BUFFERSIZE   );
    ioDump.toDump(&nConcludedPart,1);

    dumpSize = ioDump.memory();

    bufferBinary = (unsigned char*) malloc(sizeof(unsigned char)*dumpSize);
    
  }

  inline const unsigned char* readBuffer() const {return bufferBinary;}
  inline size_t memory() const {return dumpSize;}
  inline unsigned long bufferSize() const {return BUFFERSIZE;}

  inline double last() const {return lastHist;}
  inline void setLast(const double newLast) {lastHist = newLast;}
  inline unsigned long stored() const {return nStored;}
  inline unsigned long partConcluded() const {return nConcludedPart;}
  inline unsigned long partLastHist() const {return particlesLastHist;}
  
  long int store(const double nhist,
		 const unsigned kpar,
		 const pen_particleState& state);

  inline int dump(unsigned char*& pout,
		  size_t& written,
		  const size_t outputSize,
		  const unsigned verbose,
		  const bool dumpAll = false) {
    
    if(dumpAll){
      nConcludedPart += particlesLastHist;
      int ret = ioDump.dump(pout,written,outputSize,verbose);
      nConcludedPart -= particlesLastHist;
      return ret;      
    }
    else{
      return ioDump.dump(pout,written,outputSize,verbose);
    }
  }

  inline int dump(const unsigned verbose = 0,
		  const bool dumpAll = false){
    if(dumpAll){
      nConcludedPart += particlesLastHist;
      size_t written;
      int ret = ioDump.dump(bufferBinary,written,dumpSize,verbose);
      nConcludedPart -= particlesLastHist;
      return ret;      
    }
    else{
      size_t written;
      return ioDump.dump(bufferBinary,written,dumpSize,verbose);
    }
  }
  
  int write(FILE* fout, const unsigned verbose, const bool dumpAll = false);

  int clearBuffer();
  inline void clear(){nStored = 0;
    lastHist = 0.0;
    particlesLastHist = 0;
    nConcludedPart = 0.0;
  }

  ~pen_psfwriter(){
    free(bufferBinary);
    free(bufferDoubles);
    free(bufferInts);
  }
  
};

class pen_psfreader{
  
public:
  static const unsigned baseStateDoubles = pen_psfwriter::baseStateDoubles;
  static const unsigned baseStateInts = pen_psfwriter::baseStateInts;
  static const unsigned long BUFFERSIZE = pen_psfwriter::BUFFERSIZE;  
private:

  size_t dumpSize;
  
  double* bufferDoubles;
  long int*  bufferInts;
  unsigned char* bufferBinary;
  unsigned long nStored;
  unsigned long nRead;

  //Dump handler
  pen_dump ioDump;
  
public:

  pen_psfreader() : bufferDoubles(nullptr),
		    bufferInts(nullptr),
		    bufferBinary(nullptr),
		    nStored(0),
		    nRead(0)
  {

    bufferDoubles = (double*) malloc(sizeof(double)*baseStateDoubles*BUFFERSIZE);
    bufferInts    = (long int*) malloc(sizeof(long int)*baseStateInts*BUFFERSIZE);
    
    ioDump.toDump(bufferDoubles  ,baseStateDoubles*BUFFERSIZE);
    ioDump.toDump(bufferInts     ,baseStateInts*BUFFERSIZE   );
    ioDump.toDump(&nStored,1);

    dumpSize = ioDump.memory();

    bufferBinary = (unsigned char*) malloc(sizeof(unsigned char)*dumpSize);
    
  }

  inline const unsigned char* readBuffer() const {return bufferBinary;}
  inline size_t memory() const {return dumpSize;}
  inline unsigned long bufferSize() const {return BUFFERSIZE;}

  inline unsigned long stored() const {return nStored;}
  inline unsigned long nread() const {return nRead;}

  double skip(const double nHists2Skip);
  
  unsigned long get(unsigned long& dhist,
		    unsigned& kpar,
		    pen_particleState& state);
  
  inline int read(const unsigned char* const pin,
		  size_t& pos,
		  const unsigned verbose){
    clear();
    return ioDump.read(pin, pos, verbose);
  }

  int read(FILE* fin, const unsigned verbose);
  
  inline void clear(){nStored = 0; nRead = 0;}

  ~pen_psfreader(){
    free(bufferBinary);
    free(bufferDoubles);
    free(bufferInts);
  }
  
};

#endif
