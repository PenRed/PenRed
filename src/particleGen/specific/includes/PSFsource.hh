
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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



#ifndef __PSF_SPECIFIC_SAMPLER__
#define __PSF_SPECIFIC_SAMPLER__

#include <mutex>
#include "pen_phaseSpaceFile.hh"
#include "sharedFile.hh"

class psf_specificSampler : public abc_specificSampler<pen_particleState>{
  DECLARE_SAMPLER(psf_specificSampler)
  private:

  static std::vector<
    std::pair<std::string,std::shared_ptr<pen_sharedFile>>
    > sharedFiles;

  static std::mutex SFlock;

  std::shared_ptr<pen_sharedFile> pSF;

  pen_psfreader psf;
  size_t nChunks;
  size_t chunksPerPart;
  size_t offsetChunks;
  size_t remainingChunks;

  unsigned char* buffer;
  size_t bufferSize;

  unsigned NSPLIT;     //Splitting factor
  double WGHTL, WGHTU; //Weight window [WGHTL,WGHTU]
  double RWGHTL;       //Minimum weight Inverse
  unsigned splitted;   //Number of returned splits from the same particle
  unsigned requiredSplits;
  pen_KPAR lastKpar;
  pen_particleState splitState;

  bool VRR(double& WGHT, pen_rand& random) const ;
  
  public:

  const unsigned MAXsplit = 10000;

  psf_specificSampler() : abc_specificSampler<pen_particleState>(USE_NONE),
			  pSF(nullptr),
			  nChunks(0),
			  chunksPerPart(0),
			  offsetChunks(0),
			  remainingChunks(0),
			  buffer(nullptr),
			  bufferSize(0),
			  NSPLIT(1),
			  WGHTL(0.0),
			  WGHTU(1.0),
			  RWGHTL(1.0e35),
			  splitted(0),
			  requiredSplits(0),
			  lastKpar(ALWAYS_AT_END)
  {
    //Check phase space file chunk size
    if(psf.memory() > 2000000000){
      throw std::out_of_range("Phase space file library uses a binary memory buffer greater than 2000000000B, which is not compatible with this library. Smaller buffer is required.");
    }

    bufferSize = psf.memory();
    buffer = (unsigned char*) malloc(sizeof(unsigned char)*bufferSize);
    if(buffer == nullptr){
      throw std::bad_alloc();
    }
  }

  void skip(const double dhists);
  
  void sample(pen_particleState& state,
	      pen_KPAR& genKpar,
	      unsigned long long& dhist,
	      pen_rand& /*random*/);
  
  int configure(double& Emax,
		const abc_spatialSampler* /*pSpatial*/,
		const abc_directionSampler* /*pDirection*/,
		const abc_energySampler* /*pEnergy*/,
		const abc_timeSampler* /*pTime*/,
		const pen_parserSection& /*config*/,
		const unsigned /*verbose*/);

  ~psf_specificSampler(){
    if(buffer != nullptr){
      free(buffer);
      buffer = nullptr;
    }
  }
};

inline bool psf_specificSampler::VRR(double& WGHT,
				     pen_rand& random) const {
  //  This subroutine applies the Russian roulette technique. PSURV is the
  //  survival probability; when the particle survives, its weight is
  //  increased by a factor 1/PSURV.
  
  //  NOTE: PSURV must be larger than zero and less than one.

  //  Return true if the particle survive and false otherwise

  double PSURV = WGHT*RWGHTL; //Survive probability
  if(random.rand() > PSURV){
    return false;
  }
  else{
    WGHT = WGHTL;
    return true;
  }
}

#endif
