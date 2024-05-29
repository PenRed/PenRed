
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



#ifndef __PSF_MEMORY_SPECIFIC_SAMPLER__
#define __PSF_MEMORY_SPECIFIC_SAMPLER__

#include "pen_phaseSpaceFile.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace pen_psfMemort{

  struct particle{
    unsigned long dhist;
    unsigned kpar;
    pen_particleState state;
  };

  struct chunk{
    std::array<pen_psfMemort::particle,pen_psfreader::BUFFERSIZE> particles;
    size_t nPart;
  };
  
} // namespace pen_psfMemort

class psfMemory_specificSampler : public abc_specificSampler<pen_particleState>{
  DECLARE_SPECIFIC_SAMPLER(psfMemory_specificSampler, pen_particleState)
  private:

  std::vector<pen_psfMemort::chunk> chunks;
  const std::vector<pen_psfMemort::chunk>* pchunks;
  size_t actualPart; //Actual particle index inside the actual chunk
  size_t actualChunk; //Actual chunk to read from
  unsigned splitted;  //Number of returned splits from the actual particle
  unsigned requiredSplits; //Number of required splits for the actual partice
  pen_particleState splitState;

  unsigned NSPLIT;     //Splitting factor
  double WGHTL, WGHTU; //Weight window [WGHTL,WGHTU]
  double RWGHTL;       //Minimum weight Inverse
  double particleRot[9];
  double dx, dy, dz;

  unsigned npartitions;
  
  bool rotation = false;
  bool VRR(double& WGHT, pen_rand& random) const ;

  double psfEmax;
  
  public:

  const unsigned MAXsplit = 10000;

  psfMemory_specificSampler() : abc_specificSampler<pen_particleState>(USE_NONE),
				actualPart(0),
				actualChunk(0),
				splitted(0),
				requiredSplits(0),
				NSPLIT(1),
				WGHTL(0.0),
				WGHTU(1.0),
				RWGHTL(1.0e35),
				npartitions(0){}

  inline int partitions() const {return npartitions;}
  inline double emax() const {return psfEmax;}
  inline void reset() {
    splitted = 0;
    requiredSplits = 0;
    actualChunk = getThread();
    actualPart = 0;
  }
  
  void skip(const unsigned long long dhists);
  
  void sample(pen_particleState& state,
	      pen_KPAR& genKpar,
	      unsigned long long& dhist,
	      pen_rand& random);
  
  int configure(double& Emax,
		const pen_parserSection& config,
		const unsigned nthreads,
		const unsigned verbose);

  int sharedConfig(const psfMemory_specificSampler&);

};

inline bool psfMemory_specificSampler::VRR(double& WGHT,
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
