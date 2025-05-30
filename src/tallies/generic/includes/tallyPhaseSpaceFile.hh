
//
//
//    Copyright (C) 2019-2025 Universitat de València - UV
//    Copyright (C) 2019-2025 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es
//        vicente.gimenez@uv.es
//    
//


#ifndef __PEN_PHASE_SPACE_FILE_TALLY__
#define __PEN_PHASE_SPACE_FILE_TALLY__

#include <vector>
#include <memory>
#include "splittedFile.hh"
#include "pen_phaseSpaceFile.hh"

class pen_tallyPhaseSpaceFile : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_tallyPhaseSpaceFile,pen_particleState)

private:

  std::shared_ptr<pen_splittedFile> pSF;
  
  int detector;
  int material;

  double emin;
  double emax;

  bool inside;

  pen_psfwriter psf;
  std::array<bool,pen_KPAR::ALWAYS_AT_END> enabledKpars;
  
public:

  pen_tallyPhaseSpaceFile() : pen_genericTally( USE_INTERFCROSS |
						USE_MATCHANGE |
						USE_ENDSIM |
						USE_MOVE2GEO |
						USE_BEGINPART |
						USE_LASTHIST),
			      pSF(nullptr),
			      detector(-1),
			      material(-1),
			      emin(0),
			      emax(1.0e6),
			      inside(false)
  {}

  void tally_lastHist(const unsigned long long lastHist);
  
  void tally_interfCross(const unsigned long long nhist,
			 const unsigned kdet,
		       const pen_KPAR kpar,
		       const pen_particleState& state);

  void tally_matChange(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const unsigned /*prevMat*/);
  
  void tally_move2geo(const unsigned long long nhist,
		      const unsigned kdet,
		      const pen_KPAR kpar,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);

  void tally_beginPart(const unsigned long long /*nhist*/,
		       const unsigned kdet,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& /*state*/);
    
  void tally_endSim(const unsigned long long /*nhist*/);
  
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);

  inline int sharedConfig(const pen_tallyPhaseSpaceFile& tally){
    pSF = tally.pSF;

    //try to create a partition for our thread
    int err = pSF.get()->createPartition(getThread());
    if(err != SPLITTED_FILE_SUCCESS){
      printf("PhaseSpaceFile:sharedConfig: Error: Unable to create a partition for "
	     "tally %s thread %u.\n",readName().c_str(),getThread());
      printf("                 Error code: %d\n",err);
      return -1;
    }
    return 0;
  }

  inline void store(const unsigned long long nhist,
		    const unsigned kpar,
		    const pen_particleState& state){

    //Try to store input state
    if(psf.store(nhist,kpar,state) == -1){
      //Store is full, dump buffers
      psf.dump(1);

      //Write dumped data to file
      pSF.get()->write(getThread(),psf.readBuffer(),psf.memory());
      
      //Clear buffers
      if(psf.clearBuffer() != PEN_PSF_SUCCESS){
	throw std::out_of_range("pen_tallyPhaseSpaceFile:store:Error: History doesn't fit in a single array.");
      }

      //Store input state
      psf.store(nhist,kpar,state);
    }
  }
  
  inline void dump(const bool dumpAll = false){
    //Dump buffers
    psf.dump(1,dumpAll);

    //Write dumped data to file
    pSF.get()->write(getThread(),psf.readBuffer(),psf.memory());

    //Clear buffers
    if(psf.clearBuffer() != PEN_PSF_SUCCESS){
      throw std::out_of_range("pen_tallyPhaseSpaceFile:store:Error: History doesn't fit in a single array.");
    }
  }
  
  void saveData(const unsigned long long /*nhist*/) const;
  void flush();
  int sumTally(const pen_tallyPhaseSpaceFile& /*tally*/);

  ~pen_tallyPhaseSpaceFile(){
  }
};






#endif
