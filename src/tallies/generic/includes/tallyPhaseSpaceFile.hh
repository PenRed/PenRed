
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


#ifndef __PEN_PHASE_SPACE_FILE_TALLY__
#define __PEN_PHASE_SPACE_FILE_TALLY__

#include <vector>
#include <memory>
#include <mutex>
#include "splittedFile.hh"
#include "pen_phaseSpaceFile.hh"

class pen_tallyPhaseSpaceFile : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_tallyPhaseSpaceFile,pen_particleState)

private:

  static std::vector<
    std::pair<std::string,std::shared_ptr<pen_splittedFile>>
    > splittedFiles;

  static std::mutex SFlock;

  std::shared_ptr<pen_splittedFile> pSF;
  
  int detector;
  int material;

  double emin;
  double emax;

  bool inside;

  pen_psfwriter psf;
  
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
