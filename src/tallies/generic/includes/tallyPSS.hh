
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

#ifndef __PEN_PSS_TALLY__
#define __PEN_PSS_TALLY__

#include <algorithm>
#include <array>
#include "pen_constants.hh"

class pen_PSS: public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_PSS,pen_particleState)
  
  private:

  pen_genericTally<pen_particleState>* primary;
  pen_genericTally<pen_particleState>* scatter;
  pen_genericTally<pen_particleState>* multiScatter;
  
  std::array<bool,constants::MAXMAT> sourceMats;

  //Create a structure with bkp information for each stored particle
  struct bkpInfo{
    //bkp value, depending on the number of interactions
    unsigned val;
    //Store the bkp value to consider when the begin particle callback is
    //called. This must be taken into account because particles
    //"extract" energy from the material when begins their simulation.
    //Therefore, the energy must be extracted in the tally where the energy
    //has been deposited, to avoid, for example, negative energy counting.
    unsigned bkp2Extract; 

    bkpInfo() = default;
    bkpInfo(const unsigned valIn,
	    const unsigned bkp2ExtractIn) : val(valIn),
					    bkp2Extract(bkp2ExtractIn){}
  };
  
  struct particleBKP{
    std::array<bkpInfo,constants::NMS> bkp;
    unsigned int n;

    particleBKP() : n(0){}

    inline bkpInfo get(){ return bkp[--n]; }
    inline unsigned int stored() const { return n; }
    inline void store(const unsigned value){
      bkp[n++] = bkpInfo(value,value);
    }
    inline void store(const unsigned value, const unsigned bkp2Extract){
      bkp[n++] = bkpInfo(value,bkp2Extract);
    }
    inline void clear(){ n = 0; }
  };
  
  std::array<particleBKP,constants::nParTypes> bkps;
  pen_KPAR enabledKPAR;

  //Store the actual BKP from the current particle
  unsigned actualBKP;

  //Flag if a knock has been called after last call to edep
  bool lastKnock;
  
public:

  const static std::array<std::string,11> compatibleTallies;
  const static std::array<std::string,6> noCompatibleTallies;
    
  pen_PSS() : pen_genericTally(USE_BEGINSIM |
			       USE_ENDSIM |
			       USE_SAMPLEDPART |
			       USE_ENDHIST |
			       USE_MOVE2GEO |
			       USE_BEGINPART |
			       USE_ENDPART |
			       USE_JUMP |
			       USE_STEP |
			       USE_INTERFCROSS |
			       USE_MATCHANGE |
			       USE_KNOCK |
			       USE_LOCALEDEP |
			       USE_LASTHIST)

  {
    primary = nullptr;
    scatter = nullptr;
    multiScatter = nullptr;     
    std::fill(sourceMats.begin(),
	      sourceMats.end(),
	      false);
    enabledKPAR = PEN_PHOTON;
    actualBKP = 0;
    lastKnock = false;
  }
    
  void clear(void);
  
  void tally_beginSim();
  
  void tally_endSim(const unsigned long long nhist);
  
  void tally_sampledPart(const unsigned long long nhist,
			 const unsigned long long dhist,
			 const unsigned kdet,
			 const pen_KPAR kpar,
			 const pen_particleState& state);
  
  void tally_endHist(const unsigned long long nhist);

  void tally_move2geo(const unsigned long long nhist,
		      const unsigned kdet,
		      const pen_KPAR kpar,
		      const pen_particleState& state,
		      const double dsef,
		      const double dstot);
  
  void tally_beginPart(const unsigned long long nhist,
		       const unsigned kdet,
		       const pen_KPAR kpar,
		       const pen_particleState& state);
  
  void tally_endPart(const unsigned long long nhist,
		     const pen_KPAR kpar,
		     const pen_particleState& state);
  
  void tally_localEdep(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const double dE);
  
  void tally_step(const unsigned long long nhist,
		  const pen_KPAR kpar,
		  const pen_particleState& state,
		  const tally_StepData& stepData);
  
  void tally_interfCross(const unsigned long long nhist,
			 const unsigned kdet,
			 const pen_KPAR kpar,
			 const pen_particleState& state);
  
  void tally_matChange(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const unsigned prevMat);
  
  void tally_jump(const unsigned long long nhist,
		  const pen_KPAR kpar,
		  const pen_particleState& state,
		  const double ds);
  
  void tally_knock(const unsigned long long nhist,
		   const pen_KPAR kpar,
		   const pen_particleState& state,
		   const int icol);

  void tally_lastHist(const unsigned long long lasthist);
    
  int configure(const wrapper_geometry& geometry,
		const abc_material* const materials[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  void saveData(const unsigned long long nhist) const;
  void flush(void);
  int sumTally(const pen_PSS& tally);
  
  ~pen_PSS(){clear();}
  
};


#endif
