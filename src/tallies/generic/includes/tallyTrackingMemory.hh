
//
//
//    Copyright (C) 2025 Vicent Giménez Alventosa
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
//    
//


#ifndef __PEN_TRACKING_MEMORY_TALLY__
#define __PEN_TRACKING_MEMORY_TALLY__

#include "pen_constants.hh"

class pen_tallyTrackingMemory : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_tallyTrackingMemory,pen_particleState)

private:
  unsigned long long nhists;
  unsigned long long lastHist;
  bool active;
  bool enabledKpar[constants::nParTypes];

  double xlast, ylast, zlast;
  bool sampledToPrint;
  pen_particleState auxState;
  std::vector<int> tracks;
  int nStates;
public:

  static constexpr const double PRECISION = 1000.0;
  static constexpr const int STATESIZE = 4;
  static constexpr const int MAXSTATES = 500000;

  pen_tallyTrackingMemory() : pen_genericTally( USE_BEGINSIM |
						USE_BEGINPART |
						USE_SAMPLEDPART |
						USE_ENDPART |
						USE_STEP |
						USE_JUMP |
						USE_KNOCK |
						USE_ENDSIM |
						USE_LASTHIST),
			      nhists(0),
			      lastHist(0),
			      active(true),
			      sampledToPrint(false),
			      tracks(STATESIZE*MAXSTATES),
			      nStates(0)
  {}
  
  
  void tally_beginSim();

  void tally_endSim(const unsigned long long /*nhist*/);
  
  void tally_beginPart(const unsigned long long /*nhist*/,
		       const unsigned /*kdet*/,
		       const pen_KPAR kpar,
		       const pen_particleState& state); 
  void tally_sampledPart(const unsigned long long nhist,
			 const unsigned long long dhist,
			 const unsigned /*kdet*/,
			 const pen_KPAR kpar,
			 const pen_particleState& state); 
  void tally_endPart(const unsigned long long /*nhist*/,
		     const pen_KPAR kpar,
		     const pen_particleState& state);  
  void tally_step(const unsigned long long /*nhist*/,
		  const pen_KPAR /*kpar*/,
		  const pen_particleState& state,
		  const tally_StepData& stepData);
  void tally_jump(const unsigned long long /*nhist*/,
		const pen_KPAR /*kpar*/,
		const pen_particleState& state,
		const double ds);
  void tally_knock(const unsigned long long /*nhist*/,
		 const pen_KPAR /*kpar*/,
		 const pen_particleState& state,
		 const int icol);

  inline void tally_lastHist(const unsigned long long lasthistIn){
    lastHist = lasthistIn + nhists; 
    active = true;
  }
  
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  void saveData(const unsigned long long /*nhist*/) const;
  void flush();
  int sumTally(const pen_tallyTrackingMemory& /*tally*/);

  inline void encodeState(const pen_particleState& state){
    if(nStates < MAXSTATES){
      int offset = 1 + nStates*STATESIZE;
      tracks[offset] = static_cast<int>(state.E*PRECISION);
      tracks[offset+1] = static_cast<int>(state.X*PRECISION);
      tracks[offset+2] = static_cast<int>(state.Y*PRECISION);
      tracks[offset+3] = static_cast<int>(state.Z*PRECISION);
      ++nStates;
    }
  }

  inline void getResults(std::vector<int>& results) override {
    tracks[0] = nStates;
    results = std::move(tracks);
    tracks.resize(STATESIZE*MAXSTATES);
    nStates = 0;
  }
  
};






#endif
