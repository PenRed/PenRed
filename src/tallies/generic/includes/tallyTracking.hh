
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


#ifndef __PEN_TRACKING_TALLY__
#define __PEN_TRACKING_TALLY__

#include "pen_constants.hh"

class pen_tallyTracking : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_tallyTracking,pen_particleState)

private:
  FILE* fout;
  unsigned long long nhists;
  bool active;
public:

  pen_tallyTracking() : pen_genericTally( USE_BEGINSIM |
					  USE_BEGINPART |
					  USE_BEGINHIST |
					  USE_MOVE2GEO |
					  USE_ENDPART |
					  USE_STEP |
					  USE_INTERFCROSS |
					  USE_JUMP |
					  USE_KNOCK),
			fout(nullptr),
			nhists(0),
			active(true)
  {}
  
  
  void tally_beginSim();
  void tally_beginPart(const unsigned long long /*nhist*/,
		       const unsigned /*kdet*/,
		       const pen_KPAR kpar,
		       const pen_particleState& state); 
  void tally_beginHist(const unsigned long long nhist,
		       const unsigned /*kdet*/,
		       const pen_KPAR kpar,
		       const pen_particleState& state); 
  void tally_move2geo(const unsigned long long /*nhist*/,
		      const unsigned /*kdet*/,
		      const pen_KPAR /*kpar*/,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);  
  void tally_endPart(const unsigned long long /*nhist*/,
		     const pen_KPAR /*kpar*/,
		     const pen_particleState& state);  
  void tally_step(const unsigned long long /*nhist*/,
		  const pen_KPAR /*kpar*/,
		  const pen_particleState& state,
		  const tally_StepData& stepData);
  
  void tally_interfCross(const unsigned long long /*nhist*/,
			 const unsigned /*kdet*/,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& /*state*/);
  void tally_jump(const unsigned long long /*nhist*/,
		const pen_KPAR /*kpar*/,
		const pen_particleState& state,
		const double ds);
  void tally_knock(const unsigned long long /*nhist*/,
		 const pen_KPAR /*kpar*/,
		 const pen_particleState& state,
		 const int icol);
  
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  void saveData(const unsigned long long /*nhist*/) const;
  void flush();
  int sumTally(const pen_tallyTracking& /*tally*/);
};






#endif
