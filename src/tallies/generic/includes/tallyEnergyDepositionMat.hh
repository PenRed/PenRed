
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


#ifndef __PEN_ENERGY_EDEP_MAT_TALLY__
#define __PEN_ENERGY_EDEP_MAT_TALLY__

#include "pen_constants.hh"

class pen_EdepMat : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_EdepMat,pen_particleState)

private:
    double edptmp[constants::MAXMAT];
    double edep[constants::MAXMAT];
    double edep2[constants::MAXMAT];
    int nmat;
  
public:

  pen_EdepMat() : pen_genericTally( USE_LOCALEDEP |
				    USE_BEGINPART |
				    USE_BEGINHIST |
				    USE_STEP |
				    USE_ENDHIST |
				    USE_MOVE2GEO),
		  nmat(0)
  {}
  
  void tally_localEdep(const double /*nhist*/,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state,
		       const double dE);
  void tally_beginPart(const double /*nhist*/,
		       const unsigned /*kdet*/,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state);
  void tally_beginHist(const double /*nhist*/,
		       const unsigned /*kdet*/,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state);

  void tally_step(const double /*nhist*/,
		  const pen_KPAR /*kpar*/,
		  const pen_particleState& state,
		  const tally_StepData& stepData);
  
  void tally_move2geo(const double /*nhist*/,
		      const unsigned /*kdet*/,
		      const pen_KPAR /*kpar*/,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);
    
  void tally_endHist(const double /*nhist*/);

  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],     
		const pen_parserSection& config, const unsigned verbose);
  void flush();
  void saveData(const double nhist) const;
  int sumTally(const pen_EdepMat& tally);
};






#endif
