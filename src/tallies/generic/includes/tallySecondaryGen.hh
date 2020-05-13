
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


#ifndef __PEN_SECONDARY_GEN_TALLY__
#define __PEN_SECONDARY_GEN_TALLY__

#include "pen_constants.hh"

class pen_tallySecondary : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_tallySecondary,pen_particleState)

  private:
  double absSec[constants::nParTypes],
    absSec2[constants::nParTypes],
    absSecTmp[constants::nParTypes];

  double upBoundSec[constants::nParTypes],
    upBoundSec2[constants::nParTypes],
    upBoundSecTmp[constants::nParTypes];

  double downBoundSec[constants::nParTypes],
    downBoundSec2[constants::nParTypes],
    downBoundSecTmp[constants::nParTypes];

  double absPrim[constants::nParTypes],
    absPrim2[constants::nParTypes],
    absPrimTmp[constants::nParTypes];

  double upBoundPrim[constants::nParTypes],
    upBoundPrim2[constants::nParTypes],
    upBoundPrimTmp[constants::nParTypes];

  double downBoundPrim[constants::nParTypes],
    downBoundPrim2[constants::nParTypes],
    downBoundPrimTmp[constants::nParTypes];    
public:

  pen_tallySecondary() : pen_genericTally(USE_ENDPART |
					  USE_ENDHIST)
  {}

  void tally_endPart(const unsigned long long /*nhist*/,
		     const pen_KPAR kpar,
		     const pen_particleState& state);  
  
  void tally_endHist(const unsigned long long /*nhist*/);

  
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],
		const pen_parserSection& /*config*/,
		const unsigned verbose);
  
  void saveData(const unsigned long long nhist) const;
  void flush();
  int sumTally(const pen_tallySecondary& tally);
};

#endif
