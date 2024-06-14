//
//
//    Copyright (C) 2020-2021 Universitat de València - UV
//    Copyright (C) 2020-2021 Universitat Politècnica de València - UPV
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

#ifndef __PEN_RUSSIAN_ROULETTE__
#define __PEN_RUSSIAN_ROULETTE__

#include "pen_constants.hh"

class pen_VRRussianRoulette : public pen_genericVR<pen_particleState>{
  DECLARE_VR(pen_VRRussianRoulette)

  private:

  //  ----  splitting numbers, IXRSPL(IBODY).
  double DRR[pen_geoconst::NB];
  bool   LRR[pen_geoconst::NB];

  double minWght;
  double maxWght;
  
  public:

  pen_VRRussianRoulette() : pen_genericVR(VR_USE_INTERFCROSS)
  {
    for(unsigned i = 0; i < pen_geoconst::NB; ++i)
      LRR[i] = false;
    minWght = maxWght = 0.0;
  }

  int configure(const pen_parserSection& config,
		const wrapper_geometry& geometry,
		const unsigned verbose);
  
  void vr_interfCross(const unsigned long long /*nhist*/,
		      const pen_KPAR /*kpar*/,
		      const unsigned /*kdet*/,
		      pen_particleState& state,
		      std::vector<pen_particleState>& /*stack*/,
		      unsigned& /*created*/,
		      const unsigned /*available*/,
		      pen_rand& random) const;
};

#endif
