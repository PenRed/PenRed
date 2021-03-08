//
//
//    Copyright (C) 2020 Universitat de València - UV
//    Copyright (C) 2020 Universitat Politècnica de València - UPV
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

#ifndef __PEN_VR_XRAY_SPLITTING__
#define __PEN_VR_XRAY_SPLITTING__

#include "pen_constants.hh"

class pen_VRxraysplitting : public pen_genericVR<pen_state_gPol>{
  DECLARE_VR(pen_VRxraysplitting)

  private:

  //  ----  x-ray splitting numbers, IXRSPL(IBODY).
  unsigned int IXRSPL[pen_geoconst::NB];
  bool         LXRSPL[pen_geoconst::NB];
  
  public:

  pen_VRxraysplitting() : pen_genericVR(VR_USE_PARTICLESTACK)
  {
    for(unsigned i = 0; i < pen_geoconst::NB; ++i)
      LXRSPL[i] = false;
  }

  int configure(const pen_parserSection& config,
		const wrapper_geometry& geometry,
		const unsigned verbose);
  
  void vr_particleStack(const unsigned long long /*nhist*/,
			const pen_KPAR /*kpar*/,
			const unsigned /*kdet*/,
			pen_state_gPol& state,
			std::array<pen_state_gPol,constants::NMS>& stack,
			unsigned& created,
			const unsigned available,
			pen_rand& random) const;
};

#endif
