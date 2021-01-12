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

#ifndef __PEN_RADIAL_SPLITTING__
#define __PEN_RADIAL_SPLITTING__

#include "pen_constants.hh"

class pen_VRradialSplitting : public pen_genericVR<pen_particleState>{
  DECLARE_VR(pen_VRradialSplitting)

  private:

  double x0,y0,z0;
  double rmin,rmax;
  double rmin2,rmax2;
  
  double minWght;
  double maxWght;

  bool linear;
  
  public:

  pen_VRradialSplitting() : pen_genericVR(VR_USE_PARTICLESTACK),
			    minWght(0), maxWght(0), linear(false)
  {}

  int configure(const pen_parserSection& config,
		const wrapper_geometry& geometry,
		const unsigned verbose);
  
  void vr_particleStack(const unsigned long long /*nhist*/,
			const pen_KPAR /*kpar*/,
			const unsigned /*kdet*/,
			pen_particleState& state,
			std::array<pen_particleState,constants::NMS>& stack,
			unsigned& created,
			const unsigned available,
			pen_rand& /*random*/) const;
};

#endif
