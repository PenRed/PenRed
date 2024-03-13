//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
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
//    
//

#ifndef __PEN_X_RAY_FILTER__
#define __PEN_X_RAY_FILTER__

#include "x-ray-common.hh"
#include <iostream>

namespace penred{

  namespace xray{

    void createBaseFilter(double dx,
			  double dy,
			  double dz,
			  unsigned nVG,
			  std::ofstream& out,
			  const unsigned matIndex,
			  const std::string& filterName,
			  const std::string& parentName,
			  const bool numObjects);
    
  };
};


#endif
