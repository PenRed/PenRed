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

#ifndef __PEN_X_RAY_COMMON__
#define __PEN_X_RAY_COMMON__

#include "pen_simulation.hh"

namespace penred{

  namespace xray{

    namespace errors{
      enum errors{
	SUCCESS = 0,
	NEGATIVE_DISTANCE,
	INVALID_Z,
	INVALID_SIZE,
	INVALID_DISTANCE,
	NO_FILTER_PROVIDED,
	UNABLE_TO_CREATE_MATERIAL,
	UNABLE_TO_GENERATE_SPECTRUM,
	ERROR_ON_CONTEXT_INITIALIZATION,
	ERROR_ON_GEOMETRY_INITIALIZATION,
	ERROR_UNABLE_TO_OPEN_FILE,
      };
    };
    
    using detectedPart = simulation::detectedPart;

    struct preloadGeos{
      //Get anode geometry file at compile time
      static constexpr const char* anodeGeoFile = {
        #include "baseAnode.geo"
      };

      static constexpr const char* phantomCylGeoFile = {
        #include "baseCylPhantom.geo"
      };
      
    };

  } // namespace xray
} // namespace penred


#endif
