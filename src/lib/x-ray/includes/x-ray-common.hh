//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
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
	USING_RESERVED_MATERIAL,
	UNABLE_TO_GENERATE_SPECTRUM,
	ERROR_ON_CONTEXT_INITIALIZATION,
	ERROR_ON_GEOMETRY_INITIALIZATION,
	ERROR_UNABLE_TO_OPEN_FILE,
	ERROR_INVALID_FILE,
	MINIMUM_ENERGY_TOO_LOW,
	BEAM_ENERGY_TOO_HIGH,
	INVALID_ENERGY_RANGE,
	INVALID_FOCAL_SPOT,
	INVALID_LIMITING_ANGLE,
	NOTHING_TO_SIMULATE,
	INVALID_CONFIGURATION,
	UNKNOWN_ERROR,
      };

      constexpr const char* message(const int val){
	switch(val){
	case SUCCESS: return "success";
	case NEGATIVE_DISTANCE: return "negative distance";
	case INVALID_Z: return "invalid Z value";
	case INVALID_SIZE: return "invalid size";
	case INVALID_DISTANCE: return "invalid distance";
	case NO_FILTER_PROVIDED: return "no filter provided";
	case UNABLE_TO_CREATE_MATERIAL: return "unable to create material file";
	case USING_RESERVED_MATERIAL: return "user geometry uses a reserved material index";
	case UNABLE_TO_GENERATE_SPECTRUM: return "unable to generate spectrum";
	case ERROR_ON_CONTEXT_INITIALIZATION: return "error on context initialization. Check context report";
	case ERROR_ON_GEOMETRY_INITIALIZATION: return "error on geometry initialization. Check geometry report, if created, and configured elements distances";
	case ERROR_UNABLE_TO_OPEN_FILE: return "unable to open required file. Check permissions";
	case ERROR_INVALID_FILE: return "invalid file";
	case MINIMUM_ENERGY_TOO_LOW: return "minimum energy too low";
	case BEAM_ENERGY_TOO_HIGH: return "beam energy too high";
	case INVALID_ENERGY_RANGE: return "invalid energy range";
	case INVALID_FOCAL_SPOT: return "invalid focal spot. Must be in the range (0,1] cm";
	case INVALID_LIMITING_ANGLE: return "invalid limiting angle. Must be in the range [0,90) degrees";
	case NOTHING_TO_SIMULATE: return "no histories to be simulated";
	case INVALID_CONFIGURATION: return "invalid configuration format";
	case UNKNOWN_ERROR: return "unknown error";
	default: return "Unknown error";
	}
      }
    }
    
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
