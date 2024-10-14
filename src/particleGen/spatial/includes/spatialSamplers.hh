
//
//
//    Copyright (C) 2019-2023 Universitat de València - UV
//    Copyright (C) 2019-2023 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

 
#ifndef __PENRED_SPATIAL_SAMPLERS__
#define __PENRED_SPATIAL_SAMPLERS__

#include "box_spatialSampling.hh"
#include "point_spatialSampling.hh"
#include "image_spatialSampling.hh"
#include "cylinder_spatialSampling.hh"
#include "measure_spatialSampling.hh"

namespace penred{
  namespace sampler{

#ifdef _PEN_USE_DICOM_
    using typesGenericSpatial = std::tuple<box_spatialSampling,
					   point_spatialSampling,
					   image_spatialSampling,
					   cylinder_spatialSampling,
					   measure3D_spatialSampling,
					   measure2D_spatialSampling,
					   measure1D_spatialSampling>;
#else
    using typesGenericSpatial = std::tuple<box_spatialSampling,
					   point_spatialSampling,
					   cylinder_spatialSampling,
					   measure3D_spatialSampling,
					   measure2D_spatialSampling,
					   measure1D_spatialSampling>;    
#endif

    template<size_t T>
    typename std::enable_if<T >= std::tuple_size<typesGenericSpatial>::value, bool>::type
    checkRegistersSpatial(const unsigned){ return true; }
    
    template<size_t T>
    typename std::enable_if<T < std::tuple_size<typesGenericSpatial>::value, bool>::type
    checkRegistersSpatial(const unsigned verbose){
      using tupleType = typename std::tuple_element<T, typesGenericSpatial>::type;
      int val = tupleType::registerStatus();
      if(val != 0){
	if(verbose > 0){
	  printf("Warning: Spatial sampler type '%s' register failed."
		 " Error code: %d\n",
		 tupleType::___ID, val);
	}
	checkRegistersSpatial<T+1>(verbose);
	return false;
      }else{
	return checkRegistersSpatial<T+1>(verbose);
      }
    }
  }
}
    
#endif
