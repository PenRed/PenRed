
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
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//

 
#ifndef __PENRED_OBJECT_GEOMETRIES__
#define __PENRED_OBJECT_GEOMETRIES__

//Include implemented geometries include files
#include "quadric_geo.hh"
#include "dummy_geo.hh"
#include "meshBody_geo.hh"
#include "combo_geo.hh"
#include "filter_geo.hh"

namespace penred{
  namespace geometry{

    using typesObjectGeos = std::tuple<pen_quadricGeo,
				       pen_dummyGeo,
				       pen_meshBodyGeo,
				       pen_comboGeo,
				       pen_filterGeo>;

    template<size_t T>
    typename std::enable_if<T >= std::tuple_size<typesObjectGeos>::value, bool>::type
    checkRegistersObject(const unsigned){ return true; }
    
    template<size_t T>
    typename std::enable_if<T < std::tuple_size<typesObjectGeos>::value, bool>::type
    checkRegistersObject(const unsigned verbose){
      using tupleType = typename std::tuple_element<T, typesObjectGeos>::type;
      int val = tupleType::registerStatus();
      if(val != 0){
	if(verbose > 0){
	  printf("Warning: Object geometry type '%s' register failed."
		 " Error code: %d\n",
		 tupleType::___ID, val);
	}
	checkRegistersObject<T+1>(verbose);
	return false;
      }else{
	return checkRegistersObject<T+1>(verbose);
      }
    }
  }
}

#endif
