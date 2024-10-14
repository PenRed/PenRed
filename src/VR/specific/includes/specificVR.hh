
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#ifndef __PEN_SPECIFIC_VR__ 
#define __PEN_SPECIFIC_VR__ 

#include "VRxraySplitting.hh"

namespace penred{
  namespace vr{

    using typesSpecificVR = std::tuple<pen_VRxraysplitting>;    

    template<size_t T>
    typename std::enable_if<(T == std::tuple_size<typesSpecificVR>::value), bool>::type
    checkRegistersSpecific(const unsigned){ return true; }
    
    template<size_t T>
    typename std::enable_if<(T < std::tuple_size<typesSpecificVR>::value), bool>::type
      checkRegistersSpecific(const unsigned verbose){
      using tupleType = typename std::tuple_element<T, typesSpecificVR>::type;
      int val = tupleType::registerStatus();
      if(val != 0){
	if(verbose > 0){
	  printf("Warning: Specific VR type '%s' register failed."
		 " Error code: %d\n",
		 tupleType::___ID, val);
	}
	checkRegistersSpecific<T+1>(verbose);
	return false;
      }else{
	return checkRegistersSpecific<T+1>(verbose);
      }
    }
    
    template<>
    bool checkRegistered<pen_state_gPol>(const unsigned verbose);
  }
}

#endif
