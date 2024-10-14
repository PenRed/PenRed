
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


#ifndef __PENRED_ENERGY_SAMPLERS__
#define __PENRED_ENERGY_SAMPLERS__

#include "intervals_energySampling.hh"
#include "monoenergetic.hh"
#include "fileSpectrum_energySampling.hh"

namespace penred{
  namespace sampler{

    using typesGenericEnergy = std::tuple<intervals_energySampling,
					  monoenergetic,
					  fileSpectrum_energySampling>;


    template<size_t T>
    typename std::enable_if<T >= std::tuple_size<typesGenericEnergy>::value, bool>::type
    checkRegistersEnergy(const unsigned){ return true; }
    
    template<size_t T>
    typename std::enable_if<T < std::tuple_size<typesGenericEnergy>::value, bool>::type
    checkRegistersEnergy(const unsigned verbose){
      using tupleType = typename std::tuple_element<T, typesGenericEnergy>::type;
      int val = tupleType::registerStatus();
      if(val != 0){
	if(verbose > 0){
	  printf("Warning: Energy sampler type '%s' register failed."
		 " Error code: %d\n",
		 tupleType::___ID, val);
	}
	checkRegistersEnergy<T+1>(verbose);
	return false;
      }else{
	return checkRegistersEnergy<T+1>(verbose);
      }
    }
  }
}

#endif
