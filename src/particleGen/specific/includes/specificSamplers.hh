
//
//
//    Copyright (C) 2019-2022 Universitat de València - UV
//    Copyright (C) 2019-2022 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#ifndef __PENRED_SPECIFIC_SAMPLERS__
#define __PENRED_SPECIFIC_SAMPLERS__

#include "randomState.hh"
#include "PSFsource.hh"
#include "gammaPolarised.hh"
#include "CTsource.hh"
#include "pennuc.hh"
#include "brachySource.hh"
#include "memoryPSFsource.hh"

namespace penred{
  namespace sampler{

#ifdef _PEN_USE_DICOM_
    using typesSpecificCommonState = std::tuple<random_specificSampler,
						psf_specificSampler,
						ct_specificSampler,
						pennuc_specificSampler,
						brachy_specificSampler,
						psfMemory_specificSampler>;
#else
    using typesSpecificCommonState = std::tuple<random_specificSampler,
						psf_specificSampler,
						ct_specificSampler,
						pennuc_specificSampler,
						psfMemory_specificSampler>;    
#endif

    using typesSpecificGPolState = std::tuple<gammaPolarised_specificSampler>;
    

    template<size_t T>
    typename std::enable_if<T >= std::tuple_size<typesSpecificCommonState>::value, bool>::type
    checkRegistersSpecificCommon(const unsigned){ return true; }
    
    template<size_t T>
    typename std::enable_if<T < std::tuple_size<typesSpecificCommonState>::value, bool>::type
    checkRegistersSpecificCommon(const unsigned verbose){
      using tupleType = typename std::tuple_element<T, typesSpecificCommonState>::type;
      int val = tupleType::registerStatus();
      if(val != 0){
	if(verbose > 0){
	  printf("Warning: Specific sampler type '%s' register failed."
		 " Error code: %d\n",
		 tupleType::___ID, val);
	}
	checkRegistersSpecificCommon<T+1>(verbose);
	return false;
      }else{
	return checkRegistersSpecificCommon<T+1>(verbose);
      }
    }
      
      
    template<size_t T>
    typename std::enable_if<T >= std::tuple_size<typesSpecificGPolState>::value, bool>::type
    checkRegistersSpecificGPol(const unsigned){ return true; }
    
    template<size_t T>
    typename std::enable_if<T < std::tuple_size<typesSpecificGPolState>::value, bool>::type
    checkRegistersSpecificGPol(const unsigned verbose){
      using tupleType = typename std::tuple_element<T, typesSpecificGPolState>::type;
      int val = tupleType::registerStatus();
      if(val != 0){
	if(verbose > 0){
	  printf("Warning: Specific sampler type '%s' register failed."
		 " Error code: %d\n",
		 tupleType::___ID, val);
	}
	checkRegistersSpecificGPol<T+1>(verbose);
	return false;
      }else{
	return checkRegistersSpecificGPol<T+1>(verbose);
      }
    }    
  }
}

#endif
