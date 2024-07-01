 
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//
//

#ifndef __PENELOPE_COMPOS_DBS__
#define __PENELOPE_COMPOS_DBS__

#include <memory>
#include "dataBasesCommon.hh"

#include "ICRP_AF.hh"
#include "ICRP_AM.hh"

namespace penred{

  namespace dataBases{

    namespace compositions{

    constexpr const std::array<const char*, 2> dataBases{
      "ICRP_AF",
      "ICRP_AM",      
    };

    inline size_t indexDB(const std::string& name){
      for(size_t i = 0; i < dataBases.size(); ++i){
	if(name.compare(dataBases[i]) == 0){
	  return i;
	}
      }      
      return dataBases.size();
    }

    inline std::shared_ptr<materials> getDB(const unsigned index){

      switch(index){
      case 0: return std::make_shared<ICRP::AF>();
      case 1: return std::make_shared<ICRP::AM>();
      default: return std::shared_ptr<materials>();
      };      
    }

    inline std::shared_ptr<materials> getDB(const std::string& name){
      return getDB(indexDB(name));
    }

    } //namespace compositions
    
  } //namespace dataBases
  
} //namespace penred

#endif
