
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

 
#ifndef __PENRED_GEOMETRY_LIB__
#define __PENRED_GEOMETRY_LIB__

#include "geometry_classes.hh"
#include "instantiator.hh"

#define DECLARE_GEOMETRY(Class) \
  public: \
  inline int registerStatus() const { return ___register_return;} \
  inline const char* readID() const { return ___ID;}\
  private: \
  static const char* ___ID;\
  static const int ___register_return;

#define REGISTER_GEOMETRY(Class, ID) \
  const int Class::___register_return = penGeoRegister_add<Class>(static_cast<const char *>(#ID)); \
  const char* Class::___ID = static_cast<const char *>(#ID);



instantiator<wrapper_geometry>& ___wrapper_pen_geometry_register();

template <class geometryClass>
int penGeoRegister_add(const char* typeID){
  return ___wrapper_pen_geometry_register().template addSubType<geometryClass>(typeID);
}
wrapper_geometry* penGeoRegister_create(const char* typeID);

inline std::string geometryList(){
  return ___wrapper_pen_geometry_register().typesList();
}
inline std::string geometryList(std::vector<std::string>& list){
  return ___wrapper_pen_geometry_register().typesList(list);
}

#endif
