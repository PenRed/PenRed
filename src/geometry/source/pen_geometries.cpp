
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

#include "pen_geometries.hh"

//Geometry register (don't change this line)
#include "pen_geometry_register.cpp"

#include "../objects/source/pen_object_geos.cpp"
#include "../meshes/source/pen_mesh_geos.cpp"

namespace penred{
  namespace geometry{

    bool checkRegisteredMesh(const unsigned verbose){
      return checkRegistersMesh<0>(verbose);
    }

    bool checkRegisteredObj(const unsigned verbose){
      return checkRegistersObject<0>(verbose);
    }    
    
  }
}
