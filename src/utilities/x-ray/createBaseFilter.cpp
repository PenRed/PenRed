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

#include "x-ray.hh"


int main(int argc, char** argv){

  //Check arguments
  if(argc < 5){
    printf("usage: %s dx dy dz n-x-vergex-groups", argv[0]);
    return 1;
  }

  //Get arguments
  double dx = std::atof(argv[1]);
  double dy = std::atof(argv[2]);
  double dz = std::atof(argv[3]);
  unsigned nVG = std::atoi(argv[4]);

  std::ofstream out("filter.msh", std::ofstream::out);
    
  penred::xray::createBaseFilter(dx,dy,dz,nVG,out,1,"filter","void",true);

  out.close();
  
  return 0;
}
