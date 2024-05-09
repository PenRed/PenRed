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
  if(argc < 9){
    printf("usage: %s angle dx dy dz mat-index body-name "
	   "parent-name to-be-included\n",
	   argv[0]);
    return 1;
  }

  //Get arguments
  double angle = std::atof(argv[1]);
  double dx = std::atof(argv[2]);
  double dy = std::atof(argv[3]);
  double dz = std::atof(argv[4]);
  int imat = std::atoi(argv[5]);

  int toBeIncluded = 0;
  if(argc >= 5)
    toBeIncluded = std::atoi(argv[8]);

  std::ofstream out("anode.msh", std::ofstream::out);


  penred::xray::createAnode(out, angle, imat,
			    dx,dy,dz,
			    argv[6], argv[7],
			    toBeIncluded != 0);

  out.close();
  
  return 0;
}
