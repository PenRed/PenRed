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
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

#include "materialCreator.hh"
#include <string>

int main (int argc, const char** argv){

  if(argc < 3){
    printf("usage: %s material-list-id filename\n", argv[0]);
    printf("usage: %s composDB material filename\n", argv[0]);
    return 1;
  }

  if(argc == 3){
    unsigned id = std::stoul(argv[1]);
    std::string errorString;
    int err = penred::penMaterialCreator::createMat(id, argv[2], errorString);
    if(err != 0){
      printf ("%s\n", errorString.c_str());
      printf ("IRETRN =%d\n", err);
      return -1;
    }
  }else{
    std::string errorString;
    int err = penred::penMaterialCreator::createMat(argv[2], argv[1], argv[2], errorString, argv[3]);
    if(err != 0){
      printf ("%s\n", errorString.c_str());
      printf ("IRETRN =%d\n", err);
      return -1;
    }    
  }

  return 0;
}
