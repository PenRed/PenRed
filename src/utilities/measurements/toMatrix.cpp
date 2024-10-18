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

#include "math_classes.hh"
#include <fstream>

int main(int argc, char** argv){

  if(argc < 3){
    printf("usage: %s data-file output-prefix\n", argv[0]);
    return 1;
  }

  //Create a results structure with maximum dimensions
  penred::measurements::results<double, penred::measurements::maxDims> toConvert;
  
  //Read data
  std::ifstream fin(argv[1], std::ifstream::in);
  if(!fin){
    printf("Unable to open file\n");
    return 2;
  }
  
  int err = toConvert.read(fin);
  if(err != 0){
    printf("Error reading data file.\n"
	   "  Error code: %d\n"
	   "  Error message: %s\n",
	   err,
	   penred::measurements::errorToString(err));
    return 3;
  }
  fin.close();

  //Convert to matrix
  return toConvert.printMatrix(argv[2], true);
}
