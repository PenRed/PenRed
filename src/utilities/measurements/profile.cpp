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
    printf("usage: %s data-file dimension-to-profile "
	   "1st-dim-index 1st-dim-lower-bin 1st-dim-max-bin ... "
	   "nth-dim-index nth-dim-lower-bin nth-dim-max-bin\n",
	   argv[0]);
    return 1;
  }

  if((argc - 3) % 3 != 0){
    printf("usage: %s data-file dimension-to-profile "
	   "1s5-dim-index 1st-dim-lower-bin 1st-dim-max-bin ... "
	   "nth-dim-index nth-dim-lower-bin nth-dim-max-bin\n",
	   argv[0]);
    return 1;
  }

  //Create a results structure with maximum dimensions
  penred::measurements::results<double, penred::measurements::maxDims> toProfile;
  
  //Read data
  std::ifstream fin(argv[1], std::ifstream::in);
  if(!fin){
    printf("Unable to open file\n");
    return 2;
  }
  
  int err = toProfile.read(fin);
  if(err != 0){
    printf("Error reading data file.\n"
	   "  Error code: %d\n"
	   "  Error message: %s\n",
	   err,
	   penred::measurements::errorToString(err));
    return 3;
  }
  fin.close();

  //Get profile dimension
  int auxProfDim = std::atoi(argv[2]);
  if(auxProfDim < 0){
    printf("Error: Profile dimension must be greater or equal to 0\n");
  }
  if(auxProfDim > static_cast<long int>(penred::measurements::maxDims)){
    printf("Error: Profile dimension isgreater than maximum results dimension (%lu)\n",
	   static_cast<unsigned long>(penred::measurements::maxDims));
  }

  unsigned long profDim = static_cast<unsigned long>(auxProfDim);

  //Generate profile
  penred::measurements::results<double, 1> profile;

  std::vector<std::array<unsigned long, 3>> limits;
  for(int i = 3; i < argc; i += 3){

    int idim = std::atoi(argv[i]);
    if(idim < 0){
      printf("Error: Negative dimension specified for limit number %d\n", i);
      return 4;
    }

    int ilow = std::atoi(argv[i+1]);
    if(ilow < 0){
      printf("Error: Negative low bin limit specified for dimension %d\n", i);
      return 4;
    }

    int itop = std::atoi(argv[i+2]);
    if(itop < 0){
      printf("Error: Negative top bin limit specified for dimension %d\n", i);
      return 4;
    }
    
    limits.push_back({static_cast<unsigned long>(idim),
	static_cast<unsigned long>(ilow),
	static_cast<unsigned long>(itop)});
  }
  
  err = toProfile.profile1D(profDim, limits, profile);
  if(err != 0){
    printf("Error profiling data.\n"
	   "  Error code: %d\n"
	   "  Error message: %s\n",
	   err,
	   penred::measurements::errorToString(err));
    return 5;
  }
  
  FILE* fout = nullptr;
  fout = fopen("profile.dat", "w");
  profile.print(fout, 2, true, false);
  fclose(fout);
  
  return 0;
}
