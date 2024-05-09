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
#include <cstring>

int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s config-file [--numMats]\n\n %s\n",
	   argv[0],
	   pen_readerSection::stringifyObjectSection<penred::xray::readerXRayDeviceSimulate>().c_str());    
    return 1;
  }

  // ** Parse configuration

  //Parse configuration file
  pen_parserSection config;
  std::string errorLine;
  unsigned long errorLineNum;
  int err = parseFile(argv[1],config,errorLine,errorLineNum);

  //printf("Configuration:\n");
  //printf("%s\n", config.stringify().c_str());
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }

  if(argc > 2){
    if(strcmp(argv[2],"--numMats") == 0){
      // ** Check device sim
      unsigned nMats;
      err = penred::xray::checkSimDevice(config, nMats, 2);
      if(err != 0){
	printf("Error on device configuration.\n");
	return -2;
      }
      printf("Materials used by device geometry: %u\n", nMats);
      return 0;
    }
  }

  // ** Simulate device
  err = penred::xray::simDevice(config, 2);
  if(err != 0){
    printf("Error on device simulation.\n");
    return -2;
  }
  
  return 0;  
}
