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
  if(argc < 2){
    printf("usage: %s configuration-file\n\n%s\n",
	   argv[0],
	   pen_readerSection::stringifyObjectSection<penred::xray::readerXRayDeviceCreate>().c_str());
    return 1;
  }

  //Read the configuration file
  pen_parserSection config;
  std::string errorLine;
  unsigned long errorLineNum;
  int err = parseFile(argv[1], config, errorLine, errorLineNum);

  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }

  //Read configuration
  penred::xray::readerXRayDeviceCreate reader;
  err = reader.read(config, 2);
  if(err != penred::xray::readerXRayDeviceCreate::SUCCESS){
    printf("Error: Bad configuration values\n");
    return -2;
  }

  std::ofstream out("device.msh", std::ofstream::out);

  penred::xray::constructDevice(out,config,2);
  out.close();

  if(err != 0){
    printf("Error: Unable to create device geometry\n");
    return -3;
  }
  
  return 0;
}
