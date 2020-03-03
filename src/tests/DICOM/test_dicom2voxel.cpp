
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


int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s configuration-file\n",argv[0]);
    return 1;
  }

  unsigned verbose = 3;
  int err;
  
  //**************************
  // Parse configuration
  //**************************

  //Parse configuration file
  pen_parserSection config;
  err = parseFile(argv[1],config);

  printf("Configuration:\n");
  printf("%s\n", config.stringify().c_str());
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    return -1;
  }

  //**************************
  // Configure geometry
  //**************************  
  
  //Create dicom voxelized geometry
  pen_dicomGeo geometry;
  
  err = geometry.configure(config,verbose);
  if(err != 0){
    printf("Error at geometry configuration.\n");
    printf("                 Error code: %d\n",err);
    return -2;
  }

  err = geometry.printImage("image.dat");
  if(err != 0){
    printf("Error printing processed image.\n");
    printf("                 Error code: %d\n",err);
    return -3;
  }
  
  return 0;
}
