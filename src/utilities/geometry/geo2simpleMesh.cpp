
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
//    Copyright (C) 2026 Vicent Giménez Alventosa
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


#include <cstdio>
#include <vector>
#include <thread>
#include <atomic>
#include "geoHandler.hh"

int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s config-filename\n",argv[0]);
    return 0;
  }

  unsigned verbose = 3;

  //**************************
  // Parse configuration
  //**************************

  //Parse configuration file
  pen_parserSection config;
  std::string errorLine;
  long unsigned errorLineNum;
  int err = parseFile(argv[1],config,errorLine,errorLineNum);
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }
  
  //Create the geometry handler
  penred::geometry::Handler geoHandler;

  //Configure the geometry using the provided file
  penred::errors::Error error = geoHandler.configure(config, verbose, "geometry");
  if(error){
    printf("%s\n", error.stringify().c_str());
    return -5;
  }else{
    printf("Geometry configured!\n");
  }  

  //Read parameters
  int nx,ny,nz;
  double dx,dy,dz;
  double ox,oy,oz;
  int granul;
  const unsigned maxGranul = 200;
  
  //Number of voxels
  if(config.read("voxelized/nx",nx) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/nx' not specified. integer expected.\n");
    }
    return -6;
  }
  if(config.read("voxelized/ny",ny) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/ny' not specified. integer expected.\n");
    }
    return -6;
  }
  if(config.read("voxelized/nz",nz) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'nz' not specified. integer expected.\n");
    }
    return -6;
  }

  //Voxel sizes
  if(config.read("voxelized/dx",dx) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'dx' not specified. Double expected.\n");
    }
    return -7;
  }
  if(config.read("voxelized/dy",dy) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'dy' not specified. Double expected.\n");
    }
    return -7;
  }
  if(config.read("voxelized/dz",dz) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'dz' not specified. Double expected.\n");
    }
    return -7;
  }

  //Origin
  if(config.read("voxelized/ox",ox) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'ox' not specified. Double expected.\n");
    }
    return -8;
  }
  if(config.read("voxelized/oy",oy) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'oy' not specified. Double expected.\n");
    }
    return -8;
  }
  if(config.read("voxelized/oz",oz) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'oz' not specified. Double expected.\n");
    }
    return -8;
  }  

  double time;
  if(config.read("voxelized/time",time) != INTDATA_SUCCESS){
    time = 0.0;
  }  

  //Read granularity
  if(config.read("voxelized/granularity",granul) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'granularity' not specified. integer expected.\n");
    }
    return -9;
  }


  //Read smooth parameters
  unsigned smoothSteps = 50;
  int smoothStepsi;  
  if(config.read("smooth/steps",smoothStepsi) == INTDATA_SUCCESS){
    if(smoothStepsi >= 0)
      smoothSteps = static_cast<unsigned>(smoothStepsi);
  }

  float smoothFactor = 0.2;
  double smoothFactord;  
  if(config.read("smooth/factor",smoothFactord) == INTDATA_SUCCESS){
    if(smoothFactord > 0)
      smoothFactor = static_cast<float>(smoothFactord);
  }
  
  //Check values
  if( nx < 1 || ny < 1 || nz < 1 ){
    printf("Error: Number of voxels must be, as least, 1 on each axis.\n");
    printf("            nx: %d\n",nx);
    printf("            ny: %d\n",ny);
    printf("            nz: %d\n",nz);
    return -9;
  }

  if( dx < 1.0e-15 || dy < 1.0e-15 || dz < 1.0e-15 ){
    printf("Error: Voxel sizes must be greater than zero.\n");
    printf("            dx: %12.4E\n",dx);
    printf("            dy: %12.4E\n",dy);
    printf("            dz: %12.4E\n",dz);
    return -10;
  }

  if(granul <= 0 || granul > (int)maxGranul){
    printf("Error: Granularity must be, as least, 1 and less than %u\n",maxGranul);
    printf("   granularity: %d\n",granul);
    return -11;
  }

  error = geoHandler.createMesh(nx, ny, nz, dx, dy, dz, ox, oy, oz, time, "geo", granul, smoothSteps, smoothFactor);
  if(error){
    printf("Error creating mesh:\n");
    printf("%s\n", error.stringify().c_str());
    return -5;
  }else{
    printf("Mesh created!\n");
  }
    
  return 0;
}
