
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
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
//        vicente.gimenez@uv.es
//    
//


#include <cstdio>
#include <vector>
#include "geoHandler.hh"

int nominalDens(const pen_parserSection& config,
		std::array<double, constants::MAXMAT>& densities,
		unsigned verbose = 0);

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
  
  //Get nominal densities
  std::array<double, constants::MAXMAT> densities;
  if(nominalDens(config, densities,verbose) != 0){
    printf("Error reading material densities!\n");
    return -12;
  }
  
  //Create the geometry handler
  penred::geometry::Handler geoHandler;

  //Set material densities
  geoHandler.setDensities(densities);

  //Configure the geometry using the provided file
  penred::errors::Error error = geoHandler.configFromFile(argv[1], verbose, "geometry");
  if(error){
    printf("%s\n", error.stringify().c_str());
    return -5;
  }else{
    printf("Geometry configured!\n");
  }

  // Create voxelized geometry information
  //****************************************

  //Read voxelized geometry parameters
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
      printf("createGeometry: Error: field 'voxelized/nz' not specified. integer expected.\n");
    }
    return -6;
  }

  //Voxel sizes
  if(config.read("voxelized/dx",dx) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/dx' not specified. Double expected.\n");
    }
    return -7;
  }
  if(config.read("voxelized/dy",dy) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/dy' not specified. Double expected.\n");
    }
    return -7;
  }
  if(config.read("voxelized/dz",dz) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/dz' not specified. Double expected.\n");
    }
    return -7;
  }

  //Origin
  if(config.read("voxelized/ox",ox) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/ox' not specified. Double expected.\n");
    }
    return -8;
  }
  if(config.read("voxelized/oy",oy) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/oy' not specified. Double expected.\n");
    }
    return -8;
  }
  if(config.read("voxelized/oz",oz) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/oz' not specified. Double expected.\n");
    }
    return -8;
  }  

  double time;
  if(config.read("voxelized/time",oz) != INTDATA_SUCCESS){
    time = 0.0;
  }  

  //Read granularity
  if(config.read("voxelized/granularity",granul) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'voxelized/granularity' not specified. integer expected.\n");
    }
    return -9;
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

  //******************************
  // Create voxel geometry
  //******************************
  
  std::shared_ptr<pen_voxelGeo> voxelGeo;
  error = geoHandler.createVoxelized(voxelGeo, nx, ny, nz, dx, dy, dz, ox, oy, oz, time, granul);
  if(error){
    printf("%s\n", error.stringify().c_str());
    return -2;
  }

  //Print ASCII file
  voxelGeo->saveASCII("voxelGeo.ascii");
  //Dump binary file
  voxelGeo->dump2File("voxelGeo.bin");

  //Print plotable data
  voxelGeo->printImage("voxelGeo.plot");
  
  return 0;
}



int nominalDens(const pen_parserSection& config,
		std::array<double, constants::MAXMAT>& densities,
		unsigned verbose){

  //Get nominal material densities
  //*******************************

  //Initialize densities to 1.0
  for(unsigned imat = 0; imat < constants::MAXMAT; imat++){
    densities[imat] = 1.1;
  }

  //Read material densities section
  pen_parserSection matSec;
  std::vector<std::string> matNames;
  if(config.readSubsection("materials",matSec) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("Warning: No material information provided. All materials will be "
	     "considered with density 1.0\n");
    }
    return 0;
  }
  
  //Extract material names
  matSec.ls(matNames);  

  if(verbose > 1){
    printf(" Material Name |  ID  | Density (g/cm^3)\n");
  }
  
  //Iterate over all material
  for(unsigned imat = 0; imat < matNames.size(); imat++){

    double auxDens;
    int auxID;

    std::string idField = matNames[imat] + std::string("/ID");
    std::string densField = matNames[imat] + std::string("/density");

    // Material
    //**********
	
    //Read material ID
    if(matSec.read(idField.c_str(),auxID) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Unable to read material ID for material '%s'. Integer expected.\n",matNames[imat].c_str());
      }
      return -1;
    }
    //Check material ID
    if(auxID < 1 || auxID > (int)constants::MAXMAT){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Invalid ID specified for material '%s'.\n",matNames[imat].c_str());
	printf("                         ID: %d\n",auxID);
	printf("Maximum number of materials: %d\n",constants::MAXMAT);
      }
      return -2;
    }

    // Density
    //**********
	
    //Read density
    if(matSec.read(densField.c_str(),auxDens) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Unable to read density for material '%s'. Double expected.\n",matNames[imat].c_str());
      }
      return -3;
    }
    //Check density
    if(auxDens <= 0.0){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Invalid density specified for material '%s'. Must be greater than zero.\n",matNames[imat].c_str());
	printf("                        density: %12.4E g/cm^3\n",auxDens);
      }
      return -4;
    }

    //Store density values for specified materials
    densities[auxID-1] = auxDens;

    if(verbose > 1){
      printf("%15.15s  %4d  %12.5E\n",matNames[imat].c_str(),auxID,auxDens);
    }
    
  }
  return 0;
}
