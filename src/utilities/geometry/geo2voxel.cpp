
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


#include <cstdio>
#include <vector>
#include "pen_geometries.hh"

int nominalDens(const pen_parserSection& config,
		double densities[constants::MAXMAT],
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
  
  // Create origin geometry
  //************************

  //Get geometry section
  pen_parserSection geometrySection;
  err = config.readSubsection("geometry",geometrySection);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("Error: Configuration 'geometry' section doesn't exist.\n");
    }
    return -2;
  }

  //Get geometry type
  std::string geoType;
  if(geometrySection.read("type",geoType) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'geometry/type' not specified. String expected.\n");
    }
    return -3;
  }  
  
  wrapper_geometry* geometry = nullptr;
  //Create geometry
  geometry = penGeoRegister_create(geoType.c_str());
  if(geometry == nullptr){
    if(verbose > 0){
      printf("Error creating geometry instance of type '%s'\n", geoType.c_str());
    }
    return -4;
  }

  //Configure geometry  
  geometry->name.assign("geometry");    
  geometry->configure(geometrySection,verbose);    

  //Check errors
  if(geometry->configureStatus() != 0){
    if(verbose > 0)
      printf("Error: Fail on geometry configuration.\n");
    return -5;
  }

  // Create voxelized geometry
  //***************************

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

  //Get nominal densities
  double densities[constants::MAXMAT];
  if(nominalDens(config, densities,verbose) != 0){
    return -12;
  }
  
  //Create voxelized geometry
  unsigned* voxMats = nullptr;
  double*   voxDensFact = nullptr;
  
  //Create random voxelized geometry
  unsigned nvox[3] = {(unsigned)nx,(unsigned)ny,(unsigned)nz};
  size_t tvox = nx*ny*nz;
  double sizes[3] = {dx,dy,dz};
  voxMats = (unsigned*) calloc(tvox,sizeof(unsigned));
  voxDensFact = (double*)   malloc(tvox*sizeof(double));

  if(voxMats == nullptr || voxDensFact == nullptr){
    printf("Error: Unable to allocate memory for voxel data\n");
    return -12;
  }

  //Fill voxel matrix
  double subSteps[3] = {dx/(double)granul,dy/(double)granul,dz/(double)granul};
  size_t index = 0;
  const size_t nsubPoints = granul*granul*granul;
  pen_particleState state;
  for(int k = 0; k < nz; ++k)
    for(int j = 0; j < ny; ++j)
      for(int i = 0; i < nx; ++i){

	unsigned pointsPerMat[constants::MAXMAT];
	for(unsigned imat = 0; imat < constants::MAXMAT; ++imat)
	  pointsPerMat[imat] = 0;
	
	double meanDens = 0.0;

	//Set initial state position
	double xInit = (double)i*dx-subSteps[0]*0.5+ox;
	double yInit = (double)j*dy-subSteps[1]*0.5+oy;
	double zInit = (double)k*dz-subSteps[2]*0.5+oz;
	state.X = xInit;
	state.Y = yInit;
	state.Z = zInit;
	for(int k2 = 0; k2 < granul; ++k2){
	  state.Z += subSteps[2];
	  state.Y = yInit;
	  for(int j2 = 0; j2 < granul; ++j2){
	    state.Y += subSteps[1];
	    state.X = xInit;
	    for(int i2 = 0; i2 < granul; ++i2){
	      state.X += subSteps[0];

	      //Locate the state
	      geometry->locate(state);
	      ++pointsPerMat[state.MAT];

	      if(state.MAT > 0){ //Ensure non void material
		if(densities[state.MAT-1] <= 0.0){
		  printf("Error: Missing density for material with index %u \n",state.MAT);
		  free(voxMats);	  
		  free(voxDensFact);
		  return -13;
		}
		meanDens += densities[state.MAT-1];
	      }
	    }
	  }
	}

	//Set material
	unsigned maxPoints = 0;
	for(unsigned imat = 1; imat < constants::MAXMAT; ++imat){
	  //Skyp void, which is a non valid voxel material
	  if(pointsPerMat[imat] > maxPoints){
	    maxPoints = pointsPerMat[imat];
	    voxMats[index] = imat;
	  }
	}

	if(voxMats[index] == 0 || meanDens <= 0.0){
	  printf("Error: Voxels can't be set to void material nor null density!\n");
	  free(voxMats);	  
	  free(voxDensFact);
	  return -14;
	}
	
	voxDensFact[index] = (meanDens/double(nsubPoints))/densities[voxMats[index]-1];
	
	++index;
      }


  //******************************
  // Initialize voxel geometry
  //******************************
  pen_voxelGeo voxelgeo;
  err = voxelgeo.setVoxels(nvox,sizes,voxMats,voxDensFact,3);
  if(err != 0){
    printf("Error using 'setVoxels': %d\n",err);
    free(voxMats);
    free(voxDensFact);
    return -13;
  }  

  //Print ASCII file
  voxelgeo.printImage("voxelGeo.ascii");
  //Dump binary file
  voxelgeo.dump2File("voxelGeo.bin");
  
  free(voxMats);
  free(voxDensFact);
  return 0;
}



int nominalDens(const pen_parserSection& config,
		double densities[constants::MAXMAT],
		unsigned verbose){

  //Get nominal material densities
  //*******************************

  //Initialize densities to -1
  for(unsigned imat = 0; imat < constants::MAXMAT; imat++){
    densities[imat] = -1.0;
  }

  //Read material densities section
  pen_parserSection matSec;
  std::vector<std::string> matNames;
  if(config.readSubsection("materials",matSec) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: No material information provided\n");
    }
    return 4;
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
