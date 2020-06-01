
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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

 
#ifdef _PEN_USE_DICOM_

#include "DICOM_geo.hh"

int pen_dicomGeo::configure(const pen_parserSection& config, const unsigned verbose){

  int err = 0;

  // Read DICOM directory path
  //****************************
  std::string directoryPath;
  if(config.read("directory",directoryPath) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure:Error: Unable to read field 'directory'. String spected.\n");
    }
    return 1;
  }

  // Check for calibration array
  //*****************************
  pen_parserArray calibration;
  if(config.read("calibration",calibration) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("pen_dicomGeo:configure: No calibration specified, raw image data will be used for conversion.\n");
    }
  }

  // Check for default material and density
  //****************************************
  double defDens;
  int defMat;
  if(config.read("default/material",defMat) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: No default material provided ('default/material')\n");
    }
    defMat = 1;
  }

  printf("DICOM default material: %d\n",defMat);
  if(defMat < 1){
    if(verbose > 0)
      printf("pen_dicomGeo:configure: Error: Default material must be greater than zero\n");
    return 2;
  }

  if(config.read("default/density",defDens) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: No default density provided ('default/density')\n");
    }
    defDens = 0.0012;  //assign air density
  }
  printf("DICOM default density: %14.5E g/cm^3\n",defDens);
  if(defDens <= 0.0){
    if(verbose > 0)
      printf("pen_dicomGeo:configure: Error: Invalid default density.\n");
    return 3;
  }

  //Check for intensity range material assign
  //****************************************
  pen_parserSection intensityRanges;
  std::vector<std::string> intensityRangesNames;

  //Create vectors to store information
  std::vector<unsigned> intensityRangesMat;
  std::vector<double>   intensityRangesLow;
  std::vector<double>   intensityRangesTop;
  std::vector<double>   intensityRangesDens;
  if(config.readSubsection("intensity-ranges",intensityRanges) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf(" No image intensity ranges field ('intensity-ranges') provided to assign materials\n");
    }
  }else{
    //Extract material names
    intensityRanges.ls(intensityRangesNames);  
  }
  
  if(intensityRangesNames.size() > 0){
    if(verbose > 1)
      printf("\nRange name  | MAT ID | density (g/cm^3) | Intensity range\n");
    for(unsigned long i = 0; i < intensityRangesNames.size(); i++){
      //Read material assigned to this contour
      int auxMat;
      double auxIntensityLow;
      double auxIntensityTop;
      double auxIntensityDensity;
      //Create field strings
      std::string matField = intensityRangesNames[i] + std::string("/material");
      std::string intensityLowField = intensityRangesNames[i] + std::string("/low");
      std::string intensityTopField = intensityRangesNames[i] + std::string("/top");
      std::string intensityDensityField = intensityRangesNames[i] + std::string("/density");

      // Material
      //**********
	
      //Read material ID
      if(intensityRanges.read(matField,auxMat) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read material ID for intensity range '%s'. Integer expected.\n",intensityRangesNames[i].c_str());
	}
	return -1;
      }
      //Check material ID
      if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Invalid material ID for intensity range '%s'.\n",intensityRangesNames[i].c_str());
	  printf("                         ID: %d\n",auxMat);
	  printf("Maximum number of materials: %d\n",constants::MAXMAT);
	}
	return -2;
      }

      // Intensity range
      //*****************
	
      //Read low intensity limit
      if(intensityRanges.read(intensityLowField,auxIntensityLow) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read low intensity for range '%s'. Double expected.\n",intensityRangesNames[i].c_str());
	}
	return -3;
      }

      //Read top voxel intensity
      if(intensityRanges.read(intensityTopField,auxIntensityTop) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read top intensity for range '%s'. Double expected.\n",intensityRangesNames[i].c_str());
	}
	return -4;
      }
      
      //Check voxel intensities
      if(auxIntensityLow >= auxIntensityTop){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Invalid intensity range specified for range '%s'.\n",intensityRangesNames[i].c_str());
	  printf("                    low intensity: %12.4E\n",auxIntensityLow);
	  printf("                    top intensity: %12.4E\n",auxIntensityTop);
	}
	return -5;
      }

      // Density (g/cm^3)
      //******************

      //Read low voxel intensity
      if(intensityRanges.read(intensityDensityField,auxIntensityDensity) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read density for intensity range '%s'. Double expected.\n",intensityRangesNames[i].c_str());
	}
	return -3;
      }

      if(verbose > 1)
	printf("%10.10s -> %4d   %12.4E    %12.4E - %12.4E\n",
	       intensityRangesNames[i].c_str(),auxMat,
	       auxIntensityDensity,auxIntensityLow,auxIntensityTop);
      
      //Store values
      intensityRangesMat.push_back(unsigned(auxMat));
      intensityRangesLow.push_back(auxIntensityLow);
      intensityRangesTop.push_back(auxIntensityTop);
      intensityRangesDens.push_back(auxIntensityDensity);
    }
  }
  else if(verbose > 1){
    printf("\nNo image intensity ranges specified to assign voxel materials and densities.\n");
  }

  
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
  
  //Check for contour material assign
  //**********************************
  pen_parserSection contourSec;
  std::vector<std::string> contourNames;

  //Create vectors to store material
  std::vector<unsigned> contourMat;
  std::vector<double>   contourDens;
  std::vector<double>   contourPrio;
  if(config.readSubsection("contours",contourSec) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: No contour information provided\n");
    }
  }else{
    //Extract contour names
    contourSec.ls(contourNames);  
  }

  if(contourNames.size() > 0){
    if(verbose > 1)
      printf("\nContour name   Material  Density  Priority\n");
    for(unsigned i = 0; i < contourNames.size(); i++){
      //Read material assigned to this contour
      int auxMat;
      double auxDens;
      double auxPrio;
      //Create field strings
      std::string matField = contourNames[i] + std::string("/material");
      std::string densField = contourNames[i] + std::string("/density");
      std::string prioField = contourNames[i] + std::string("/priority");

      // Material
      //**********
	
      //Read material ID
      if(contourSec.read(matField.c_str(),auxMat) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read material ID for specified contour '%s'. Integer expected.\n",contourNames[i].c_str());
	}
	return -1;
      }
      //Check material ID
      if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Invalid material ID assigned to contour '%s'.\n",contourNames[i].c_str());
	  printf("                         ID: %d\n",auxMat);
	  printf("Maximum number of materials: %d\n",constants::MAXMAT);
	}
	return -2;
      }

      // Density
      //**********
	
      //Read density
      if(contourSec.read(densField.c_str(),auxDens) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read material density for specified contour '%s'. Double expected.\n",contourNames[i].c_str());
	}
	return -3;
      }
      //Check density
      if(auxDens <= 0.0){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Invalid material density specified for contour '%s'. Must be greater than zero.\n",contourNames[i].c_str());
	  printf("                        density: %12.4E\n",auxDens);
	}
	return -4;
      }

      // Priority
      //**********
	
      //Read density
      if(contourSec.read(prioField.c_str(),auxPrio) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read priority for contour '%s'. Double expected.\n",contourNames[i].c_str());
	}
	return -5;
      }
      
      if(verbose > 1)
	printf("%13.13s -> %4d %12.4E %12.4E\n",contourNames[i].c_str(),auxMat,auxDens,auxPrio);

      contourMat.push_back(unsigned(auxMat));
      contourDens.push_back(auxDens);
      contourPrio.push_back(auxPrio);
    }
  }
  else if(verbose > 1){
    printf("\nNo contours specified to assign material and density.\n");

    //Check if we can assign density using calibration
    if(calibration.size() == 0 && intensityRangesNames.size() == 0){
      printf("pen_dicomGeo:configure: Error: No contour information nor calibration nor intensity ranges provided to assign density to voxels.\n");
      return 4;
    }
  }
  

  //Check for density range material assign
  //****************************************
  pen_parserSection ranges;
  std::vector<std::string> rangesNames;

  //Create vectors to store information
  std::vector<unsigned> rangesMat;
  std::vector<double>   rangesDensLow;
  std::vector<double>   rangesDensTop;
  if(config.readSubsection("ranges",ranges) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: No density ranges provided to assign materials\n");
    }
  }else{
    //Extract material names
    ranges.ls(rangesNames);  
  }
  
  if(rangesNames.size() > 0){
    if(calibration.size() == 0 && verbose > 0){
      printf("\npen_dicomGeo:configure: Warning: Density ranges can't be used to assign materials without calibration. Ranges will be ignored.\n");
    }
    if(verbose > 1)
      printf("\nRange name  | MAT ID | density range (g/cm^3)\n");
    for(unsigned i = 0; i < rangesNames.size(); i++){
      //Read material assigned to this contour
      int auxMat;
      double auxDensLow;
      double auxDensTop;
      //Create field strings
      std::string matField = rangesNames[i] + std::string("/material");
      std::string densLowField = rangesNames[i] + std::string("/density-low");
      std::string densTopField = rangesNames[i] + std::string("/density-top");

      // Material
      //**********
	
      //Read material ID
      if(ranges.read(matField.c_str(),auxMat) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read material ID for density range '%s'. Integer expected.\n",rangesNames[i].c_str());
	}
	return -1;
      }
      //Check material ID
      if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Invalid material ID for density range '%s'.\n",rangesNames[i].c_str());
	  printf("                         ID: %d\n",auxMat);
	  printf("Maximum number of materials: %d\n",constants::MAXMAT);
	}
	return -2;
      }

      // Density range
      //***************
	
      //Read low density
      if(ranges.read(densLowField.c_str(),auxDensLow) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read low density for range '%s'. Double expected.\n",rangesNames[i].c_str());
	}
	return -3;
      }

      //Read top density
      if(ranges.read(densTopField.c_str(),auxDensTop) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read top density for range '%s'. Double expected.\n",rangesNames[i].c_str());
	}
	return -4;
      }
      
      //Check densities
      if(auxDensLow <= 0.0 || auxDensLow >= auxDensTop){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Invalid density range specified for range '%s'.\n",rangesNames[i].c_str());
	  printf("                    low density: %12.4E\n",auxDensLow);
	  printf("                    top density: %12.4E\n",auxDensTop);
	}
	return -5;
      }
      if(verbose > 1)
	printf("%10.10s -> %4d   %12.4E - %12.4E\n",rangesNames[i].c_str(),auxMat,auxDensLow,auxDensTop);

      //Store values
      rangesMat.push_back(unsigned(auxMat));
      rangesDensLow.push_back(auxDensLow);
      rangesDensTop.push_back(auxDensTop);
      
    }
  }
  else if(verbose > 1){
    printf("\nNo density ranges specified to assign materials.\n");
  }
  
  //*******************
  // Try to load DICOM
  //*******************
  err = dicom.loadDicom(directoryPath.c_str(),verbose);
  if(err != PEN_DICOM_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: Error loading DICOM '%s'\n",directoryPath.c_str());
      printf("                 Error code: %d\n",err);
    }
    return 5;
  }

  // Get contour indexes and assign priorities
  //********************************************
  std::vector<int> contourIndexes;
  contourIndexes.resize(dicom.nContours());
  //Initialize contour indexes
  for(unsigned icont = 0; icont < dicom.nContours(); icont++){
    contourIndexes[icont] = -1;
  }
  
  //Get dicom contour index for each configuration contour
  for(unsigned icont = 0; icont < contourNames.size(); icont++){
    unsigned index;
    index = dicom.getContourIndex(contourNames[icont].c_str());
    if(index >= dicom.nContours()){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Contour '%s' doesn't exist in specified DICOM (%s)\n",contourNames[icont].c_str(),directoryPath.c_str());
      }
      return 6;
    }
    //Store configuration index
    contourIndexes[index] = icont;
    //Set priority
    dicom.setContourPriority(index,contourPrio[icont]);
  }

  // Assign contour to each voxel
  //*******************************  
  err = dicom.assignContours();
  if(err != PEN_DICOM_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: Error on DICOM contour assign '%s'\n",directoryPath.c_str());
      printf("                 Error code: %d\n",err);
    }
    return 7;
  }

  if(std::numeric_limits<unsigned>::max() < dicom.getNZ()){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: Error: DICOM '%s' is too large\n"
	     ,directoryPath.c_str());
      printf("                Maximum z planes: %u\n"
	     "                   Read z planes: %lu\n",
	     std::numeric_limits<unsigned>::max(),dicom.getNZ());
    }
    return 8;
  }

  unsigned nvox[3] = {static_cast<unsigned>(dicom.getNX()),
		      static_cast<unsigned>(dicom.getNY()),
		      static_cast<unsigned>(dicom.getNZ())};
  double dvox[3] = {dicom.getDX(),dicom.getDY(),dicom.getDZ()};

  unsigned* mats = nullptr;
  double*   dens = nullptr;

  unsigned long tnvox = dicom.getNVox();
  
  mats = (unsigned*) malloc(sizeof(unsigned)*tnvox);
  dens = (double*)   malloc(sizeof(double)*tnvox);

  if(mats == nullptr || dens == nullptr){
    throw std::bad_alloc();
  }

  //Set materials and densities to default values
  for(unsigned long ivox = 0; ivox < tnvox; ivox++){mats[ivox] = defMat;}
  for(unsigned long ivox = 0; ivox < tnvox; ivox++){dens[ivox] = defDens;}
  
  //Get contour and image information from loaded DICOM
  const double* image = dicom.readImage();
  const int* contours = dicom.readContour();

  // Intensity ranges assign
  //*************************
  
  //Check if intensity ranges has been provided
  if(intensityRangesTop.size() > 0){
    if(verbose > 1)
      printf("Using voxel intensities to assign densities and materials\n");

    for(unsigned long ivox = 0; ivox < tnvox; ++ivox){
      for(unsigned long irange = 0; irange < intensityRangesTop.size(); ++irange)
	if(image[ivox] >= intensityRangesLow[irange] && image[ivox] < intensityRangesTop[irange]){
	  //Voxel in range, assign material and density
	  dens[ivox] = intensityRangesDens[irange];
	  mats[ivox] = intensityRangesMat[irange];
	  break;
	}
    }
  }
  else if(verbose > 1){
    printf("No voxel intensity ranges used to set voxel densities and materials.\n");
  }

  // Density ranges assign
  //************************
  
  //Check for calibrated HU - g/cm^3 convertion
  if(calibration.size() > 0){
    if(verbose > 1)
      printf("Using calibration to assign densities\n");

    //Extract calibration
    std::vector<double> calibrationVect;
    calibrationVect.resize(calibration.size());
    for(unsigned long i = 0; i < calibration.size(); i++){
      double aux;
      err = calibration.read(aux,i);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error on calibration coefficient "
		 "%lu. Number expected.\n",i);
	  printf("                        Error code: %d\n",err);
	}
      }
      calibrationVect[i] = aux;
    }

    //Convert intensity to density
    for(unsigned long ivox = 0; ivox < tnvox; ivox++){
      dens[ivox] = 0.0;
      for(unsigned long i = calibrationVect.size()-1; i > 0; i--){
	dens[ivox] = image[ivox]*(calibrationVect[i] + dens[ivox]);
      }
      dens[ivox] += calibrationVect[0];
    }

    //Assign materials using density ranges
    for(unsigned long ivox = 0; ivox < tnvox; ivox++){
      for(unsigned long irange = 0; irange < rangesMat.size(); irange++){
	if(dens[ivox] >= rangesDensLow[irange] &&
	   dens[ivox] < rangesDensTop[irange]){
	  mats[ivox] = rangesMat[irange];
	  break;
	}
      }
    }
  }
  else{
    printf("No calibration found.\n");
  }

  // Contours assign
  //******************
  
  //Assign materials and indexes using contours
  if(contourNames.size() > 0){
    for(unsigned long ivox = 0; ivox < tnvox; ivox++){
      if(contours[ivox] >= 0){
	int index = contourIndexes[contours[ivox]];
	if(index >= 0){
	  mats[ivox] = contourMat[index];
	  dens[ivox] = contourDens[index];
	}
      }
    }
  }

  //Convert densities to density factor: density(voxel)/density(material)
  for(unsigned long ivox = 0; ivox < tnvox; ivox++){
    if(densities[mats[ivox]-1] < 0.0){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Nominal density not provided for material index %d, which is used in the provided DICOM image.\n",mats[ivox]);
      }
      return 8;
    }
    //Calculate density factor
    dens[ivox] /= densities[mats[ivox]-1];
  }
  
  //Create voxelized geometry
  err = setVoxels(nvox,dvox,mats,dens,verbose);
  if(err != 0){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: Error: Unable to create voxel geometry from DICOM with provided configuration.\n");
      printf("                        Error code: %d\n",err);
    }
    return 9;
  }

  // Check print image option
  //***************************
  bool toASCII = true;
  if(config.read("print-ASCII",toASCII) != INTDATA_SUCCESS){
    toASCII = false;
  }

  if(toASCII){
    printImage("dicomASCII.rep");
    dicom.printSeeds("dicomSeeds.dat");
    if(dicom.nContours() > 0){
      dicom.printContours("dicomContours.dat");
      dicom.printContourVox("dicomContourMask.dat");
    }
  }
  
  
  return 0;
}

int pen_dicomGeo::printImage(const char* filename) const{

  if(filename == nullptr)
    return -1;
  
  //Create a file to store contours data
  FILE* OutVox = nullptr;
  OutVox = fopen(filename,"w");
  if(OutVox == nullptr){
    return -2;
  }

  fprintf(OutVox,"# \n");
  fprintf(OutVox,"# Voxel geometry file\n");
  fprintf(OutVox,"# Nº of voxels (nx,ny,nz):\n");  
  fprintf(OutVox,"# %5u %5u %5u\n",nx,ny,nz);
  fprintf(OutVox,"# Voxel sizes (dx,dy,dz):\n");  
  fprintf(OutVox,"# %8.5E %8.5E %8.5E\n",dx,dy,dz);
  fprintf(OutVox,"# Voxel data:\n");
  fprintf(OutVox,"#    X(cm)   |    Y(cm)   | MAT | density(g/cm^3)\n");

  //Iterate over Z planes
  for(unsigned k = 0; k < nz; k++){
    unsigned long indexZ = nxy*static_cast<unsigned long>(k);
    fprintf(OutVox,"# Index Z = %4d\n",k);

    //Iterate over rows
    for(unsigned j = 0; j < ny; j++){
      unsigned long indexYZ = indexZ +
	static_cast<unsigned long>(j)*static_cast<unsigned long>(nx);
      fprintf(OutVox,"# Index Y = %4d\n",j);

      //Iterate over columns
      for(unsigned i = 0; i < nx; i++){
	unsigned long ivoxel = indexYZ + static_cast<unsigned long>(i);

	//Save voxel X Y and intensity
	fprintf(OutVox," %12.5E %12.5E %4u   %12.5E\n",
		i*dx, j*dy, mesh[ivoxel].MATER,
		densities[mesh[ivoxel].MATER-1]*mesh[ivoxel].densityFact);
      }
      
    }
    //Set a space between planes
    fprintf(OutVox,"\n\n\n");    
  }

  return 0;
}


REGISTER_GEOMETRY(pen_dicomGeo,DICOM)

#endif
