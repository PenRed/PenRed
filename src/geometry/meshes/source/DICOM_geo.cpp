
//
//
//    Copyright (C) 2019-2022 Universitat de València - UV
//    Copyright (C) 2019-2022 Universitat Politècnica de València - UPV
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
      printf("pen_dicomGeo:configure:Error: Unable to read field 'directory'. String expected.\n");
    }
    configStatus = 1;
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
    configStatus = 2;
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
    configStatus = 3;
    return 3;
  }

  //Check for intensity range material assign
  //****************************************

  std::vector<intensityRange> intensityRanges;
  if(readIntensityRanges(config, intensityRanges, verbose) != 0){
    if(verbose > 0)
      printf("pen_dicomGeo:configure: Error reading intensity ranges.\n");
    return 4;
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
    configStatus = 4;
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
      configStatus = -1;
      return -1;
    }
    //Check material ID
    if(auxID < 1 || auxID > (int)constants::MAXMAT){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Invalid ID specified for material '%s'.\n",matNames[imat].c_str());
	printf("                         ID: %d\n",auxID);
	printf("Maximum number of materials: %d\n",constants::MAXMAT);
      }
      configStatus = -2;
      return -2;
    }

    // Density
    //**********
	
    //Read density
    if(matSec.read(densField.c_str(),auxDens) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Unable to read density for material '%s'. Double expected.\n",matNames[imat].c_str());
      }
      configStatus = -3;
      return -3;
    }
    //Check density
    if(auxDens <= 0.0){
      if(verbose > 0){
	printf("pen_dicomGeo:configure: Error: Invalid density specified for material '%s'. Must be greater than zero.\n",matNames[imat].c_str());
	printf("                        density: %12.4E g/cm^3\n",auxDens);
      }
      configStatus = -4;
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

  //Create vectors to store contour assign information
  std::vector<contourAssign> contourAssigns;
  
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
	auxMat = 0;
      }else{
	//Check material ID
	if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	  if(verbose > 0){
	    printf("pen_dicomGeo:configure: Error: Invalid material ID assigned to contour '%s'.\n",contourNames[i].c_str());
	    printf("                         ID: %d\n",auxMat);
	    printf("Maximum number of materials: %d\n",constants::MAXMAT);
	  }
	  configStatus = -2;
	  return -2;
	}
      }

      // Density
      //**********
	
      //Read density
      if(contourSec.read(densField.c_str(),auxDens) != INTDATA_SUCCESS){
	auxDens = -1.0;
      }else{
	//Check density
	if(auxDens <= 0.0){
	  if(verbose > 0){
	    printf("pen_dicomGeo:configure: Error: Invalid material density specified for contour '%s'. Must be greater than zero.\n",contourNames[i].c_str());
	    printf("                        density: %12.4E\n",auxDens);
	  }
	  configStatus = -4;
	  return -4;
	}
      }
      // Priority
      //**********
	
      //Read priority
      if(contourSec.read(prioField.c_str(),auxPrio) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:configure: Error: Unable to read priority for contour '%s'. Double expected.\n",contourNames[i].c_str());
	}
	configStatus = -5;
	return -5;
      }
      
      if(verbose > 1){
	//Print contour name
	printf("%13.13s -> ",contourNames[i].c_str());

	//Print contour material
	if(auxMat == 0)
	  printf("none ");
	else
	  printf("%4d ",auxMat);

	//Print contour density
	if(std::signbit(auxDens))
	  printf("none ");
	else
	  printf("%12.4E ",auxDens);

	//Print priority
	printf("%12.4E\n",auxPrio);	
      }

      contourAssigns.emplace_back();
      contourAssign& contAssigns = contourAssigns.back();

      contAssigns.defaultMat = static_cast<unsigned>(auxMat);
      contAssigns.defaultDens = auxDens;
      contAssigns.priority = auxPrio;

      //Load intensity and density ranges
      if(readIntensityRanges(config, contAssigns.intensityRanges, verbose) != 0){
	if(verbose > 0)
	  printf("pen_dicomGeo:configure: Error reading intensity"
		 " ranges in contour %s.\n", contourNames[i].c_str());
	configStatus = -5;
	return -5;
      }

      if(readDensityRanges(config, contAssigns.densityRanges, verbose) != 0){
	if(verbose > 0)
	  printf("pen_dicomGeo:configure: Error reading density"
		 " ranges in contour %s.\n", contourNames[i].c_str());
	configStatus = -5;
	return -5;
      }

    }
  }
  else if(verbose > 1){
    printf("\nNo contours specified to assign material and density.\n");

    //Check if we can assign density using calibration
    if(calibration.size() == 0 && intensityRanges.size() == 0){
      printf("pen_dicomGeo:configure: Error: No contour information nor calibration nor intensity ranges provided to assign density to voxels.\n");
      configStatus = 4;
      return 4;
    }
  }
  

  //Check for density range material assign
  //****************************************

  std::vector<densityRange> densityRanges;
  if(readDensityRanges(config, densityRanges, verbose) != 0){
    if(verbose > 0)
      printf("pen_dicomGeo:configure: Error reading density ranges.\n");
    return 4;
  }

  if(densityRanges.size() > 0 && calibration.size() == 0 && verbose > 0){
    printf("\npen_dicomGeo:configure: Warning: Density ranges can't be "
	   "used to assign materials without calibration. Density ranges "
	   "will be ignored.\n");
  }
  
  //Try to read print-ASCII output dir
  std::string OutputDirPath;
  if(config.read("outputdir",OutputDirPath) != INTDATA_SUCCESS){
      if(verbose > 0){
  printf("pen_dicomGeo: configure: Warning: unable to read field geometry/outputdir. Assumed default output dir path ./\n");
      }
  }

  if(verbose > 1){
      printf("geometry print-ASCII outputdir: '%s'\n\n",OutputDirPath.c_str());
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
    configStatus = 5;
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
      configStatus = 6;
      return 6;
    }
    //Store configuration index
    contourIndexes[index] = icont;
    //Set priority
    dicom.setContourPriority(index,contourAssigns[icont].priority);
  }

  // Assign contour to each voxel
  //*******************************  
  err = dicom.assignContours();
  if(err != PEN_DICOM_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: Error on DICOM contour assign '%s'\n",directoryPath.c_str());
      printf("                 Error code: %d\n",err);
    }
    configStatus = 7;
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
    configStatus = 8;
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
  if(intensityRanges.size() > 0){
    if(verbose > 1)
      printf("Using voxel intensities to assign densities and materials\n");

    for(unsigned long ivox = 0; ivox < tnvox; ++ivox){
      for(const intensityRange& range : intensityRanges){
	if(range.inner(image[ivox])){
	  //Voxel in range, assign material and density
	  dens[ivox] = range.dens;
	  mats[ivox] = range.mat;
	  break;	  
	}
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
      for(const densityRange& range : densityRanges){
	if(range.inner(dens[ivox])){
	  //Voxel in range, assign material
	  mats[ivox] = range.mat;
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

	  //Get contour assign structure
	  const contourAssign& cAssign = contourAssigns[index];

	  //Flag if the voxel has been assigned by contour data
	  bool assigned = false;

	  //Try to assign material and density with
	  //specific intensity ranges for this contour
	  for(const intensityRange& range : cAssign.intensityRanges){
	    if(range.inner(image[ivox])){
	      //Voxel in range, assign material and density
	      dens[ivox] = range.dens;
	      mats[ivox] = range.mat;
	      assigned = true;
	      break;
	    }
	  }

	  //Try to assign material and density with
	  //specific density ranges for this contour
	  for(const densityRange& range : cAssign.densityRanges){
	    if(range.inner(dens[ivox])){
	      //Voxel in range, assign material
	      mats[ivox] = range.mat;
	      assigned = true;
	      break;	  
	    }	
	  }

	  //If the voxel material and density has not been assigned
	  //with contour data, use the default values for the contour
	  if(!assigned){
	    unsigned cmat = cAssign.defaultMat;
	    if(cmat > 0)
	      mats[ivox] = cmat;
	    double cdens = cAssign.defaultDens;
	    if(!std::signbit(cdens))
	      dens[ivox] = cdens;
	  }
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
      configStatus = 8;
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
    configStatus = 9;
    return 9;
  }

  // Read enclosure information
  ///////////////////////////////
  double dr;
  if(config.read("enclosure-margin",dr) != INTDATA_SUCCESS){
      if(verbose > 0){
        printf("pen_voxelGeo:configure:Error: Enclosure margin value"
        " 'enclosure-margin' not found. "
        "Double expected.\n");
      }
      err++;
  }else{
      if(dr < 1.0e-6){
          if(verbose > 0){
            printf("pen_voxelGeo:configure:Error: Enclosure "
            "margin value must be greater than zero\n");
          }
          err++;
      }
      
    enclosureMargin = dr;
    
    //Precalculate enclosure high margins (+x,+y,+z)
    enclosureXlimit = Mdx + dr;
    enclosureYlimit = Mdy + dr;
    enclosureZlimit = Mdz + dr;

    //Precalculate the enclosure limits x,y,z moving it to the origin (0,0,0)
    enclosureXlimit0 = enclosureXlimit + dr;
    enclosureYlimit0 = enclosureYlimit + dr;
    enclosureZlimit0 = enclosureZlimit + dr;
  }  
  
  
  
   // Material enclosure
   //**********
	
    //Read material ID
    int auxMat;
    if(config.read("enclosure-material",auxMat) != INTDATA_SUCCESS){
      if(verbose > 0){
        printf("pen_dicomGeo:configure: Error: Unable to read enclosure material ID for material. Integer expecteed");
      }
      err++;
    }
    //Check material ID
    if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
      if(verbose > 0){
        printf("pen_dicomGeo:configure: Error: Invalid ID specified for enclosure material.");
        printf("                         ID: %d\n",auxMat);
        printf("Maximum number of materials: %d\n",constants::MAXMAT);
      }
      err++;
    }
    else
        enclosureMat=auxMat;
    //Create a interface between the enclosure and the mesh
    //assigning a detector to the enclosure
    KDET[0] = 1;  
  
  // Check print image option
  //***************************
  bool toASCII = true;
  if(config.read("print-ASCII",toASCII) != INTDATA_SUCCESS){
    toASCII = false;
  }

  bool mhdMasks = true;
  if(config.read("mhd-masks",mhdMasks) != INTDATA_SUCCESS){
    mhdMasks = false;
  }  

  if(toASCII){

    std::string finalFilename;
    finalFilename = OutputDirPath + std::string("dicomASCII.rep");
    printImage(finalFilename.c_str());
    finalFilename = OutputDirPath + std::string("dicomSeeds.dat");
    dicom.printSeeds(finalFilename.c_str());
    if(dicom.nContours() > 0){
      finalFilename = OutputDirPath + std::string("dicomContours.dat");
      dicom.printContours(finalFilename.c_str());
      finalFilename = OutputDirPath + std::string("dicomContourMask.dat");
      dicom.printContourVox(finalFilename.c_str());
      finalFilename = OutputDirPath + std::string("roi");
      printContourMasks(finalFilename.c_str());
    }
  }

  if(mhdMasks){

    unsigned nVoxelsAxis[3] = {
      static_cast<unsigned>(dicom.getNX()),
      static_cast<unsigned>(dicom.getNY()),
      static_cast<unsigned>(dicom.getNZ())};
    
    float elementSizes[3] = {
      static_cast<float>(dicom.getDX()),
      static_cast<float>(dicom.getDY()),
      static_cast<float>(dicom.getDZ())};
    double origin[3];
    getOffset(origin);
    
    size_t cont = 0;
    for(const std::vector<unsigned char>& mask : dicom.readContourMasks()){

      //Get contour name and remove white spaces
      std::string cname = dicom.contour(cont).name;
      cname.erase(std::remove(cname.begin(),cname.end(),' '),cname.end());

      std::string filename =
	OutputDirPath + std::string("roi_") + cname;

      std::function<std::uint8_t(unsigned long long, size_t)> f =
	[=, &mask](unsigned long long,
		   size_t i) -> std::uint8_t{
	  
	  return static_cast<std::uint8_t>(mask[i]);
	  
	};
      
      pen_imageExporter exporter(f);

      exporter.baseName = filename;
      exporter.setDimensions(3,nVoxelsAxis,elementSizes);
      exporter.setOrigin(origin);

      exporter.exportImage(1,pen_imageExporter::formatTypes::MHD);
      ++cont;
    }
  }

  if(toASCII || mhdMasks){
    std::string finalFilename = OutputDirPath + std::string("roi");
    printContourMaskSummary(finalFilename.c_str());
  }
  
  configStatus = 0;
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

int pen_dicomGeo::printContourMasks(const char* filename) const{

  if(filename == nullptr)
    return PEN_DICOM_ERROR_NULL_FILENAME;

  double voxVol = dicom.getVoxVol();

  std::vector<std::vector<unsigned char>> contourMasks;
  contourMasks = dicom.readContourMasks();

  for(size_t imask = 0; imask < contourMasks.size(); ++imask){

    //Get mask
    const std::vector<unsigned char>& mask = contourMasks[imask];
    
    pen_contour contour;
    contour = dicom.contour(imask);

    //Create a file to store contour mask
    std::string sfilename(filename);
    sfilename.append("_");
    //remove white spaces
    std::string auxstr(contour.name.c_str());
    auxstr.erase(std::remove(auxstr.begin(),auxstr.end(),' '),auxstr.end());
    //sfilename.append(contours[imask].name.c_str());
    sfilename.append(auxstr.c_str());
    std::string sfilenameORIG=sfilename;
    sfilename.append(".dat");
    
    FILE* OutMask = nullptr;
    OutMask = fopen(sfilename.c_str(),"w");
    if(OutMask == nullptr){
      return PEN_DICOM_ERROR_CREATING_FILE;
    }
    
    unsigned long nInner=std::accumulate(mask.begin(),mask.end(),static_cast<unsigned long>(0));
    
    //Count number of contour points
    unsigned long nPoints = 0;

    for(unsigned j = 0; j < contour.NPlanes(); j++)
      {
        nPoints += contour.nPoints(j);
      }

    long int nbin = dicom.getNVox();
    double contMass = 0.0;
    for(long int i = 0; i < nbin; ++i){
      if(mask[i] == 1){
        contMass += densities[mesh[i].MATER-1]*mesh[i].densityFact;
      }
    }    

    fprintf(OutMask,"# PenRed MASK DATA\n");
    fprintf(OutMask,"# Inner voxels are flagged with a '1' in the mask.\n"
      "# Instead, outer voxels are flagged with a '0'.\n");
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# Contour Name:\n");    
    fprintf(OutMask,"#    %s\n",contour.name.c_str());
    fprintf(OutMask,"#\n");    
    fprintf(OutMask,"# Number of voxels:\n");
    fprintf(OutMask,"#    %lu\n",static_cast<unsigned long>(mask.size()));
    fprintf(OutMask,"# Number of inner voxels:\n");
    fprintf(OutMask,"#    %lu\n", nInner);    
    fprintf(OutMask,"# Number contour points:\n");
    fprintf(OutMask,"#    %lu\n", nPoints);    
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# Mass (g):\n");
    fprintf(OutMask,"#    %E\n",
      voxVol*contMass);
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# Volume (cm^3):\n");
    fprintf(OutMask,"#    %E\n",
      voxVol*static_cast<double>(nInner));
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# Mean density (g/cm^3):\n");
    fprintf(OutMask,"#    %E\n",
      contMass/static_cast<double>(nInner));
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# n voxel |  mask \n");
    
    for(size_t i = 0; i < mask.size(); ++i){
      fprintf(OutMask,"  %7lu     %u\n",i,mask[i]);
    }
    
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# End of mask data\n");
    fclose(OutMask);  

  }
    
  return 0;
}

int pen_dicomGeo::printContourMaskSummary(const char* filename) const{

  if(filename == nullptr)
    return PEN_DICOM_ERROR_NULL_FILENAME;

  std::string sfilenameSum(filename);
  sfilenameSum.append("-summary.dat");
    
  FILE* OutSumMask = nullptr;
  OutSumMask = fopen(sfilenameSum.c_str(),"w");
  if(OutSumMask == nullptr){
    return PEN_DICOM_ERROR_CREATING_FILE;
  }

  fprintf(OutSumMask,"# PenRed MASK SUMMARY\n");
  fprintf(OutSumMask,"#\n");
  fprintf(OutSumMask,"# Number of contours:\n");
  fprintf(OutSumMask,"#    %lu\n",dicom.nContours());
  fprintf(OutSumMask,"#\n");
  fprintf(OutSumMask,"# Note that overlapping/contour priority is NOT taking into account.\n");
  fprintf(OutSumMask,"#\n");
  fprintf(OutSumMask,"# Name, number of voxels in contour, contour mass (g), contour volume (cm^3), average density (g/cm^3), number of contour points\n");
  fprintf(OutSumMask,"#\n");

  double voxVol = dicom.getVoxVol();

  std::vector<std::vector<unsigned char>> contourMasks;
  contourMasks = dicom.readContourMasks();

  for(size_t imask = 0; imask < contourMasks.size(); ++imask){

    //Get mask
    const std::vector<unsigned char>& mask = contourMasks[imask];
     
    unsigned long nInner=std::accumulate(mask.begin(),mask.end(),static_cast<unsigned long>(0));
        
    //Count number of contour points
    unsigned long nPoints = 0;
    pen_contour contour;
    contour = dicom.contour(imask);

    for(unsigned j = 0; j < contour.NPlanes(); j++)
      {
        nPoints += contour.nPoints(j);
      }

    long int nbin = dicom.getNVox();
    double contMass = 0.0;
    for(long int i = 0; i < nbin; ++i){
      if(mask[i] == 1){
        contMass += densities[mesh[i].MATER-1]*mesh[i].densityFact;
      }
    }    

    fprintf(OutSumMask,"#    Contour %li: %s %lu %.5E %.5E %.5E %ld\n",
              imask, contour.name.c_str(), nInner,voxVol*contMass,
              voxVol*static_cast<double>(nInner), contMass/static_cast<double>(nInner),
              nPoints);
  }
  fprintf(OutSumMask,"#\n");
  fclose(OutSumMask);  
  
  return 0;
}

int readIntensityRanges(const pen_parserSection& config,
			std::vector<intensityRange>& data,
			const unsigned verbose){

  pen_parserSection intensityRanges;
  std::vector<std::string> intensityRangesNames;
  
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
	  printf("pen_dicomGeo:readIntensityRanges: Error: Unable to read material ID for intensity range '%s'. Integer expected.\n",intensityRangesNames[i].c_str());
	}
	return -1;
      }
      //Check material ID
      if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_dicomGeo:readIntensityRanges: Error: Invalid material ID for intensity range '%s'.\n",intensityRangesNames[i].c_str());
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
	  printf("pen_dicomGeo:readIntensityRanges: Error: Unable to read low intensity for range '%s'. Double expected.\n",intensityRangesNames[i].c_str());
	}
	return -3;
      }

      //Read top voxel intensity
      if(intensityRanges.read(intensityTopField,auxIntensityTop) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:readIntensityRanges: Error: Unable to read top intensity for range '%s'. Double expected.\n",intensityRangesNames[i].c_str());
	}
	return -4;
      }
      
      //Check voxel intensities
      if(auxIntensityLow >= auxIntensityTop){
	if(verbose > 0){
	  printf("pen_dicomGeo:readIntensityRanges: Error: Invalid intensity range specified for range '%s'.\n",intensityRangesNames[i].c_str());
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
	  printf("pen_dicomGeo:readIntensityRanges: Error: Unable to read density for intensity range '%s'. Double expected.\n",intensityRangesNames[i].c_str());
	}
	return -3;
      }

      if(verbose > 1)
	printf("%10.10s -> %4d   %12.4E    %12.4E - %12.4E\n",
	       intensityRangesNames[i].c_str(),auxMat,
	       auxIntensityDensity,auxIntensityLow,auxIntensityTop);
      
      //Store values
      data.emplace_back(static_cast<unsigned>(auxMat),
			auxIntensityLow,
			auxIntensityTop,
			auxIntensityDensity);
    }
  }
  else if(verbose > 1){
    printf("\nNo image intensity ranges specified to assign voxel materials and densities.\n");
  }  

  return 0;
}

int readDensityRanges(const pen_parserSection& config,
		      std::vector<densityRange>& data,
		      const unsigned verbose){


  pen_parserSection ranges;
  std::vector<std::string> rangesNames;

  if(config.readSubsection("ranges",ranges) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:readDensityRanges: No density ranges provided to assign materials\n");
    }
  }else{
    //Extract material names
    ranges.ls(rangesNames);  
  }
  
  if(rangesNames.size() > 0){

    if(verbose > 1)
      printf("\nRange name  | MAT ID | density range (g/cm^3)\n");
    for(unsigned i = 0; i < rangesNames.size(); i++){
      //Read material assigned to this range
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
	  printf("pen_dicomGeo:readDensityRanges: Error: Unable to read material ID for density range '%s'. Integer expected.\n",rangesNames[i].c_str());
	}
	return -1;
      }
      //Check material ID
      if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_dicomGeo:readDensityRanges: Error: Invalid material ID for density range '%s'.\n",rangesNames[i].c_str());
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
	  printf("pen_dicomGeo:readDensityRanges: Error: Unable to read low density for range '%s'. Double expected.\n",rangesNames[i].c_str());
	}
	return -3;
      }

      //Read top density
      if(ranges.read(densTopField.c_str(),auxDensTop) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_dicomGeo:readDensityRanges: Error: Unable to read top density for range '%s'. Double expected.\n",rangesNames[i].c_str());
	}
	return -4;
      }
      
      //Check densities
      if(auxDensLow <= 0.0 || auxDensLow >= auxDensTop){
	if(verbose > 0){
	  printf("pen_dicomGeo:readDensityRanges: Error: Invalid density range specified for range '%s'.\n",rangesNames[i].c_str());
	  printf("                    low density: %12.4E\n",auxDensLow);
	  printf("                    top density: %12.4E\n",auxDensTop);
	}
	return -5;
      }
      if(verbose > 1)
	printf("%10.10s -> %4d   %12.4E - %12.4E\n",rangesNames[i].c_str(),auxMat,auxDensLow,auxDensTop);

      //Store values
      data.emplace_back(static_cast<unsigned>(auxMat),
			auxDensLow,
			auxDensTop);      
    }
  }
  else if(verbose > 1){
    printf("\nNo density ranges specified to assign materials.\n");
  }

  return 0;
}

REGISTER_GEOMETRY(pen_dicomGeo,DICOM)

#endif
