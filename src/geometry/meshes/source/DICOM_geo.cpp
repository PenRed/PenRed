
//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2024-2026 Vicent Giménez Alventosa
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

 
#ifdef _PEN_USE_DICOM_

#include "DICOM_geo.hh"

penred::errors::Error pen_dicomGeo::specificConfigure(const pen_parserSection& config,
						      const unsigned verbose){

  penred::errors::SpecificError<pen_dicomGeo> error;

  // Read DICOM directory path
  //****************************
  std::string directoryPath;
  if(config.read("directory",directoryPath) != INTDATA_SUCCESS){
    error.code = MISSING_CONFIG_PARAMETER;
    error.description = "pen_dicomGeo:configure:Error: Unable to read "
      "field 'directory'. String expected.";
    return error;
  }

  // Check for calibration array
  //*****************************
  pen_parserArray calibration;
  if(config.read("calibration",calibration) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("pen_dicomGeo:configure: No calibration specified, raw image "
	     "data will be used for conversion.\n");
    }
  }

  // Check for default material and density
  //****************************************
  double defDens;
  int defMat;
  if(config.read("default/material",defMat) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("pen_dicomGeo:configure: Warning: No default material "
	     "provided ('default/material'). Material 1 will be used by default\n");
    }
    defMat = 1;
  }

    if(verbose > 1){
      printf("DICOM default material: %d\n",defMat);
    }
  if(defMat < 1){
    error.code = BAD_VALUE;
    error.description = "pen_dicomGeo:configure: Error: Default material must be greater than zero.";
    return error;
  }

  if(config.read("default/density",defDens) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("pen_dicomGeo:configure: Warning: No default density "
	     "provided ('default/density'). Air density will be used by default\n");
    }
    defDens = 0.0012;  //assign air density
  }
  
    if(verbose > 1){
      printf("DICOM default density: %14.5E g/cm^3\n",defDens);
    }
  if(defDens <= 0.0){
    error.code = BAD_VALUE;
    error.description = "pen_dicomGeo:configure: Error: Invalid default density.";
    return error;
  }

  //Check for intensity range material assign
  //****************************************

  std::vector<intensityRange> intensityRanges;
  penred::errors::Error errorIR;
  errorIR = readIntensityRanges(config, intensityRanges, verbose);
  if(errorIR){
    error.code = INTENSITY_RANGE_CONFIG_ERROR;
    error.description = "pen_dicomGeo:configure: Error configuring global intensity ranges.";
    error.setTrace(errorIR);
    return error;
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
    error.code = MISSING_CONFIG_PARAMETER;
    error.description = "pen_dicomGeo:configure: No material information provided.";
    return error;
  }
  
  //Extract material names
  matSec.ls(matNames);  

  if(verbose > 1){
    printf(" Material Name |  ID  | Density (g/cm^3)\n");
  }
  
  //Iterate over all materials
  for(unsigned imat = 0; imat < matNames.size(); imat++){

    double auxDens;
    int auxID;

    std::string idField = matNames[imat] + std::string("/ID");
    std::string densField = matNames[imat] + std::string("/density");

    // Material
    //**********
	
    //Read material ID
    if(matSec.read(idField.c_str(),auxID) != INTDATA_SUCCESS){
      error.code = MISSING_CONFIG_PARAMETER;
      error.description = "pen_dicomGeo:configure: Error: Unable to read material ID for material '";
      error.description += matNames[imat];
      error.description += "'. Integer expected.";
      return error;
    }
    //Check material ID
    if(auxID < 1 || auxID > (int)constants::MAXMAT){
      error.code = BAD_VALUE;
      error.description = "pen_dicomGeo:configure: Error: Invalid ID specified for material '";
      error.description += matNames[imat];
      error.description += "'. Maximum number of materials: " + std::to_string(constants::MAXMAT);
      return error;
    }

    // Density
    //**********
	
    //Read density
    if(matSec.read(densField.c_str(),auxDens) != INTDATA_SUCCESS){
      error.code = MISSING_CONFIG_PARAMETER;
      error.description = "pen_dicomGeo:configure: Error: Unable to read density for material '";
      error.description += matNames[imat];
      error.description += "'.";
      return error;
    }
    //Check density
    if(auxDens <= 0.0){
      error.code = BAD_VALUE;
      error.description = "pen_dicomGeo:configure: Error: Invalid density specified for material '";
      error.description += matNames[imat];
      error.description += "'. Must be greater than zero.";
      return error;
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
    if(verbose > 1){
      printf(" + No contour information provided\n");
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
	  error.code = BAD_VALUE;
	  error.description = "pen_dicomGeo:configure: Error: "
	    "Invalid material ID assigned to contour '" + contourNames[i] + "'.";
	  error.description += "The maximum number of materials is ";
	  error.description += std::to_string(constants::MAXMAT);
	  return error;
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
	  error.code = BAD_VALUE;
	  error.description = "pen_dicomGeo:configure: Error: "
	    "Invalid material density specified for contour '" + contourNames[i] + "'.";
	  error.description += " Must be greater than zero.";
	  return error;
	}
      }
      // Priority
      //**********
	
      //Read priority
      if(contourSec.read(prioField.c_str(),auxPrio) != INTDATA_SUCCESS){
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:configure: Error: Unable to read priority for contour '";
	error.description += contourNames[i];
	error.description += "'.";
	return error;
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

      pen_parserSection contourNameSec;      
      if(contourSec.readSubsection(contourNames[i],contourNameSec) !=
	 INTDATA_SUCCESS){
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:configure: Error: Unable to read contour '";
	error.description += contourNames[i];
	error.description += "' section.";
	return error;
      }

      //Check if intensity ranges have been specified for this contour
      if(contourNameSec.isSection("intensity-ranges")){

	if(verbose > 1){
	  printf("  - Intensity ranges defined in this contour:\n");
	}
	
	//Load intensity and density ranges
	penred::errors::Error errorIRC;
	errorIRC = readIntensityRanges(contourNameSec,
                                   contAssigns.intensityRanges,
                                   verbose);
	if(errorIRC){
	  error.code = INTENSITY_RANGE_CONFIG_ERROR;
	  error.description = "pen_dicomGeo:configure: Error configuring intensity ranges "
	    "in contour " + contourNames[i];
	  error.setTrace(errorIRC);
	  return error;
	}
      }

      //Check if density ranges have been specified for this contour
      if(contourNameSec.isSection("ranges")){

	if(verbose > 1){
	  printf("  - Density ranges defined in this contour:\n");
	}
	
	penred::errors::Error errorDRC;
	errorDRC = readDensityRanges(contourNameSec,
				     contAssigns.densityRanges,
				     verbose);
	if(errorDRC){
	  error.code = DENSITY_RANGE_CONFIG_ERROR;
	  error.description = "pen_dicomGeo:configure: Error configuring density ranges "
	    "in contour " + contourNames[i];
	  error.setTrace(errorDRC);
	  return error;
	}
      }

    }
  }
  else if(verbose > 1){
    printf("\nNo contours specified to assign material and density.\n");

    //Check if we can assign density using calibration
    if(calibration.size() == 0 && intensityRanges.size() == 0){
      error.code = MISSING_CONFIG_PARAMETER;
      error.description = "pen_dicomGeo:configure: Error: No contour information"
	" nor calibration nor intensity ranges provided to assign density to voxels.";
      return error;
    }
  }
  

  //Check for density range material assign
  //****************************************

  penred::errors::Error errorDR;
  std::vector<densityRange> densityRanges;
  errorDR = readDensityRanges(config, densityRanges, verbose);
  if(errorDR){
    error.code = DENSITY_RANGE_CONFIG_ERROR;
    error.description = "pen_dicomGeo:configure: Error configuring global density ranges.";
    error.setTrace(errorDR);
    return error;
  }

  if(densityRanges.size() > 0 && calibration.size() == 0 && verbose > 1){
    printf("\npen_dicomGeo:configure: Warning: Density ranges can't be "
	   "used to assign materials without calibration. Density ranges "
	   "will be ignored.\n");
  }
  
  //Try to read print-ASCII output dir
  std::string OutputDirPath;
  if(config.read("outputdir",OutputDirPath) == INTDATA_SUCCESS){
    if(verbose > 1){
      printf("geometry print-ASCII outputdir: '%s'\n\n",OutputDirPath.c_str());
    }
  }
   
  //Check for segmentation constrains
  //****************************************
  std::vector<segmentConstraints> constraints;
  penred::errors::Error errorSC;
  errorSC = readSegmentConstraints(config, constraints, verbose);
  if(errorSC){
    error.code = MISSING_CONFIG_PARAMETER;
    error.description = "pen_dicomGeo:configure: Error reading segmentation constraints.";
    error.setTrace(errorSC);
    return error;
  }
  
  //*******************
  // Try to load DICOM
  //*******************
  penred::errors::Error errDicom = dicom.loadDicom(directoryPath.c_str(),verbose);
  if(errDicom){
    error.code = DICOM_LOAD_ERROR;
    error.description = "pen_dicomGeo:configure: Error loading DICOM from " + directoryPath;
    error.setTrace(errDicom);
    return error;
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

    //Convert contour name to lower case
    std::string contourNameLowerCase(contourNames[icont]);
    std::transform(contourNameLowerCase.begin(),
		   contourNameLowerCase.end(),
		   contourNameLowerCase.begin(), ::tolower);
	      
    index = dicom.getContourIndex(contourNameLowerCase.c_str());
    if(index >= dicom.nContours()){
      error.code = DICOM_CONTOUR_ERROR;
      error.description = "pen_dicomGeo:configure:Error: Contour " + contourNames[icont] +
	" not found within DICOM stored at " + directoryPath;
      return error;
    }
    //Store configuration index
    contourIndexes[index] = icont;
    //Set priority
    dicom.setContourPriority(index,contourAssigns[icont].priority);
  }

  // Assign contour to each voxel
  //*******************************  
  errDicom = dicom.assignContours();
  if(errDicom){
    error.code = DICOM_CONTOUR_ERROR;
    error.description = "pen_dicomGeo:configure: Error on DICOM contour assign.";
    error.setTrace(errDicom);
    return error;
  }

  if(std::numeric_limits<unsigned>::max() < dicom.getNZ()){
    error.code = LIMIT_EXCEEDED;
    error.description = "pen_dicomGeo:configure:Error: Number of Z planes (" +
      std::to_string(dicom.getNZ()) + ") larger than maximum (" +
      std::to_string(std::numeric_limits<unsigned>::max()) + ")";
    return error;
  }

  unsigned nvox[3] = {static_cast<unsigned>(dicom.getNX()),
		      static_cast<unsigned>(dicom.getNY()),
		      static_cast<unsigned>(dicom.getNZ())};
  double dvox[3] = {dicom.getDX(),dicom.getDY(),dicom.getDZ()};

  unsigned long tnvox = dicom.getNVox();

  std::vector<unsigned> mats(tnvox);
  std::vector<double>   dens(tnvox);

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
      int err = calibration.read(aux,i);
      if(err != INTDATA_SUCCESS){
	error.code = BAD_VALUE;
	error.description = "pen_dicomGeo:configure: Error on calibration coefficient " +
	  std::to_string(i) + ". Number expected.";
	return error;
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
    if(verbose > 1)
      printf("No calibration found.\n");
  }

  // Segmentation constraints
  //***************************

  if(constraints.size() > 1 && verbose > 1)
    printf(" + Apply constraints:\n");
    
  
  bool empty = false;
  double voxVol = dvox[0]*dvox[1]*dvox[2];
  unsigned long nvoxZ = nvox[0]*nvox[1];
  for(const segmentConstraints& c : constraints){

    if(verbose > 1)
      printf("    - %s\n", c.stringify().c_str());

    //Check min volume and number of clusters
    if(c.minVolume > 0.0 ||
       c.maxVolume > 0.0 ||
       c.maxClusters > 0){

      //Store voxel clusters
      std::vector<std::pair<double, std::set<unsigned long>>> clusters;
      
      //Test constraints
      for(unsigned long ivox = 0; ivox < tnvox; ivox++){
	if(mats[ivox] == c.mat){
	  //Constraint material found

	  //Check if has already been processed
	  bool found = false;
	  for(const auto& cluster : clusters){
	    if(cluster.second.find(ivox) != cluster.second.cend()){
	      found = true;
	      break;
	    }
	  }
	  if(found)
	    continue;

	  //Create a new cluster
	  clusters.emplace_back();
	  std::set<unsigned long>& includedVoxels = clusters.back().second;

	  //Create a set with non processed voxels
	  std::unordered_set<unsigned long> nonProcessed;
	  nonProcessed.insert(ivox);

	  //Process all voxels in the cluster
	  while(nonProcessed.size() > 0){

	    //Get first element
	    unsigned long testVox = *nonProcessed.cbegin();

	    //Remove the element in the to be processed list
	    nonProcessed.erase(nonProcessed.cbegin());

	    //Add it to the processed list
	    includedVoxels.insert(testVox);
	    
	    unsigned ix = testVox % nvox[0];
	    unsigned iy = (testVox % nvoxZ) / nvox[0];
	    unsigned iz = testVox / nvoxZ;

	    unsigned xlow = ix > 0 ? ix - 1 : ix;
	    unsigned xtop = std::min(nvox[0]-1, ix+1);

	    unsigned ylow = iy > 0 ? iy - 1 : iy;
	    unsigned ytop = std::min(nvox[1]-1, iy+1);

	    unsigned zlow = iz > 0 ? iz - 1 : iz;
	    unsigned ztop = std::min(nvox[2]-1, iz+1);

	    for(size_t kk = zlow; kk <= ztop; ++kk){
	      for(size_t jj = ylow; jj <= ytop; ++jj){
		for(size_t ii = xlow; ii <= xtop; ++ii){
		  const unsigned long toAdd = kk*nvoxZ + jj*nvox[0] + ii;
		  if(toAdd != testVox &&
		     mats[toAdd] == c.mat &&
		     includedVoxels.count(toAdd) == 0)
		    nonProcessed.insert(toAdd);
		}
	      }
	    }
	    
	  }

	  const double clusterVol = static_cast<double>(includedVoxels.size())*voxVol;
	  clusters.back().first = clusterVol;

	  if(verbose > 1)
	    printf("      * Cluster found with material %u at voxel index %lu "
		   "with volume %E cm**3.\n",
		   c.mat,
		   static_cast<unsigned long>(ivox),
		   clusterVol);
  
	  //Check if the cluster volume is in the limit
	  if(c.minVolume > clusterVol ||
	     (c.maxVolume > 0.0 && c.maxVolume < clusterVol)){

	    if(verbose > 1)
	      printf("         Volume constraints not fulfilled, remove it.\n");
	    
	    //Remove this cluster
	    for(const unsigned long& index : includedVoxels){
	      mats[index] = 0;
	    }
	    empty = true;
	    clusters.pop_back();
	    continue;
	  }

	  //Sort clusters
	  std::sort(clusters.begin(), clusters.end(),
		    std::greater<std::pair<double, std::set<unsigned long>>>());

	  //Check maximum number of clusters
	  if(c.maxClusters > 0 &&
	     clusters.size() > c.maxClusters){
	    //Remove smallest cluster

	    if(verbose > 1)
	      printf("         Maximum number of clusters reached, "
		     "remove the smallest one with volume %E cm**3.\n",
		     clusters.back().first);
	    
	    for(const unsigned long& index : clusters.back().second){
	      mats[index] = 0;
	    }
	    clusters.pop_back();
	    empty = true;
	  }
	}
      }
    }
    
  }

  // Fill empty voxels due constraints
  while(empty){
    empty = false;
    std::vector<voxelInfo> toSet;
    for(unsigned long ivox = 0; ivox < tnvox; ivox++){

      if(mats[ivox] == 0){
	//Set material according to neighbour voxels
	unsigned ix = ivox % nvox[0];
	unsigned iy = (ivox % nvoxZ) / nvox[0];
	unsigned iz = ivox / nvoxZ;
	
	unsigned xlow = ix > 1 ? ix - 2 : ( ix > 0 ? ix - 1 : ix );
	unsigned xtop = std::min(nvox[0]-1, ix+2);

	unsigned ylow = iy > 1 ? iy - 2 : ( iy > 0 ? iy - 1 : iy );
	unsigned ytop = std::min(nvox[1]-1, iy+2);

	unsigned zlow = iz > 1 ? iz - 2 : ( iz > 0 ? iz - 1 : iz );
	unsigned ztop = std::min(nvox[2]-1, iz+2);

	std::array<unsigned, constants::MAXMAT> matWeights;
	std::fill(matWeights.begin(), matWeights.end(), 0);
	double meanDens = 0.0;
	unsigned nonVoid = 0;
	for(size_t kk = zlow; kk <= ztop; ++kk){
	  for(size_t jj = ylow; jj <= ytop; ++jj){
	    for(size_t ii = xlow; ii <= xtop; ++ii){
	      const unsigned long index = kk*nvoxZ + jj*nvox[0] + ii;
	      const unsigned mat = mats[index];

	      if(mat > 0){
		matWeights[mat] += 1;
		meanDens += dens[index];
		nonVoid += 1;
	      }
	    }
	  }
	}

	//Get material with maximum weight
	if(nonVoid == 0){
	  empty = true;
	}
	else{
	  meanDens /= static_cast<double>(nonVoid);
	  toSet.emplace_back(ivox,
			     std::distance(matWeights.cbegin(),
					   std::max_element(matWeights.cbegin(),
							    matWeights.cend())),
			     meanDens);
	}
      }
      
    }

    //Apply changes for this iteration
    if(calibration.size() == 0){
      for(const voxelInfo& v : toSet){
	mats[v.index] = v.mat;
	dens[v.index] = v.dens;
      }
    }
    else{
      for(const voxelInfo& v : toSet){
	mats[v.index] = v.mat;
      }
    }
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
      error.code = BAD_VALUE;
      error.description = "pen_dicomGeo:configure: Error: Nominal density not "
	"provided for material index " + std::to_string(mats[ivox]) + ", but used in the "
	"provided DICOM image";
      return error;
    }
    //Calculate density factor
    dens[ivox] /= densities[mats[ivox]-1];
  }
  
  //Create voxelized geometry
  penred::errors::Error errVox = setVoxels(nvox,dvox,mats.data(),dens.data());
  if(errVox){
    error.code = VOXEL_ASSIGN_ERROR;
    error.description = "pen_dicomGeo:configure:Error: Unable to create voxel "
      "geometry from DICOM with the provided configuration";
    error.setTrace(errVox);
    return error;
  }

  // Read enclosure information
  ///////////////////////////////
  double dr;
  if(config.read("enclosure-margin",dr) != INTDATA_SUCCESS){
    error.code = MISSING_CONFIG_PARAMETER;
    error.description = "pen_voxelGeo:configure:Error: Enclosure margin value"
      " 'enclosure-margin' not found.";
    return error;
  }else{
    if(dr < 1.0e-6){
      error.code = BAD_VALUE;
      error.description = "pen_voxelGeo:configure:Error: Enclosure "
	"margin value must be greater than zero.";      
      return error;
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
      error.code = MISSING_CONFIG_PARAMETER;
      error.description = "pen_voxelGeo:configure:Error: Unable to read enclosure "
	"material ID. Integer expected.";
      return error;
    }
    //Check material ID
    if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
      error.code = BAD_VALUE;
      error.description = "pen_voxelGeo:configure:Error: Enclosure material must be greater than 0 "
	"and lesser than " + std::to_string(constants::MAXMAT);
      return error;
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
  
  return error;
}

penred::errors::Error pen_dicomGeo::printImage(const char* filename) const{
  
  penred::errors::SpecificError<pen_dicomGeo> error;
  
  if(filename == nullptr){
    error.code = NULL_FILENAME;
    error.description = "pen_dicomGeo:printImage:Error: No filename provided.";
    return error;
  }
  
  //Create a file to store contours data
  FILE* OutVox = nullptr;
  OutVox = fopen(filename,"w");
  if(OutVox == nullptr){
    error.code = UNABLE_TO_CREATE_FILE;
    error.description = "pen_dicomGeo:printImage:Error: Unable to create output file.";
    return error;
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

  fclose(OutVox);
  
  return error;
}

penred::errors::Error pen_dicomGeo::printContourMasks(const char* filename) const{

  penred::errors::SpecificError<pen_dicomGeo> error;
  
  if(filename == nullptr){
    error.code = NULL_FILENAME;
    error.description = "pen_dicomGeo:printContourMasks:Error: No filename provided.";
    return error;
  }

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
      error.code = UNABLE_TO_CREATE_FILE;
      error.description = "pen_dicomGeo:printContourMasks:Error: Unable to create output file.";
      return error;
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
    
  return error;
}

penred::errors::Error pen_dicomGeo::printContourMaskSummary(const char* filename) const{

  penred::errors::SpecificError<pen_dicomGeo> error;
  
  if(filename == nullptr){
    error.code = NULL_FILENAME;
    error.description = "pen_dicomGeo:printContourMaskSummary:Error: No filename provided.";
    return error;
  }

  std::string sfilenameSum(filename);
  sfilenameSum.append("-summary.dat");
    
  FILE* OutSumMask = nullptr;
  OutSumMask = fopen(sfilenameSum.c_str(),"w");
  if(OutSumMask == nullptr){
    error.code = UNABLE_TO_CREATE_FILE;
    error.description = "pen_dicomGeo:printContourMaskSummary:Error: Unable to create output file.";
    return error;
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
  
  return error;
}

penred::errors::Error pen_dicomGeo::readIntensityRanges(const pen_parserSection& config,
                                                        std::vector<intensityRange>& data,
                                                        const unsigned verbose){

  penred::errors::SpecificError<pen_dicomGeo> error;

  pen_parserSection intensityRanges;
  std::vector<std::string> intensityRangesNames;
  
  if(config.readSubsection("intensity-ranges",intensityRanges) != INTDATA_SUCCESS){
    if(verbose > 1){
      penred::logs::logger::printf(penred::logs::CONFIGURATION,
				   " No image intensity ranges field ('intensity-ranges') "
				   "provided to assign materials\n");
    }
    return error;
  }else{
    //Extract material names
    intensityRanges.ls(intensityRangesNames);  
  }

  if(intensityRangesNames.size() > 0){
    if(verbose > 1)
      penred::logs::logger::printf(penred::logs::CONFIGURATION,
				   "\nRange name  | MAT ID | density (g/cm^3) | Intensity range\n");
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
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:readIntensityRanges: Error: "
	  "Unable to read material ID for intensity range '" + intensityRangesNames[i] +
	  "'. Integer expected.";
	return error;
      }
      //Check material ID
      if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	error.code = BAD_VALUE;
	error.description = "pen_dicomGeo:readIntensityRanges: Error: "
	  "Invalid material ID for intensity range '" + intensityRangesNames[i] + "'.";
	error.description += "The maximum number of materials is ";
	error.description += std::to_string(constants::MAXMAT);
	return error;
      }

      // Intensity range
      //*****************
	
      //Read low intensity limit
      if(intensityRanges.read(intensityLowField,auxIntensityLow) != INTDATA_SUCCESS){
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:readIntensityRanges: Error: "
	  "Unable to read low intensity for intensity range '" + intensityRangesNames[i];
	return error;
      }

      //Read top voxel intensity
      if(intensityRanges.read(intensityTopField,auxIntensityTop) != INTDATA_SUCCESS){
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:readIntensityRanges: Error: "
	  "Unable to read top intensity for intensity range '" + intensityRangesNames[i];
	return error;
      }
      
      //Check voxel intensities
      if(auxIntensityLow >= auxIntensityTop){
	error.code = BAD_VALUE;
	error.description = "pen_dicomGeo:readIntensityRanges: Error: "
	  "Invalid intensity range specified for range '" + intensityRangesNames[i] + "': ";
	error.description += "[" + std::to_string(auxIntensityLow) + ",";
	error.description += std::to_string(auxIntensityTop) + ").";
	return error;
      }

      // Density (g/cm^3)
      //******************

      //Read low voxel intensity
      if(intensityRanges.read(intensityDensityField,auxIntensityDensity) != INTDATA_SUCCESS){
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:readIntensityRanges: Error: "
	  "Unable to read density for intensity range '" + intensityRangesNames[i];
	return error;
      }

      if(verbose > 1)
	penred::logs::logger::printf(penred::logs::CONFIGURATION,
				     "%10.10s -> %4d   %12.4E    %12.4E - %12.4E\n",
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
    penred::logs::logger::printf(penred::logs::CONFIGURATION,
				 "\nNo image intensity ranges specified to "
				 "assign voxel materials and densities.\n");
  }  

  return error;
}

penred::errors::Error pen_dicomGeo::readDensityRanges(const pen_parserSection& config,
                                                      std::vector<densityRange>& data,
                                                      const unsigned verbose){

  penred::errors::SpecificError<pen_dicomGeo> error;
  
  pen_parserSection ranges;
  std::vector<std::string> rangesNames;

  if(config.readSubsection("ranges",ranges) != INTDATA_SUCCESS){
    if(verbose > 1){
      penred::logs::logger::printf(penred::logs::CONFIGURATION,
				   "pen_dicomGeo:readDensityRanges: "
				   "No density ranges field ('ranges') provided "
				   "to assign materials\n");
    }
    return error;
  }else{
    //Extract material names
    ranges.ls(rangesNames);  
  }
  
  if(rangesNames.size() > 0){

    if(verbose > 1)
      penred::logs::logger::printf(penred::logs::CONFIGURATION,
				   "\nRange name  | MAT ID | density range (g/cm^3)\n");
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
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:readDensityRanges: Error: "
	  "Unable to read material ID for density range '" + rangesNames[i] +
	  "'. Integer expected.";
	return error;
      }
      //Check material ID
      if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	error.code = BAD_VALUE;
	error.description = "pen_dicomGeo:readDensityRanges: Error: "
	  "Invalid material ID for density range '" + rangesNames[i] + "'.";
	error.description += "The maximum number of materials is ";
	error.description += std::to_string(constants::MAXMAT);
	return error;
      }

      // Density range
      //***************
	
      //Read low density
      if(ranges.read(densLowField.c_str(),auxDensLow) != INTDATA_SUCCESS){
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:readDensityRanges: Error: "
	  "Unable to read low density for range '" + rangesNames[i];
	return error;
      }

      //Read top density
      if(ranges.read(densTopField.c_str(),auxDensTop) != INTDATA_SUCCESS){
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:readDensityRanges: Error: "
	  "Unable to read top density for range '" + rangesNames[i];
	return error;
      }
      
      //Check densities
      if(auxDensLow <= 0.0 || auxDensLow >= auxDensTop){
	error.code = BAD_VALUE;
	error.description = "pen_dicomGeo:readDensityRanges: Error: "
	  "Invalid density range specified for range '" + rangesNames[i] + "': ";
	error.description += "[" + std::to_string(auxDensLow) + ",";
	error.description += std::to_string(auxDensTop) + "].";
	return error;
      }
      if(verbose > 1)
	penred::logs::logger::printf(penred::logs::CONFIGURATION,
				     "%10.10s -> %4d   %12.4E - %12.4E\n",
				     rangesNames[i].c_str(),auxMat,auxDensLow,auxDensTop);

      //Store values
      data.emplace_back(static_cast<unsigned>(auxMat),
			auxDensLow,
			auxDensTop);      
    }
  }
  else if(verbose > 1){
    penred::logs::logger::printf(penred::logs::CONFIGURATION,
				 "\nNo density ranges specified to assign materials.\n");
  }

  return error;
}

penred::errors::Error pen_dicomGeo::readSegmentConstraints(const pen_parserSection& config,
                                                           std::vector<segmentConstraints>& data,
                                                           const unsigned verbose){

  penred::errors::SpecificError<pen_dicomGeo> error;

  pen_parserSection constraints;
  std::vector<std::string> constraintsNames;

  if(config.readSubsection("constraints",constraints) != INTDATA_SUCCESS){
    if(verbose > 2){
      penred::logs::logger::printf(penred::logs::CONFIGURATION,
				   "pen_dicomGeo:readSegmentConstraints: No constraints "
				   "field ('constraints') provided to apply on segmentation\n");
    }
    return error;
  }else{
    //Extract constraints names
    constraints.ls(constraintsNames);  
  }
  
  if(constraintsNames.size() > 0){

    if(verbose > 1)
      penred::logs::logger::printf(penred::logs::CONFIGURATION,
				   "\nConstraint name  | MAT ID | min volume (cm^3)"
				   " | max volume (cm^3) | max clusters\n");
    for(unsigned i = 0; i < constraintsNames.size(); i++){
      //Read material assigned to this range
      int auxMat;
      double auxMinVol;
      double auxMaxVol;
      int auxMaxClusters;
      //Create field strings
      std::string matField = constraintsNames[i] + std::string("/material");
      std::string minVolField = constraintsNames[i] + std::string("/min-volume");
      std::string maxVolField = constraintsNames[i] + std::string("/max-volume");
      std::string maxClustersField = constraintsNames[i] + std::string("/max-clusters");

      // Material
      //**********
	
      //Read material ID
      if(constraints.read(matField.c_str(),auxMat) != INTDATA_SUCCESS){
	error.code = MISSING_CONFIG_PARAMETER;
	error.description = "pen_dicomGeo:readSegmentConstraints: Error: "
	  "Unable to read material ID for constraint '" + constraintsNames[i] +
	  "'. Integer expected.";
	return error;
      }
      //Check material ID
      if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
	error.code = BAD_VALUE;
	error.description = "pen_dicomGeo:readSegmentConstraints: Error: "
	  "Invalid material ID for constraint '" + constraintsNames[i] + "'.";
	error.description += "The maximum number of materials is ";
	error.description += std::to_string(constants::MAXMAT);
	return error;
      }

      // Volume constraints
      //*********************
	
      //Read min volume
      if(constraints.read(minVolField.c_str(),auxMinVol) != INTDATA_SUCCESS){
	auxMinVol = -2.0;
      }
      if(constraints.read(maxVolField.c_str(),auxMaxVol) != INTDATA_SUCCESS){
	auxMaxVol = -1.0;
      }

      if(auxMinVol > 0.0 &&
	 auxMaxVol > 0.0 &&
	 auxMinVol >= auxMaxVol){
	error.code = BAD_VALUE;
	error.description = "pen_dicomGeo:readSegmentConstraints: Error in "
	  "constraint '" + constraintsNames[i] + "': Minimum volume greater or equal to maximum ";
	error.description += "[" + std::to_string(auxMinVol) + ",";
	error.description += std::to_string(auxMaxVol) + ").";
	return error;
      }

      // Cluster constraints
      //*********************
      if(constraints.read(maxClustersField.c_str(),auxMaxClusters) != INTDATA_SUCCESS){
	auxMaxClusters = 0;
      }

      if(verbose > 1){
	penred::logs::logger::printf(penred::logs::CONFIGURATION,
				     "%10.10s -> %4d   %12.4E - %12.4E  %d\n",
				     constraintsNames[i].c_str(),
				     auxMat,
				     auxMinVol <= 0.0 ? 0.0 : auxMinVol,
				     auxMaxVol <= 0.0 ? 1.0e35 : auxMaxVol,
				     auxMaxClusters <= 0 ? 100000 : auxMaxClusters);
      }

      //Store values
      data.emplace_back(static_cast<unsigned>(auxMat),
			auxMinVol,
			auxMaxVol,
			auxMaxClusters <= 0 ? 0 : static_cast<unsigned>(auxMaxClusters));
    }
  }
  else if(verbose > 1){
    penred::logs::logger::printf(penred::logs::CONFIGURATION,
				 "\nNo constraints specified for segmentation.\n");
  }

  return error;
}

REGISTER_GEOMETRY(pen_dicomGeo,DICOM)

#endif
