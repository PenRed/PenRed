
//
//
//    Copyright (C) 2019-2023 Universitat de València - UV
//    Copyright (C) 2019-2023 Universitat Politècnica de València - UPV
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


#include "pen_dicom.hh"

//-------------------------
// Penelope dicom contour
//-------------------------

pen_contour::pen_contour(){

  name.assign("NULL");
  //Set number of planes
  nPlanes = 0;
  nPlanePoints = nullptr;
  contourPoints = nullptr;
  priority = 0.0;
}

pen_contour::pen_contour(const pen_contour& c){

  //Initialize
  nPlanes = 0;
  nPlanePoints = nullptr;
  contourPoints = nullptr;

  //Set number of planes
  setPlanes(c.NPlanes());
  //Copy points of each plane
  for(unsigned i = 0; i < nPlanes; i++){
    unsigned long npoints = c.nPoints(i);
    const double* ppoints = c.readPoints(i);
    setPoints(i,npoints,ppoints);      
  }

  name.assign(c.name);
  priority = c.priority;
}

void pen_contour::setPlanes(const unsigned np){

  //Clear previous configuration
  clear();

  //Requires at least one plane
  if(np < 1)
    return;
  
  //Allocate memory for "np" planes and init to zero
  nPlanePoints = (unsigned long *) calloc(np,sizeof(unsigned long));
  if(nPlanePoints == nullptr){
    throw std::bad_alloc();
  }

  //Allocate memory for contour points arrays
  contourPoints = (double**) malloc(sizeof(double*)*np);
  if(contourPoints == nullptr){
    throw std::bad_alloc();
  }

  //Init point arrays to nullptr
  for(unsigned i = 0; i < np; i++)
    contourPoints[i] = nullptr;
    
  //Set number of planes
  nPlanes = np;
}

void pen_contour::setPoints(const unsigned nplane,
			    const unsigned long npoints,
			    const double* points){

    
  //Check if plane number is in range
  if(nplane >= nPlanes)
    throw std::out_of_range("pen_contour:setPoints: Invalid plane number");
  
  //Check if the array must be resized
  if(nPlanePoints[nplane] != npoints ||
     contourPoints[nplane] == nullptr){

    //Clear contour points array
    if(contourPoints[nplane] != nullptr)
      free(contourPoints[nplane]);
    contourPoints[nplane] = nullptr;
    nPlanePoints[nplane] = npoints;
    
    //Allocate memory for new points
    if(npoints > 0){
      contourPoints[nplane] = (double*) malloc(sizeof(double)*npoints*3);
      if(contourPoints[nplane] == nullptr)
	throw std::bad_alloc();
    }
  }

  //Copy memory
  if(npoints > 0 && points != nullptr){
    memcpy(contourPoints[nplane],points,sizeof(double)*npoints*3);
  }
  
}

void pen_contour::updatePoints(const unsigned nplane,
			       const unsigned long init,
			       const unsigned long npoints,
			       const double* points){
  if(npoints == 0)
    return;
  //Check if plane number is in range
  if(nplane >= nPlanes)
    throw std::out_of_range("pen_contour:updatePoints: Invalid plane number");

  //Check last point to update
  if(init+npoints > nPlanePoints[nplane])
    throw std::out_of_range("pen_contour:updatePoints: Points to update out of range");

  //Update points
  memcpy(&contourPoints[nplane][3*init],points,sizeof(double)*npoints*3);
  
}

void pen_contour::convertPoints(const double factor){
  for(unsigned i = 0; i < nPlanes; i++){
    unsigned long ncomponents = 3l*nPlanePoints[i];
    for(unsigned long j = 0; j < ncomponents; j++)
      contourPoints[i][j] *= factor;
  }
}

bool pen_contour::hasPoints() const {
  for(unsigned i = 0; i < nPlanes; i++)
    if(nPoints(i) > 0l)
      return true;
  return false;
}

pen_contour& pen_contour::operator=(const pen_contour& c){

  if(this != &c){
    //Set number of planes
    setPlanes(c.NPlanes());
    //Copy points of each plane
    for(unsigned i = 0; i < nPlanes; i++){
      unsigned long npoints = c.nPoints(i);
      const double* ppoints = c.readPoints(i);
      setPoints(i,npoints,ppoints);      
    }

    name.assign(c.name);
    priority = c.priority;
    
  }
  return *this;
}

void pen_contour::clear(){
  for(unsigned i = 0; i < nPlanes; i++){
    if(contourPoints[i] != nullptr){
      free(contourPoints[i]);
      contourPoints[i] = nullptr;
    }
  }
  
  if(nPlanePoints != nullptr){
    free(nPlanePoints);
    nPlanePoints = nullptr;
  }
  if(contourPoints != nullptr){
    free(contourPoints);
    contourPoints = nullptr;
  }
  
  nPlanes = 0;
}

pen_contour::~pen_contour(){
  clear();
}

pen_seed::pen_seed(){
  positions = nullptr;
  weights = nullptr;
  distances = nullptr;
  directions = nullptr;
  controlPoints = 0;
  size = 0;
  time = 0.0;
}

pen_seed::pen_seed(const pen_seed& c){

  //Initialize
  positions = nullptr;
  weights = nullptr;
  distances = nullptr;
  directions = nullptr;
  controlPoints = 0;
  size = 0;
  time = 0.0;
  
  //Get number of control points
  unsigned long nCP = c.nPoints();
  setCP(nCP);

  //Get control positions pointer
  const double* ppositions = c.readPositions();
  //Set control points positions
  setPositions(ppositions,0,nCP);

  //Get weights
  const double* pweights = c.readWeights();
  //Set weights
  setWeights(pweights,0,nCP);

  //Get distances
  const double* pdist = c.readDistances();
  //Set distances
  setDistances(pdist,0,nCP);

  //Get diractions
  const double* pdir = c.readDirections();
  //Set directions
  setDirections(pdir,0,nCP);
    
  time = c.time;
  ID = c.ID;  
}

pen_seed& pen_seed::operator=(const pen_seed& c){

  if(this != &c){
    //Clear previous configuration
    clear();
    //Get number of control points
    unsigned long nCP = c.nPoints();
    setCP(nCP);

    //Get control positions pointer
    const double* ppositions = c.readPositions();
    //Set control points positions
    setPositions(ppositions,0,nCP);

    //Get weights
    const double* pweights = c.readWeights();
    //Set weights
    setWeights(pweights,0,nCP);

    //Get distances
    const double* pdist = c.readDistances();
    //Set distances
    setDistances(pdist,0,nCP);

    //Get diractions
    const double* pdir = c.readDirections();
    //Set directions
    setDirections(pdir,0,nCP);
    
    time = c.time;
    ID = c.ID;
  }
  
  return *this;
}

void pen_seed::clear(){
  if(positions != nullptr)
    free(positions);
  if(weights != nullptr)
    free(weights);
  if(distances != nullptr)
    free(distances);
  if(directions != nullptr)
    free(directions);

  positions = nullptr;
  weights = nullptr;
  distances = nullptr;
  directions = nullptr;
  controlPoints = 0;
  time = 0.0;
  size = 0;
}

void pen_seed::setCP(unsigned long n){

  //Check if n CPs fits in the existing arrays
  if(n <= size){
    controlPoints = n;
    return;
  }
  
  //Clear previous configuration
  clear();

  //Allocate memory for all arrays
  if(n > 0){
    positions = (double*) calloc(n*3,sizeof(double));
    weights = (double*) calloc(n,sizeof(double));
    distances = (double*) calloc(n,sizeof(double));
    directions = (double*) calloc(n*3,sizeof(double));
    
    if(positions == nullptr ||
       weights == nullptr ||
       distances == nullptr ||
       directions == nullptr){
      throw std::bad_alloc();
    }

    controlPoints = n;
    size = n;
  }
  
}

void pen_seed::setPositions(const double* pos, const unsigned long ni){
  setPositions(pos,ni,1);
}
void pen_seed::setPositions(const double* pos,
			    const unsigned long ni,
			    const unsigned long npos){

  if(npos == 0)
    return;
  if(positions == nullptr)
    throw std::out_of_range("pen_seed:setPositions: No memory allocated for control points positions");  
  //Check final element number
  if(ni+npos > controlPoints)
    throw std::out_of_range("pen_seed:setPositions: Points to update out of range");

  //Copy elements
  memcpy(&positions[3*ni],pos,3*npos*sizeof(double));
}

void pen_seed::setWeights(const double wght, const unsigned long ni){
  setWeights(&wght,ni,1);
}
void pen_seed::setWeights(const double* wght,
			  const unsigned long ni,
			  const unsigned long npos){

  if(npos == 0)
    return;
  if(weights == nullptr)
    throw std::out_of_range("pen_seed:setWeights: No memory allocated for control points weights");  
  
  //Check final element number
  if(ni+npos > controlPoints)
    throw std::out_of_range("pen_seed:setWeights: Points to update out of range");

  //Copy elements
  memcpy(&weights[ni],wght,npos*sizeof(double));
}

void pen_seed::setDistances(const double ds, const unsigned long ni){
  setDistances(&ds,ni,1);
}
void pen_seed::setDistances(const double* ds,
			    const unsigned long ni,
			    const unsigned long npos){

  if(npos == 0)
    return;  
  if(distances == nullptr)
    throw std::out_of_range("pen_seed:setDistances: No memory allocated for control points distances");
  
  //Check final element number
  if(ni+npos > controlPoints)
    throw std::out_of_range("pen_seed:setDistances: Points to update out of range");

  //Copy elements
  memcpy(&distances[ni],ds,npos*sizeof(double));
}

void pen_seed::setDirections(const double* dir, const unsigned long ni){
  setDirections(dir,ni,1);
}
void pen_seed::setDirections(const double* dir,
			    const unsigned long ni,
			    const unsigned long npos){

  if(npos == 0)
    return;
  if(directions == nullptr)
    throw std::out_of_range("pen_seed:setDirections: No memory allocated for control points directions");
  
  //Check final element number
  if(ni+npos > controlPoints)
    throw std::out_of_range("pen_seed:setDistances: Points to update out of range");

  //Copy elements
  memcpy(&directions[3*ni],dir,3*npos*sizeof(double));
}


pen_seed::~pen_seed(){
  clear();
}


void pen_ctData::load(DcmDataset* metainfo){

  // Acquisition type
  const char* acTypeChar;
  OFCondition status = metainfo->findAndGetString(DCM_AcquisitionType,
				      acTypeChar);
  if(status.bad()){
    acquisitionTypeString.assign("UNKNOWN");
    acquisitionType = 5;
  }else{
    acquisitionTypeString.assign(acTypeChar);
    if(acquisitionTypeString.compare("SEQUENCED") == 0)
      acquisitionType = 0;
    else if(acquisitionTypeString.compare("SPIRAL") == 0)
      acquisitionType = 1;
    else if(acquisitionTypeString.compare("CONSTANT_ANGLE") == 0)
      acquisitionType = 2;
    else if(acquisitionTypeString.compare("STATIONARY") == 0)
      acquisitionType = 3;
    else if(acquisitionTypeString.compare("FREE") == 0)
      acquisitionType = 4;
    else
      acquisitionType = 5;
  }

  // Revolution time
  status = metainfo->findAndGetFloat64(DCM_RevolutionTime,
				       revolutionTime);
  if(status.bad()){
    revolutionTime = -1;
  }

  // Total Collimation Width
  status = metainfo->findAndGetFloat64(DCM_TotalCollimationWidth,
				       totalCollimationWidth);
  if(status.bad()){
    totalCollimationWidth = -1;
  }else{
    //Convert to cm
    totalCollimationWidth /= 10.0;
  }

  // Table speed
  status = metainfo->findAndGetFloat64(DCM_TableSpeed,
				       tableSpeed);
  if(status.bad()){
    tableSpeed = -1;
  }else{
    //Convert to cm/s
    tableSpeed /= 10.0;
  }

  // Table Feed Per Rotation
  status = metainfo->findAndGetFloat64(DCM_TableFeedPerRotation,
				       tableFeedPerRotation);
  if(status.bad()){
    tableFeedPerRotation = -1;
  }else{
    //Convert to cm
    tableFeedPerRotation /= 10.0;
  }
	
  // Spiral Pitch Factor
  status = metainfo->findAndGetFloat64(DCM_SpiralPitchFactor,
				       spiralPitchFactor);
  if(status.bad()){
    spiralPitchFactor = -1;
  }

  // Data Collection Diameter
  status = metainfo->findAndGetFloat64(DCM_DataCollectionDiameter,
				       dataCollectionDiameter);
  if(status.bad()){
    dataCollectionDiameter = -1;
  }else{
    //Convert to cm
    dataCollectionDiameter /= 10.0;
  }

  // Filter Type
  const char* filterTypeChar;
  status = metainfo->findAndGetString(DCM_FilterType,
				      filterTypeChar);
  if(status.bad()){
    filterTypeString.assign("NONE");
    filterType = 0;
  }else{
    filterTypeString.assign(filterTypeChar);
    if(filterTypeString.compare("NONE") == 0)
      filterType = 0;
    else if(filterTypeString.compare("STRIP") == 0)
      filterType = 1;
    else if(filterTypeString.compare("WDGE") == 0)
      filterType = 2;
    else if(filterTypeString.compare("BUTTERFLY") == 0)
      filterType = 3;
    else if(filterTypeString.compare("MULTIPLE") == 0)
      filterType = 4;
    else if(filterTypeString.compare("FLAT") == 0)
      filterType = 5;
    else
      filterType = 6;
  }

  // Filter Material
  for(unsigned i = 0; i < maxFilters; ++i){
    const char* filterMatChar;
    status = metainfo->findAndGetString(DCM_FilterMaterial,
					filterMatChar, i);
    if(!status.bad()){
      filterMaterial.push_back(std::string(filterMatChar));
    }else{
      break;
    }
  }

  // Focal Spots
  for(unsigned i = 0; i < maxFilters; ++i){
    double auxFocalSpot;
    status = metainfo->findAndGetFloat64(DCM_FocalSpots,auxFocalSpot,i);
    if(status.bad()){
      break;
    }else{
      //Convert to cm
      auxFocalSpot /= 10.0;
      focalSpots.push_back(auxFocalSpot);
    }
  }

  // CTDIvol
  status = metainfo->findAndGetFloat64(DCM_CTDIvol,
				       CTDIvol);
  if(status.bad()){
    CTDIvol = -1;
  }

  // Reconstruction Diameter
  status = metainfo->findAndGetFloat64(DCM_ReconstructionDiameter,
				       reconstructionDiameter);
  if(status.bad()){
    reconstructionDiameter = -1;
  }else{
    //Convert to cm
    reconstructionDiameter /= 10.0;
  }

  // Distance Source To Detector
  status = metainfo->findAndGetFloat64(DCM_DistanceSourceToDetector,
				       distanceSourceToDetector);
  if(status.bad()){
    distanceSourceToDetector = -1;
  }else{
    //Convert to cm
    distanceSourceToDetector /= 10.0;
  }

  // Distance Source To Data Collection Center
  status = metainfo->findAndGetFloat64(DCM_DistanceSourceToDataCollectionCenter,
				       distanceSourceToDataCollectionCenter);
  if(status.bad()){
    distanceSourceToDataCollectionCenter = -1;
  }else{
    //Convert to cm
    distanceSourceToDataCollectionCenter /= 10.0;
  }
  
  // Rotation Direction
  const char* rotationDirectionChar;
  status = metainfo->findAndGetString(DCM_RotationDirection,
				      rotationDirectionChar);
  if(status.bad()){
    rotationDirectionString.assign("UNKNOWN");
    rotationDirection = 2;
  }else{
    rotationDirectionString.assign(rotationDirectionChar);
    if(rotationDirectionString.compare("CW") == 0)
      rotationDirection = 0;
    else if(rotationDirectionString.compare("CC") == 0)
      rotationDirection = 1;
    else
      rotationDirection = 2;
  }

  // Exposure Time
  Sint32 auxSint32;
  status = metainfo->findAndGetSint32(DCM_ExposureTime,
				      auxSint32);
  if(status.bad()){
    exposureTime = -1;
  }else{
    exposureTime = static_cast<long int>(auxSint32);
  }

  // XRay Tube Current
  status = metainfo->findAndGetSint32(DCM_XRayTubeCurrent,
				      auxSint32);
  if(status.bad()){
    xRayTubeCurrent = -1;
  }else{
    xRayTubeCurrent = static_cast<long int>(auxSint32);
  }

  // Exposure
  status = metainfo->findAndGetSint32(DCM_Exposure,
				      auxSint32);
  if(status.bad()){
    exposure = -1;
  }else{
    exposure = static_cast<long int>(auxSint32);
  }

  // Exposure In uAs
  status = metainfo->findAndGetSint32(DCM_ExposureInuAs,
				      auxSint32);
  if(status.bad()){
    exposureInuAs = -1;
  }else{
    exposureInuAs = static_cast<long int>(auxSint32);
  }

  // Generator Power
  status = metainfo->findAndGetSint32(DCM_GeneratorPower,
				      auxSint32);
  if(status.bad()){
    generatorPowerInt = -1;
  }else{
    generatorPowerInt = static_cast<long int>(auxSint32);
  }


  // Data Collection Center Patient
  status =
    metainfo->findAndGetFloat64(DCM_DataCollectionCenterPatient,
				dataCollectionCenterPatient.x,0);
  if(status.bad()){
    dataCollectionCenterPatient.x = 0.0;
    dataCollectionCenterPatient.y = 0.0;
    dataCollectionCenterPatient.z = 0.0;
  }else{
    status =
      metainfo->findAndGetFloat64(DCM_DataCollectionCenterPatient,
				  dataCollectionCenterPatient.y,1);
    status =
      metainfo->findAndGetFloat64(DCM_DataCollectionCenterPatient,
				  dataCollectionCenterPatient.z,2);


    dataCollectionCenterPatient.x /= 10.0;
    dataCollectionCenterPatient.y /= 10.0;
    dataCollectionCenterPatient.z /= 10.0;
  }

  // Reconstruction Target Center Patient
  status =
    metainfo->findAndGetFloat64(DCM_ReconstructionTargetCenterPatient,
				reconstructionTargetCenterPatient.x,0);
  if(status.bad()){
    reconstructionTargetCenterPatient.x = 0.0;
    reconstructionTargetCenterPatient.y = 0.0;
    reconstructionTargetCenterPatient.z = 0.0;
  }else{
    status =
      metainfo->findAndGetFloat64(DCM_ReconstructionTargetCenterPatient,
				  reconstructionTargetCenterPatient.y,1);
    status =
      metainfo->findAndGetFloat64(DCM_ReconstructionTargetCenterPatient,
				  reconstructionTargetCenterPatient.z,2);


    reconstructionTargetCenterPatient.x /= 10.0;
    reconstructionTargetCenterPatient.y /= 10.0;
    reconstructionTargetCenterPatient.z /= 10.0;
  }

  // Isocenter Position
  status =
    metainfo->findAndGetFloat64(DCM_IsocenterPosition,
				isocenterPosition.x,0);
  if(status.bad()){
    isocenterPosition.x = 0.0;
    isocenterPosition.y = 0.0;
    isocenterPosition.z = 0.0;
  }else{
    status =
      metainfo->findAndGetFloat64(DCM_IsocenterPosition,
				  isocenterPosition.y,1);
    status =
      metainfo->findAndGetFloat64(DCM_IsocenterPosition,
				  isocenterPosition.z,2);


    isocenterPosition.x /= 10.0;
    isocenterPosition.y /= 10.0;
    isocenterPosition.z /= 10.0;
  }

  // Mulienergy CT Acquisition
  const char* usedMultiEnergyCT;
  status = metainfo->findAndGetString(DCM_MultienergyCTAcquisition,
				      usedMultiEnergyCT);
  if(status.bad()){
    multiEnergyAcquisition = false;
  }else{
    if(std::string(usedMultiEnergyCT).compare("YES") == 0)
      multiEnergyAcquisition = true;
    else
      multiEnergyAcquisition = false;
  }

  // KVP
  status = metainfo->findAndGetFloat64(DCM_KVP,kvp);
  if(status.bad()){
    kvp = -1;
  }

  // Single Collimation Width
  status = metainfo->findAndGetFloat64(DCM_SingleCollimationWidth,singleCollimationWidth);
  if(status.bad()){
    singleCollimationWidth = -1;
  }else{
    //Convert to cm
    singleCollimationWidth /= 10.0;
  }

  // Total Collimation Width
  status = metainfo->findAndGetFloat64(DCM_TotalCollimationWidth,totalCollimationWidth);
  if(status.bad()){
    totalCollimationWidth = -1;
  }else{
    //Convert to cm
    totalCollimationWidth /= 10.0;
  }
}

pen_dicom::pen_dicom()
{
  //Set null image modality
  imageModality.assign("NULL");
  //Set default origin
  dicomOrigin[0] = 0.0;
  dicomOrigin[1] = 0.0;
  dicomOrigin[2] = 0.0;
  originalOrigin[0] = dicomOrigin[0];
  originalOrigin[1] = dicomOrigin[1];
  originalOrigin[2] = dicomOrigin[2];

  //Set default dicom voxel dimensions
  dvox_x = 1.0;
  dvox_y = 1.0;
  dvox_z = 1.0;
  voxVol = 1.0;
  //Set default number of voxels
  nvox_x = 0;
  nvox_y = 0;
  nvox_z = 0;
  tnvox = 0;
  //Set dicom dimension
  dimDicom[0] = 0.0;
  dimDicom[1] = 0.0;
  dimDicom[2] = 0.0;

  ////////////////
  //Contour vars//
  ////////////////

  //Save container contour for each voxel (-1 for non contour)
  voxelContour = nullptr;   

  ////////////////
  // Image vars //
  ////////////////

  dicomImage = NULL;    //Store dicom image  
}

int pen_dicom::loadDicom(const char* dirName,
			 const unsigned verbose,
			 const bool onlyMetadata)
{
  //Load and process dicom file
  
  const double minvoxSide=1.0E-6;

  if(dirName == nullptr)
    return PEN_DICOM_FOLDER_NOT_SPECIFIED;

  //Remove '/' at the end of dirname if exists
  std::string dirString(dirName);
  if(dirString[dirString.length()-1] == '/'){
    dirString.pop_back();
  }
  
  //Clear previous configuration
  clear();

  //Open the directory
  DIR *auxDir;
  auxDir=opendir(dirString.c_str());
  if(auxDir==NULL)	//check directory existence
    {
      if(verbose > 0)
	printf("\npen_dicom:loadDicom: Error: Folder not found (%s)\n", dirString.c_str());
      clear();
      return PEN_DICOM_FOLDER_NOT_FOUND;
    }
  if(verbose > 1)
    printf("\nLoading dicom from folder '%s'\n", dirString.c_str());
  
  // use different var to check files in the specified folder
  dirent* data;
  //Stores filename of loaded dicom
  std::string AuxChar;
  //Stores filename of dicom with contours data
  std::string Dicomrtss("");
  //Stores filename of dicom with seeds data
  std::string Dicomrtplan("");
  //Create a vector to store pairs Z-plane/npixels
  //for each dicom image
  std::vector<std::pair<double,unsigned long>> pairs_ZNpix;
  //Store image dicom filenames
  std::vector<std::string> filenamesDicom;
  //Store bottom plane dicom filename position
  std::string bottomDicomFilename;
  bool firstImageDicom = true;
  //Store if we expect only one dicom image file
  bool onlyOne = false;
  //Stores minimum z origin found
  double zmin = 1e20;

  //Read one file in each iteration
  while((data = readdir(auxDir)) != nullptr)           
    {
      //Save filename (d_name field)
      std::string filename(data->d_name);

      //skip hide files (.filename)
      if(filename[0] == '.'){ continue;}

      //append path to filename
      filename = dirName + std::string("/") + filename;      

      //load dicom
      DcmFileFormat fileformat;
      OFCondition status = fileformat.loadFile(filename.c_str()); 
      if(status.good())
	{
	  //Get metadata
	  DcmDataset* metainfo = fileformat.getDataset();
	  //Read dicom modality
	  OFString  auxOfString;
	  status = metainfo->findAndGetOFString(DCM_Modality, auxOfString);
	  if(status.good()) //check modality dataelement
	    {
	      //store dicom image modality
	      std::string dicomModality(auxOfString.c_str()); 
	      auxOfString.clear();
	      if(checkImgModality(dicomModality.c_str()))
		{
		  if(firstImageDicom)
		    {
		      //store dicom image modality
		      firstImageDicom = false;
		      imageModality.assign(dicomModality);
		    }
		  else
		    {
		      //Check if image modalities match
		      if(imageModality.compare(dicomModality) != 0)
			{
			  if(verbose > 0){
			    printf("pen_dicom:loadDicom:Error: multiple image "
				   "modality detected: %s %s\n",
				   imageModality.c_str(), dicomModality.c_str());
			  }
			  metainfo->clear();
			  fileformat.clear();
			  clear();
			  return PEN_DICOM_MULTIPLE_MODALITIES;
			}
		    }
		  
		  //Read Z position of the first plane of this dicom.
		  //If the dicom only have a single image and doesn't contain the
		  //"ImagePositionPatient" field, we asume a (0,0,0) origin. 
		  double AuxZ = 0.0;
		  //check element existence
		  if(metainfo->findAndGetFloat64(DCM_ImagePositionPatient,
						 AuxZ, 2).bad()){
		    //This file doesn't contain "ImagePositionPatient", check if
		    //it's the first and unique dicom image file 
		    AuxZ = 0.0;
		    if(filenamesDicom.size() > 0){
		      if(verbose > 0){
			printf("pen_dicom:loadDicom:Error: Can't access "
			       "'ImagePositionPatient' data from dicom:\n "
			       "'%s'\n",filename.c_str());
		      }
		      metainfo->clear();
		      fileformat.clear();
		      clear();
		      return PEN_DICOM_BAD_READ_IMAGE_POSITION;
		    }
		    onlyOne = true;
		  }
		  else if(onlyOne){
		    if(verbose > 0){
		      printf("pen_dicom:loadDicom:Error: Loading multiple dicom "
			     "image files when previous dicom ('%s') has not "
			     "'ImagePositionPatient' field.\n",
			     bottomDicomFilename.c_str());
		    }
		    metainfo->clear();
		    fileformat.clear();
		    clear();
		    return PEN_DICOM_BAD_READ_IMAGE_POSITION;		    
		  }
		  //read width, height, pixel number, number
		  //of frames and bytes per pixel		      
		  short unsigned int width, height;
		  // - Columns
		  status = metainfo->findAndGetUint16(DCM_Columns,width);
		  if(status.bad())
		    {
		      if(verbose > 0){
			printf("pen_dicom:loadDicom:Error: Can't extract "
			       "'Columns' fom\n   %s\n",filename.c_str());
			printf("   Error: %s\n",status.text());
		      }
		      metainfo->clear();
		      fileformat.clear();
		      clear();
		      return PEN_DICOM_BAD_READ_COLUMNS;
		    }
		  // - Rows
		  status = metainfo->findAndGetUint16(DCM_Rows,height);
		  if(status.bad())
		    {
		      if(verbose > 0){
			printf("loadDicom:Error: can't extract "
			       "'Rows' fom\n   %s\n",filename.c_str());
			printf("   Error: %s\n",status.text());
		      }
		      metainfo->clear();
		      fileformat.clear();
		      clear();
		      return PEN_DICOM_BAD_READ_ROWS;
		    }
		  unsigned long nframes;  //Planes in the image
		  Sint32 auxFrames;
		  // - Frames
		  status = metainfo->findAndGetSint32(DCM_NumberOfFrames,
						      auxFrames);
		  if(status.bad())
		    {
		      //Assume single frame
		      auxFrames = 1;
		    }
		  nframes = static_cast<unsigned long>(abs(auxFrames));
		      
		  // - Store number of pixels in this file
		  unsigned long npixelsFile =
		    static_cast<unsigned long>(width)*
		    static_cast<unsigned long>(height)*nframes;
		  pairs_ZNpix.push_back(std::make_pair(AuxZ,npixelsFile));

		  //Check actual lowest z plane
		  if(AuxZ < zmin)
		    {
		      zmin = AuxZ;
		      bottomDicomFilename.assign(filename);
		    }
		      
		  // - Add the number of frames
		  nvox_z += nframes;

		  //store filename
		  filenamesDicom.push_back(filename);
		    
		}
	      else if(dicomModality.compare("RTSTRUCT") == 0 ||
		      dicomModality.compare("rtstruct") == 0 )
		{
		  //This modality store contours data.
		  //Check the existence of element "ROI Contour Sequence"
		  if(metainfo->tagExists(DCM_ROIContourSequence))
		    {
		      Dicomrtss.assign(filename);
		    }
		}
	      else if(dicomModality.compare("RTPLAN") == 0 ||
		      dicomModality.compare("rtplan") == 0)
		{
		  //This modality store seeds data.
		  //Check the existente of element "Application Setup Sequence"
		  if(metainfo->tagExists(DCM_ApplicationSetupSequence))
		    {
		      Dicomrtplan.assign(filename);
		    }
		}
	      else
		{
		  if(verbose > 0)
		    printf("pen_dicom:loadDicom:Warning: Invalid modality (%s) "
			   "found in dicom:\n   %s\n",
			   dicomModality.c_str(), filename.c_str());
		}
	      metainfo->clear();
	    }
	  else
	    {
	      if(verbose > 0){
		printf("pen_dicom:loadDicom:warning: Unable to read modality from dicom:\n   %s\n", filename.c_str());
		printf("   Error: %s\n",status.text());
	      }
	      metainfo->clear();
	    }

	}
      else
	{
	  if(verbose > 1){
	    printf("%s is not a DICOM, skiping.\n", filename.c_str());
	  }
	}
      fileformat.clear();
    }

  //Close directory
  closedir(auxDir);
  
  if(filenamesDicom.size() == 0)
    {
      if(verbose > 0)
	printf("No dicom files found\n");
      clear();
      return PEN_DICOM_NO_DICOM_FOUND;	   
    }
	  
  if(verbose > 1)
    printf("\n************  %lu DICOMs to load\n",filenamesDicom.size());


  //***********************************
  // Try to read dicom contours data 
  //***********************************
  
  DcmFileFormat fileformat;
  OFCondition status = fileformat.loadFile(Dicomrtss.c_str()); //load dicom
  DcmDataset* metainfo_dfrtss = fileformat.getDataset();
  
  if(status.good())
    {
      if(verbose > 1){
	printf("\n*** Loading DICOMs contours\n");
      }
      //Control if there are more contours to load
      bool remainingContours=true;
      //Store next contour to read
      long int dicomContour = -1;
      while(remainingContours)
	{
	  dicomContour++;
	  //Take contour number "dicomContour" from "Structure Set Roi Sequence"
	  DcmItem *Item_Contour = nullptr;
	  if(metainfo_dfrtss->findAndGetSequenceItem(DCM_StructureSetROISequence,
						     Item_Contour,
						     dicomContour).good())
	    {
	      //Create a new empty contour
	      pen_contour dummyCont;
	      contours.push_back(dummyCont);
	      nVoxContour.push_back(0);
	      pen_contour& contourRef = contours[contours.size()-1];

	      //Take contour name
	      const char* auxconstChar;
	      Item_Contour->findAndGetString(DCM_ROIName,auxconstChar);
	      contourRef.name.assign(auxconstChar);
	      

	      //Convert to lower case
	      std::transform(contourRef.name.begin(),
			     contourRef.name.end(),
			     contourRef.name.begin(), ::tolower);

	      if(verbose > 1){
		printf(" Contour %ld: %s\n",
		       dicomContour,contourRef.name.c_str());
	      }
	      
	      Item_Contour = nullptr;
	      
	      //Read member number "dicomContour" from "Roi Contour Sequence"
	      if(metainfo_dfrtss->findAndGetSequenceItem(DCM_ROIContourSequence,
							 Item_Contour,
							 dicomContour).good())
		{
		  //read points from "Contour Sequence", ordered by planes.
		  //first, count number of planes
		  DcmItem *Item_CPlane = nullptr;
		  long int contPlanes = 0;
		  bool validContour = true;
		  while(Item_Contour->findAndGetSequenceItem(DCM_ContourSequence,
							     Item_CPlane,
							     contPlanes).good())
		    {
		      contPlanes++;
		      //Check if this contour plane geometry is "CLOSED_PLANAR"
		      const char* geoType;
		      Item_CPlane->findAndGetString(DCM_ContourGeometricType,
						    geoType);
		      if(strcmp(geoType,"CLOSED_PLANAR") != 0)
			{
			  if(verbose > 0)
			    printf("pen_dicom:loadDicom:Warning: Invalid "
				   "contour (%s) type. Countours must be"
				   " 'CLOSED_PLANAR'\n",
				   contourRef.name.c_str());
			  validContour = false;
			  break;
			}
		      Item_CPlane = nullptr;
		    }

		  Item_CPlane = nullptr;

		  //Check if is a valid contour
		  if(!validContour){
		    contours.pop_back();
		    nVoxContour.pop_back();
		    continue;
		  }
		  
		  if(verbose > 1)
		    printf("  Number of Planes: %4ld\n",contPlanes);

		  //Allocate memory fot these planes
		  if(contPlanes>0)
		    {
		      contourRef.setPlanes((unsigned)contPlanes);
		      contPlanes=0;
		      while(Item_Contour->
			    findAndGetSequenceItem(DCM_ContourSequence,
						   Item_CPlane,
						   contPlanes).good()){
			//Count the number of points in 'ContourData'
			//All contour planes are 'CLOSED_PLANAR'
			contPlanes++;
			  
			unsigned long nPoints = 0;
			double dummyPoint;
			while(Item_CPlane->findAndGetFloat64(DCM_ContourData,
							     dummyPoint,
							     3*nPoints).good()){
			  nPoints++;
			}
			  
			if(nPoints <= 0)
			  {
			    if(verbose > 0)
			      printf("pen_dicom:loadDicom:Warning: Plane %ld "
				     "of contour '%s' is empty.'\n",
				     contPlanes-1,contourRef.name.c_str());
			    Item_CPlane = nullptr;
			    continue;
			  }

			//Allocate memory for contour plane points
			contourRef.setPoints(contPlanes-1,nPoints);

			//Store points
			for(unsigned long contPoint = 0; contPoint < nPoints; ++contPoint){
			  double p[3];
			  unsigned long index = 3l*contPoint; 
			  Item_CPlane->findAndGetFloat64(DCM_ContourData,p[0],
							 index  );
			  Item_CPlane->findAndGetFloat64(DCM_ContourData,p[1],
							 index+1);
			  Item_CPlane->findAndGetFloat64(DCM_ContourData,p[2],
							 index+2);

			  //Convert from mm to cm
			  p[0] *= 0.1;
			  p[1] *= 0.1;
			  p[2] *= 0.1;
			    
			  contourRef.updatePoints(contPlanes-1,contPoint,1,p);
			}
			  
			Item_CPlane = nullptr;		      
		      }

		      if(contourRef.hasPoints())
			{
			  if(verbose > 1)
			    printf("Contour %s loaded\n",
				   contourRef.name.c_str());
			}
		      else
			{
			  if(verbose > 1)
			    printf("Contour %s has not any point\n",
				   contourRef.name.c_str());
			}
		    }
		  else
		    {
		      if(verbose > 1)
			printf("Contour %s has not any plane\n",
			       contourRef.name.c_str());
		    }

		  Item_Contour = nullptr;
		}
	      else
		{
		  printf("loadDicom:warning: can't read element "
			 "'Roi Contour Sequence' from dicom:\n   %s\n",
			 Dicomrtss.c_str());
		  contours.pop_back();
		  nVoxContour.pop_back();
		}
	    }
	  else
	    {
	      remainingContours=false;
	    }
	}
      metainfo_dfrtss->clear();
      fileformat.clear();      
    }
  else
    {
      if(verbose > 0)
	printf("\n*** Missing contours DICOM file (rtss.dcm)\n");
    }
      	      
  /////////////////////////////////
  // Check dicom with seeds data //
  /////////////////////////////////
  
  //Try to open dicom with seeds data
  status = fileformat.loadFile(Dicomrtplan.c_str()); //load dicom
  DcmDataset* metainfo_dfrtplan = fileformat.getDataset();
  
  if(status.good()) 
    {
      if(verbose > 1)
	printf("\n*** Loading seeds planning\n");

      //control when all seed types was loaded
      bool remainingSeedTypes=true;
      long int nSeeds = 0;
      while(remainingSeedTypes)
	{
	  //read number "dicomPenelope::numSeedTypes" of sequence
	  //"Application Setup Sequence". Check existence of seed
	  //type number "dicomPenelope::numSeedTypes"
	  DcmItem *Item_Seeds = 0;
	  if(metainfo_dfrtplan->
	     findAndGetSequenceItem(DCM_ApplicationSetupSequence,
				    Item_Seeds,nSeeds).good())
	    {
	      //Create new seed
	      size_t nseed = seeds.size();
	      pen_seed dummySeed;
	      seeds.push_back(dummySeed);
	      seeds[nseed].ID = nSeeds;
	      
	      long int channelNumber=0;
	      //read channel number (cateters) from "Channel Sequence"
	      DcmItem *Item_Channel = nullptr;
	      unsigned long nPoints = 0;
	      for(;;)
		{
		  //read member number "channelNumber" from "Channel Sequence"
		  if(Item_Seeds->findAndGetSequenceItem(DCM_ChannelSequence,
							Item_Channel,
							channelNumber).good()) //check existence of channel "channelNumber"
		    {
		      channelNumber++;
		      //Read the number of checkpoints in this channel
		      Sint32 nCheckPoints = 0;
		      status = Item_Channel->
			findAndGetSint32(DCM_NumberOfControlPoints,nCheckPoints);

		      if(status.bad())
			{
			  if(verbose > 0){
			    printf("loadDicom:Error: Can't obtain "
				   "'NumberOfControlPoints' for channel "
				   "%ld.\n",channelNumber+1);
			    printf("                 Error: %s\n",status.text());
			  }
			  //Remove new seed
			  seeds.pop_back();
			  break;
			}
		      
		      if(nCheckPoints <= 0)
			{
			  if(verbose > 0){
			    printf("loadDicom:Error: Expected positive value "
				   "for 'NumberOfControlPoints' in channel"
				   " %ld.\n",channelNumber+1);
			    printf("                 NumberOfControlPoints: "
				   "%ld\n",static_cast<long int>(nCheckPoints));
			  }
			  //Remove new seed
			  seeds.pop_back();
			  break;			  
			}

		      nPoints += nCheckPoints;		      
		      Item_Channel = 0; 
		    }
		  else
		    {
		      break;
		    }
		}

	      //Check if new seed has been removed
	      if(nseed == seeds.size()){
		if(verbose > 0)
		  printf("loadDicom:warning: Ignoring seed\n");
		continue;
	      }
	      
	      //Now we know the number of positions. Set number of control
	      //points to seed.
	      seeds[nseed].setCP(nPoints);

	      //Global seed checkpoint conter
	      unsigned long seedCont = 0;
	      //Channel checkpoint conter
	      long int seedChannelCont = 0;
	      seeds[nseed].time = 0.0;
	      while(seedChannelCont<channelNumber)
		{						    
		  //read member number "seedChannelCont" from "Channel Sequence"
		  if(Item_Seeds->findAndGetSequenceItem(DCM_ChannelSequence,
							Item_Channel,
							seedChannelCont).good())
		    {
		      //Get the 'Final Cumulative Time Weight' and 'Channel Total Time'
		      double finalCumulativeTimeWeight;
		      double ChannelTotalTime;

		      status = Item_Channel->
			findAndGetFloat64(DCM_FinalCumulativeTimeWeight,
					  finalCumulativeTimeWeight);
		      if(status.bad())
			{
			  if(verbose > 0){
			    printf("loadDicom:Error: Can't get "
				   "'Final Cumulative Time Weight' from "
				   "channel %ld.\n",seedChannelCont+1);
			    printf("                Error: %s\n",status.text());
			  }
			  seeds.pop_back();
			  break;
			}

		      status = Item_Channel->
			findAndGetFloat64(DCM_ChannelTotalTime,ChannelTotalTime);
		      if(status.bad())
			{
			  if(verbose > 0){
			    printf("loadDicom:Error: Can't get 'Final "
				   "Cumulative Time Weight' from channel "
				   "%ld.\n",seedChannelCont+1);
			    printf("                Error: %s\n",status.text());
			  }
			  seeds.pop_back();
			  break;
			}
		      
		      seeds[nseed].time += ChannelTotalTime;

		      double previousTime = -1.0; 
		      double previousRelativePosition = -1.0;
		      
		      DcmItem *Item_ControlPoint = 0;
		      bool remainingCheckPoints=true;
		      //Check point conter for specific channel
		      long int contCheckPoints=0;
		      while(remainingCheckPoints)
			{
			  //read member number "contCheckPoints" from
			  //"Brachy Control Point Sequence"
			  if(Item_Channel->
			     findAndGetSequenceItem(DCM_BrachyControlPointSequence,
						    Item_ControlPoint,
						    contCheckPoints).good())
			    {

			      //Store 3D position
			      double AuxPosicio[3];
			      status = Item_ControlPoint->
				findAndGetFloat64(DCM_ControlPoint3DPosition,
						  AuxPosicio[0],0);
			      if(status.bad())
				{
				  if(verbose > 0){
				    printf("loadDicom:Error: Can't get "
					   "'ControlPoint3DPosition' from "
					   "checkpoint %ld in channel %ld of "
					   "seed %lu in dicom:\n   %s\n",
					   contCheckPoints,seedChannelCont,
					   nseed+1,Dicomrtplan.c_str());
				    printf("   Error: %s\n",status.text());
				  }
				  seeds.pop_back();
				  break;
				}
			      Item_ControlPoint->
				findAndGetFloat64(DCM_ControlPoint3DPosition,
						  AuxPosicio[1],1);
			      Item_ControlPoint->
				findAndGetFloat64(DCM_ControlPoint3DPosition,
						  AuxPosicio[2],2);

			      double auxPos[3];
			      
			      auxPos[0]=AuxPosicio[0]*0.1; //convert to cm
			      auxPos[1]=AuxPosicio[1]*0.1; //convert to cm
			      auxPos[2]=AuxPosicio[2]*0.1; //convert to cm

			      //Store position
			      seeds[nseed].setPositions(auxPos,seedCont);

			      //Check if this is the first checkpoint
			      if(previousTime == -1 &&
				 previousRelativePosition == -1)
				{
				  Item_ControlPoint->
				    findAndGetFloat64(DCM_CumulativeTimeWeight,
						      previousTime);
				  Item_ControlPoint->
				    findAndGetFloat64(DCM_ControlPointRelativePosition,
						      previousRelativePosition);
                    
                  //Find the next channel point on a different position
                  DcmItem* Item_ControlPoint_AUX = nullptr;
                  long int contCheckPointsAUX = 1;
                  double nextPos[3] = {auxPos[0],auxPos[1],auxPos[2]+1.0};
                  for(;;){
                    if(Item_Channel->findAndGetSequenceItem(DCM_BrachyControlPointSequence,Item_ControlPoint_AUX,contCheckPointsAUX++).good()){
                            double relativePositionAUX = -1;
                            Item_ControlPoint_AUX->findAndGetFloat64(DCM_ControlPointRelativePosition,relativePositionAUX);
                            if(previousRelativePosition != relativePositionAUX){
                                    
                                Item_ControlPoint_AUX->findAndGetFloat64(DCM_ControlPoint3DPosition,nextPos[0],0);
                                Item_ControlPoint_AUX->findAndGetFloat64(DCM_ControlPoint3DPosition,nextPos[1],1);
                                Item_ControlPoint_AUX->findAndGetFloat64(DCM_ControlPoint3DPosition,nextPos[2],2);
                                
                                nextPos[0] *= 0.1; //to cm
                                nextPos[1] *= 0.1; //to cm
                                nextPos[2] *= 0.1; //to cm
                                
                                break;
                            }
                    }else{
                        //No other points in different positions
                        break;
                    }
                  }
                    
				  //Channel init
                  double u[3];
                  u[0] = nextPos[0] - auxPos[0];
                  u[1] = nextPos[1] - auxPos[1];
                  u[2] = nextPos[2] - auxPos[2];

                  double d = sqrt(pow(u[0],2) + pow(u[1],2) + pow(u[2],2));
                      
                  //Normalize direction
                  u[0] /= d;
                  u[1] /= d;
                  u[2] /= d;
				  
				  seeds[nseed].setDistances(-1.0,seedCont); //First point has not distance
				  seeds[nseed].setDirections(u,seedCont);
				  seeds[nseed].setWeights(0.0,seedCont); //First point has null weight
				}
			      else
				{
				  double actualTime;
				  double actualRelativePosition;
				  Item_ControlPoint->
				    findAndGetFloat64(DCM_CumulativeTimeWeight,
						      actualTime);
				  Item_ControlPoint->
				    findAndGetFloat64(DCM_ControlPointRelativePosition,
						      actualRelativePosition);
				  
				  if(actualTime < previousTime)
				    {
				      remainingCheckPoints = false;
				      Item_ControlPoint = 0;
				      //no sense.
				      if(verbose > 0){
					printf("loadDicom:Warning: in seed type "
					       "%lu, channel %ld checkpoint %ld "
					       "has cumulative time weight "
					       "smaller than previous one.\n",
					       seeds.size()+1,seedChannelCont+1,
					       contCheckPoints);
					printf("                   Skipping all "
					       "next checkpoints.\n");
				      }
				      break;
				    }
				  
				  if(previousTime != actualTime &&
				     previousRelativePosition == actualRelativePosition)
				    {
				      //If times are different but relative positions are equal,
				      //a seed remain in this position for a while.

				      //Seeds remain at the same position
				      seeds[nseed].setDistances(0.0,seedCont);
				      //Take direction of previous checkpoint
				      double auxDirect[3];
				      seeds[nseed].getDirection(auxDirect,seedCont-1);
				      seeds[nseed].setDirections(auxDirect,seedCont);
				      //Store time weight at this position
				      double time = ChannelTotalTime*(actualTime-previousTime)/
					finalCumulativeTimeWeight;
				      
				      seeds[nseed].setWeights(time,seedCont);
				    }
				  else if(previousTime != actualTime &&
					  previousRelativePosition != actualRelativePosition)
				    {
				      //Seed moves between consecutive points with not null
				      //time interval

				      //Get previous position
				      double prevPos[3];
				      seeds[nseed].getPos(prevPos,seedCont-1);
				      
				      double u[3];
				      u[0] = auxPos[0] - prevPos[0];
				      u[1] = auxPos[1] - prevPos[1];
				      u[2] = auxPos[2] - prevPos[2];

				      double d = sqrt(pow(u[0],2) + pow(u[1],2) + pow(u[2],2));
                      
                      //Normalize direction
                      u[0] /= d;
                      u[1] /= d;
                      u[2] /= d;
				      
				      //Seed moves to next point
				      seeds[nseed].setDistances(d,seedCont);
				      //Store direction
				      seeds[nseed].setDirections(u,seedCont);
				      //Store time weight at this step
				      double time = ChannelTotalTime*(actualTime-previousTime)/
					finalCumulativeTimeWeight;
				      
				      seeds[nseed].setWeights(time,seedCont);
				    }
				  else if(previousTime == actualTime &&
					  previousRelativePosition != actualRelativePosition)
				    {
				      //Seed moves between consecutive points without time interval

				      //Get previous position
				      double prevPos[3];
				      seeds[nseed].getPos(prevPos,seedCont-1);

				      
				      double u[3];
				      u[0] = auxPos[0] - prevPos[0];
				      u[1] = auxPos[1] - prevPos[1];
				      u[2] = auxPos[2] - prevPos[2];

				      double d = sqrt( pow(u[0],2) + pow(u[1],2) + pow(u[2],2));

                      //Normalize direction
                      u[0] /= d;
                      u[1] /= d;
                      u[2] /= d;
                      
				      //Seed moves to next point
				      seeds[nseed].setDistances(d,seedCont);
				      //Store direction
				      seeds[nseed].setDirections(u,seedCont);
				      //Null time spended to travel the step
				      seeds[nseed].setWeights(0.0,seedCont);
				    }
				  else if(previousTime == actualTime &&
					  previousRelativePosition == actualRelativePosition)
				    {
				      //no sense.
				      if(verbose > 0){
					printf("pen_dicom:loadDicom:Warning: in "
					       "seed type %lu, channel %ld "
					       "checkpoints %ld and %ld are the "
					       "same checkpoint.\n",
					       nseed+1,seedChannelCont+1,
					       contCheckPoints-1,
					       contCheckPoints);
				      }
				      seedCont--;
				    }
				  previousTime = actualTime;
				  previousRelativePosition = actualRelativePosition;

				  //Check if this is the last point
				  if(actualTime - finalCumulativeTimeWeight >= -1.0e-12)
				    {
				      remainingCheckPoints = false;
				    }
				}
			      seedCont++;
			      Item_ControlPoint = 0;
			      contCheckPoints++;
			    }
			  else
			    {
			      remainingCheckPoints=false;
			    }			  
			}
		      Item_Channel = 0;

		      //Check if new seed has been removed
		      if(nseed == seeds.size()){
			break;
		      }
		    }
		  seedChannelCont++;
		}

	      //Check if new seed has been removed
	      if(nseed == seeds.size()){
		if(verbose > 0)
		  printf("loadDicom:warning: Ignoring seed\n");
		continue;
	      }

	      //Update number of control points
	      seeds[nseed].setCP(seedCont);

	      //Renormalize positions/movements time weights
	      double sumW = seeds[nseed].normWeights();

	      if(verbose > 1)
		printf("%ld control points of seed type %lu loaded. "
		       "Time weights sum: %12.4E/%12.4E\n",
		       seeds[nseed].nPoints(),seeds[nseed].ID+1,
		       sumW,seeds[nseed].time);
	      Item_Seeds = 0;
	      nSeeds++;
	    }
	  else
	    {
	      remainingSeedTypes=false;
	    }					    
	}
      metainfo_dfrtplan->clear();
      fileformat.clear();      
    }
  else
    {
      printf("\n*** Missing seeds planing DICOM file (rtplan.dcm)\n");
    }
  
  //Seeds read.

  //Sort dicom Z values
  std::sort(pairs_ZNpix.begin(),pairs_ZNpix.end());
  
  /////////////////////////
  //  Read image data  ////
  /////////////////////////

  printf("\n");

  //Read data from fisrt dicom (with zmin)
  status = fileformat.loadFile(bottomDicomFilename.c_str()); //load dicom

  DcmDataset* metainfo_bottomDicom = fileformat.getDataset();
  if(status.bad())
    {
      if(verbose > 0)
	printf("pen_dicom:loadDicom:Error: can't re-read fist dicom.");
      metainfo_bottomDicom->clear();
      fileformat.clear();      
      clear();
      return PEN_DICOM_BAD_READ;
    }
  else
    {
      if(imageModality.compare("CT") == 0){

	if(verbose > 1)
	  printf("************  Reading dicom CT properties\n\n");

	ctData.load(metainfo_bottomDicom);
	
      }
      
      if(verbose > 1)
	printf("************  Reading dicom voxels properties\n\n");
      //Extract number of pixels in this dicom image file
      //read width, height, pixel number, number of frames and bytes per pixel
      short unsigned int width, height; 
      status = metainfo_bottomDicom->findAndGetUint16(DCM_Columns,width);
      if(status.bad())
	{
	  if(verbose > 0){
	    printf("pen_dicom:loadDicom:Error: can't extract 'Columns' fom\n   %s\n",bottomDicomFilename.c_str());
	    printf("   Error: %s\n",status.text());
	  }
	  metainfo_bottomDicom->clear();
	  fileformat.clear();      
	  clear();
	  return PEN_DICOM_BAD_READ_COLUMNS;
	}
      status = metainfo_bottomDicom->findAndGetUint16(DCM_Rows,height);
      if(status.bad())
	{
	  if(verbose > 0){
	    printf("pen_dicom:loadDicom:Error: can't extract 'Rows' fom\n   %s\n",bottomDicomFilename.c_str());
	    printf("   Error: %s\n",status.text());
	  }
	  metainfo_bottomDicom->clear();
	  fileformat.clear();      
	  clear();
	  return PEN_DICOM_BAD_READ_ROWS;
	}
      nvox_x = width;
      nvox_y = height;
      nvox_xy = nvox_x*nvox_y;

      // - Frames
      Sint32 bottomFrames;
      status = metainfo_bottomDicom->findAndGetSint32(DCM_NumberOfFrames,
						      bottomFrames);
      if(status.bad())
	{
	  //Assume single frame
	  bottomFrames = 1;
	}
      bottomFrames = abs(bottomFrames);

      //Get z spacing between first and second DICOM
      double zspacing = -1.0;
      if(pairs_ZNpix.size() > 1){
	zspacing = (pairs_ZNpix[1].first - pairs_ZNpix[0].first)/
	  static_cast<double>(bottomFrames);
      }
      
      ///////////////////////////////////////////////////
      //take the origin of the lower plane of the image//
      double AuxOrigenDicom[3];
      status = metainfo_bottomDicom->findAndGetFloat64(DCM_ImagePositionPatient,
						       AuxOrigenDicom[0],0);
      if(status.bad())
	{
	  if(verbose > 1){
	    printf("pen_dicom:loadDicom:Warning: can't extract "
		   "'Image Position Patient' fom\n   %s\n",
		   bottomDicomFilename.c_str());
	    printf("     Origin will be set to (0,0,0)\n");
	  }
	  AuxOrigenDicom[0] = 0.0;
	  AuxOrigenDicom[1] = 0.0;
	  AuxOrigenDicom[2] = 0.0;
	  //metainfo_bottomDicom->clear();
	  //fileformat.clear();      
	  //clear();
	  //return PEN_DICOM_BAD_READ_ORIGIN;
	}
      else{
	metainfo_bottomDicom->findAndGetFloat64(DCM_ImagePositionPatient,
						AuxOrigenDicom[1],1); 
	metainfo_bottomDicom->findAndGetFloat64(DCM_ImagePositionPatient,
						AuxOrigenDicom[2],2); 
      }
      
      dicomOrigin[0]=AuxOrigenDicom[0]*0.1; //transform to cm
      dicomOrigin[1]=AuxOrigenDicom[1]*0.1; //transform to cm
      dicomOrigin[2]=AuxOrigenDicom[2]*0.1; //transform to cm

      originalOrigin[0] = dicomOrigin[0];
      originalOrigin[1] = dicomOrigin[1];
      originalOrigin[2] = dicomOrigin[2];
      
      if(verbose > 1){
	printf("Dicom origin x,y,z(cm):\n");
	printf(" %10.5e %10.5e %10.5e \n",
	       dicomOrigin[0],dicomOrigin[1],dicomOrigin[2]);
      }
      
      /////////////////////////////////////////////////
      //read pixel size
      double DimYPixel;
      double DimXPixel;
      double DimZPixel;
      status = metainfo_bottomDicom->findAndGetFloat64(DCM_PixelSpacing,
						       DimYPixel,0);
      if(status.bad())
	{
	  if(verbose > 0){
	    printf("pen_dicom:loadDicom:Error: can't extract "
		   "'PixelSpacing' from\n   %s\n",
		   bottomDicomFilename.c_str());
	    printf("   Error: %s\n",status.text());
	  }
	  metainfo_bottomDicom->clear();
	  fileformat.clear();      
	  clear();
	  return PEN_DICOM_BAD_READ_PIXEL_SPACING;
	}
      status = metainfo_bottomDicom->findAndGetFloat64(DCM_PixelSpacing,
						       DimXPixel,1);
      if(status.bad())
	{
	  if(verbose > 0){
	    printf("pen_dicom:loadDicom:Error: can't extract "
		   "'PixelSpacing' from\n   %s\n",
		   bottomDicomFilename.c_str());
	    printf("   Error: %s\n",status.text());
	  }
	  metainfo_bottomDicom->clear();
	  fileformat.clear();      
	  clear();
	  return PEN_DICOM_BAD_READ_PIXEL_SPACING;
	}
      status = metainfo_bottomDicom->findAndGetFloat64(DCM_SliceThickness,
						       DimZPixel);
      if(status.bad())
	{
	  if(std::signbit(zspacing)){
	    if(verbose > 0){
	      printf("pen_dicom:loadDicom:Error: can't extract "
		     "'SliceThickness' from\n   %s\n",
		     bottomDicomFilename.c_str());
	      printf("   Error: %s\n",status.text());
	    }
	    metainfo_bottomDicom->clear();
	    fileformat.clear();      
	    clear();
	    return PEN_DICOM_BAD_READ_SLICE_THICKNESS;
	  }
	  else{
	    if(verbose > 1){
	      printf("pen_dicom:loadDicom:Warning: can't extract "
		     "'SliceThickness' from\n   %s\n",
		     bottomDicomFilename.c_str());
	      printf("   Error: %s\n",status.text());
	      printf(" The calculated z spacing between first and second "
		     "DICOM will be used.\n");
	    }
	    DimZPixel = zspacing;
	  }
	}

      //Check if spacings coincide
      if(!std::signbit(zspacing)){
	if(fabs(zspacing-DimZPixel) > 1.0e-8){
	  if(verbose > 1){
	    printf("pen_dicom:loadDicom:Warning: Read slice thickness and "
		   "the calculated spacing between first and "
		   "second DICOM mismatch.\n"
		   "             slice thickness: %12.4E mm\n"
		   "            calculated space: %12.4E mm\n"
		   " The calculated z spacing between first and second "
		   "DICOM will be used.\n"		   
		   ,DimZPixel,zspacing);
	  }
	  DimZPixel = zspacing;
	}
      }

      //Obtain image orientation
      const double* imageOrientation;
      const double defaultOrientation[6] = {1.0,0.0,0.0,
					    0.0,1.0,0.0};
      
      status = metainfo_bottomDicom->findAndGetFloat64Array(DCM_ImageOrientationPatient,
							    imageOrientation);
      if(status.bad())
	{
	  //If image orientation information is not stored in dicom file
	  //use default image orientation
	  imageOrientation = defaultOrientation;
	}

      if(verbose > 1){
	printf("Image orientation:\n");
	printf(" x:(%4.2f,%4.2f,%4.2f) y:(%4.2f,%4.2f,%4.2f)\n",
	       imageOrientation[0],imageOrientation[1],imageOrientation[2],
	       imageOrientation[3],imageOrientation[4],imageOrientation[5]);
      }
	
      dvox_x=DimXPixel/10.0; //transform to cm
      dvox_y=DimYPixel/10.0; //transform to cm
      dvox_z=DimZPixel/10.0; //transform to cm
      dimDicom[0] = static_cast<double>(nvox_x)*dvox_x;
      dimDicom[1] = static_cast<double>(nvox_y)*dvox_y;
      dimDicom[2] = static_cast<double>(nvox_z)*dvox_z;

      //Move dicom origin
      dicomOrigin[0] -= 0.5*dvox_x;
      dicomOrigin[1] -= 0.5*dvox_y;
      dicomOrigin[2] -= 0.5*dvox_z;
      
      if(verbose > 1){
	printf("Voxel dimensions in x,y,z (cm):\n");
	printf(" %12.5E %12.5E %12.5E\n",dvox_x,dvox_y,dvox_z);
	printf("Move origin to first voxel down left bottom corner x,y,z(cm):\n");
	printf(" %10.5e %10.5e %10.5e \n",dicomOrigin[0],dicomOrigin[1],dicomOrigin[2]);fflush(stdout);
      }
	
      double tmpaux=(dvox_x < dvox_y ? dvox_x : dvox_y);
      double tmp2aux=(tmpaux < dvox_z ? tmpaux : dvox_z);
      if (tmp2aux < minvoxSide)
	{
	  if(verbose > 0){
	    printf("pen_dicom:loadDicom:ERROR: voxel side too small, tracking algorithm\n");
	    printf("  requires voxel sides to be larger than (cm):%12.5E\n",minvoxSide);
	  }
	  metainfo_bottomDicom->clear();
	  fileformat.clear();      
	  clear();
	  return PEN_DICOM_SMALL_VOXELS;
	}					      
      
      voxVol = dvox_x*dvox_y*dvox_z;
      if(verbose > 1){
	printf("Voxels volume (cm^3):\n");
	printf("%12.5E\n",voxVol);fflush(stdout);
      }

      //Transform contour points and seed positions  
      transformContoursAndSeeds(imageOrientation,verbose); //revisar

    }
        
  metainfo_bottomDicom->clear();
  fileformat.clear();

  tnvox = nvox_xy*nvox_z;

  if(onlyMetadata){
    if(verbose > 1){
      printf(" Loaded DICOM metada\n");fflush(stdout);
    }    
    return PEN_DICOM_SUCCESS;
  }
	  
  // Allocate arrays:
  dicomImage   = nullptr;
  voxelContour = nullptr;
  dicomImage   = (double*) malloc(sizeof(double)*tnvox);
  voxelContour = (int*)    malloc(sizeof(int)*tnvox);
      
  if (dicomImage == nullptr || voxelContour == nullptr)
    {
      if(verbose > 0)
	printf("pen_dicom:loadDicom: Error allocating memory.\n");
      clear();
      return PEN_DICOM_BAD_ALLOCATION;
    }
  

  if(verbose > 0)
    printf("\n");
      
  //Check if spacing is consistent between all DICOMs
  double dvox_zmm = dvox_z*10.0;
  if(verbose > 2){
      printf("\nSlicePosition-SliceThickness consistency check\n");
      printf("----------------------------------------------\n\n");
      printf(" Index Z       z[iz]            z[iz+1]        SliceSpacing        dvox_z\n");
      printf(" =======  ================  ================  ===============  =============== \n");
  }
  for(size_t k = 0; k < pairs_ZNpix.size()-1; ++k){
    double nframes = static_cast<double>(pairs_ZNpix[k].second)/static_cast<double>(nvox_xy);
    double spacing = (pairs_ZNpix[k+1].first - pairs_ZNpix[k].first)/nframes;
    if(verbose > 2){
      printf("  %3ld     %15.9E  %15.9E  %15.9E  %15.9E\n",k,pairs_ZNpix[k].first,pairs_ZNpix[k+1].first,spacing,dvox_zmm);
      fflush(stdout);
    }
    if(fabs(spacing-dvox_zmm)/dvox_zmm > 5.0e-5){
      clear();
      if(verbose > 0)
	printf("pen_dicom:loadDicom: Error: Spacing between images mismatch.\n");
      return PEN_DICOM_ERROR_SPACING_MISMATCH;
    }
  }
  
  unsigned long loadedDicoms = 0;

  //Read one file data in each iteration
  for(size_t k = 0; k < filenamesDicom.size(); k++)           
    {
      AuxChar.assign(filenamesDicom[k]); // extract filename

      if(verbose > 1)
	printf("Processing dicomfile %s\n",AuxChar.c_str());
      
      long int imagePos = -1; //Save image position
      unsigned long firstPixelPos = 0; //stores first pixel to be filled

      //check if it's a dicom image
      status = fileformat.loadFile(AuxChar.c_str()); //load dicom
      DcmDataset* metainfo_Dicom;
      
      if(status.good())
	{
	  metainfo_Dicom = fileformat.getDataset();

	  //Check image position,"ImagePositionPatient",
	  //and read Z position of this file.
	  //If we have only one file to load, missing "ImagePositionPatient" field
	  //will be interpreted as (0.0,0.0,0.0) origin
	  double AuxZPos = 0.0;
	  status = metainfo_Dicom->findAndGetFloat64(DCM_ImagePositionPatient,AuxZPos,2);
	  if(status.good() || filenamesDicom.size() == 1)
	    {
	      if(status.bad())
		AuxZPos = 0.0;
	      //Search the position of this file
	      imagePos = -1;
	      for(long unsigned i = 0; i < pairs_ZNpix.size(); i++)
		if(AuxZPos == pairs_ZNpix[i].first){
		  imagePos = i;
		  break;
		}
		      
	      if(imagePos==-1) //Not found
		{
		  if(verbose > 0){
		    printf("pen_dicom:loadDicom: Error: Missing dicom file %s\n", AuxChar.c_str());
		  }
		  metainfo_Dicom->clear();
		  fileformat.clear();
		  clear();
		  return PEN_DICOM_ERROR_REOPENING_DICOM;
		}
	    }
	  else
	    {
	      //this dicom has not image data and there are more
	      //than single image file to load
	      if(verbose > 1){
		printf("pen_dicom:loadDicom: Error: Dicom file %s doesn't contain "
		       "Z plane position information.\n",AuxChar.c_str());
	      }	      
	      metainfo_Dicom->clear();
	      fileformat.clear();
	      continue;
	    }
		  	   	    
	}
      else
	{
	  //Dicom is corrupted
	  if(verbose > 1){
	    printf("pen_dicom:loadDicom: Error: Unable to read previous "
		   "read dicom file %s.\n",AuxChar.c_str());
	  }
	  metainfo_bottomDicom->clear();
	  fileformat.clear();
	  clear();
	  return PEN_DICOM_ERROR_REOPENING_DICOM;
	}
      if(verbose > 2)
	printf("Dicom loaded\n");
      
      
      //Calculate first pixel position for this dicom
      for(int long k2 = 0; k2 < imagePos; k2++)
	{
	  firstPixelPos += pairs_ZNpix[k2].second;
	}
      
      //get dicom image
      //DicomImage *image = new DicomImage(AuxChar.c_str());
      DicomImage *image = new DicomImage(AuxChar.c_str(),CIF_IgnoreModalityTransformation);
      if(image->getStatus() == EIS_Normal && image != 0)
	{
	  if(!(image->isMonochrome())) //Check if this image uses monochromatic format ("SamplesPerPixel" = 1)
	    {
	      if(verbose > 0)
		printf("pen_dicom:load_dicom: Error: Image format is not grayscale:\n   %s\n",
		       AuxChar.c_str());
	      delete image;
	      clear();
	      return PEN_DICOM_NON_MONOCHROME_IMAGE;
	    }
	  //read width, height
	  short unsigned int width, height; 
	  status = metainfo_Dicom->findAndGetUint16(DCM_Columns,width);
	  if(status.bad())
	    {
	      if(verbose > 0){
		printf("pen_dicom:loadDicom:Error: Can't extract 'Columns' fom\n   %s\n",
		       AuxChar.c_str());
		printf("   Error: %s\n",status.text());
	      }
	      delete image;
	      clear();
	      return PEN_DICOM_BAD_READ_COLUMNS;
	    }
      if(verbose > 2)
	printf("Width: %d\n",width);
      
	  status = metainfo_Dicom->findAndGetUint16(DCM_Rows,height);
	  if(status.bad())
	    {
	      if(verbose > 0){
		printf("pen_dicom:loadDicom: Error: can't extract 'Rows' fom\n   %s\n",
		       AuxChar.c_str());
		printf("   Error: %s\n",status.text());
	      }
	      delete image;
	      clear();
	      return PEN_DICOM_BAD_READ_ROWS;
	    }
	  
	  if(verbose > 2)
	    printf("Height: %d\n",height);

	  //Check dimensions
	  if(nvox_x != width || nvox_y != height)
	    {
	      if(verbose > 0){
		printf("pen_dicom:load_dicom: Error: DICOMs image size mismatched\n");
	      }
	      delete image;
	      clear();
	      return PEN_DICOM_MISMATCH_DIMENSIONS;
	    }

	  //Get planes in image file
	  unsigned long nframes;
	  Sint32 auxFrames;
	  status = metainfo_Dicom->findAndGetSint32(DCM_NumberOfFrames,auxFrames);
	  if(status.bad())
	    {
	      //Assume only 1 frame
	      auxFrames = 1;
	    }
	  nframes = static_cast<unsigned long>(abs(auxFrames));
	  
	  if(verbose > 2)
	    printf("Nframes: %ld\n",nframes);

	  unsigned long totalDicomVoxels =
	    static_cast<unsigned long>(width)*
	    static_cast<unsigned long>(height)*nframes;
	  if(totalDicomVoxels != pairs_ZNpix[imagePos].second){
	    if(verbose > 0){
	      printf("pen_dicom:loadDicom: Error: Number of voxels mismatch into "
		     "dicom file %s\n",AuxChar.c_str());
	    }
	    delete image;
	    clear();
	    return PEN_DICOM_NVOXELS_MISMATCH;
	  }

	  if(verbose > 2){
	    printf("TotalDicomVoxels: %ld\n",totalDicomVoxels);
	    printf("Reading Intercept and Slope ... \n");fflush(stdout);
	  }
	  //RescaleIntercept and RescaleSlope, are used to transform each pixel
	  //value with the formula:  RescaleSlope*pixelValue+RescaleIntercept
	  double RescaleIntercept;
	  status = metainfo_Dicom->findAndGetFloat64(DCM_RescaleIntercept,RescaleIntercept);
	  double RescaleSlope;
	  status = metainfo_Dicom->findAndGetFloat64(DCM_RescaleSlope,RescaleSlope);
	  if(RescaleSlope==0.0){RescaleSlope=1.0;}
	  
	  if(verbose > 2){
	    printf("RescaleIntercept: %15.8E\n RescaleSlope: %15.8E\n",
		   RescaleIntercept,RescaleSlope);
	    fflush(stdout);
	  }
    //RescaleIntercept=0.0;
    //RescaleSlope=1.0;
    //printf("but RescaleIntercept: %15.8E\n RescaleSlope: %15.8E will be used instead.\n",RescaleIntercept,RescaleSlope);fflush(stdout);

	  //read if pixel data is signed (1) or unsigned (0)
	  short unsigned int PixelRep;
	  status = metainfo_Dicom->findAndGetUint16(DCM_PixelRepresentation,PixelRep);
	  if(status.bad())
	    {
	      if(verbose > 0){
		printf("pen_dicom:loadDicom:Error: can't extract 'PixelRepresentation' from\n   %s\n",AuxChar.c_str());
		printf("   Error: %s\n",status.text());
	      }
	      delete image;
	      clear();
	      return PEN_DICOM_BAD_READ_PIXEL_REPRESENTATION;
	    }
	  
	  if(verbose > 2){
	    printf("Pixel representation: %u\n",PixelRep);
	    printf("Getting data ...\n");
	    fflush(stdout);
	  }
	  
	  //load pixels
	  const DiPixel *inter = image->getInterData();
	  if(inter == nullptr)
	    {
	      if(verbose > 0){
		printf("pen_dicom:load_dicom: Error: Can't extract pixel data"
		       " form dicom:\n   %s\n",AuxChar.c_str());
	      }
	      delete image;
	      clear();
	      return PEN_DICOM_BAD_READ_PIXEL_DATA;
	    }

	  //Check number of pixels
	  if(inter->getCount() != totalDicomVoxels){
	    if(verbose > 0){
	      printf("pen_dicom:load_dicom: Error: Image pixel count and"
		     " expected number of pixels mismatch"
		     " form dicom:\n   %s\n"
		     "     Read: %lu\n"
		     " Expected: %lu\n",
		     AuxChar.c_str(),inter->getCount(),totalDicomVoxels);
	    }
	    delete image;
	    clear();
	    return PEN_DICOM_BAD_READ_PIXEL_DATA;	    
	  }

	  //Get pixel representation
	  EP_Representation rep = inter->getRepresentation();

	  if(rep == EPR_Uint8){
      printf("Pixel representation: Uint8\n");fflush(stdout);
	    const Uint8* pixeldata =
	      static_cast<const Uint8*>(inter->getData());
	    for(unsigned long int ii = 0; ii < totalDicomVoxels; ii++)
	      dicomImage[firstPixelPos + ii] =
		//static_cast<double>(pixeldata[ii]);
    RescaleIntercept+RescaleSlope*static_cast<double>(pixeldata[ii]);
	  }
	  else if(rep == EPR_Sint8){
      printf("Pixel representation: Sint8\n");fflush(stdout);
	    const Sint8* pixeldata =
	      static_cast<const Sint8*>(inter->getData());
	    for(unsigned long int ii = 0; ii < totalDicomVoxels; ii++)
	      dicomImage[firstPixelPos + ii] =
		//static_cast<double>(pixeldata[ii]);	    
    RescaleIntercept+RescaleSlope*static_cast<double>(pixeldata[ii]);
	  }
	  else if(rep == EPR_Uint16){
      printf("Pixel representation: Uint16\n");fflush(stdout);
	    const Uint16* pixeldata =
	      static_cast<const Uint16*>(inter->getData());
	    for(unsigned long int ii = 0; ii < totalDicomVoxels; ii++)
	      dicomImage[firstPixelPos + ii] =
		//static_cast<double>(pixeldata[ii]);	    	    
    RescaleIntercept+RescaleSlope*static_cast<double>(pixeldata[ii]);
	  }
	  else if(rep == EPR_Sint16){
      printf("Pixel representation: Sint16\n");fflush(stdout);
	    const Sint16* pixeldata =
	      static_cast<const Sint16*>(inter->getData());
	    for(unsigned long int ii = 0; ii < totalDicomVoxels; ii++)
	      dicomImage[firstPixelPos + ii] =
		//static_cast<double>(pixeldata[ii]);	    
    RescaleIntercept+RescaleSlope*static_cast<double>(pixeldata[ii]);
	  }
	  else if(rep == EPR_Uint32){
      printf("Pixel representation: Uint32\n");fflush(stdout);
	    const Uint32* pixeldata =
	      static_cast<const Uint32*>(inter->getData());
	    for(unsigned long int ii = 0; ii < totalDicomVoxels; ii++)
	      dicomImage[firstPixelPos + ii] =
		//static_cast<double>(pixeldata[ii]);	    	    
    RescaleIntercept+RescaleSlope*static_cast<double>(pixeldata[ii]);
	  }
	  else if(rep == EPR_Sint32){
      printf("Pixel representation: Sint32\n");fflush(stdout);
	    const Sint32* pixeldata =
	      static_cast<const Sint32*>(inter->getData());
	    for(unsigned long int ii = 0; ii < totalDicomVoxels; ii++)
	      dicomImage[firstPixelPos + ii] =
		//static_cast<double>(pixeldata[ii]);	    	    
    RescaleIntercept+RescaleSlope*static_cast<double>(pixeldata[ii]);
	  }
	  else
	    {
	      if(verbose > 0){
		printf("pen_dicom:load_dicom: Error: Invalid pixel "
		       "representation (%d) form dicom:\n   %s\n",
		       rep,AuxChar.c_str());
	      }
	      delete image;
	      clear();
	      return PEN_DICOM_BAD_PIXEL_REPRESENTATION;
	    } 
	}
      else
	{
	  if(verbose > 0){
	    printf("pen_dicom:load_dicom: Error: Can't open image "
		   "from dicom:\n   %s\n",AuxChar.c_str());
	  }
	  delete image;
	  clear();
	  return PEN_DICOM_BAD_IMAGE_OPEN;
	}
      delete image;
      
      loadedDicoms++;
      metainfo_Dicom->clear();
      fileformat.clear();
    }

  if(verbose > 1){
    printf(" Loaded image DICOMs: %lu\n",loadedDicoms);fflush(stdout);
  }
  return PEN_DICOM_SUCCESS;
}

int pen_dicom::assignContours(){

  //////////////////////////////
  // Check voxels in contours //
  //////////////////////////////

  if(voxelContour == nullptr){
    return PEN_DICOM_NO_DICOM_LOADED;
  }
  
  if( tnvox > 0 )
    {

      //To know if a voxel is in a contour the algorithm is the solution 1 (2D)
      //described in the notes on "Determining if a point lies on the interior
      //of a polygon" written by Paul Bourke 1997.

      //Fill pairs priority/index with the information of each contour
      std::vector<std::pair<float,unsigned>> pairs;
      for(unsigned icontour = 0; icontour < contours.size(); icontour++){
	pairs.push_back(std::make_pair(contours[icontour].priority,icontour));
      }

      //Sort pairs by priority
      std::sort(pairs.begin(), pairs.end());

      int* voxelContour2 = nullptr;
      voxelContour2 = (int*) malloc(sizeof(int)*tnvox);

      if(voxelContour2 == nullptr){
	return PEN_DICOM_BAD_ALLOCATION;
      }

      //Clear voxels contour
      for(unsigned long i = 0; i < tnvox; i++){
	voxelContour[i] = -1;
      }
      //Clear auxiliar array
      memcpy(voxelContour2,voxelContour,sizeof(int)*tnvox);

      //Initialize individual contour masks
      contourMasks.clear();
      contourMasks.resize(contours.size());

      for(auto& contMask : contourMasks){
	contMask.resize(tnvox,0);
      }
      
      //Now, assign a contour to each voxel
      for(unsigned long ipair = 0; ipair < contours.size(); ipair++){

	//Get the corresponding contour index
	unsigned icontour = pairs[ipair].second;
	//Get contour reference
	pen_contour& refCont = contours[icontour];

	//Get individual contour mask for this contour
	std::vector<unsigned char>& contourMask = contourMasks[icontour];
	
	//Iterate over all contour planes
	for(unsigned nplane = 0; nplane < refCont.NPlanes(); nplane++){

	  //Need at least 3 points in the same plane to create a closed contour
	  if(refCont.nPoints(nplane) < 3)
	    continue;
	  
	  //Calculate Z plane
	  double p0[3];
	  refCont.getPoint(p0,nplane,0);
	  double zPlanef = p0[2]/dvox_z;
	  long int zPlane = p0[2]/dvox_z;
	  
	  //Ensure that the contour Z coordinate is not on the voxel Z edge
	  if(zPlanef - static_cast<double>(zPlane) < 1.0e-6){
	    //The contour is on the voxel Z edge. Assign it to lower voxel
	    zPlane = static_cast<long int>(zPlanef - dvox_z*0.5);
	  }
      
	  if(zPlane < 0 || zPlane >= (int)nvox_z)
	    continue;

	  long int zIndex = nvox_xy*zPlane;
	  //Iterate over all rows (y axis)
	  for(unsigned long i = 0; i < nvox_y; i++)
	    {
	      //Calculate row index in this plane
	      unsigned long rowIndex = i*nvox_x;
	      //Calculate central point of row 'i'
	      double yCenter = ((double)i+0.5)*dvox_y;
              for(unsigned long j = 0; j < nvox_x; j++)
                {
                  unsigned long index = zIndex + rowIndex + j;
                  //Calculate central point of column 'j'
                  double xCenter = ((double)j+0.5)*dvox_x;
                  //Iterate over all contour plane segments
                  int counter = 0;
                  for(unsigned long npoint = 0; npoint < refCont.nPoints(nplane); npoint++)
                    {
                      //Calculate next point index
                      unsigned long nextPoint = npoint+1;
                      if(npoint == refCont.nPoints(nplane)-1){
                        //Final point must be connected with initial one
                        nextPoint = 0; 
                      }
                      //Get actual and next point
                      double p1[3], p2[3];
                      refCont.getPoint(p1,nplane,npoint);
                      refCont.getPoint(p2,nplane,nextPoint);
                      
                      //Calculate the vector between these 2 points on z plane
                      double u[2];
                      u[0] = p2[0] - p1[0];
                      u[1] = p2[1] - p1[1];

                      if(yCenter > std::min(p1[1],p2[1])){
                        if(yCenter <= std::max(p1[1],p2[1])){
                          if(xCenter <= std::max(p1[0],p2[0])){
                            if(p1[1] != p2[1]){
                              double xinters = (yCenter-p1[1])*u[0]/u[1]+p1[0];
                              if(p1[0] == p2[0] || xCenter <= xinters){
                                counter++;
                              }
                            }
                          }
                        }
                      }
                    }

                  if(counter % 2 != 0)
                    {

                      if(voxelContour[index] == (int)icontour)
                        {
                          //This can happen because contour 'icontour' have 2 or more
                          //sequences at the same plane. If this sequences
                          //intersect, the program will consider that this sequences
                          //create a in-hole contour. So, we must remove this
                          //voxel from this contour and restore previous one.
                          voxelContour[index] = voxelContour2[index];
                          contourMask[index] = 0;
                        }
                      else
                        {
                          voxelContour[index] = icontour;
                          contourMask[index] = 1;
                        }
                    }
                }
             }
          }   
               
        //Copy actual voxel contour data to auxiliar array
        memcpy(voxelContour2,voxelContour,sizeof(int)*tnvox);
      }

      free(voxelContour2);	
      
      //Calculate number of voxels in each contour
      //Note that contourMasks does not take into account the contour priority; i.e., no overlap is considered.
      //On the contrary, voxelContour, and thus also nVoxContour, do take overlapping into account.
      for(unsigned long i = 0; i < tnvox; i++)
	if(voxelContour[i] > -1)
	  nVoxContour[voxelContour[i]]++;
    }

  return PEN_DICOM_SUCCESS;
}


void pen_dicom::clear(){

  //Contours
  contours.clear();
  nVoxContour.clear();
  
  if(voxelContour != nullptr){
    free(voxelContour);
    voxelContour = nullptr;
  }

  //Seeds
  seeds.clear();

  //Image
  if(dicomImage != nullptr){
    free(dicomImage);
    dicomImage = nullptr;
  }

  //Set null image modality
  imageModality.assign("NULL");
  //Set default origin
  dicomOrigin[0] = 0.0;
  dicomOrigin[1] = 0.0;
  dicomOrigin[2] = 0.0;
  originalOrigin[0] = dicomOrigin[0];
  originalOrigin[1] = dicomOrigin[1];
  originalOrigin[2] = dicomOrigin[2];
  
  //Set default dicom voxel dimensions
  dvox_x = 1.0;
  dvox_y = 1.0;
  dvox_z = 1.0;
  voxVol = 1.0;
  //Set default number of voxels
  nvox_x = 0;
  nvox_y = 0;
  nvox_z = 0;
  tnvox = 0;
  //Set dicom dimension
  dimDicom[0] = 0.0;
  dimDicom[1] = 0.0;
  dimDicom[2] = 0.0;

}

bool pen_dicom::checkImgModality(const char* Modality)
{
  //Check if image modality is avaible.
  
  if( strcmp(Modality, "US") == 0 || strcmp(Modality, "CT") == 0  || strcmp(Modality, "PT") == 0  || strcmp(Modality, "NM") == 0)
    {
      return true;
    }
  else
    {
      return false;
    }
}

int pen_dicom::transformContoursAndSeeds(const double* imageOrientation,
					  const unsigned verbose)
{
  //Transform positions of contour points and seeds
  //according to specified image orientation stored in 'imageOrientation'
  //and dicom origin
  const double defaultOrientation[6] = {1.0,0.0,0.0,
                                        0.0,1.0,0.0};

  bool defaultOr = true;
  for(int i = 0; i < 6; i++)
    if(defaultOrientation[i] != imageOrientation[i])
      {
	defaultOr = false;
	break;
      }
  
  if(!defaultOr)  //Rotation needed
    {

      //We modelize the rotation with three rotations: Rx(Phi)Rz(Theta)Rx(Gamma)

      //First Rx(Gamma) only moves Y axis. Next, Rx(Phi)Rz(Theta) move X axis to vector specified in dicom (imageOrientation[0],imageOrientation[1],imageOrientation[2]) 
      //We need to calculate inverse rotation Rx^-1(Phi)Rz^-1(Theta)Rx^-1(Gamma) = Rx(-Gamma)Rz(-Theta)Rx(-Phi)

      if(verbose > 1)
	printf(" *Rotation needed*\n");
      //Calculate rotation to convert x' from image orientation to (1.0,0.0,0.0)
      double modX = sqrt(pow(imageOrientation[0],2) + pow(imageOrientation[1],2) + pow(imageOrientation[2],2)); //mod(Vx,Vy,Vz)
      double modXyz = sqrt(pow(imageOrientation[1],2) + pow(imageOrientation[2],2));   //mod(0,Vy,Vz)

      double cosT;
      double cosP;
      
      if(modX == 0.0)
	{
	  if(verbose > 0)
	    printf("Oriantation-transform: Error: Invalid image orientation with mod equal to 0:\n");
	  return PEN_DICOM_INVALID_ORIENTATION;
	}
      if(modXyz == 0)
	{
	  // No y and z components.
	  cosP = 1.0;  // 0º
	}
      else
	{
	  cosP = imageOrientation[1]/modXyz;   //  Vy/mod(0,Vy,Vz)
	}
      
      cosT = imageOrientation[0]/modX;

      double sinT = sqrt(1.0 - cosT*cosT);
      double sinP = sqrt(1.0 - cosP*cosP);

      //Construct inverse rotation matrix Rz^-1 Rx^-1
      double RzInv[9] ={
			cosT,sinT,0.0,
			-sinT,cosT,0.0,
			0.0, 0.0,1.0
      };

      double RxPInv[9] ={
			 1.0,  0.0, 0.0,	
			 0.0, cosP,sinP,
			 0.0,-sinP,cosP	 
      };

      double RXZinv[9];

      matmul3(RzInv,RxPInv,RXZinv);

      //Rinv transform X' to (1,0,0). Now, we need to calculate gamma to transform Rinv*Y' to (0,1,0)

      double Ytrans[3];
      double auxYOrientation[3] = {imageOrientation[3],imageOrientation[4],imageOrientation[5]};
      matvect3(RXZinv,auxYOrientation,Ytrans);

      if(fabs(Ytrans[0]) > 1.0e-3)
	{
	  if(verbose > 0)
	    printf("loadDicom: Warning: After transformation Y axis have no zero 'x' component (%e)\n",Ytrans[0]);
	}

      double modYyz = sqrt(pow(Ytrans[1],2) + pow(Ytrans[2],2));
      double cosG = Ytrans[1]/modYyz;
      double sinG = sqrt(1.0 - cosG*cosG);

      double RxGInv[9] ={
			 1.0,  0.0, 0.0,	
			 0.0, cosG,sinG,
			 0.0,-sinG,cosG	 
      };

      double Rinv[9];

      matmul3(RxGInv,RXZinv,Rinv);  //Obtain final inverted rotation

      if(verbose > 1){
	printf("Rotation matrix:\n");
	printf(" %6.4e %6.4e %6.4e\n",Rinv[0],Rinv[1],Rinv[2]);
	printf(" %6.4e %6.4e %6.4e\n",Rinv[3],Rinv[4],Rinv[5]);
	printf(" %6.4e %6.4e %6.4e\n",Rinv[6],Rinv[7],Rinv[8]);
      }
      
      //Apply rotation to dicom origin
      matvect3(Rinv,dicomOrigin);

      //Because we set the (0,0,0) on the
      //corner of first voxel we must translate the origin
      double trans[3] = {dicomOrigin[0],
			 dicomOrigin[1],
			 dicomOrigin[2]};

	if(verbose > 1){
			printf("New dicom origin:\n");
			printf(" %6.4e %6.4e %6.4e\n",
			       dicomOrigin[0],dicomOrigin[1],dicomOrigin[2]);
	}      
      //Apply rotation and translation to contours points and seeds positions

      printf(" *Rotate and translate contour points* \n");
      for(size_t j = 0; j<contours.size(); j++)
	{
	  for(unsigned i = 0; i < contours[j].NPlanes(); i++)
	    {
	      for(unsigned long k = 0; k < contours[j].nPoints(i); k++)
		{
		  double p[3];
		  contours[j].getPoint(p,i,k);
		  matvect3(Rinv,p);
		  p[0] -= trans[0];
		  p[1] -= trans[1];
		  p[2] -= trans[2];
		  contours[j].updatePoints(i,k,1,p);	  
		}
	    }
	}

      if(verbose > 1)
	printf(" *Rotate and translate seed positions and directions* \n");
      for(size_t j = 0; j < seeds.size(); j++)
	{
	  for(unsigned long i = 0; i < seeds[j].nPoints(); i++) 
	    {
	      double p[3];
	      seeds[j].getPos(p,i);
	      matvect3(Rinv,p);
	      p[0] -= trans[0];
	      p[1] -= trans[1];
	      p[2] -= trans[2];
	      seeds[j].setPositions(p,i);
	      double dir[3];
	      seeds[j].getDirection(dir,i);
	      matvect3(Rinv,dir);
	      seeds[j].setDirections(dir,i);
	      
	    }
	}
      
    }
  else //Only translate
    {  

      if(verbose > 1)
	printf(" *Only origin translation is required*\n");

      //Because we set the (0,0,0) on the
      //corner of first voxel we must translate the origin
      double trans[3] = {dicomOrigin[0],
			 dicomOrigin[1],
			 dicomOrigin[2]};
    
	if(verbose > 1)
	  printf(" *Translate contour points* \n");
      for(size_t j = 0; j<contours.size(); j++)
	{
	  for(unsigned i = 0; i < contours[j].NPlanes(); i++)
	    {
	      for(unsigned long k = 0; k < contours[j].nPoints(i); k++)
		{
		  double p[3];
		  contours[j].getPoint(p,i,k);
		  p[0] -= trans[0];
		  p[1] -= trans[1];
		  p[2] -= trans[2];


		  contours[j].updatePoints(i,k,1,p);	  
		}
	    }
	}

      if(verbose > 1)
	printf(" *Translate seed positions* \n");
      for(size_t j = 0; j < seeds.size(); j++)
	{
	  for(unsigned long i = 0; i < seeds[j].nPoints(); i++) 
	    {
	      double p[3];
	      seeds[j].getPos(p,i);
	      p[0] -= trans[0];
	      p[1] -= trans[1];
	      p[2] -= trans[2];
	      seeds[j].setPositions(p,i);
	    }
	}    
    }
  return PEN_DICOM_SUCCESS;
}

int pen_dicom::printContourMasks(const char* filename) const{

  if(filename == nullptr)
    return PEN_DICOM_ERROR_NULL_FILENAME;

  for(size_t imask = 0; imask < contourMasks.size(); ++imask){

    //Get mask
    const std::vector<unsigned char>& mask = contourMasks[imask];
    
    //Create a file to store contour mask
    std::string sfilename(filename);
    sfilename.append("_");
    //remove white spaces
    std::string auxstr(contours[imask].name.c_str());
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
    
    fprintf(OutMask,"# PenRed MASK DATA\n");
    fprintf(OutMask,"# Inner voxels are flagged with a '1' in the mask.\n"
	    "# Instead, outer voxels are flagged with a '0'.\n");
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# Contour Name:\n");    
    fprintf(OutMask,"#    %s\n",contours[imask].name.c_str());
    fprintf(OutMask,"#\n");    
    fprintf(OutMask,"# Number of voxels:\n");
    fprintf(OutMask,"#    %lu\n",static_cast<unsigned long>(mask.size()));
    fprintf(OutMask,"# Number of inner voxels:\n");
    fprintf(OutMask,"#    %lu\n", nInner);    
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# Volume (cm^3):\n");
    fprintf(OutMask,"#    %E\n",
	    voxVol*static_cast<double>(nInner));
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# n voxel |  mask \n");
    
    for(size_t i = 0; i < mask.size(); ++i){
      fprintf(OutMask,"  %7lu     %u\n",i,mask[i]);
    }
    
    fprintf(OutMask,"#\n");
    fprintf(OutMask,"# End of mask data\n");
    fclose(OutMask);  

  }
    
  return PEN_DICOM_SUCCESS;
}

int pen_dicom::printContourMasksMHD(const char* filename) const{

  if(filename == nullptr)
    return PEN_DICOM_ERROR_NULL_FILENAME;

  FILE* OutMask = nullptr;

  for(size_t imask = 0; imask < contourMasks.size(); ++imask){

    //Get mask
    const std::vector<unsigned char>& mask = contourMasks[imask];
    
    //Create a file to store contour mask
    std::string sfilename(filename);
    sfilename.append("_");
    //remove white spaces
    std::string auxstr(contours[imask].name.c_str());
    auxstr.erase(std::remove(auxstr.begin(),auxstr.end(),' '),auxstr.end());
    //sfilename.append(contours[imask].name.c_str());
    sfilename.append(auxstr.c_str());
        
      //**************************
      //*  RAW IMAGE OF THE MASK *
      //**************************

      // Write header .mhd
  
      std::string sfilenameORIG=sfilename;
      sfilename.append(".mhd");
      OutMask = fopen(sfilename.c_str(),"w");
      if(OutMask == nullptr){
        printf("****************************************************\n");
        printf("pen_dicom: Error: Cannot open output mhd image mask\n");
        printf("****************************************************\n");
        return PEN_DICOM_ERROR_CREATING_FILE;
      }

      fprintf(OutMask,"ObjectType = Image\n");
      fprintf(OutMask,"NDims = 3\n");
      fprintf(OutMask,"BinaryData = True\n");
      fprintf(OutMask,"BinaryDataByteOrderMSB = False\n");
      fprintf(OutMask,"CompressedData = False\n");
      fprintf(OutMask,"TransformMatrix = 1 0 0 0 1 0 0 0 1\n");
      fprintf(OutMask,"Offset = %6.2f %6.2f %6.2f\n",(dicomOrigin[0]+dvox_x/2.0)*10.,(dicomOrigin[1]+dvox_y/2.0)*10.,(dicomOrigin[2]+dvox_z/2.0)*10.);
      fprintf(OutMask,"CenterOfRotation = 0 0 0\n");
      fprintf(OutMask,"AnatomicalOrientation = RAI\n");
      fprintf(OutMask,"ElementSpacing = %10.6f %10.6f %10.6f\n",dvox_x*10.,dvox_y*10.,dvox_z*10.);
      fprintf(OutMask,"DimSize = %lu %lu %lu\n",nvox_x,nvox_y,nvox_z);
      fprintf(OutMask,"ElementType = MET_UCHAR\n");

      sfilename = sfilenameORIG;
      std::size_t found = sfilename.find_last_of("/\\");
      fprintf(OutMask,"ElementDataFile = %-s\n",((sfilename.substr(found+1)).append(".raw")).c_str());
      fclose(OutMask);

      // Write data .raw

      sfilename = sfilenameORIG;

      sfilename.append(".raw");
      OutMask = fopen(sfilename.c_str(),"wb");
      if(OutMask == nullptr){
        printf("****************************************************\n");
        printf("pen_dicom: Error: Cannot open output raw image mask\n");
        printf("****************************************************\n");
        return PEN_DICOM_ERROR_CREATING_FILE;
      }
     

      //Write data
    
      for(size_t i = 0; i < mask.size(); ++i){
          unsigned char q = static_cast<unsigned char>(mask[i]);
          fwrite(&q,sizeof(q),1,OutMask);
      }

      fclose(OutMask);  
    }
  
  return PEN_DICOM_SUCCESS;
}

int pen_dicom::printContours(const char* filename) const{

  if(filename == nullptr)
    return PEN_DICOM_ERROR_NULL_FILENAME;
  
  //Create a file to store contours data
  FILE* OutContorns = nullptr;
  OutContorns = fopen(filename,"w");
  if(OutContorns == nullptr){
    return PEN_DICOM_ERROR_CREATING_FILE;
  }
		    
  fprintf(OutContorns,"# PenRed CONTOUR DATA\n");
  fprintf(OutContorns,"#\n");
  fprintf(OutContorns,"# Number of contours:\n");
  fprintf(OutContorns,"#    %lu\n",contours.size());
  fprintf(OutContorns,"#\n");
  fprintf(OutContorns,"# Note that overlapping/contour priority is taking into account.\n");
  fprintf(OutContorns,"#\n");
  fprintf(OutContorns,"# Name, number of voxels in contour, volume of each contour (cm^3), number of contour points\n");
  fprintf(OutContorns,"#\n");
  for(size_t i = 0; i < contours.size(); i++)
    {
      //Count number of contour points
      unsigned long nPoints = 0;
      for(unsigned j = 0; j < contours[i].NPlanes(); j++)
	{
	  nPoints += contours[i].nPoints(j);
	}
      
      fprintf(OutContorns,"#    Contour %li: %s %lu %.5E %ld\n",
	      i, contours[i].name.c_str(), nVoxContour[i],
	      nVoxContour[i]*voxVol,nPoints);
    }
  fprintf(OutContorns,"#\n");
  fprintf(OutContorns,"# contour points data\n");
  fprintf(OutContorns,"# contour number | sequence number | Z Position | X Position | Y Position |\n");

  for(size_t i = 0; i < contours.size(); i++)
    {
      for(unsigned j = 0; j < contours[i].NPlanes(); j++)
	{
	  unsigned long planePoints = contours[i].nPoints(j);
	  for(unsigned long k = 0; k < planePoints; k++)
	    {
	      //Get point
	      double p[3];
	      contours[i].getPoint(p,j,k);
	      fprintf(OutContorns," %7lu          %8u          %.5E  %.5E  %.5E\n",
		      i,j+1,p[2],p[0],p[1]);
	    }
	}
    }
  fprintf(OutContorns,"#\n");
  fprintf(OutContorns,"# End of contour data\n");
  fclose(OutContorns);  

  return PEN_DICOM_SUCCESS;
}

int pen_dicom::printSeeds(const char* filename) const{

  if(filename == nullptr)
    return PEN_DICOM_ERROR_NULL_FILENAME;
  
  //Create a file to store contours data
  FILE* OutSeeds = nullptr;
  OutSeeds = fopen(filename,"w");
  if(OutSeeds == nullptr){
    return PEN_DICOM_ERROR_CREATING_FILE;
  }
  
  fprintf(OutSeeds,"# \n");
  fprintf(OutSeeds,"# Number of seed types:\n");
  fprintf(OutSeeds,"# %lu \n",seeds.size());
  fprintf(OutSeeds,"# Seeds positions/movements (x,y,z) weight:\n");
  fprintf(OutSeeds,"# \n");

  for(unsigned long i = 0; i < seeds.size(); i++)
    {
      fprintf(OutSeeds,"#   Type %ld: %ld control points (only points with no "
	      "null weight will be printed)\n", i+1,seeds[i].nPoints());
      fprintf(OutSeeds,"#\n");
      
      //Stores mass center of seeds positions and steps
      double CM[3] = {0.0,0.0,0.0};  
      long int contChannel = 0;
      unsigned long nCP = seeds[i].nPoints();
      for(unsigned long j = 0; j < nCP; j++)
	{
	  double wght = seeds[i].getWeight(j);
	  double pos[3];
	  seeds[i].getPos(pos,j);
	  double d = seeds[i].getDistance(j);
	  double dir[3];
	  seeds[i].getDirection(dir,j);
	  if(d == -1.0)
	    {
	      contChannel++;
	      fprintf(OutSeeds,"#     Channel %ld, first position: "
		      "(%9.5e,%9.5e,%9.5e)\n",
		      contChannel,pos[0],pos[1],pos[2]);
	    }
	  else
	    {
	      //Skip 0 weight positions/movements
	      if(wght == 0.0)
		{
		  continue;
		}

	      double center[3];  //Store the center of the step traveled by the seed
	      double posPrev[3]; //Previous position
	      if(j > 0)
		seeds[i].getPos(posPrev,j-1);
	      else
		seeds[i].getPos(posPrev,0);
	      
	      //Check if the seed stay on this position
	      //or is moving to next checkpoint
	      if(d > 0.0)
		{
		  fprintf(OutSeeds,"#         move to: ");
		  center[0] = posPrev[0]+d*dir[0];
		  center[1] = posPrev[1]+d*dir[1];
		  center[2] = posPrev[2]+d*dir[2];
		}
	      else
		{
		  fprintf(OutSeeds,"#         stay on: ");
		  center[0] = pos[0];
		  center[1] = pos[1];
		  center[2] = pos[2];
		}
	      fprintf(OutSeeds,"%9.5e %9.5e %9.5e weight: %9.5e\n",
		      pos[0],pos[1],pos[2],wght);
	      CM[0] += center[0]*wght;
	      CM[1] += center[1]*wght;
	      CM[2] += center[2]*wght;
	    }
	}
      fprintf(OutSeeds,"#   Positions mass center: %9.5e %9.5e %9.5e\n",CM[0],CM[1],CM[2]);
      fprintf(OutSeeds,"#\n");
    }
  fclose(OutSeeds);

  return PEN_DICOM_SUCCESS;
}

int pen_dicom::printImage(const char* filename) const{

  if(filename == nullptr)
    return PEN_DICOM_ERROR_NULL_FILENAME;
  
  //Create a file to store contours data
  FILE* OutImage = nullptr;
  OutImage = fopen(filename,"w");
  if(OutImage == nullptr){
    return PEN_DICOM_ERROR_CREATING_FILE;
  }

  fprintf(OutImage,"# \n");
  fprintf(OutImage,"# DICOM processed with PenRed\n");
  fprintf(OutImage,"# Nº of voxels (nx,ny,nz):\n");  
  fprintf(OutImage,"# %5lu %5lu %5lu\n",nvox_x,nvox_y,nvox_z);
  fprintf(OutImage,"# %8.5E %8.5E %8.5E\n",dvox_x,dvox_y,dvox_z);
  fprintf(OutImage,"# DICOM image data:\n");

  //Iterate over Z planes
  for(unsigned long k = 0; k < nvox_z; k++){
    unsigned long indexZ = nvox_xy*k;
    fprintf(OutImage,"# Index Z = %4lu\n",k);

    //Iterate over rows
    for(unsigned long j = 0; j < nvox_y; j++){
      unsigned long indexYZ = indexZ + j*nvox_x;
      fprintf(OutImage,"# Index Y = %4lu\n",j);

      //Iterate over columns
      for(unsigned long i = 0; i < nvox_x; i++){
	unsigned long ivoxel = indexYZ + i;

	//Save voxel X Y and intensity
	fprintf(OutImage,"%12.5E %12.5E %12.5E\n", i*dvox_x, j*dvox_y, dicomImage[ivoxel]);
      }
      
    }
    //Set a space between planes
    fprintf(OutImage,"\n\n\n");    
  }

  return PEN_DICOM_SUCCESS;
}

int pen_dicom::printContourVox(const char* filename) const{

  if(filename == nullptr)
    return PEN_DICOM_ERROR_NULL_FILENAME;
  
  //Create a file to store contours data
  FILE* OutImage = nullptr;
  OutImage = fopen(filename,"w");
  if(OutImage == nullptr){
    return PEN_DICOM_ERROR_CREATING_FILE;
  }

  fprintf(OutImage,"# PenRed CONTOUR VOXEL DATA\n");
  fprintf(OutImage,"#\n");
  fprintf(OutImage,"# Number of contours:\n");
  fprintf(OutImage,"#    %lu\n",contours.size());
  fprintf(OutImage,"#\n");
  fprintf(OutImage,"# Note that overlapping/contour priority is taking into account.\n");
  fprintf(OutImage,"#\n");
  fprintf(OutImage,"# Name, number of voxels in contour, volume of each contour (cm^3), number of contour points\n");
  fprintf(OutImage,"#\n");
  for(size_t i = 0; i < contours.size(); i++)
    {
      //Count number of contour points
      unsigned long nPoints = 0;
      for(unsigned j = 0; j < contours[i].NPlanes(); j++)
        {
          nPoints += contours[i].nPoints(j);
        }
      
      fprintf(OutImage,"#    Contour %li: %s %lu %.5E %ld\n",
              i, contours[i].name.c_str(), nVoxContour[i],
              nVoxContour[i]*voxVol,nPoints);
    }

  fprintf(OutImage,"# \n");
  fprintf(OutImage,"# Nº of voxels (nx,ny,nz):\n");  
  fprintf(OutImage,"# %5lu %5lu %5lu\n",nvox_x,nvox_y,nvox_z);
  fprintf(OutImage,"# %8.5E %8.5E %8.5E\n",dvox_x,dvox_y,dvox_z);
  fprintf(OutImage,"# DICOM contour mask image:\n");

  //Iterate over Z planes
  for(unsigned long k = 0; k < nvox_z; k++){
    unsigned long indexZ = nvox_xy*k;
    fprintf(OutImage,"# Index Z = %4lu\n",k);

    //Iterate over rows
    for(unsigned long j = 0; j < nvox_y; j++){
      unsigned long indexYZ = indexZ + j*nvox_x;
      fprintf(OutImage,"# Index Y = %4lu\n",j);

      //Iterate over columns
      for(unsigned long i = 0; i < nvox_x; i++){
	unsigned long ivoxel = indexYZ + i;

	//Save voxel X Y and intensity
	fprintf(OutImage,"%12.5E %12.5E %3d\n",
		i*dvox_x, j*dvox_y, voxelContour[ivoxel]);
      }
      
    }
    //Set a space between planes
    fprintf(OutImage,"\n\n\n");    
  }

  return PEN_DICOM_SUCCESS;
}


pen_dicom::~pen_dicom(){
  clear();
}
//----------------------
// Auxiliar functions
//----------------------

void matvect3(const double* A, const double* V, double* C)
{
  //*******************************************************************
  //*    Computes the product of 3x3 matrix with 3x1 vector.          *
  //*******************************************************************

  for(int i = 0; i < 3; i++) //Matrix row
    C[i] = A[i*3]*V[0] + A[i*3+1]*V[1] + A[i*3+2]*V[2];  
}

void matvect3(const double* A, double* V)
{
  //*******************************************************************
  //*    Computes the product of 3x3 matrix with 3x1 vector in place. *
  //*******************************************************************

  double auxV[3];
  
  for(int i = 0; i < 3; i++) //Matrix row
    auxV[i] = A[i*3]*V[0] + A[i*3+1]*V[1] + A[i*3+2]*V[2];

  memcpy(V,auxV,sizeof(double)*3);
}

void matmul3(double* A, double* B, double* C)
{
  //*******************************************************************
  //*    Computes the product of the two input matrices 3x3.          *
  //*******************************************************************

  int i,j,k;

  //Init C to zero
  for(i = 0; i < 9; i++)
    C[i] = 0.0;

  for(i = 0; i < 3; i++)    
    for(k = 0; k < 3; k++)	
      for(j = 0; j < 3; j++)	
	C[i*3+j] += A[i*3+k]*B[k*3+j];	
}
