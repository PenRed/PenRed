
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


#ifndef __PEN_DICOM__
#define __PEN_DICOM__

#include <vector>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <dirent.h>
#include <cmath>
#include <numeric>

/* make sure OS specific configuration is included first */
#include "dcmtk/config/osconfig.h"    
#include "dcmtk/dcmdata/dctk.h" 
#include "dcmtk/dcmdata/dcpxitem.h"
#include "dcmtk/dcmimgle/dcmimage.h"

#include "math_classes.hh"

enum pen_dicom_status{
		      PEN_DICOM_SUCCESS = 0,
		      PEN_DICOM_FOLDER_NOT_SPECIFIED,
		      PEN_DICOM_FOLDER_NOT_FOUND,
		      PEN_DICOM_MULTIPLE_MODALITIES,
		      PEN_DICOM_BAD_IMAGE_OPEN,
		      PEN_DICOM_BAD_READ,
		      PEN_DICOM_BAD_READ_COLUMNS,
		      PEN_DICOM_BAD_READ_ROWS,
		      PEN_DICOM_BAD_READ_IMAGE_POSITION,
		      PEN_DICOM_BAD_READ_PIXEL_SPACING,
		      PEN_DICOM_BAD_READ_SLICE_THICKNESS,
		      PEN_DICOM_NO_DICOM_FOUND,
		      PEN_DICOM_SMALL_VOXELS,
		      PEN_DICOM_BAD_ALLOCATION,
		      PEN_DICOM_NON_MONOCHROME_IMAGE,
		      PEN_DICOM_MISMATCH_DIMENSIONS,
		      PEN_DICOM_BAD_READ_PIXEL_REPRESENTATION,
		      PEN_DICOM_BAD_PIXEL_REPRESENTATION,
		      PEN_DICOM_BAD_READ_PIXEL_DATA,
		      PEN_DICOM_BAD_PIXEL_CONVERSION,
		      PEN_DICOM_MISMATCH_EXPECTED_READ,
		      PEN_DICOM_BAD_READ_ORIGIN,
		      PEN_DICOM_NO_DICOM_LOADED,
		      PEN_DICOM_ERROR_REOPENING_DICOM,
		      PEN_DICOM_NVOXELS_MISMATCH,
		      PEN_DICOM_INVALID_ORIENTATION,
		      PEN_DICOM_ERROR_CREATING_FILE,
		      PEN_DICOM_ERROR_NULL_FILENAME,
		      PEN_DICOM_ERROR_SPACING_MISMATCH,
		      PEN_DICOM_IMPOSSIBLE_ERROR
};

struct pen_contour{

private:
  //Store the nummber of Zplanes
  unsigned nPlanes; 
  //Store number of points of each plane
  unsigned long* nPlanePoints;
  //Store contours points (x,y,z) of each plane of
  //[nplane][x1,y1,z1,x2,y2,z2...]  
  double** contourPoints;
  
public:
  
  //Store contour name
  std::string name;
  //Store contour priority
  float priority;

  pen_contour();
  pen_contour(const pen_contour& c);
  void setPlanes(const unsigned np);
  void setPoints(const unsigned nplane,
		 const unsigned long npoints,
		 const double* points = nullptr);
  void updatePoints(const unsigned nplane,
		    const unsigned long init,
		    const unsigned long npoints,
		    const double* points);  
  void convertPoints(const double factor);
  bool hasPoints() const;
  
  inline unsigned NPlanes() const {return nPlanes;}
  inline unsigned long nPoints(const unsigned np) const {
    if(np >= nPlanes)
      throw std::out_of_range("pen_contour:nPoints: Out of range");
    return nPlanePoints[np];
  }
  inline void getPoint(double* p,
		       const unsigned nplane,
		       const unsigned long npoint) const {
    if(nplane >= nPlanes){
      char err[300];
      sprintf(err,"pen_contour:getPoint: Plane out of range (%u/%u)",nplane,nPlanes);
      throw std::out_of_range(err);
    }

    if(npoint >= nPlanePoints[nplane]){
      char err[300];
      sprintf(err,"pen_contour:getPoint: Point out of range (%lu/%lu)",
	      npoint,nPlanePoints[nplane]);
      throw std::out_of_range(err);
    }

    unsigned long np3 = npoint*3;
    p[0] = contourPoints[nplane][np3];
    p[1] = contourPoints[nplane][np3+1];
    p[2] = contourPoints[nplane][np3+2];
    
  }
  
  inline const double* readPoints(const unsigned nplane) const{
    //Check if plane number is in range
    if(nplane >= nPlanes){
      char err[300];
      sprintf(err,"pen_contour:readPoints: Plane out of range (%u/%u)",nplane,nPlanes);
      throw std::out_of_range(err);
    }

    return contourPoints[nplane];
  }

  pen_contour& operator=(const pen_contour& c);
  
  void clear();
  ~pen_contour();
  
};

struct pen_seed{

private:
  //Store control points positions (x1,y1,z1,x2,y2,z2...)
  double* positions;
  //Store control points weights
  double* weights;
  //Store distances between succesive control points
  double* distances;
  //Store directions on each control point
  double* directions;
  //Store number of control points
  unsigned long controlPoints;
  //Store array size
  unsigned long size;

public:
  //Store the sum of sequence time events
  double time;
  long int ID;
    
  pen_seed();
  pen_seed(const pen_seed& c);

  pen_seed& operator=(const pen_seed& c);
  
  void clear();
  void setCP(unsigned long n);
  void setPositions(const double* pos, const unsigned long ni);
  void setPositions(const double* pos,
		    const unsigned long ni,
		    const unsigned long npos);
  void setWeights(const double wght, const unsigned long ni);
  void setWeights(const double* wght,
		  const unsigned long ni,
		  const unsigned long npos);
  void setDistances(const double ds, const unsigned long ni);
  void setDistances(const double* ds,
		    const unsigned long ni,
		    const unsigned long npos);
  void setDirections(const double* dir, const unsigned long ni);
  void setDirections(const double* dir,
		     const unsigned long ni,
		     const unsigned long npos);

  inline double normWeights(){
    double sum = 0.0;
    for(unsigned long i = 0; i < controlPoints; i++)
      if(weights[i] > 0.0)
	sum += weights[i];
    for(unsigned long i = 0; i < controlPoints; i++)
      if(weights[i] > 0.0)
	weights[i] /= sum;
    return sum;
  }

  inline unsigned long nPoints() const {return controlPoints;}
  inline void getPos(double* pos, const unsigned long n) const {
    if(n >= controlPoints)
      throw std::out_of_range("pen_seedSequence:getPos: Out of range");
    unsigned long n3 = 3*n;
    pos[0] = positions[n3];
    pos[1] = positions[n3+1];
    pos[2] = positions[n3+2];
  }
  inline double getWeight(const unsigned long n) const {
    if(n >= controlPoints)
      throw std::out_of_range("pen_seedSequence:getPos: Out of range");
    return weights[n];
  }
  inline double getDistance(const unsigned long n) const {
    if(n >= controlPoints)
      throw std::out_of_range("pen_seedSequence:getPos: Out of range");
    return distances[n];
  }
  inline void getDirection(double* dir, const unsigned long n) const {
    if(n >= controlPoints)
      throw std::out_of_range("pen_seedSequence:getPos: Out of range");
    unsigned long n3 = 3*n;
    dir[0] = directions[n3];
    dir[1] = directions[n3+1];
    dir[2] = directions[n3+2];
  }
  inline const double* readPositions() const{
    return positions;
  }
  inline const double* readWeights() const{
    return weights;
  }
  inline const double* readDistances() const{
    return distances;
  }
  inline const double* readDirections() const{
    return directions;
  }
  
  
  ~pen_seed();
};

struct pen_ctData{
  // 0-> SEQUENCED, 1-> SPIRAL, 2-> CONSTANT_ANGLE, 3-> STATIONARY, 4->FREE, 5-> UNKNOWN
  std::string acquisitionTypeString;
  unsigned acquisitionType;

  double revolutionTime; //In seconds
  double singleCollimationWidth; //In cm
  double totalCollimationWidth; //In cm
  double tableSpeed; //In cm/s
  double tableFeedPerRotation; //In cm
  double spiralPitchFactor;
  double dataCollectionDiameter; //In cm

  static const unsigned maxFilters = 1000;
  // 0->NONE, 1->STRIP, 2->WDGE, 3->BUTTERFLY, 4->MULTIPLE, 5->FLAT, 6->OTHERS
  std::string filterTypeString;
  unsigned filterType;
  std::vector<std::string> filterMaterial;
  std::vector<double> focalSpots; //In cm
  

  double CTDIvol; //In mGy
  double reconstructionDiameter; //In cm
  double distanceSourceToDetector; //In cm
  double distanceSourceToDataCollectionCenter; //In cm

  // 0-> CW (clockwise (decreasing angle))
  // 1-> CC (counter-clockwise (increasing angle))
  // 2-> UNKNOWN
  std::string rotationDirectionString;
  unsigned rotationDirection;

  long int exposureTime; //In msec
  long int xRayTubeCurrent; //In mA
  long int exposure; //In mAs
  long int exposureInuAs; //In muAs
  long int generatorPowerInt; //In kW

  vector3D<double> dataCollectionCenterPatient; //In cm
  vector3D<double> reconstructionTargetCenterPatient; //In cm
  vector3D<double> isocenterPosition; //In cm

  bool multiEnergyAcquisition; //If the acquisition has been done with multiple energies

  double kvp;
  
  void load(DcmDataset* metainfo);
};

class pen_dicom{

private:

  std::string dicomDirname;
  std::string imageModality;
  
  double dicomOrigin[3];
  //Voxel dimensions
  double dvox_x, dvox_y, dvox_z;
  double voxVol;
  //Number of voxels
  unsigned long nvox_x, nvox_y, nvox_z;
  unsigned long nvox_xy, tnvox;
  double dimDicom[3];
        
  ////////////////
  //Contour vars//
  ////////////////

  //Store contours information
  std::vector<pen_contour> contours;

  //Store number of voxels in each contour
  std::vector<unsigned long> nVoxContour;
    
  //Save container contour for each voxel (-1 for non contour)
  int* voxelContour;
  std::vector<std::vector<unsigned char>> contourMasks;

  ////////////////
  // Seeds vars //
  ////////////////

  std::vector<pen_seed> seeds;

  ////////////////
  //  CT  vars  //
  ////////////////

  pen_ctData ctData;

public:
  inline bool isCT() const {return imageModality.compare("CT") == 0;}
  inline std::string ctReadAcquisitionType() const {return ctData.acquisitionTypeString;}
  inline unsigned ctReadAcquisitionTypeIndex() const {return ctData.acquisitionType;}
  inline double ctReadRevolutionTime() const {return ctData.revolutionTime;}
  inline double ctReadTotalCollimationWidth() const {return ctData.totalCollimationWidth;}
  inline double ctReadTableSpeed() const {return ctData.tableSpeed;}
  inline double ctReadTableFeedPerRotation() const {return ctData.tableFeedPerRotation;}
  inline double ctReadSpiralPitchFactor() const {return ctData.spiralPitchFactor;}
  inline double ctReadDataCollectionDiameter() const {return ctData.dataCollectionDiameter;}
  inline std::string ctReadFilterType() const {return ctData.filterTypeString;}
  inline double ctReadFilterTypeIndex() const {return ctData.filterType;}
  inline size_t ctReadNFilterMaterial() const {return ctData.filterMaterial.size();}
  inline std::string ctReadFilterMaterial(const unsigned i) const {return ctData.filterMaterial[i];}
  
  inline size_t ctReadNFocalSpots() const {return ctData.focalSpots.size();}
  inline double ctReadFocalSpots(const unsigned i) const {return ctData.focalSpots[i];}
  
  inline double ctReadCTDIvol() const {return ctData.CTDIvol;}
  inline double ctReadReconstructionDiameter() const {return ctData.reconstructionDiameter;}
  inline double ctReadDistanceSourceToDetector() const {return ctData.distanceSourceToDetector;}
  inline std::string ctReadRotationDirection() const {return ctData.rotationDirectionString;}
  inline unsigned ctReadRotationDirectionIndex() const {return ctData.rotationDirection;}
  inline long int ctReadExposureTime() const {return ctData.exposureTime;}
  inline long int ctReadXRayTubeCurrent() const {return ctData.xRayTubeCurrent;}
  inline long int ctReadExposure() const {return ctData.exposure;}
  inline long int ctReadExposureInuAs() const {return ctData.exposureInuAs;}
  inline long int ctReadGeneratorPowerInt() const {return ctData.generatorPowerInt;}
  inline vector3D<double> ctReadDataCollectionCenterPatient() const {
    return ctData.dataCollectionCenterPatient;
  }
  inline vector3D<double> ctReadReconstructionTargetCenterPatient() const {
    return ctData.reconstructionTargetCenterPatient;
  }
  inline vector3D<double> ctReadIsocenterPosition() const {
    return ctData.isocenterPosition;
  }
  inline bool ctReadMultipleEnergyUsed() const {return ctData.multiEnergyAcquisition;}
  inline double ctReadKVP() const {return ctData.kvp;}
  inline double ctReadSingleCollimationWidth() const {return ctData.singleCollimationWidth;}
protected:
  
  ////////////////
  // Image vars //
  ////////////////

  double* dicomImage;    //Store dicom image

  double originalOrigin[3];
  
  int transformContoursAndSeeds(const double* imageOrientation,
				 const unsigned verbose);    
protected:
public:

  pen_dicom();
  int loadDicom(const char* dirName,
		const unsigned verbose,
		const bool onlyMetadata = false);  
  int assignContours();
  
  
  static bool checkImgModality(const char* Modality);
  inline void setContourPriority(const unsigned icontour,
				 const float priority){
    if(icontour >= contours.size()){
      throw std::out_of_range("pen_dicom:setContourPriority: Out of range");
    }
    contours[icontour].priority = priority;
  }
  inline void setContourName(const unsigned icontour,
			     const char* name){
    if(icontour >= contours.size()){
      throw std::out_of_range("pen_dicom:setContourPriority: Out of range");
    }
    contours[icontour].name.assign(name);
  }
  
  inline unsigned getContourIndex(const char* name) const {
    for(unsigned long i = 0; i < contours.size(); i++)
      if(contours[i].name.compare(name) == 0)
	return i;
    return contours.size();
  }
  inline float getContourPriority(const unsigned icontour) const {
    if(icontour >= contours.size()){
      throw std::out_of_range("pen_dicom:setContourPriority: Out of range");
    }
    return contours[icontour].priority;
  }
  
  inline const double* readImage() const {return dicomImage;}
  inline const int* readContour() const {return voxelContour;}
  inline const std::vector<std::vector<unsigned char>>& readContourMasks() const{
    return contourMasks;
  }

  inline unsigned long nContours() const {return contours.size();}
  inline unsigned long nSeeds() const {return seeds.size();}
  
  inline unsigned long getNX() const {return nvox_x;}
  inline unsigned long getNY() const {return nvox_y;}
  inline unsigned long getNZ() const {return nvox_z;}
  inline unsigned long getNVox() const {return tnvox;}

  inline double getDX() const {return dvox_x;}
  inline double getDY() const {return dvox_y;}
  inline double getDZ() const {return dvox_z;}
  inline double getVoxVol() const {return voxVol;}

  inline double getOriginX() const {return dicomOrigin[0];}
  inline double getOriginY() const {return dicomOrigin[1];}
  inline double getOriginZ() const {return dicomOrigin[2];}

  inline const pen_contour& contour(unsigned long icont) const {
      return contours[icont];
  }

  inline const pen_seed& seed(unsigned long icont) const {
    return seeds[icont];
  }

  int printContours(const char* filename) const;
  int printContourMasks(const char* filename) const;
  int printContourMasksMHD(const char* filename) const;
  int printSeeds(const char* filename) const;
  int printImage(const char* filename) const;
  int printContourVox(const char* filename) const;
  void clear();

  ~pen_dicom();
};


//----------------------
// Auxiliar functions
//----------------------

void matvect3(const double* A, const double* V, double* C);

void matvect3(const double* A, double* V);

void matmul3(double* A, double* B, double* C);

#endif
