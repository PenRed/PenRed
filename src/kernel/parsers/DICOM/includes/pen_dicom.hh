
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


#ifndef __PEN_DICOM__
#define __PEN_DICOM__

#include <vector>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <dirent.h>
#include <cmath>

/* make sure OS specific configuration is included first */
#include "dcmtk/config/osconfig.h"    
#include "dcmtk/dcmdata/dctk.h" 
#include "dcmtk/dcmdata/dcpxitem.h"
#include "dcmtk/dcmimgle/dcmimage.h"

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
		      PEN_DICOM_IMPOSSIBLE_ERROR
};

struct pen_contour{

private:
  //Store the nummber of Zplanes
  unsigned nPlanes; 
  //Store number of points of each plane
  unsigned* nPlanePoints;
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
		 const unsigned npoints,
		 const double* points = nullptr);
  void updatePoints(const unsigned nplane,
		    const unsigned init,
		    const unsigned npoints,
		    const double* points);  
  void convertPoints(const double factor);
  bool hasPoints() const;
  
  inline unsigned NPlanes() const {return nPlanes;}
  inline unsigned nPoints(const unsigned np) const {
    if(np >= nPlanes)
      throw std::out_of_range("pen_contour:nPoints: Out of range");
    return nPlanePoints[np];
  }
  inline void getPoint(double* p,
		       const unsigned nplane,
		       const unsigned npoint) const {
    if(nplane >= nPlanes){
      char err[300];
      sprintf(err,"pen_contour:getPoint: Plane out of range (%u/%u)",nplane,nPlanes);
      throw std::out_of_range(err);
    }

    if(npoint >= nPlanePoints[nplane]){
      char err[300];
      sprintf(err,"pen_contour:getPoint: Point out of range (%u/%u)",npoint,nPlanePoints[nplane]);
      throw std::out_of_range(err);
    }

    unsigned np3 = npoint*3;
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
  unsigned controlPoints;
  //Store array size
  unsigned size;

public:
  //Store the sum of sequence time events
  double time;
  int ID;
    
  pen_seed();
  pen_seed(const pen_seed& c);

  pen_seed& operator=(const pen_seed& c);
  
  void clear();
  void setCP(unsigned n);
  void setPositions(const double* pos, const unsigned ni);
  void setPositions(const double* pos, const unsigned ni, const unsigned npos);
  void setWeights(const double wght, const unsigned ni);
  void setWeights(const double* wght, const unsigned ni, const unsigned npos);
  void setDistances(const double ds, const unsigned ni);
  void setDistances(const double* ds, const unsigned ni, const unsigned npos);
  void setDirections(const double* dir, const unsigned ni);
  void setDirections(const double* dir, const unsigned ni, const unsigned npos);

  inline double normWeights(){
    double sum = 0.0;
    for(unsigned i = 0; i < controlPoints; i++)
      if(weights[i] > 0.0)
	sum += weights[i];
    for(unsigned i = 0; i < controlPoints; i++)
      if(weights[i] > 0.0)
	weights[i] /= sum;
    return sum;
  }

  inline unsigned nPoints() const {return controlPoints;}
  inline void getPos(double* pos, const unsigned n) const {
    if(n >= controlPoints)
      throw std::out_of_range("pen_seedSequence:getPos: Out of range");
    unsigned n3 = 3*n;
    pos[0] = positions[n3];
    pos[1] = positions[n3+1];
    pos[2] = positions[n3+2];
  }
  inline double getWeight(const unsigned n) const {
    if(n >= controlPoints)
      throw std::out_of_range("pen_seedSequence:getPos: Out of range");
    return weights[n];
  }
  inline double getDistance(const unsigned n) const {
    if(n >= controlPoints)
      throw std::out_of_range("pen_seedSequence:getPos: Out of range");
    return distances[n];
  }
  inline void getDirection(double* dir, const unsigned n) const {
    if(n >= controlPoints)
      throw std::out_of_range("pen_seedSequence:getPos: Out of range");
    unsigned n3 = 3*n;
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


class pen_dicom{

private:

  std::string dicomDirname;
  std::string imageModality;
  
  double dicomOrigin[3];
  //Voxel dimensions
  double dvox_x, dvox_y, dvox_z;
  double voxVol;
  //Number of voxels
  unsigned nvox_x, nvox_y, nvox_z;
  unsigned nvox_xy, tnvox;
  double dimDicom[3];
        
  ////////////////
  //Contour vars//
  ////////////////

  //Store contours information
  std::vector<pen_contour> contours;

  //Store number of voxels in each contour
  std::vector<unsigned> nVoxContour;
    
  //Save container contour for each voxel (-1 for non contour)
  int* voxelContour;
  ////////////////
  // Seeds vars //
  ////////////////

  std::vector<pen_seed> seeds;  
  
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
		const unsigned verbose);  
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
    for(unsigned i = 0; i < contours.size(); i++)
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

  inline unsigned nContours() const {return contours.size();}
  inline unsigned nSeeds() const {return seeds.size();}
  
  inline unsigned getNX() const {return nvox_x;}
  inline unsigned getNY() const {return nvox_y;}
  inline unsigned getNZ() const {return nvox_z;}
  inline unsigned getNVox() const {return tnvox;}

  inline double getDX() const {return dvox_x;}
  inline double getDY() const {return dvox_y;}
  inline double getDZ() const {return dvox_z;}
  inline double getVoxVol() const {return voxVol;}

  int printContours(const char* filename) const;
  int printSeeds(const char* filename) const;
  int printImage(const char* filename) const;
  int printContourVox(const char* filename) const;
  void clear();

  ~pen_dicom();
};


//----------------------
// Auxiliar functions
//----------------------



bool checkEndian();

inline Uint32 U32bitSwap(Uint32 x)
{
  return (((x>>24) & 0x000000ff) | ((x>>8) & 0x0000ff00) | ((x<<8) & 0x00ff0000) | ((x<<24) & 0xff000000));
}

template <class store> int bytes2int(const store* chunks, unsigned int chunkSize, unsigned int outputSize, unsigned int bitsStored, bool endianness, unsigned int highBit, Uint32& output, unsigned* usedChunks)
{

  // chunks -> Array with data to join
  // chunkSize -> size in bites of each element in chunks
  // outputSize -> size in bits of output integer (maximum 32 bits)
  // bitsStored -> number of bits used to construct the number
  // endianness ->  Specifies how the data in chunks is stored: 0 if is little-endian or 1 for big-endian
  // output -> output integer
  // usedChunks -> counter where the number of chunks used will be added

  //   1  2  3  4  5  
  //  [0][1][2][3][4]    big-endian, la primera posició (0) va més a l'esquerra. L'agafarem la primera

  //   1  2  3  4  5  
  //  [4][3][2][1][0]    little-endian, l'última posició (4) va més a l'esquerra. L'agafarem la primera
     
  if(outputSize > 32)
    return -1;

  if(outputSize < chunkSize)
    return -2;

  if(outputSize == 0 || chunkSize == 0)
    return -3;
  
  unsigned int chunks2add = outputSize/chunkSize;
  output = 0;

  if(endianness == 0)
    {
      for(int i = chunks2add-1; i >= 0; i--)
	{
	  output = (output << chunkSize) | chunks[i];
	}
    }
  else
    {
      for(unsigned i = 0; i < chunks2add; i++)
	{
	  output = (output << chunkSize) | chunks[i];
	}
    }
      
  //Check if stored endianness coincide with
  //machine ones.
  bool localEndianness = checkEndian();
  if(localEndianness != endianness)
    {
      //Convert stored number endianess
      output = U32bitSwap(output);
    }

  //Remove unused bits
  Uint32 unusedBits = outputSize-bitsStored;

  if(unusedBits > 0)
    {
      if(highBit == bitsStored-1)
	{
	  if(localEndianness == 0) //Little endian
	    {
	      //Remove x: [xxxxxfffffffff]
	      output = (output << unusedBits);
	      output = (output >> unusedBits);
	    }
	  else //Big endian
	    {
	      //Remove x: [fffffffffxxxxx]
	      output = (output >> unusedBits);	  
	      output = (output << unusedBits);	  
	    }
	}
      else if(highBit == outputSize-1)
	{
	  if(localEndianness == 0) //Little endian
	    {
	      //Remove x: [fffffffffxxxxx]
	      output = (output >> unusedBits);
	      output = (output << unusedBits);
	    }
	  else //Big endian
	    {
	      //Remove x: [xxxxxfffffffff]
	      output = (output << unusedBits);
	      output = (output >> unusedBits);
	    }
	}
      else{
	// No valid high bit value
	return -4;
      }
    }
  
  if(usedChunks != 0)
    *usedChunks += chunks2add;
  
  return 0;
}

template <typename type> int searchArray(int n, type* v, type element)
{
  //search "element" in array "v" of size "n" and return its position if exists.
  //If doesn't exist, return -1

  for(int i = 0; i < n; i++)
    {
      if(element == v[i])
	{
	  return i;
	}
    }
  return -1;
  
}

inline char* stringToLow(char* text)
{
  int c = 0;
  while(text[c]){text[c]=tolower(text[c]); c++;} //convert to lowecase
  return text;
}

void matvect3(const double* A, const double* V, double* C);

void matvect3(const double* A, double* V);

void matmul3(double* A, double* B, double* C);

#endif
