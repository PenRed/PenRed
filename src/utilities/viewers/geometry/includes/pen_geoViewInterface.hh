
//
//
//    Copyright (C) 2021-2022 Universitat de València - UV
//    Copyright (C) 2021-2022 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicent Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Olver Gil)
//

#ifndef __PEN_GEO_VIEWER_INTERFACE__
#define __PEN_GEO_VIEWER_INTERFACE__

#ifndef __PEN_GEO_VIEWER__

#include <cstdio>
#include <vector>

struct geoError{
  double from[3];
  double to[3];
  unsigned int iIBODY, iMAT; //Initial
  unsigned int eIBODY, eMAT; //Expected
  unsigned int fIBODY, fMAT; //Final
  bool endInVoid;
};
#else
struct geoError;
#endif

class pen_geoViewInterface{

private:

public:
  
  virtual double z2dir(const double u,
		       const double v,
		       const double w,
		       const double omega,
		       double rotation[9],
		       const double phiAux,
		       const double threshold) const = 0;
  
  virtual float z2dirf(const float u,
		       const float v,
		       const float w,
		       const float omega,
		       float rotation[9],
		       const float phiAux,
		       const float threshold) const = 0;
  
  virtual int init(const char* filename,
		   const unsigned verbose = 5) = 0;
  
  virtual void testX(std::vector<geoError>& errors,
		     const float x, const float y, const float z,
		     const float dy, const float dz,
		     const unsigned ny, const unsigned nz) const = 0;
  
  virtual void renderX(unsigned char* renderMat,unsigned int* renderBody,
		       const float x, const float y, const float z,
		       const float dy, const float dz,
		       const unsigned ny, const unsigned nz,
		       const unsigned nthreads = 1) const = 0;

  virtual void renderXtoLeft(unsigned char* renderMat,
			     unsigned int* renderBody,
			     const unsigned nPixels,				
			     const float x, const float y, const float z,
			     const float dy, const float dz,
			     const unsigned ny, const unsigned nz) const = 0;

  virtual void renderXtoRight(unsigned char* renderMat,
			     unsigned int* renderBody,
			     const unsigned nPixels,				
			     const float x, const float y, const float z,
			     const float dy, const float dz,
			     const unsigned ny, const unsigned nz) const = 0;

  virtual void renderXtoUp(unsigned char* renderMat,
			     unsigned int* renderBody,
			     const unsigned nPixels,				
			     const float x, const float y, const float z,
			     const float dy, const float dz,
			     const unsigned ny, const unsigned nz) const = 0;

  virtual void renderXtoDown(unsigned char* renderMat,
			     unsigned int* renderBody,
			     const unsigned nPixels,				
			     const float x, const float y, const float z,
			     const float dy, const float dz,
			     const unsigned ny, const unsigned nz) const = 0;
  
  virtual void testY(std::vector<geoError>& errors,
		     const float x, const float y, const float z,
		     const float dx, const float dz,
		     const unsigned nx, const unsigned nz) const = 0;
  
  virtual void renderY(unsigned char* renderMat,unsigned int* renderBody,
		       const float x, const float y, const float z,
		       const float dx, const float dz,
		       const unsigned nx, const unsigned nz,
		       const unsigned nthreads = 1) const = 0;

  virtual void renderYtoLeft(unsigned char* renderMat,
			     unsigned int* renderBody,
			     const unsigned nPixels,				
			     const float x, const float y, const float z,
			     const float dx, const float dz,
			     const unsigned nx, const unsigned nz) const = 0;

  virtual void renderYtoRight(unsigned char* renderMat,
			      unsigned int* renderBody,
			      const unsigned nPixels,
			      const float x, const float y, const float z,
			      const float dx, const float dz,
			      const unsigned nx, const unsigned nz) const = 0;

  virtual void renderYtoUp(unsigned char* renderMat,
			   unsigned int* renderBody,
			   const unsigned nPixels,				
			   const float x, const float y, const float z,
			   const float dx, const float dz,
			   const unsigned nx, const unsigned nz) const = 0;
  
  virtual void renderYtoDown(unsigned char* renderMat,
			     unsigned int* renderBody,
			     const unsigned nPixels,
			     const float x, const float y, const float z,
			     const float dx, const float dz,
			     const unsigned nx, const unsigned nz) const = 0;
  
  virtual void testZ(std::vector<geoError>& errors,
		     const float x, const float y, const float z,
		     const float dx, const float dy,
		     const unsigned nx, const unsigned ny) const = 0;
  
  virtual void renderZ(unsigned char* renderMat,unsigned int* renderBody,
		       const float x, const float y, const float z,
		       const float dx, const float dy,
		       const unsigned nx, const unsigned ny,
		       const unsigned nthreads = 1) const = 0;

  virtual void renderZtoLeft(unsigned char* renderMat,
			     unsigned int* renderBody,
			     const unsigned nPixels,				
			     const float x, const float y, const float z,
			     const float dx, const float dy,
			     const unsigned nx, const unsigned ny) const = 0;

  virtual void renderZtoRight(unsigned char* renderMat,
			      unsigned int* renderBody,
			      const unsigned nPixels,
			      const float x, const float y, const float z,
			      const float dx, const float dy,
			      const unsigned nx, const unsigned ny) const = 0;

  virtual void renderZtoUp(unsigned char* renderMat,
			   unsigned int* renderBody,
			   const unsigned nPixels,				
			   const float x, const float y, const float z,
			   const float dx, const float dy,
			   const unsigned nx, const unsigned ny) const = 0;
  
  virtual void renderZtoDown(unsigned char* renderMat,
			     unsigned int* renderBody,
			     const unsigned nPixels,
			     const float x, const float y, const float z,
			     const float dx, const float dy,
			     const unsigned nx, const unsigned ny) const = 0;
  
  virtual int render3Dortho(unsigned char* renderMat,unsigned int* renderBody,
			    const float x, const float y, const float z,
			    const float u, const float v, const float w,
			    const float roll, float& phi,
			    float* distances,
			    float& minDistance, float& maxDistance,
			    const float threshold = 1.0e-4) const = 0;
  
  virtual int render3Dortho(unsigned char* renderMat,unsigned int* renderBody,
			    const float x, const float y, const float z,
			    const float u, const float v, const float w,
			    const float roll, float& phi,
			    const float dx, const float dy,
			    const unsigned nx, const unsigned ny,
			    float* distances,
			    float& minDistance, float& maxDistance,
			    const float threshold = 1.0e-4) const = 0;
  
  virtual void set3DResolution(const unsigned nx, const unsigned ny,
			       const float dx, const float dy,
			       const float perspective) = 0;


  virtual int render3D(unsigned char* renderMat,unsigned int* renderBody,
		       const float x, const float y, const float z,
		       const float u, const float v, const float w,
		       const float roll, float& phi,
		       float* distances,
		       float& minDistance, float& maxDistance,
		       const float threshold = 1.0e-4) const = 0;

  virtual ~pen_geoViewInterface(){};
};

#endif
