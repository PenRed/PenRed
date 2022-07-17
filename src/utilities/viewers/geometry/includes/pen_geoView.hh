
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

#ifndef __PEN_GEO_VIEWER__
#define __PEN_GEO_VIEWER__

#include <vector>
#include <thread>
#include "pen_geoViewInterface.hh"
#include "pen_geometries.hh"
#include "pen_matrix.hh"

struct geoError{
  double from[3];
  double to[3];
  unsigned int iIBODY, iMAT; //Initial
  unsigned int eIBODY, eMAT; //Expected
  unsigned int fIBODY, fMAT; //Final
  bool endInVoid;
};

class pen_geoView : public pen_geoViewInterface{

private:
  wrapper_geometry* geometry;

  float perspective3D;
  unsigned nx3D, ny3D;
  unsigned long nxy3D;
  float dx3D,dy3D;
  float ox3D, oy3D;

  std::vector<float> perspectiveVect;
  
  bool checkStep(pen_particleState& state,
		 double toTravel,
		 geoError& err) const;
  
public:

  pen_geoView() : geometry(nullptr){
    set3DResolution(512,512,0.1,0.1,0.3490658503988659);
  }

  inline double z2dir(const double u,
		      const double v,
		      const double w,
		      const double omega,
		      double rotation[9],
		      const double phiAux,
		      const double threshold) const {
    return rollAlign(u,v,w,omega,rotation,phiAux,threshold);
  }
  
  inline float z2dirf(const float u,
		      const float v,
		      const float w,
		      const float omega,
		      float rotation[9],
		      const float phiAux,
		      const float threshold) const {
    return rollAlignf(u,v,w,omega,rotation,phiAux,threshold);
  }
  
  inline int init(const char* filename,
		  const unsigned verbose = 5){

    //Parse configuration file
    pen_parserSection config;
    std::string errorLine;
    unsigned long errorLineNum;
    int err = parseFile(filename,config,errorLine,errorLineNum);
  
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_geoView: init: Error parsing configuration.\n");
	printf("pen_geoView: init: Error code: %d\n",err);
	printf("pen_geoView: init: Error message: %s\n",pen_parserError(err));
	printf("pen_geoView: init: Error located at line %lu, at text: %s\n",
	       errorLineNum,errorLine.c_str());
	fflush(stdout);
      }
      return -1;
    }

    return init(config,verbose);
  }

  int init(const pen_parserSection& configIn,
	   const unsigned verbose = 5);

  void testX(std::vector<geoError>& errors,
	     const float x, const float y, const float z,
	     const float dy, const float dz,
	     const unsigned ny, const unsigned nz) const;
  
  void renderX(unsigned char* renderMat,unsigned int* renderBody,
	       const float x, const float y, const float z,
	       const float dy, const float dz,
	       const unsigned ny, const unsigned nz,
	       const unsigned nthreads = 1) const;

  void renderXtoLeft(unsigned char* renderMat,
		     unsigned int* renderBody,
		     const unsigned nPixels,				
		     const float x, const float y, const float z,
		     const float dy, const float dz,
		     const unsigned ny, const unsigned nz) const;

  void renderXtoRight(unsigned char* renderMat,
		      unsigned int* renderBody,
		      const unsigned nPixels,
		      const float x, const float y, const float z,
		      const float dy, const float dz,
		      const unsigned ny, const unsigned nz) const;

  void renderXtoUp(unsigned char* renderMat,
		     unsigned int* renderBody,
		     const unsigned nPixels,				
		     const float x, const float y, const float z,
		     const float dy, const float dz,
		     const unsigned ny, const unsigned nz) const;

  void renderXtoDown(unsigned char* renderMat,
		      unsigned int* renderBody,
		      const unsigned nPixels,
		      const float x, const float y, const float z,
		      const float dy, const float dz,
		      const unsigned ny, const unsigned nz) const;
  
  void testY(std::vector<geoError>& errors,
	     const float x, const float y, const float z,
	     const float dx, const float dz,
	     const unsigned nx, const unsigned nz) const;
  
  void renderY(unsigned char* renderMat,unsigned int* renderBody,
	       const float x, const float y, const float z,
	       const float dx, const float dz,
	       const unsigned nx, const unsigned nz,
	       const unsigned nthreads = 1) const;

  void renderYtoLeft(unsigned char* renderMat,
		     unsigned int* renderBody,
		     const unsigned nPixels,				
		     const float x, const float y, const float z,
		     const float dx, const float dz,
		     const unsigned nx, const unsigned nz) const;

  void renderYtoRight(unsigned char* renderMat,
		      unsigned int* renderBody,
		      const unsigned nPixels,
		      const float x, const float y, const float z,
		      const float dx, const float dz,
		      const unsigned nx, const unsigned nz) const;

  void renderYtoUp(unsigned char* renderMat,
		     unsigned int* renderBody,
		     const unsigned nPixels,				
		     const float x, const float y, const float z,
		     const float dx, const float dz,
		     const unsigned nx, const unsigned nz) const;

  void renderYtoDown(unsigned char* renderMat,
		      unsigned int* renderBody,
		      const unsigned nPixels,
		      const float x, const float y, const float z,
		      const float dx, const float dz,
		      const unsigned nx, const unsigned nz) const;  
  
  void testZ(std::vector<geoError>& errors,
	     const float x, const float y, const float z,
	     const float dx, const float dy,
	     const unsigned nx, const unsigned ny) const;
  
  void renderZ(unsigned char* renderMat,unsigned int* renderBody,
	       const float x, const float y, const float z,
	       const float dx, const float dy,
	       const unsigned nx, const unsigned ny,
	       const unsigned nthreads = 1) const;

  void renderZtoLeft(unsigned char* renderMat,
		     unsigned int* renderBody,
		     const unsigned nPixels,				
		     const float x, const float y, const float z,
		     const float dx, const float dy,
		     const unsigned nx, const unsigned ny) const;

  void renderZtoRight(unsigned char* renderMat,
		      unsigned int* renderBody,
		      const unsigned nPixels,
		      const float x, const float y, const float z,
		      const float dx, const float dy,
		      const unsigned nx, const unsigned ny) const;

  void renderZtoUp(unsigned char* renderMat,
		     unsigned int* renderBody,
		     const unsigned nPixels,				
		     const float x, const float y, const float z,
		     const float dx, const float dy,
		     const unsigned nx, const unsigned ny) const;

  void renderZtoDown(unsigned char* renderMat,
		      unsigned int* renderBody,
		      const unsigned nPixels,
		      const float x, const float y, const float z,
		      const float dx, const float dy,
		      const unsigned nx, const unsigned ny) const;  
  
  inline int render3Dortho(unsigned char* renderMat,unsigned int* renderBody,
			   const float x, const float y, const float z,
			   const float u, const float v, const float w,
			   const float roll, float& phi,
			   float* distances,
			   float& minDistance, float& maxDistance,
			   const float threshold = 1.0e-4) const{
    return render3Dortho(renderMat,renderBody,
			 x,y,z,u,v,w,roll,phi,dx3D,dy3D,nx3D,ny3D,
			 distances,minDistance,maxDistance,threshold);
  }
  
  int render3Dortho(unsigned char* renderMat,unsigned int* renderBody,
		    const float x, const float y, const float z,
		    const float u, const float v, const float w,
		    const float roll, float& phi,
		    const float dx, const float dy,
		    const unsigned nx, const unsigned ny,
		    float* distances,
		    float& minDistance, float& maxDistance,
		    const float threshold = 1.0e-4) const;
  
  void set3DResolution(const unsigned nx, const unsigned ny,
		       const float dx, const float dy,
		       const float perspective);


  int render3D(unsigned char* renderMat,unsigned int* renderBody,
	       const float x, const float y, const float z,
	       const float u, const float v, const float w,
	       const float roll, float& phi, float* distances,
	       float& minDistance, float& maxDistance,
	       const float threshold = 1.0e-4) const;
  
  ~pen_geoView(){
    if(geometry != nullptr)
      delete geometry;
  }
};

#endif
