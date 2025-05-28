
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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


#ifndef __PEN_GEOMETRY_CLASES_
#define __PEN_GEOMETRY_CLASES_

#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <cmath>

#include "pen_classes.hh"

template<class body> class pen_geometry;

namespace penred{
  namespace geometry{
    
    //Check registered types
    bool checkRegisteredObj(const unsigned verbose);
    bool checkRegisteredMesh(const unsigned verbose);
  }
}

//-------------------
// Base body struct
//-------------------

struct pen_baseBody{
  
  unsigned int MATER;
  unsigned int KDET;
  double DSMAX;
  double localEABS[constants::nParTypes];

  pen_baseBody() : MATER(0), KDET(0), DSMAX(1.0e35){
    //Set body energy cutoffs to "infinite" by default
    for(unsigned i = 0; i < constants::nParTypes; i++){
      localEABS[i] = 1.0e-15;
    }
  }
};

//---------------------------
// Base mesh element struct
//---------------------------

struct pen_baseMesh{
  
  unsigned int MATER;
  
  pen_baseMesh() : MATER(0){}
};

template<class body>
class abc_geometry : public wrapper_geometry{

public:
  
  //static const unsigned int NB  = pen_geoconst::NB; //Maximum number of bodies
  static const unsigned int NB = pen_geoconst::NB;
  static const bool checkBody = true;
  
protected:
  //Bodies stack
  unsigned int NBODYS;
  body bodies[NB];

public:
  
  abc_geometry() : wrapper_geometry(),
		   NBODYS(0){
    //Check if used body is a subtype of pen_baseBody
    if(checkBody){
      body testBody;
      pen_baseBody* pbody = dynamic_cast<pen_baseBody*>(&testBody);
      if(pbody == nullptr){
	printf("Warning: abc_geometry: Used bodies are not a subtype of 'pen_baseBody.\n");
	printf("          Used body type can be incompatible with other PenRed classes, use it by your own responsability.\n");
	printf("          To avoid this warning change the value of 'checkBody' flag in abc_geometry class.\n");
      }
    }

  //Check registered types to ensure static library linking of the register variable    
    if(!penred::geometry::checkRegisteredObj(1)){
      printf("Warning: Some geometry object types are not "
	     "properly registered\n");
    }
    
  }

  inline const char* getType() const {return "BODIES";}
  void usedMat(bool used[constants::MAXMAT+1]) const {

    //Set all materials to unused
    for(unsigned i = 0; i <= constants::MAXMAT; i++)
      used[i] = false;
    
    for(unsigned i = 0; i < NBODYS; i++){
      unsigned matIndex = bodies[i].MATER;
      if(matIndex <= constants::MAXMAT)
	if(!used[matIndex])
	  used[matIndex] = true;
    }
  }
  inline double getEabs(const unsigned ibody, const unsigned kpar) const {
    return bodies[ibody].localEABS[kpar];
  }
  
  inline unsigned getMat(const unsigned ibody) const{
    if(ibody >= NBODYS)
      return 0;
    return bodies[ibody].MATER;
  }
  inline double getDSMAX(const unsigned ibody) const{
    if(ibody >= NBODYS)
      return 1.0e35;
    return bodies[ibody].DSMAX;
  }
  inline unsigned getDET(const unsigned ibody) const{
    if(ibody >= NBODYS)
      return 0;
    return bodies[ibody].KDET;
  }  
  
  inline unsigned long getElements() const{
    return NBODYS;
  }

  inline unsigned getBodies() const{
    return NBODYS;
  }
  
  inline const body& readBody(unsigned int kb) const {
    if(kb >= NB)
      {
	char error[300];
	sprintf(error,"getBody: Ordered body (%u) out of range (%d)",kb,NB);
	throw std::out_of_range(error);	
      }
    return bodies[kb];
  }
  int setBodyEabs(unsigned bodyIndex, const double bodyEABS[constants::nParTypes]){
    if(bodyIndex >= NBODYS)
      return 1;
    //Get body reference
    body& bodyRef = bodies[bodyIndex];
    for(unsigned j = 0; j < constants::nParTypes; j++){
      if(bodyEABS[j] > 0.0){
	// Set local absorption energy for specified body
	bodyRef.localEABS[j] = bodyEABS[j];
      }
    }
    return 0;
  }
  
  virtual ~abc_geometry(){}
  
};


template<class meshElement>
class abc_mesh : public wrapper_geometry{

public:
  enum pen_meshStates{
		      PEN_MESH_SUCCESS = 0,
		      PEN_MESH_UNINITIALIZED,
		      PEN_MESH_INITIALIZED,
		      PEN_MESH_ALLOCATION_FAILED
  };
  
private:
  static const bool checkMesh = true;
  pen_meshStates meshStatus;

protected:
  meshElement* mesh;
  unsigned long meshDim;
  unsigned long nElements;

  double DSMAX[constants::MAXMAT+1];  //+1 for the enclosure
  unsigned KDET[constants::MAXMAT+1]; //+1 for the enclosure
  unsigned nBodies;
  
  unsigned enclosureMat;
  
  void clearMesh(){
    if(mesh != nullptr){
      delete[] mesh;
      mesh = nullptr;
    }
    meshDim = 0;
    nElements = 0;
    nBodies = 0;
    meshStatus = PEN_MESH_UNINITIALIZED;
  }

  void resizeMesh(unsigned dim){
    //Clear mesh
    clearMesh();

    //Allocate new mesh
    mesh = new meshElement[dim];
    if(mesh == nullptr)
      meshStatus = PEN_MESH_ALLOCATION_FAILED;
    else{
      meshDim = dim;
      meshStatus = PEN_MESH_INITIALIZED;
    }
  }  
  
public:
  
  abc_mesh() : wrapper_geometry(),
	       meshStatus(PEN_MESH_UNINITIALIZED),
	       mesh(nullptr),
	       meshDim(0),
	       nElements(0),
	       nBodies(0){
    //Check if used body is a subtype of pen_baseBody
    if(checkMesh){
      meshElement testMesh;
      pen_baseMesh* ptest = dynamic_cast<pen_baseMesh*>(&testMesh);
      if(ptest == nullptr){
	printf("Warning: abc_mesh: Used mesh elements are not a subtype of 'pen_baseMesh.\n");
	printf("          Used mesh type can be incompatible with other PenRed classes, use it by your own responsability.\n");
	printf("          To avoid this warning change the value of 'checkMesh' flag in abc_mesh class.\n");
      }
    }

    //Check registered types to ensure static library linking of the register variable    
    if(!penred::geometry::checkRegisteredMesh(1)){
      printf("Warning: Some geometry mesh types are not "
	     "properly registered\n");
    }
    
    for(unsigned i = 0; i <= constants::MAXMAT; i++){
      DSMAX[i] = 1.0e35;
      KDET[i] = 0;
    }
  }

  inline const char* getType() const {return "MESH";}
  
  inline pen_meshStates getStatus() const {return meshStatus;}
  inline unsigned long getElements() const {
    return nElements;
  }

  inline unsigned getBodies() const{
    return nBodies;
  }

  inline unsigned getMat(const unsigned ibody) const{
    if(ibody == 0)
        return enclosureMat;
    if(ibody >= nBodies)
      return 0;
    return ibody;
  }
  
  inline double getDSMAX(const unsigned ibody) const{
    if(ibody >= nBodies)
      return 1.0e35;
    return DSMAX[ibody];
  }  
  inline unsigned getDET(const unsigned ibody) const{
    if(ibody >= nBodies)
      return 0;
    return KDET[ibody];
  }  
  
  void usedMat(bool used[constants::MAXMAT+1]) const {

    //Set all materials to unused
    for(unsigned i = 0; i <= constants::MAXMAT; i++)
      used[i] = false;
    
    for(unsigned i = 0; i < nElements; i++){
      unsigned matIndex = mesh[i].MATER;
      if(matIndex <= constants::MAXMAT)
	if(!used[matIndex])
	  used[matIndex] = true;
    }
    
    //Add enclosure material as used
    used[enclosureMat] = true;
  }
  virtual double getEabs(const unsigned /*index*/, const unsigned /*kpar*/) const{
    return 1.0e-15;
  };

  inline const meshElement& readElement(const unsigned long index) const{
    if(index >= nElements)
      {
	char error[300];
	sprintf(error,"readElement: Ordered element (%lu) out of range (%lu)",index,nElements);
	throw std::out_of_range(error);	
      }
    return mesh[index];    
  }
  
  inline const meshElement* readMesh() const{
      return mesh;
  }
  
  virtual ~abc_mesh(){clearMesh();}
};

#endif
