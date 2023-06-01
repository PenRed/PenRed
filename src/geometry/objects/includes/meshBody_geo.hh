// 
//
//
//    Copyright (C) 2022-2023 Universitat de València - UV
//    Copyright (C) 2022-2023 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es    
//        vicente.gimenez@uv.es
//

 
#ifndef __PEN_MESH_BODY_GEOMETRY__
#define __PEN_MESH_BODY_GEOMETRY__

#include "geometry_classes.hh"
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <atomic>

enum pen_meshBodyErr{
    PEN_MESHBODY_GEO_SUCCESS = 0,
    PEN_MESHBODY_GEO_INPUT_SECTION,
    PEN_MESHBODY_GEO_BAD_READ_DSMAX,
    PEN_MESHBODY_GEO_INVALID_DSMAX,
    PEN_MESHBODY_GEO_BAD_READ_KDET,
    PEN_MESHBODY_GEO_INVALID_KDET,
    PEN_MESHBODY_GEO_BAD_READ_REGIONSIZE,
    PEN_MESHBODY_GEO_INVALID_REGIONSIZE,
    PEN_MESHBODY_GEO_UNKNOWN_PARTICLE,
    PEN_MESHBODY_GEO_BAD_READ_EABS,
    PEN_MESHBODY_GEO_INVALID_EABS,
    PEN_MESHBODY_GEO_UNDEF_BODY_LABEL,
    PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT,
    PEN_MESHBODY_GEO_INVALID_NBODIES,
    PEN_MESHBODY_GEO_INVALID_MAT,
    PEN_MESHBODY_BAD_MEMORY_ALLOCATION,
    PEN_MESHBODY_MULTIPLE_WORLDS,
    PEN_MESHBODY_WORLD_NOT_FOUND,
    PEN_MESHBODY_GEO_INVALID_VERTEX_INDEX,
    PEN_MESHBODY_GEO_TRIANGLES_OUT_OF_REGIONS,
    PEN_MESHBODY_GEO_LOST_TRIANGLES,
    PEN_MESHBODY_GEO_BODY_INTERSECTIONS_FOUND,
    PEN_MESHBODY_GEO_CL_VECTOR_SIZES_MISMATCH,
    PEN_MESHBODY_GEO_CL_DEVICE_NOT_FOUND,
    PEN_MESHBODY_GEO_CL_MALLOC_FAIL,
    PEN_MESHBODY_GEO_CL_BUFFER_CREATION_FAIL,
    PEN_MESHBODY_GEO_CL_BUFFER_WRITE_FAIL,
    PEN_MESHBODY_GEO_CL_BUFFER_PACK_FAIL,
    PEN_MESHBODY_GEO_CL_PROGRAM_CREATION_FAIL,
    PEN_MESHBODY_GEO_CL_PROGRAM_BUILD_FAIL,
    PEN_MESHBODY_GEO_CL_KERNEL_CREATION_FAIL,
    PEN_MESHBODY_GEO_CL_KERNEL_PACK_FAIL,
    PEN_MESHBODY_GEO_CL_NO_CONFIGURED_DEVICES,    
};

class pen_meshBodyGeo;

struct meshBodyTriangle : public triangle<double>{

  static constexpr double crossThreshold = 1.0e-8;
  
  typedef vector3D<double> v3D;

  v3D edge1;
  v3D edge2;
  
  v3D c;
  double r;
  double r2;
    
  void fill(const v3D v1in, const v3D v2in, const v3D v3in);
    
  void refresh();
    
  inline bool canCross(const v3D pos,
		       const v3D dir,
		       const double maxDS) const {
    //Function to calculate de distance between the center of the
    //triangle and the ray. If this distance is shorter than the sphere
    //radious, the ray cross the sphere and can also cross the triangle.
    //Check Line-Sphere intersection algorithms for normalized direction vector

    // oc = pos - c
    // proj = dir * oc
    // L = proj**2 - ( ||oc||**2 - r**2 )
    // d2aux = ( ||oc||**2 - r**2 )
    // d = -proj +- sqrt(L)

    //If the point is inside the sphere, is not required further checks, i.e.
    //
    //  if(d2aux < 0) return true
    //
    //If the previous condition is not satisfied, then L
    //must be lesser than the square of the projection
    //   
    //  proj*proj > L  -> proj > sqrt(L)
    //
    //Therefore, if the projection is positive, the two points where the
    //sphere is cutted are behind the origin point. Thus,
    //
    //  if(proj > 0) return false
    //
    //Finally, check the sign of L and the distance
    
    const v3D oc = pos - c;
    double d2aux = oc.mod2();
    d2aux -= r2;
    
    if(d2aux < crossThreshold){
      return true;
    }
    
    const double proj = dir*oc;

    if(!std::signbit(proj))
      return false;
    
    const double sqrtArg = proj*proj - d2aux;

    if(sqrtArg > -crossThreshold ){
      double ds2Bound = -proj;
      if(sqrtArg > crossThreshold){
	ds2Bound -= sqrt(sqrtArg);
      }

      if(ds2Bound > maxDS)
	return false;
      
      return true;
    
    }else{
      return false;
    }
  }
    
  bool intersect(const v3D pos,
		 const v3D dir,
		 double& t,
		 const bool back) const;

  inline v3D readV1() const {return v1;}
  inline v3D readV2() const {return v2;}
  inline v3D readV3() const {return v3;}

  inline v3D minv() const {

    unsigned zi = minzi();
    if(zi == 1)
      return v1;
    if(zi == 2)
      return v2;
    if(zi == 3)
      return v3;

   unsigned yi = minyi();
    if(yi == 1)
      return v1;
    if(yi == 2)
      return v2;
    if(yi == 3)
      return v3;    
    
   unsigned xi = minxi();
    if(xi == 1)
      return v1;
    if(xi == 2)
      return v2;
    if(xi == 3)
      return v3;    

    return v1;
  }
};

struct pen_meshBody : public pen_baseBody{

  typedef container<meshBodyTriangle,double> triangleRegion;
  typedef container<triangleRegion,double> superRegion;
  
  static constexpr double crossThreshold = meshBodyTriangle::crossThreshold;
  
  typedef vector3D<double> v3D;
    
  char PALIAS[100]; //Parent alias name
  char BALIAS[100];

  //Mean number of triangles per region
  unsigned long nTriangles;
  unsigned long meanTrianglesRegion;
  unsigned long meanRegionsSuperRegion;

  //Body box boundary
  box<double> boundary;

  //Body regions
  static const unsigned MAX_SUP_REGIONS = 500;
  static const unsigned MAX_REGIONS = 10000;
  std::vector<superRegion> regions;

  //Sister bodies with overlap
  std::array<unsigned,300> overlapedBodies;
  //Number of syster overlaps
  unsigned nOverlap;

  
  //Flags if the body can overlap its parent
  bool canOverlapParent;

  //Daughters bodies array
  std::array<unsigned,300> daughters;
  //Number of daughter bodies
  unsigned nDaughters;

  //Parent body index
  unsigned parent;
  
  pen_meshBody() : nTriangles(0), meanTrianglesRegion(0),
		   nOverlap(0), nDaughters(0) {}
    
  bool inside(const v3D pos) const;
    
  bool cross(const v3D pos, 
	     const v3D dir, 
	     double& ds,
	     const bool back,
	     const double maxDs = 1.0e35) const;

  inline size_t nSupRegions() const {return regions.size();}
  
  inline void addDaughter(const unsigned i){
    if(nDaughters < 300){
      daughters[nDaughters++] = i;
    }else{
      printf("meshBody_geo: Error: maximum number of daughters is 300\n");
    }
  }
    
  inline void printDaughters(const int indent, const pen_meshBody* bodies) const{
        
    printf("%*s|-%s\n",indent,"",BALIAS);
    for(size_t i = 0; i < nDaughters; ++i){
      bodies[daughters[i]].printDaughters(indent+2,bodies);
    }
  }

  inline double xmin() const{return boundary.minx();}
  inline double ymin() const{return boundary.miny();}
  inline double zmin() const{return boundary.minz();}

  inline double xmax() const{return boundary.maxx();}
  inline double ymax() const{return boundary.maxy();}
  inline double zmax() const{return boundary.maxz();}

};

class pen_meshBodyGeo : public abc_geometry<pen_meshBody>{
  DECLARE_GEOMETRY(pen_meshBodyGeo)

  private:

  unsigned iworld;
  bool worldFound;

public:
      
  static constexpr double threshold = meshBodyTriangle::crossThreshold;
      
  typedef vector3D<double> v3D;    
      
  pen_meshBodyGeo() : iworld(0), worldFound(false) {
    configStatus = 0;
  }
  
  int configure(const pen_parserSection& config, const unsigned verbose) override;
  int GEOMESH(FILE* in, const unsigned verbose);
  
  void locate(pen_particleState&) const final override;
  void step(pen_particleState&,
	    double,
	    double &,
	    double &,
	    int &) const final override;
        
  bool canOverlapParent(const unsigned) const ;

  bool canOverlap(const unsigned, const unsigned) const;
  
  void checkCross(const unsigned iparent);

  inline unsigned getIBody(const char* alias) const override{
      
    for(unsigned i = 0; i < getBodies(); ++i){
      if(strcmp(bodies[i].BALIAS,alias) == 0)
	return i;
    }
    return getBodies();
  }
    
  inline void move(const double ds, pen_particleState& state) const{
    state.X += ds*state.U;
    state.Y += ds*state.V;
    state.Z += ds*state.W;
  }

  inline void move(const double ds, const v3D dir, v3D& pos) const{
    pos.x += ds*dir.x;
    pos.y += ds*dir.y;
    pos.z += ds*dir.z;
  }  
  
  
  
  inline bool solveOverlapsFlat(const double travel,
				const v3D& pos,
				const v3D& dir,
				const unsigned ibody, 
				unsigned& nextBody) const {

    // travel   -> Traveled distance in cm
    // pos      -> Position vector (x,y,z)
    // dir      -> Normalized direction (u,v,w)
    // ibody    -> Actual body index
    // nextBody -> Body index going to. Will be updated if an overlap
    //             exists.
            

    //Get body reference
    const pen_meshBody& body = bodies[ibody];
    
    //Check all overlaps at the same level
    for(unsigned iover = 0; iover < body.nOverlap; ++iover){
        
      //Get body index
      const unsigned overIndex = body.overlapedBodies[iover];
      //Get the reference of the possible overlaping body
      const pen_meshBody& overBody = bodies[overIndex];
        
      //Check if this body is crossed in this direction
      double dsOverlap;
      if(overBody.cross(pos, dir, dsOverlap, false)){
            
	//Is crossed, check if is crossed before the travel finish
	if(dsOverlap - travel < threshold){
	  //The cross is before the travel finish. Thus, the actual
	  //next body is the crossed one. However, the cross with its
	  //daughters must be checked
	  nextBody = overIndex;
	  solveOverlapsDown(travel,pos,dir,overIndex,nextBody);
	  return true;
	}
      }
    }
    
    //No overlaps
    return false;
  }
  
  inline bool solveOverlapsUp(const double travel,
			      const v3D& pos,
			      const v3D& dir,
			      const unsigned ibody, 
			      unsigned& nextBody) const {
    
    // travel   -> Traveled distance in cm
    // pos      -> Position vector (x,y,z)
    // dir      -> Normalized direction (u,v,w)
    // ibody    -> Actual body index
    // nextBody -> Body index going to. Will be updated if an overlap
    //             exists. In this function, next body should be the
    //             parent of ibody
                                
    
    //Get body reference
    const pen_meshBody& body = bodies[ibody];
    
    //Check if can overlap with the parent
    if(body.canOverlapParent){
      //Can overlap with the parent, get parent reference
      const pen_meshBody& parent = bodies[body.parent];
        
      //Check if the parent is crossed in this direction
      double dsOverlap;
      if(parent.cross(pos, dir, dsOverlap, true)){
            
	//The parent is crossed, check if the parent
	//cross is located before the travel finish
	if(dsOverlap - travel < threshold){
                
	  //ibody is overlaping with its parent within the travel 
	    //distance. Therefore, the parent is crossed and the
	    //particle escapes to the parent of the ibody parent. 
	    //So, update next body to grandparent
	    nextBody = bodies[body.parent].parent;
                
	  //However, overlaps of the parent must be also checked
	  solveOverlapsUp(travel,pos,dir,body.parent,nextBody);
	  return true;
	}
      }
    }
    
    //The parent is not overlpaed, check their sisters
    return solveOverlapsFlat(travel,pos,dir,ibody,nextBody);
  }

  inline bool solveOverlapsDown(const double travel,
				const v3D& pos,
				const v3D& dir,
				const unsigned ibody, 
				unsigned& nextBody) const {
    
    // travel   -> Traveled distance in cm
    // pos      -> Position vector (x,y,z)
    // dir      -> Normalized direction (u,v,w)
    // ibody    -> Actual body index
    // nextBody -> Body index going to. Will be updated if an overlap
    //             exists.                                  
      
    //Get body reference
    const pen_meshBody& body = bodies[ibody];
    
    //Iterate over all daughters to check possible overlaps
    for(unsigned idaught = 0; idaught < body.nDaughters; ++idaught){
        
      //Get daugther index and reference
      const unsigned daugthIndex = body.daughters[idaught];
      const pen_meshBody& daugth = bodies[daugthIndex];
        
      //Check if this daugther can overlap with ibody
      if(daugth.canOverlapParent){
	
	//Can overlap, check if is crossed in this direction
	double dsOverlap;
	if(daugth.cross(pos, dir, dsOverlap, false)){
	  
	  //Is crossed, check if the cross is before the travel end
	  if(dsOverlap - travel < threshold){
	    //Daugther is crossed, thus, the next body must be
	    //updated to this daugher. However, Their daughters
	    //must be checked also
	    
	    nextBody = daugthIndex;
	    solveOverlapsDown(travel,pos,dir,daugthIndex,nextBody);
	    return true;
	  }
	}
      }
    }
    
    //No overlapes
    return false;
  }
  
};

#endif
