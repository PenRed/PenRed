// 
//
//
//    Copyright (C) 2022-2023 Universitat de València - UV
//    Copyright (C) 2022-2023 Universitat Politècnica de València - UPV
//    Copyright (C) 2025-2026 Vicent Giménez Alventosa
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
#include <map>
#include <memory>
#include <fstream>

class pen_meshBodyGeo;

struct meshBodyTriangle : public triangle<double>{

  static constexpr const double crossThreshold = 1.0e-8;
  
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

namespace pen_meshTransform{
  
  struct base{

    typedef vector3D<double> v3D;

    v3D cm;

    constexpr base() : cm(0.0,0.0,0.0) {}
    
    virtual void apply(v3D& v) const = 0;
    virtual std::string stringify() const = 0;

    inline virtual ~base(){}
  };

  struct trans : base{

    v3D dir;
    double ds;

    constexpr trans(const v3D& dirIn, const double dsIn) : dir(dirIn), ds(dsIn) {}
    
    inline void apply(v3D& v) const override {
      v += dir*ds;
    }
    inline std::string stringify() const override {
      return std::string("Translation with direction ") +
	dir.stringify() +
	std::string(" and distance ") + std::to_string(ds);
    }
  
  };

  struct transX : base{

    double ds;

    constexpr transX(const double dsIn) : ds(dsIn) {}
  
    inline void apply(v3D& v) const override {
      v.x += ds;
    }
    inline std::string stringify() const override {
      return std::string("Translation on X axis, a distance ") +
	std::to_string(ds);
    }    
  
  };

  struct transY : base{

    double ds;

    constexpr transY(const double dsIn) : ds(dsIn) {}
    
    inline void apply(v3D& v) const override {
      v.y += ds;
    }
    inline std::string stringify() const override {
      return std::string("Translation on Y axis, a distance ") +
	std::to_string(ds);
    }    
  
  };

  struct transZ : base{

    double ds;

    constexpr transZ(const double dsIn) : ds(dsIn) {}    
  
    inline void apply(v3D& v) const override {
      v.z += ds;
    }
    inline std::string stringify() const override {
      return std::string("Translation on Z axis, a distance ") +
	std::to_string(ds);
    }    
  
  };

  struct scale : base{

    double factor;

    scale(const double factorIn) : factor(factorIn) {}
    
    inline void apply(v3D& v) const override {
      v = (v-cm)*factor + cm;
    }

    inline std::string stringify() const override {
      return std::string("Scale with factor ") +
	std::to_string(factor);
    }    
  
  };

  struct scaleX : base{

    double factor;

    scaleX(const double factorIn) : factor(factorIn) {}    
  
    inline void apply(v3D& v) const override {
      v.x = (v.x-cm.x)*factor + cm.x;
    }
    inline std::string stringify() const override {
      return std::string("Scale on X axis, with factor ") +
	std::to_string(factor);
    }
  
  };

  struct scaleY : base{

    double factor;

    scaleY(const double factorIn) : factor(factorIn) {}    
  
    inline void apply(v3D& v) const override {
      v.y = (v.y-cm.y)*factor + cm.y;
    }
    inline std::string stringify() const override {
      return std::string("Scale on Y axis, with factor ") +
	std::to_string(factor);
    }
  
  };

  struct scaleZ : base{

    double factor;

    scaleZ(const double factorIn) : factor(factorIn) {}    
  
    inline void apply(v3D& v) const override {
      v.z = (v.z-cm.z)*factor + cm.z;
    }
    inline std::string stringify() const override {
      return std::string("Scale on Z axis, with factor ") +
	std::to_string(factor);
    }
  
  };  
  
  struct scaleXY : base{

    double factor;

    scaleXY(const double factorIn) : factor(factorIn) {}    
  
    inline void apply(v3D& v) const override {
      v.x = (v.x-cm.x)*factor + cm.x;
      v.y = (v.y-cm.y)*factor + cm.y;
    }
    inline std::string stringify() const override {
      return std::string("Scale on XY plane, with factor ") +
	std::to_string(factor);
    }
  
  };

  struct scaleXZ : base{

    double factor;
    
    scaleXZ(const double factorIn) : factor(factorIn) {}
    
    inline void apply(v3D& v) const override {
      v.x = (v.x-cm.x)*factor + cm.x;
      v.z = (v.z-cm.z)*factor + cm.z;
    }
    inline std::string stringify() const override {
      return std::string("Scale on XZ plane, with factor ") +
	std::to_string(factor);
    }    
    
  
  };

  struct scaleYZ : base{

    double factor;

    scaleYZ(const double factorIn) : factor(factorIn) {}
    
    inline void apply(v3D& v) const override{
      v.y = (v.y-cm.y)*factor + cm.y;
      v.z = (v.z-cm.z)*factor + cm.z;
    }
    inline std::string stringify() const override {
      return std::string("Scale on YZ plane, with factor ") +
	std::to_string(factor);
    }    
    
  
  };


  struct group{

    typedef vector3D<double> v3D;

  private:
    std::vector<std::unique_ptr<base>> transforms;
  public:
  
    //Stores the transformations to be applied to a vertex group
    std::string name;

    inline void resize(const size_t newSize){ transforms.resize(newSize); }
    inline size_t size() const { return transforms.size(); }

    std::unique_ptr<base>&       operator[](std::size_t i)       { return transforms[i]; }
    const std::unique_ptr<base>& operator[](std::size_t i) const { return transforms[i]; }

    inline void apply(const std::vector<unsigned>& vgroup, std::vector<v3D>& vertex) const{

      //Calculate CM of the vertex group
      v3D cm(0.0,0.0,0.0);
      for(size_t iv = 0; iv < vgroup.size(); ++iv){
	cm += vertex[vgroup[iv]];
      }
      cm /= static_cast<double>(vgroup.size());

      //Apply the transformations
      for(size_t it = 0; it < transforms.size(); ++it){

	//Skip null transformations
	if(transforms[it] == nullptr)
	  continue;
      
	//Set center of mass
	transforms[it]->cm = cm;
	for(size_t iv = 0; iv < vgroup.size(); ++iv){
	  transforms[it]->apply(vertex[vgroup[iv]]);
	}
      }
    }

    //Add transform functions
    inline void addTranslation(const v3D& dir, const double ds){
      transforms.emplace_back(std::unique_ptr<base>(new trans(dir,ds)));
    }
    inline void addTranslationX(const double ds){
      transforms.emplace_back(std::unique_ptr<base>(new transX(ds)));
    }
    inline void addTranslationY(const double ds){
      transforms.emplace_back(std::unique_ptr<base>(new transY(ds)));
    }
    inline void addTranslationZ(const double ds){
      transforms.emplace_back(std::unique_ptr<base>(new transZ(ds)));
    }
    inline void addScale(const double f){
      transforms.emplace_back(std::unique_ptr<base>(new scale(f)));
    }
    inline void addScaleXY(const double f){
      transforms.emplace_back(std::unique_ptr<base>(new scaleXY(f)));
    }
    inline void addScaleXZ(const double f){
      transforms.emplace_back(std::unique_ptr<base>(new scaleXZ(f)));
    }
    inline void addScaleYZ(const double f){
      transforms.emplace_back(std::unique_ptr<base>(new scaleYZ(f)));
    }

    //Set transform functions
    inline int setTranslation(const size_t i, v3D dir, const double ds){
      if(i >= transforms.size() ) return 1;
      dir.normalize();
      transforms[i] =
	std::unique_ptr<base>(new trans(dir,ds));
      return 0;
    }
    inline int setTranslationX(const size_t i, const double ds){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new transX(ds));
      return 0;
    }
    inline int setTranslationY(const size_t i, const double ds){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new transY(ds));
      return 0;
    }
    inline int setTranslationZ(const size_t i, const double ds){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new transZ(ds));
      return 0;
    }
    inline int setScale(const size_t i, const double f){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new scale(f));
      return 0;
    }
    inline int setScaleX(const size_t i, const double f){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new scaleX(f));
      return 0;
    }
    inline int setScaleY(const size_t i, const double f){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new scaleY(f));
      return 0;
    }
    inline int setScaleZ(const size_t i, const double f){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new scaleZ(f));
      return 0;
    }    
    inline int setScaleXY(const size_t i, const double f){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new scaleXY(f));
      return 0;
    }
    inline int setScaleXZ(const size_t i, const double f){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new scaleXZ(f));
      return 0;
    }
    inline int setScaleYZ(const size_t i, const double f){
      if(i >= transforms.size() ) return 1;
      transforms[i] =
	std::unique_ptr<base>(new scaleYZ(f));
      return 0;
    }

    inline std::string stringify() const{
      std::string ret;
      for(auto&& b : transforms){
	ret += b->stringify() + "\n";
      }
      return ret;
    }
  
  };
  
} // namespace pen_meshTransform

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
  std::vector<unsigned> overlapedBodies;
  
  //Flags if the body can overlap its parent
  bool canOverlapParent;

  //Daughters bodies vector
  std::vector<unsigned> daughters;

  //Parent body index
  unsigned parent;
  
  pen_meshBody() : nTriangles(0), meanTrianglesRegion(0),
		   meanRegionsSuperRegion(0),
		   canOverlapParent(false) {}
    
  bool inside(const v3D pos) const;
    
  bool cross(const v3D pos, 
	     const v3D dir, 
	     double& ds,
	     const bool back,
	     const double maxDs = 1.0e35) const;

  inline size_t nSupRegions() const {return regions.size();}
  
  inline void addDaughter(const unsigned i){
    daughters.push_back(i);
  }
  
  inline void printDaughters(const int indent, const pen_meshBody* bodies) const{
        
    printf("%*s|-%s\n",indent,"",BALIAS);
    for(size_t i = 0; i < daughters.size(); ++i){
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

  enum errors{
    SUCCESS = 0,
    SECTION_READ_FAIL,
    LIMIT_REACHED,
    MISSING_PARAMETER,
    BAD_VALUE,
    LOW_ON_MEMORY,
    UNKNOWN_PARTICLE,
    UNKNOWN_BODY_LABEL,
    UNEXPECTED_LINE_FORMAT,
    INVALID_FILE,
    MULTIPLE_WORLDS,
    WORLD_NOT_FOUND,
    INVALID_TYPE,
    VG_NOT_FOUND,
    LOST_TRIANGLES,
    BODY_INTERSECTIONS_FOUND,
    VECTOR_SIZES_MISMATCH,
    CL_DEVICE_NOT_FOUND,
    CL_MALLOC_FAIL,
    CL_BUFFER_CREATION_FAIL,
    CL_BUFFER_WRITE_FAIL,
    CL_BUFFER_PACK_FAIL,
    CL_PROGRAM_CREATION_FAIL,
    CL_PROGRAM_BUILD_FAIL,
    CL_KERNEL_CREATION_FAIL,
    CL_KERNEL_PACK_FAIL,
    CL_NO_CONFIGURED_DEVICES,    
    UNEXPECTED_ERROR,
  };
  
  static constexpr const char* errorMessage(const int val) noexcept {
    switch(val){
    case SUCCESS: return "Success";
    case SECTION_READ_FAIL: return "Unable to read section";
    case LIMIT_REACHED: return "Limit reached";
    case MISSING_PARAMETER: return "Missing parameter";
    case BAD_VALUE: return "Invalid value";
    case LOW_ON_MEMORY: return "Low on memory";
    case UNKNOWN_PARTICLE: return "Unknown particle";
    case UNKNOWN_BODY_LABEL: return "Unknown body label";
    case UNEXPECTED_LINE_FORMAT: return "Unexpected line format";
    case INVALID_FILE: return "Invalid file";
    case MULTIPLE_WORLDS: return "Multiple worlds defined";
    case WORLD_NOT_FOUND: return "No world found";
    case INVALID_TYPE: return "Invalid type";
    case VG_NOT_FOUND: return "Vertex group not found";
    case LOST_TRIANGLES: return "Lost triangles";
    case BODY_INTERSECTIONS_FOUND: return "Intersecting bodies found";
    case VECTOR_SIZES_MISMATCH: return "Vector sizes mismatch";
    case CL_DEVICE_NOT_FOUND: return "Opencl device not found";
    case CL_MALLOC_FAIL: return "Opencl malloc failed";
    case CL_BUFFER_CREATION_FAIL: return "Opencl buffer creation failed";
    case CL_BUFFER_WRITE_FAIL: return "Opencl buffer write failed";
    case CL_BUFFER_PACK_FAIL: return "Opencl buffer pack failed";
    case CL_PROGRAM_CREATION_FAIL: return "Opencl program creation failed";
    case CL_PROGRAM_BUILD_FAIL: return "Opencl program build failed";
    case CL_KERNEL_CREATION_FAIL: return "Opencl kernel creation failed";
    case CL_KERNEL_PACK_FAIL: return "Opencl kernel pack failed";
    case CL_NO_CONFIGURED_DEVICES: return "No opencl devices configured";
    case UNEXPECTED_ERROR: return "Unexpected error";
    default: return "Unknown error";
    }
  }  

  //A preload geometry file to use instead of filename during configuration
  std::string preloadGeo;
      
  static constexpr double threshold = meshBodyTriangle::crossThreshold;
      
  typedef vector3D<double> v3D;    
      
  pen_meshBodyGeo() : iworld(0), worldFound(false) {}  
  
  penred::errors::Error specificConfigure(const pen_parserSection& config, const unsigned verbose) override;
  penred::errors::Error GEOMESH(std::istream& in,
				std::map<std::string, std::vector<pen_meshTransform::group>>& transMap,
				const unsigned verbose);
  static int meshGetLine(std::vector<std::ifstream>& included,
			 std::istream& root,
			 std::string&line,
			 unsigned long& nRead);
  
  
  void locateLocal(pen_particleState&) const final override;
  void stepLocal(pen_particleState&,
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

  inline std::string getBodyName(const unsigned ibody) const override{

    if(ibody < getBodies()){
      return std::string(bodies[ibody].BALIAS);
    }else{
      return std::string("NONE");
    }
    
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
				v3D& pos,
				v3D& dir,
				const double t,
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

    if(body.overlapedBodies.size() == 0)
      return false;
    
    //Check all overlaps at the same level
    for(unsigned iover = 0; iover < body.overlapedBodies.size(); ++iover){
        
      //Get body index
      const unsigned overIndex = body.overlapedBodies[iover];
      //Get the reference of the possible overlaping body
      const pen_meshBody& overBody = bodies[overIndex];

      //Apply the overBody inverse transform if needed.
      //This will convert the position and direction from parent local
      //coordinates to overBody local coordinates
      v3D oPos = pos;
      v3D oDir = dir;
      if(overBody.inAnimation(t)){
	overBody.readAnimation().applyInv(t, oPos, oDir);
      }
      
      //Check if this body is crossed in this direction
      double dsOverlap;
      if(overBody.cross(oPos, oDir, dsOverlap, false)){
            
	//Is crossed, check if is crossed before the travel finish
	if(dsOverlap - travel < threshold){
	  //The cross is before the travel finish. Thus, the actual
	  //next body is the crossed one. However, the cross with its
	  //daughters must be checked
	  nextBody = overIndex;
	  solveOverlapsDown(travel,oPos,oDir,t,overIndex,nextBody);

	  //Save position and direction in final body's local coordinates
	  pos = oPos;
	  dir = oDir;
	  return true;
	}
      }
    }
    
    //No overlaps
    return false;
  }
  
  inline bool solveOverlapsUp(const double travel,
			      v3D& pos,
			      v3D& dir,
			      const double t,
			      const unsigned ibody, 
			      unsigned& nextBody) const {

    // Solves the overlaps through the "ibody" boundaries, checking parent and sisters.
    // The position and direction vector are suposed to be in the "ibody" parent's non
    // transformed coordinates
    //
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
	
      //Check if the parent is crossed in this direction.
      double dsOverlap;
      if(parent.cross(pos, dir, dsOverlap, true)){
      
	//The parent is crossed, check if the parent
	//cross is located before the travel finish
	if(dsOverlap - travel < threshold){
                
	  //ibody is overlaping with its parent within the travel 
	  //distance. Therefore, the parent is crossed and the
	  //particle escapes to the ibody parent's parent. 
	  //So, update next body to grandparent
	  nextBody = parent.parent;
                
	  //However, overlaps of the parent must be also checked.
	  //Convert position and direction to grandparent non transformed
	  //coordinates
	  if(parent.inAnimation(t)){
	    parent.readAnimation().apply(t, pos, dir);
	  }
	  solveOverlapsUp(travel,pos,dir,t,body.parent,nextBody);
	  return true;
	}
      }
    }
    
    //The parent is not overlpaed, check their sisters
    return solveOverlapsFlat(travel,pos,dir,t,ibody,nextBody);
  }

  inline bool solveOverlapsDown(const double travel,
				v3D& pos,
				v3D& dir,
				const double t,
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
    for(unsigned idaught = 0; idaught < body.daughters.size(); ++idaught){
        
      //Get daught index and reference
      const unsigned daughtIndex = body.daughters[idaught];
      const pen_meshBody& daught = bodies[daughtIndex];

      //Check if this daught can overlap with ibody
      if(daught.canOverlapParent){

	//Apply the daughter inverse transform if needed.
	//This will convert the position and direction to daughter local coordinates
	v3D dPos = pos;
	v3D dDir = dir;
	if(daught.inAnimation(t)){
	  daught.readAnimation().applyInv(t, dPos, dDir);
	}
      
	//Check if is crossed in this direction
	double dsOverlap;
	if(daught.cross(dPos, dDir, dsOverlap, false)){
	  
	  //Is crossed, check if the cross is before the travel end
	  if(dsOverlap - travel < threshold){
	    //Daughter is crossed, thus, the next body must be
	    //updated to this daugher. However, Their daughters
	    //must be checked also
	    
	    nextBody = daughtIndex;
	    solveOverlapsDown(travel,dPos,dDir,t,daughtIndex,nextBody);
	    
	    //Save position and direction in final body's local coordinates
	    pos = dPos;
	    dir = dDir;
	    return true;
	  }
	}
      }
    }
    
    //No overlapes
    return false;
  }

  inline void composeTransform(const unsigned ibody,
			       v3D& pos,
			       v3D& dir,
			       const double t) const {

    const pen_meshBody& body = bodies[ibody];

    //Apply local transformation if required
    if(body.inAnimation(t)){
      body.readAnimation().apply(t, pos, dir);
    }
    
    if(ibody != iworld){
      //Apply parent transform
      composeTransform(body.parent, pos, dir, t);
    }
  }

  inline void composeTransform(const unsigned ibody,
			       v3D& pos,
			       const double t) const {

    const pen_meshBody& body = bodies[ibody];

    //Apply local transformation if required
    if(body.inAnimation(t)){
      body.readAnimation().apply(t, pos);
    }
    
    if(ibody != iworld){
      //Apply parent transform
      composeTransform(body.parent, pos, t);
    }
  }  
  
  inline void composeInvTransform(const unsigned ibody,
				  v3D& pos,
				  v3D& dir,
				  const double t) const {

    const pen_meshBody& body = bodies[ibody];
    if(ibody != iworld){
      //Apply parent inverse transform
      composeInvTransform(body.parent, pos, dir, t);
    }
    //Compose current body transformation, if needed
    if(body.inAnimation(t)){
      body.readAnimation().applyInv(t, pos, dir);
    }
  }

  inline void composeInvTransform(const unsigned ibody,
				  v3D& pos,
				  const double t) const {

    const pen_meshBody& body = bodies[ibody];
    if(ibody != iworld){
      //Apply parent inverse transform
      composeInvTransform(body.parent, pos, t);
    }
    //Compose current body transformation, if needed
    if(body.inAnimation(t)){
      body.readAnimation().applyInv(t, pos);
    }
  }

  inline bool isTransformable() const override { return true; }
  
};

template<>
constexpr const char* penred::errors::errorMessage<pen_meshBodyGeo>(const int val) noexcept {
  return pen_meshBodyGeo::errorMessage(val);
}

#endif
