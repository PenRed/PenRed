
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
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//


#ifndef __PEN_CLASSES__
#define __PEN_CLASSES__

#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <type_traits>
#include <array>
#include <functional>
#include <algorithm>
#include <numeric>

#include "../states/pen_baseState.hh"
#include "pen_constants.hh"
#include "../errors/pen_errors.hh"
#include "../rands/includes/pen_random.hh"
#include "../parsers/internalData/includes/pen_data.hh"

template<class particleType, class contextType, class materialType> class abc_interaction;
template<class stateType, class contextType> class pen_particle;
template<class stateType> class pen_particleStack;
class wrapper_geometry;

//--------------------------------
// Auxiliar structs and classes
//--------------------------------

template<class E, class T> struct container;

template<class T>
struct vector3D{
  T x,y,z;
    
  vector3D() {}
  vector3D(T xin, T yin, T zin) : x(xin), y(yin), z(zin){}
    
  inline vector3D operator+(const vector3D& a) const{
    return vector3D(x + a.x, y + a.y, z + a.z);
  }

  inline vector3D operator+(const T a) const{
    return vector3D(x + a, y + a, z + a);
  }
  
  inline vector3D operator-(const vector3D& a) const{
    return vector3D(x - a.x, y - a.y, z - a.z);
  }

  inline vector3D operator-(const T a) const{
    return vector3D(x - a, y - a, z - a);
  }
  
  inline T operator*(const vector3D& a) const{
    return x*a.x + y*a.y + z*a.z;
  }

  inline vector3D operator*(const T a) const{
    return vector3D(a*x, a*y, a*z);
  }

  inline vector3D operator/(const T a) const{
    return vector3D(x/a, y/a, z/a);
  }
    
  inline vector3D operator^(const vector3D& a) const{
    return vector3D(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
  }
    
  inline void add(const vector3D& a){
    x += a.x;
    y += a.y;
    z += a.z;
  }

  inline void add(const T a){
    x += a;
    y += a;
    z += a;
  }
    
  inline void subs(const vector3D& a){
    x -= a.x;
    y -= a.y;
    z -= a.z;
  }

  inline void subs(const T a){
    x -= a;
    y -= a;
    z -= a;
  }

  inline void pow(const unsigned a){
    x = pow(x,a);
    y = pow(y,a);
    z = pow(z,a);
  }  

  inline void crossProd(const vector3D& a){
    const double auxX = y*a.z - z*a.y;
    const double auxY = z*a.x - x*a.z;
    z = x*a.y - y*a.x;
    x = auxX;
    y = auxY;
  }
    
  inline T mod2() const{
    return x*x + y*y + z*z;
  }
    
  inline T mod() const{
    return sqrt(mod2());
  }

  inline T dist(const vector3D& a) const{
    return (a - *this).mod();
  }

  inline vector3D dir(const vector3D& a) const{
    vector3D normDir = (a - *this);
    normDir.normalize();
    return normDir;
  }
  
  inline void normalize(){
    T norm = mod();
    x /= norm;
    y /= norm;
    z /= norm;
  }

  inline std::string stringify() const{
    char str[60];
    sprintf(str,"(%.5E,%.5E,%.5E)",
	    static_cast<double>(x),
	    static_cast<double>(y),
	    static_cast<double>(z));
    return std::string(str);
  }
  
};

template<class T>
struct triangle{
  
  vector3D<T> v1;
  vector3D<T> v2;
  vector3D<T> v3;

  inline T minx() const {
    return std::min(v1.x, std::min(v2.x, v3.x));
  }

  inline T miny() const {
    return std::min(v1.y, std::min(v2.y, v3.y));
  }

  inline T minz() const {
    return std::min(v1.z, std::min(v2.z, v3.z));
  }

  inline T maxx() const {
    return std::max(v1.x, std::max(v2.x, v3.x));
  }

  inline T maxy() const {
    return std::max(v1.y, std::max(v2.y, v3.y));
  }

  inline T maxz() const {
    return std::max(v1.z, std::max(v2.z, v3.z));
  }

  inline T dx() const{
    return maxx()-minx();
  }
  inline T dy() const{
    return maxy()-miny();
  }
  inline T dz() const{
    return maxz()-minz();
  }
  
  inline unsigned minxi() const{
    if(v1.x < v2.x && v1.x < v3.x)
      return 1;
    if(v2.x < v3.x)
      return 2;
    if(v3.x < v2.x)
      return 3;
    return 0;
  }
  
  inline unsigned minyi() const{
    if(v1.y < v2.y && v1.y < v3.y)
      return 1;
    if(v2.y < v3.y)
      return 2;
    if(v3.y < v2.y)
      return 3;
    return 0;
  }  
  
  inline unsigned minzi() const{
    if(v1.z < v2.z && v1.z < v3.z)
      return 1;
    if(v2.z < v3.z)
      return 2;
    if(v3.z < v2.z)
      return 3;
    return 0;
  }

  inline vector3D<T> centroid() const{
    return (v1 + v2 + v3)/static_cast<T>(3);
  }
  
  inline std::string stringify() const{
    return v1.stringify() + " " + v2.stringify() + " " + v3.stringify();
  }
};

template<class T>
struct box{

private:
  vector3D<T> min;
  vector3D<T> max;

  vector3D<T> d;  
  
public:

  box() : min(0,0,0),
	  max(1,1,1),
	  d(1,1,1)
  {}
  box(T xmin, T ymin, T zmin, T xmax, T ymax, T zmax) : min(xmin,ymin,zmin),
							max(xmax,ymax,zmax),
							d(xmax-xmin,ymax-ymin,zmax-zmin)
  {}

  box(const triangle<T>& t){
    min.x = t.minx();
    min.y = t.miny();
    min.z = t.minz();

    max.x = t.maxx();
    max.y = t.maxy();
    max.z = t.maxz();

    d.x = max.x - min.x;
    d.y = max.y - min.y;
    d.z = max.z - min.z;
  }
  
  inline bool in(const vector3D<T> &p, const T threshold) const{
    return p.x - min.x > -threshold && max.x - p.x > -threshold &&
      p.y - min.y > -threshold && max.y - p.y > -threshold &&
      p.z - min.z > -threshold && max.z - p.z > -threshold;
  }


  inline bool in(const triangle<T> &t, const T threshold) const{
    return in(t.centroid(),threshold);
  }
  /*
  inline bool in(const triangle<T> &t, const T threshold) const{

    double ds;
    bool goIn;
    vector3D<T> dir;
    double dsMax;

    if(in(t.v1,threshold) || in(t.v2,threshold) || in(t.v3,threshold))
      return true;
    
    // * v1 to v2

    // Get dir
    dir = t.v2-t.v1;
    dsMax = dir.mod();
    dir.normalize();

    //Check cross
    goIn = toIn(t.v1, dir, ds);
    if(goIn && ds - dsMax < threshold)
      return true;
    
    // * v1 to v3

    // Get dir
    dir = t.v3-t.v1;
    dsMax = dir.mod();
    dir.normalize();

    //Check cross
    goIn = toIn(t.v1, dir, ds);
    if(goIn && ds - dsMax < threshold)
      return true;

    // * v2 to v3

    // Get dir
    dir = t.v3-t.v2;
    dsMax = dir.mod();
    dir.normalize();

    //Check cross
    goIn = toIn(t.v2, dir, ds);
    if(goIn && ds - dsMax < threshold)
      return true;
    
    return false;
  }
  */
  
  inline bool in(const box<T> &b, const T threshold) const{
    return overlap(b, threshold) > static_cast<T>(0);
  }
  
  inline bool commonBoundary(const box<T> b, const T threshold) const{
    //Check if both boxes share a boundary plane
    return fabs(max.x - b.max.x) < threshold  ||
      fabs(min.x - b.min.x) < threshold ||
      fabs(max.y - b.max.y) < threshold ||
      fabs(min.y - b.min.y) < threshold ||
      fabs(max.z - b.max.z) < threshold ||
      fabs(min.z - b.min.z) < threshold;
  }

  inline T overlapx(const box<T> b) const{
    //Calculate the boxes overlap volume
    T T0 = static_cast<T>(0);
    return std::max(T0, std::min(max.x, b.max.x) - std::max(min.x, b.min.x));
  }
  inline T overlapy(const box<T> b) const{
    //Calculate the boxes overlap volume
    T T0 = static_cast<T>(0);
    return std::max(T0, std::min(max.y, b.max.y) - std::max(min.y, b.min.y));
  }
  inline T overlapz(const box<T> b) const{
    //Calculate the boxes overlap volume
    T T0 = static_cast<T>(0);
    return std::max(T0, std::min(max.z, b.max.z) - std::max(min.z, b.min.z));
  }
  
  inline T overlapVol(const box<T> b) const{
    //Calculate the boxes overlap volume
    T T0 = static_cast<T>(0);
    T xoverlap = std::max(T0, std::min(max.x, b.max.x) - std::max(min.x, b.min.x));
    T yoverlap = std::max(T0, std::min(max.y, b.max.y) - std::max(min.y, b.min.y));
    T zoverlap = std::max(T0, std::min(max.z, b.max.z) - std::max(min.z, b.min.z));
      
    return xoverlap*yoverlap*zoverlap;
  }

  inline T overlap(const box<T>& b, const T threshold) const{
    //Calculate the boxes overlap volume
    const T T0 = static_cast<T>(0);
    const T T1 = static_cast<T>(1);    
    T xoverlap = T1;
    T yoverlap = T1;
    T zoverlap = T1;

    // X
    if(b.dx() > T0 && dx() > T0){ // Both are 3D
      xoverlap = std::max(T0, std::min(max.x, b.max.x) - std::max(min.x, b.min.x));
    }
    else if(dx() > T0){ // Local is not 2D
      //Check if b is out of local bounds
      if(b.min.x - min.x < -threshold || b.max.x - max.x > threshold){
	xoverlap = T0;
      }
    }
    else if(b.dx() > T0){ // b is not 2D
      //Check if local is out of b bounds
      if(min.x - b.min.x < -threshold || max.x - b.max.x > threshold){
	xoverlap = T0;
      }
    }
    else if(fabs(minx()-b.minx()) > threshold){ // Both are 2D, check if are close enough
      return T0;
    }
      
    // Y
    if(b.dy() > T0 && dy() > T0){ // Both are 3D
      yoverlap = std::max(T0, std::min(max.y, b.max.y) - std::max(min.y, b.min.y));
    }
    else if(dy() > T0){ // Local is not 2D
      //Check if b is out of local bounds
      if(b.min.y - min.y < -threshold || b.max.y - max.y > threshold){
	yoverlap = T0;
      }
    }
    else if(b.dy() > T0){ // b is not 2D
      //Check if local is out of b bounds
      if(min.y - b.min.y < -threshold || max.y - b.max.y > threshold){
	yoverlap = T0;
      }
    }
    else if(fabs(miny()-b.miny()) > threshold){ // Both are 2D, check if are close enough
      return T0;
    }
    
    // Z
    if(b.dz() > T0 && dz() > T0){ // Both are 3D
      zoverlap = std::max(T0, std::min(max.z, b.max.z) - std::max(min.z, b.min.z));
    }
    else if(dz() > T0){ // Local is not 2D
      //Check if b is out of local bounds
      if(b.min.z - min.z < -threshold || b.max.z - max.z > threshold){
	zoverlap = T0;
      }
    }
    else if(b.dz() > T0){ // b is not 2D
      //Check if local is out of b bounds
      if(min.z - b.min.z < -threshold || max.z - b.max.z > threshold){
	zoverlap = T0;
      }
    }
    else if(fabs(minz()-b.minz()) > threshold){ // Both are 2D, check if are close enough
      return T0;
    }
      
    return xoverlap*yoverlap*zoverlap;
  }
  
  inline double overlapfact(const box<T>& b, const T threshold) const{

    T o = overlap(b,threshold);
    T maxAV = std::max(volume(),area());
    T maxAVb = std::max(b.volume(),b.area());
    return static_cast<double>(o)/static_cast<double>(std::min(maxAV,maxAVb));
    
  }
  
  inline T volume() const{ return d.x*d.y*d.z; }

  inline T areax() const{ return d.y*d.z; }
  inline T areay() const{ return d.x*d.z; }
  inline T areaz() const{ return d.x*d.y; }
  inline T area() const{ return std::max(areax(),std::max(areay(),areaz())); }

  inline T minx() const {return min.x;}
  inline T miny() const {return min.y;}
  inline T minz() const {return min.z;}

  inline T maxx() const {return max.x;}
  inline T maxy() const {return max.y;}
  inline T maxz() const {return max.z;}
  
  inline T dx() const {return d.x;}
  inline T dy() const {return d.y;}
  inline T dz() const {return d.z;}

  inline void set(const T xmin, const T ymin, const T zmin,
		  const T xmax, const T ymax, const T zmax){
    min = vector3D<double>(xmin,ymin,zmin);
    max = vector3D<double>(xmax,ymax,zmax);
    d = max-min;
  }

  inline void set(const triangle<T>& t){
    min.x = t.minx();
    min.y = t.miny();
    min.z = t.minz();
    
    max.x = t.maxx();
    max.y = t.maxy();
    max.z = t.maxz();
    
    d.x = max.x - min.x;
    d.y = max.y - min.y;
    d.z = max.z - min.z;
  }
  
  inline void set(const box newBox){
    min = newBox.min;
    max = newBox.max;
    d = max-min;
  }

  inline void setX(const T xmin, const T xmax){
    min.x = xmin; max.x = xmax;
    d.x = max.x-min.x;
  }

  inline void setY(const T ymin, const T ymax){
    min.y = ymin; max.y = ymax;
    d.y = max.y-min.y;
  }

  inline void setZ(const T zmin, const T zmax){
    min.z = zmin; max.z = zmax;
    d.z = max.z-min.z;
  }
  
  inline bool toIn(const double x, 
		   const double y, 
		   const double z,
		   const double u,
		   const double v, 
		   const double w,
		   double& ds) const{
    
    //Check if the particle reaches the mesh in
    //each axis
    double dsxIn, dsyIn, dszIn;
    double dsxOut, dsyOut, dszOut;
    if(moveIn(x - min.x, u, dsxIn, dsxOut, d.x) &&
       moveIn(y - min.y, v, dsyIn, dsyOut, d.y) &&
       moveIn(z - min.z, w, dszIn, dszOut, d.z)){
      


      //Move the particle the required distance to reach
      //geometry mesh on all axis.
      double ds2allIn = std::max(std::max(dsxIn,dsyIn),dszIn);
      //Check if the particle go out of the mesh when this distance is traveled
      if(ds2allIn > 0.0){
        double ds2someOut = std::min(std::min(dsxOut,dsyOut),dszOut);
        if(ds2allIn >= ds2someOut){
          return false;
        }
      }

      ds = ds2allIn;

      return true;
    }
    return false;

  }
  
  inline bool toIn(const vector3D<T> p,
		   const vector3D<T> dir,
		   double& ds) const{
    return toIn(p.x, p.y, p.z, dir.x, dir.y, dir.z, ds);
  }

  static inline bool moveIn(const double pos,
                            const double dir,
                            double& ds2in,
                            double& ds2out,
                            const double max){

    //Check if a particle at position 'pos' and direction 'dir'
    //is in the geometry limits or will enter to it.  
      
    // This function asumes that geometry is in the interval [0,max)
    // Also, stores in 'ds2in' and 'ds2out' the distance to enter and exit the region
    
    const double inf = 1.0e35;
    const double eps = 1.0e-8;

    if(std::signbit(pos)){ //The position is negative

      if(std::signbit(dir)){ //The direction is negative too, can't reach the mesh
        return false;
      }

      if(dir == 0.0){
	//Particle is not moving on this axis, so never will reaches the mesh
	return false;
      }

      //The direction is positive
      ds2in  = eps-pos/dir; //Distance to travel to reach the origin (0)
      ds2out = (max-pos)/dir-eps; //Distance to travel to reach the "max"

      return true;

    }else{

      //The position is positive
      if(pos >= max){ //It is outside the region

        if(!std::signbit(dir)){ //The direction is positive too, can't reach the mesh
          return false;
        }

        if(dir == 0.0){
	  //Particle is not moving on this axis, so never reaches the mesh
	  return false;
        }

        //The direction is negative
        ds2in  = (max-pos)/dir+eps; //Distance to travel to reach the "max"
        ds2out = pos/(-dir)-eps; //Distance to travel to reach the origin (0)

        return true;

      } else {
        //It is inside the region
        ds2in = 0.0;

        if(dir == 0.0){
	  //Particle is not moving on this axis, so never will escape the limits
	  ds2out = inf;
        }else if(std::signbit(dir)){ //Negative direction
          ds2out = pos/(-dir)-eps;
        }else{ //Positive direction
          ds2out = (max-pos)/dir-eps;
        }

        return true;
      }

    }

    /*
      if(std::signbit(pos) || pos >= max){
      // position component is out of geometry,
      // check if is moving to it
      if(dir == 0.0){
      //Particle is not moving on this axis, so never will reaches the mesh
      return false;
      }

      if(std::signbit(dir) == std::signbit(pos)){
      //Particle is moving away from geometry
      return false;	
      }

      //Particle will reach the geometry
      if(std::signbit(dir)){
      //Is moving on negative direction
      ds2in  = (max-pos)/dir+eps;
      ds2out = pos/fabs(dir)-eps;
      }
      else{
      //Particle is moving on positive direction
      ds2in  = fabs(pos)/dir+eps;
      ds2out = (max-pos)/dir-eps;
      }
      return true;
      }
      ds2in = 0.0;
      if(dir == 0.0){
      ds2out = inf;
      return true;
      }

      if(std::signbit(dir)){
      //Negative direction
      ds2out = pos/fabs(dir)-eps;
      }
      else{
      //Positive direction
      ds2out = (max-pos)/dir-eps;
      }
      return true; //Particle is in geometry limits
    */
  }

  inline void enlarge(const vector3D<T> p){
    if(p.x < min.x) {min.x = p.x; d.x = max.x - min.x;}
    else if(p.x > max.x) {max.x = p.x; d.x = max.x - min.x;}

    if(p.y < min.y) {min.y = p.y; d.y = max.y - min.y;}
    else if(p.y > max.y) {max.y = p.y; d.y = max.y - min.y;}

    if(p.z < min.z) {min.z = p.z; d.z = max.z - min.z;}
    else if(p.z > max.z) {max.z = p.z; d.z = max.z - min.z;}
  }

  inline void enlarge(const triangle<T> t){
    enlarge(t.v1);
    enlarge(t.v2);
    enlarge(t.v3);
  }

  inline void enlarge(const box<T> b){
    enlarge(b.min);
    enlarge(b.max);
  }

  inline void enlarge(const T e){
    min.subs(e);
    max.add(e);
    d = max - min;
  }  
  
  inline std::string stringify() const{
    return min.stringify() + " - " + max.stringify();
  }
};

template<class E, class T>
struct container : public box<T>{

  std::vector<E> elements;

  container(T xmin, T ymin, T zmin,
	    T xmax, T ymax, T zmax) : box<T>(xmin,ymin,zmin,
					     xmax,ymax,zmax)
  {}
  container(const E& e) : box<T>(e)
  {}
  container(const box<T>& b) : box<T>(b)
  {}

  inline size_t nElements() const {return elements.size();}

  static size_t split(const unsigned n,
		      const unsigned index2Split,
		      const unsigned meanElements,
		      const T threshold,
		      std::vector<container<E,T>>& regions){

    const T sizeX = regions[index2Split].box<T>::dx()/static_cast<T>(n);
    const T sizeY = regions[index2Split].box<T>::dy()/static_cast<T>(n);
    const T sizeZ = regions[index2Split].box<T>::dz()/static_cast<T>(n);

    const T ox = regions[index2Split].box<T>::minx();
    const T oy = regions[index2Split].box<T>::miny();
    const T oz = regions[index2Split].box<T>::minz();
    
    std::vector<E> elementsCopy(regions[index2Split].elements);

    size_t origReg = regions.size(); 

    //Saves the number of performed splits
    unsigned nSplit = 0;
    //Saves remaining elements to store
    size_t remaining = regions[index2Split].nElements();

    //Create subregions
    for(unsigned k = 0; k < n; ++k){
      T zmin = oz + static_cast<T>(k)*sizeZ;
      T zmax = oz + static_cast<T>(k+1)*sizeZ;
      for(unsigned j = 0; j < n; ++j){
	T ymin = oy + static_cast<T>(j)*sizeY;
	T ymax = oy + static_cast<T>(j+1)*sizeY;
	for(unsigned i = 0; i < n; ++i){
	  T xmin = ox + static_cast<T>(i)*sizeX;
	  T xmax = ox + static_cast<T>(i+1)*sizeX;

	  //Initiate a region with init boundaries
	  container region(xmin, ymin, zmin, xmax, ymax, zmax);
	  //Create a box to store extended boundaries
	  box<T> finalRegionBox;

	  //Iterate over all unassigned elements
	  size_t next = 0;
	  while(next < remaining){

	    //Check if the element is inside this region
	    if(region.in(elementsCopy[next],threshold)){

	      //Enlarge the final region if required
	      if(region.nElements() == 0){
		finalRegionBox.set(elementsCopy[next]);
	      }else{
		finalRegionBox.enlarge(elementsCopy[next]);
	      }
	      
	      //Save this element in the region
	      region.elements.push_back(elementsCopy[next]);
	      elementsCopy[next] = elementsCopy[--remaining];

	    }else{
	      ++next;
	    }
	  }

	  //Check if there is some remaining element
	  if(remaining == 0){
	    //Check if the region has been split
	    if(nSplit == 0){
	      //It is not possible to split this region
	      return 0;
	    }
	  }
	  
	  //Once all elements in region have been assigned,
	  //set the final region box
	  region.set(finalRegionBox);
	  //Enlarge the region to avoid rounelementsding errors
	  region.enlarge(threshold);
	  
	  //Ensure that the region contains at least one element
	  if(region.nElements() > 0){

	    //Increase the number of divisions, although the
	    //new region could be merged with another one
	    ++nSplit;

	    //Check if the region has a low number of elements
	    if(region.nElements() < meanElements){
	      //Check if can be combined with some previous region to get a
	      //number of trinagles closer to the expected

	      double maxOverlapFact = 0.0;
	      unsigned maxOverlapReg = regions.size();
	      for(unsigned ireg = 0; ireg < regions.size(); ++ireg){		
		if(index2Split == ireg)
		  continue;
		//Check if this region can include more elements
		if(regions[ireg].nElements() + region.nElements() > 2.0*meanElements)
		  continue;
		//Check if both regions overlap
		double overlapFact = region.overlapfact(regions[ireg],threshold);
		if(overlapFact > maxOverlapFact){
		  maxOverlapFact = overlapFact;
		  maxOverlapReg = ireg;
		}
	      }
	      
	      //Check if some overlaping region has been found
	      if(maxOverlapReg < regions.size()){
		
		//Merge both regions
		container<E,T>& region2merge = regions[maxOverlapReg];

		//Enlarge the region to fit both regions
		region2merge.enlarge(finalRegionBox);

		//Reserve capacity for new elements
		const size_t newCap =
		  region2merge.nElements() + region.nElements();
		region2merge.elements.reserve(newCap);

		//Copy elements to the max overlaping region
		for(const E& t : region.elements){
		  region2merge.elements.push_back(t);
		}
	      }
	      else{
		//Save a new region
		regions.push_back(region);
	      }
	    }else{

	      //Check if this region intersects with another one
	      double maxOverlapFact = 0.0;
	      unsigned maxOverlapReg = regions.size();
	      for(unsigned ireg = 0; ireg < regions.size(); ++ireg){		
		if(index2Split == ireg)
		  continue;
		//Check if this region can include more elements
		if(regions[ireg].nElements() + region.nElements() > 2.0*meanElements)
		  continue;
		//Check if both regions overlap
		double overlapFact = region.overlapfact(regions[ireg],threshold);
		if(overlapFact > maxOverlapFact){
		  maxOverlapFact = overlapFact;
		  maxOverlapReg = ireg;
		}
	      }
	      
	      if(maxOverlapFact > 0.8){
		
		//Merge both regions
		container<E,T>& region2merge = regions[maxOverlapReg];

		//Enlarge the region to fit both regions
		region2merge.enlarge(finalRegionBox);

		//Reserve capacity for new elements
		const size_t newCap =
		  region2merge.nElements() + region.nElements();
		region2merge.elements.reserve(newCap);

		//Copy elements to the max overlaping region
		for(const E& t : region.elements){
		  region2merge.elements.push_back(t);
		}
	      }else{
		//Save a new region
		regions.push_back(region);
	      }
	    }
	  }	  
	}
      }
    }

    if(remaining != 0){

      printf("\n");
      for(size_t ir = origReg; ir < regions.size(); ++ir){
	printf("\n");
	printf("Region: %s\n", regions[ir].stringify().c_str());	
	for(size_t i = 0; i < remaining; ++i){
	  printf("%s: %s\n",
		 regions[index2Split].elements[i].stringify().c_str(),
		 regions[ir].in(elementsCopy[i],threshold) ? "in" : "out");
	}
      }
      printf("\n");

      
      char errmsg[200];
      sprintf(errmsg,"container:split: Error: %lu/%lu elements "
	      "out of defined regions.\n",
	      remaining,regions[index2Split].elements.size());      
      throw std::runtime_error(errmsg); 
    }

    if(nSplit > 0){
      //The container has been split, remove the original one to avoid duplicating elements
      regions.erase(regions.begin() + index2Split);
    }

    return nSplit;
  }


  static void splitUntil(const unsigned meanElements,
			 const T threshold,
			 std::vector<container<E,T>>& containers,
			 const size_t limit) {

    size_t initElements = 0;
    for(const auto& continer: containers){
      initElements += continer.nElements();
    }
    
    //Maximum achieved containers
    size_t maxContainers = 1;
    //Number of iterations with a number of containers
    //under the maximum achieved (maxContainers)
    size_t nIterUnderMaxContainers = 0;
    while(containers.size() < limit){

      bool splitDone = false;
      size_t iCont = 0;
      size_t nCont = containers.size();
      while(iCont < nCont){
	if(containers[iCont].nElements() > 1.5*meanElements){
	  
	  //Divide the container region and add them
	  size_t nDivisions = split(2,iCont,meanElements,
				    threshold,
				    containers);
	  if(nDivisions == 0){
	    nDivisions = split(5,iCont,meanElements,
			       threshold,
			       containers);
	  }

	  if(nDivisions > 0){
	    splitDone = true;
	    --nCont;	    
	  }else{	    
	    ++iCont;
	  }
	}else{	    
	  ++iCont;
	}
      }
      //Check if some split has been done
      if(!splitDone){
	//No more containers required, break the loop
	break;
      }

      //Check if the number of containers has reached a new maximum
      if(containers.size() > maxContainers){
	maxContainers = containers.size();
	nIterUnderMaxContainers = 0;
      }else{
	//Check if the number of iterations under the maximum number
	//of containers is greater than 100
	if(nIterUnderMaxContainers >= 100){
	  //Is in a infinite loop, break it
	  break;
	}
	//Increase the iteration counter
	++nIterUnderMaxContainers;
      }
    }

    //Check the integrity of the split process
    size_t finalElements = 0;
    for(const auto& continer: containers){
      finalElements += continer.nElements();
    }
    if(finalElements != initElements){
      printf("splitUntil: Error: Elements lost on container split: "
	     "      Expected elements : %lu\n"
	     "      Final elements    : %lu\n"
	     " Please, report this issue.\n",
	     initElements, finalElements);
      fflush(stdout);
      throw std::range_error("Lost elements on container split");      
    }

      
  }

};

//-------------------
// Materials
//-------------------

class abc_material{

private:
  
protected:
  
  bool initialized;
  
public:

  //  ----  Material densities and its reciprocals.
  double DEN;    //Density
  double RDEN;   //1/Density
  
  //  ----  Absorption energies, EABS(KPAR,MAT).
  double EABS[constants::nParTypes];  
  
  abc_material() : initialized(false), DEN(1.0), RDEN(1.0){
    for(unsigned i = 0; i < constants::nParTypes; i++){
      EABS[i] = 50.0; //eV
    }
  }
  inline bool initDone() const { return initialized;}
  inline double getEABS(const unsigned kpar) const {return EABS[kpar];}
  inline void setEABS(const unsigned kpar, const double eabs){
    if(eabs > 0.0){
      EABS[kpar] = eabs;
    }
  }
  inline void setDens(const double dens){
    if(dens > 1.0e-17){
      DEN = dens;
      RDEN = 1.0/dens;
    }
  }

  inline double readDens()const {return DEN;}
  inline double readIDens()const {return RDEN;}
  
  virtual ~abc_material(){};
  
};

//-------------------
// Geometry
//-------------------

class wrapper_geometry{

protected:
  int configStatus;
public:

  std::string name;
    
  wrapper_geometry() : configStatus(0),
		       name("unnamed") {}

  inline int configureStatus() const {return configStatus;}
  inline virtual const char* getType() const {return "UNKNOWN";}
  
  virtual void locate(pen_particleState& state) const = 0;
  virtual void step(pen_particleState& state, double DS, double &DSEF, double &DSTOT, int &NCROSS) const = 0;
  virtual int configure(const pen_parserSection& config, const unsigned verbose) = 0;
  virtual void usedMat(bool[constants::MAXMAT+1]) const = 0;
  virtual double getEabs(const unsigned ibody, const unsigned kpar) const = 0;
  virtual double getDSMAX(const unsigned ibody) const = 0;
  virtual unsigned getDET(const unsigned ibody) const = 0;
  virtual unsigned getMat(const unsigned ibody) const = 0;
  virtual unsigned long getElements() const = 0;
  virtual unsigned getBodies() const = 0;
  virtual unsigned getIBody(const char* elementName) const = 0;
  virtual void getOffset(double* offset) const { offset[0] = 0.0; offset[1] = 0.0; offset[2] = 0.0; }
  virtual ~wrapper_geometry(){}
  
};

//-------------------
// Variance reduction
//-------------------

template<class stateType>
class abc_VR{

public:
  virtual void run_particleStack(const unsigned long long nhist,
				 const pen_KPAR kpar,
				 const unsigned kdet,
				 stateType& state,
				 std::array<stateType,constants::NMS>& stack,
				 unsigned& created,
				 const unsigned available,
				 pen_rand& random) const = 0;

  virtual void run_matChange(const unsigned long long nhist,
			     const pen_KPAR kpar,
			     const unsigned prevMat,			    
			     stateType& state,
			     std::array<stateType,constants::NMS>& stack,
			     unsigned& created,
			     const unsigned available,
			     pen_rand& random) const = 0;

  virtual void run_interfCross(const unsigned long long nhist,
			       const pen_KPAR kpar,
			       const unsigned kdet,
			       stateType& state,
			       std::array<stateType,constants::NMS>& stack,
			       unsigned& created,
			       const unsigned available,
			       pen_rand& random) const = 0;

  virtual ~abc_VR(){}
};

//-------------------
// Interactions
//-------------------

template<class particleType, class contextType, class materialType>
class abc_interaction{

 protected:
  const int ID;

 public:

  abc_interaction() : ID(-1) {}
  abc_interaction(const int id) : ID(id) {}
  virtual void init(const contextType&){};
  virtual double iMeanFreePath(const materialType&, const particleType&) const = 0;
  virtual int interact(const contextType&, const materialType&, particleType&, double&, pen_rand&) const = 0;
  virtual int interactF(const contextType&, const materialType&, particleType&, double&, pen_rand&) const = 0;  
  inline int getID() const {return ID;}
  virtual ~abc_interaction(){}
};

//-------------------
// Particles
//-------------------

template<class stateType, class contextType, class materialType>
class abc_particle{

private:
  std::array<pen_particleState,constants::NMS> genericStates;
  std::array<stateType,constants::NMS> specificStates;
protected:

  //Particle energy at previous JUMP
  double ELAST1;

  //  **** Particle interaction constants
  
  double P[constants::MAXINTERACTIONS], ST;

  // ****  Energy grid variables for the next knock call.
  unsigned int KE;
  double XEL, XE, XEK;

  // **** Variance reduction variables

  double P0[constants::MAXINTERACTIONS];
  bool   LFORC[constants::MAXINTERACTIONS];

  //     KSOFTI ... take value 1 (0) if soft energy loss is (not) required
  int KSOFTI;

  //Particle state (This variable will be used to store particles in the stack)
  stateType state;

  //Particle current material pointer
  const materialType* pmat;
  
  //Particle last movement variables
  double dsef;   //Distance traveled on same material (Non void zones)
  double dstot;  //Total traveled distance (Included void zones)
  int ncross; //Is non zero if some interface has been crossed

  //Particle previous position variables
  unsigned MATL;
  unsigned IBODYL;
  double XL;
  double YL;
  double ZL;

  //Particle local parameters
  double EABS;
  double DSMAXbody;
  unsigned KDET;

  //Generic VR cluster
  const abc_VR<pen_particleState>* genericVR;
  //Specific VR cluster
  const abc_VR<stateType>* specificVR;
  //Reference to particle stack
  pen_particleStack<stateType>& stack;
  
public:
  
  const contextType& context;
  const pen_KPAR kpar;
  const unsigned int interactions;

  //Energy deposited at material on particle annihilation at rest (eV)
  const double annihilationEDep;
  //  ---- 

  abc_particle(const contextType& contextIn,
	       const pen_KPAR KPAR,
	       const unsigned int nInt,
	       const double annihilationEDepIn,
	       pen_particleStack<stateType>& stackIn) : KSOFTI(0),
							pmat(nullptr),
							dsef(0.0),
							dstot(0.0),
							ncross(0),
							MATL(0),
							IBODYL(0),
							XL(0.0),
							YL(0.0),
							ZL(0.0),
							EABS(0.0),
							DSMAXbody(1.0e35),
							KDET(0),
							genericVR(nullptr),
							specificVR(nullptr),
							stack(stackIn),
							context(contextIn),
							kpar(KPAR),
							interactions(nInt),
							annihilationEDep(annihilationEDepIn)
  {
    //Check if kpar is in range [0,nParTypes)
    if(KPAR < 0 || KPAR >= constants::nParTypes)
      {
	char error[300];
	sprintf(error,"Particle kpar (%d) out of range [0,%d).\n  Define an identifier for this particle at file 'kernel/particles/includes/pen_particles_ID.hh'",KPAR,constants::nParTypes);
	throw std::out_of_range(error);      
      }

    if(interactions > constants::MAXINTERACTIONS){
	char error[300];
	sprintf(error,"Number of interactions of particle with kpar (%d) out of range [0,%d)\n  Increment constant 'MAXINTERACTIONS' in constants namespace",KPAR,constants::MAXINTERACTIONS);
	throw std::out_of_range(error);      
    }

    //Check if context material is compatible
    if(!contextIn.template compatible<materialType>()){
      char error[300];
      sprintf(error,"Incompatible material type at particle instantation (kpar = %d).",KPAR);
      throw std::invalid_argument(error);
    }
  }

  virtual void START() = 0;
  virtual void JUMP(double &DS, pen_rand& penRand, const double DSMAX) = 0;
  virtual void JUMPF(double &DS, pen_rand& penRand, const double DSMAX) = 0;
  virtual void KNOCK(double &DE, int &ICOL, pen_rand& penRand) = 0;
  virtual void KNOCKF(double &DE, int &ICOL, pen_rand& penRand) = 0;

  virtual void softEloss(double& X,
			 double& Y,
			 double& Z,
			 double& DE,
			 pen_rand& penRand) = 0;
  
  virtual void dpage() = 0;
  virtual void page0() = 0;

  inline unsigned getDET() const {return KDET;}
  inline double getDSMAX() const {return DSMAXbody;}
  
  inline double getEABS() const {return EABS;}
  
  inline double DSef() const{return dsef;}
  inline double DStot() const{return dstot;}
  inline int NCross() const{return ncross;}

  inline void updateMat(){
    if(state.MAT > 0)
      pmat = context.template readMaterialPointer<materialType>(state.MAT-1);
    else
      pmat = nullptr;
  }
  inline void setMat(const unsigned imat){    
    state.MAT = imat;
    updateMat();
  }

  inline void updateBody(){
    if(state.IBODY < context.readGeometry()->getBodies()){
      EABS      = context.getEABS(state.IBODY,kpar);
      DSMAXbody = context.readGeometry()->getDSMAX(state.IBODY);
      KDET      = context.readGeometry()->getDET(state.IBODY);
    }
    else{
      EABS = 0.0;
      DSMAXbody = 1.0e35;
      KDET = 0;
    }
  }
  
  
  inline void setStep(const double dsefIn,
		      const double dstotIn,
		      const int ncrossIn){
    dsef = dsefIn;
    dstot = dstotIn;
    ncross = ncrossIn;
  }
  
  inline unsigned lastMat() const {return MATL;}
  inline unsigned lastBody() const {return IBODYL;}
  inline void lastPos(double& X, double& Y, double& Z){
    X = XL; Y = YL; Z = ZL;
  }

  inline void setLastMat(const unsigned lastMatIn) {MATL = lastMatIn;}
  inline void setLastBody(const unsigned lastBodyIn) {IBODYL = lastBodyIn;}
  inline void setLastPos(const double X, const double Y, const double Z) {
    XL = X; YL = Y; ZL = Z;
  }
  inline void saveLastPos(){
    XL = state.X; YL = state.Y; ZL = state.Z;
  }
  
  inline void jumpVolume(){

    //Store last position information
    IBODYL = state.IBODY;
    MATL = state.MAT;
    
    XL = state.X;
    YL = state.Y;
    ZL = state.Z;

    //Jump until particle crosses some interface or scapes from the geometry
    context.readGeometry()->step(state,1.0e30,dsef,dstot,ncross);
    //Calculate new particle age 
    if(state.LAGE){dpage();}

    if(MATL != state.MAT)
      updateMat();
    if(IBODYL != state.IBODY)
      updateBody();
  }
  
  inline void move(const double ds,
		   double& de,
		   double& softX,
		   double& softY,
		   double& softZ,
		   pen_rand& penRand){

    //Save actual position
    IBODYL = state.IBODY;
    MATL = state.MAT;
    
    XL = state.X;
    YL = state.Y;
    ZL = state.Z;

    //Move the particle
    context.readGeometry()->step(state, ds, dsef, dstot, ncross);
    //Calculate new particle age 
    if(state.LAGE){dpage();}

    //Check if material needs to be updated
    if(MATL != state.MAT)
      updateMat();
    if(IBODYL != state.IBODY)
      updateBody();
    
    //Check if soft energy deposition is required
    if(reqSoftELoss() == 1){
      //Calculate soft energy loss
      softEloss(softX,softY,softZ,de,penRand);
    }else{
      de = 0.0;
    }
  }
  
  virtual void annihilate(pen_rand&){}
  

  inline double xel(){return XEL;}
  inline double xe(){return XE;}
  inline double xek(){return XEK;}
  inline unsigned ke(){return KE;}

  inline void getGrid(unsigned& ke,
		      double& xel,
		      double& xe,
		      double& xek){
    ke = KE;
    xel = XEL; xe = XE; xek = XEK;
  }
  
  inline int reqSoftELoss(){return KSOFTI;}
  
  inline pen_KPAR getKpar(){
    return kpar;
  }

  inline void setBaseState(const pen_particleState& newState){
    state = newState;
  }
  
  inline void setState(const stateType& newState){
    state = newState;
  }

  inline pen_particleState& getBaseState(){
    return state;
  }
  
  inline stateType& getState(){
    return state;
  }

  inline const pen_particleState& readBaseState(){
    return state;
  }
  
  inline const stateType& readState() const{
    return state;
  }
  inline const contextType& readContext(){return context;}

  inline void registerGenericVR(const abc_VR<pen_particleState>& vrIn){
    genericVR = &vrIn;
  }

  inline void registerSpecificVR(const abc_VR<stateType>& vrIn){
    specificVR = &vrIn;
  }

  void vr_particleStack(const unsigned long long nhist,
			pen_rand& random,
			const unsigned verbose);  

  double vr_matChange(const unsigned long long nhist,
		      pen_rand& random,
		      const unsigned verbose);

  double vr_interfCross(const unsigned long long nhist,
			pen_rand& random,
			const unsigned verbose);
  
  void baseClear(){
    KSOFTI = 0;
    pmat = nullptr;
    dsef = 0.0;
    dstot = 0.0;
    ncross = 0;
    MATL = 0;
    IBODYL = 0;
    XL = 0.0;
    YL = 0.0;
    ZL = 0.0;
    state.reset();
  }
  
  virtual ~abc_particle(){}
  
  //CJUMP0
  //CJUMP1
  //CEGRID (only XEL, XE, XEK, and KE)  
};

//-------------------
// Particle stacks
//-------------------

class abc_particleStack{

 protected:
  unsigned int NSEC;
  
 public:
  abc_particleStack() : NSEC(0) {}

  inline unsigned int getNSec() const {return NSEC;}
  inline void cleans(){NSEC = 0;}

  virtual pen_particleState readBaseState(const unsigned i) const = 0;
};

template<class stateType> class pen_particleStack : public abc_particleStack{

 protected:
  stateType states[constants::NMS];
 public:
  
  pen_particleStack() : abc_particleStack() {}

  void store(const stateType& state)
  {
    if(NSEC < constants::NMS)
      {
	states[NSEC] = state;
	NSEC++;
      }
    else
      {
	printf("pen_particleStack:store:Warning: Stack full\n");
	//Stack is full remove particle with less energy
	unsigned int lessEpos = 0;
	double minE = 1.0e35;
	for(unsigned int i = 0; i < NSEC; i++)
	  {
	    if(states[i].E < minE)
	      {
		lessEpos = i;
		minE = states[i].E;
	      }
	  }
	//Store new particle in less energy position
	if(minE < state.E)
	  {
	    states[lessEpos] = state;
	  }
	//Create warning
	penError(ERR_store_FULL_STACK);
      }
  }
  
  inline unsigned get(stateType& state){
    if(NSEC > 0){
      NSEC--;
      state = states[NSEC];
      return NSEC;
    }
    return 0;
  }

  inline stateType readState(const unsigned i) const {
    return states[i];
  }

  inline pen_particleState readBaseState(const unsigned i) const{
    return pen_particleState(readState(i));
  }
  
  //CERSEC
  //SECST
  
};

//-------------------
// Context
//-------------------

template <class baseMat>
class abc_context{

private:

  //Array with cutoff energies
  double* maxEABS;

  unsigned geoBodies; // Controls geometry number of bodies

  // Geometry
  const wrapper_geometry* geometry;

  bool matsSet; //Controls if materials has been set

  // Materials
  baseMat* materials[constants::MAXMAT];
  unsigned nMats;

  void clearMats(){
    for(unsigned i = 0; i < constants::MAXMAT; i++){
      if(materials[i] != nullptr){
	delete materials[i];
	materials[i] = nullptr;
      }
    }
    matsSet = false;
    nMats = 0;
  }
  void clearGeo(){
    if(maxEABS != nullptr){
      delete [] maxEABS;
      maxEABS = nullptr;
    }
    geoBodies = 0;
    geometry = nullptr;
  }
  void clear(){
    clearMats();
    clearGeo();
  }
  
protected:

  
public:

  const static unsigned int NBV = pen_geoconst::NB;

  abc_context() : maxEABS(nullptr),
		  geoBodies(0),
		  geometry(nullptr),
		  matsSet(false),
		  nMats(0){
    
    //Set material pointers to null
    for(unsigned int i = 0; i < constants::MAXMAT; i++)
      materials[i] = nullptr;
  }

  //Function 'range' is intended to return the range in the specified
  //material for a particle with the specified energy and type
  virtual double range(const double E, const pen_KPAR kpar, const unsigned M) const = 0;

  virtual double avncol(const double E,
			const pen_KPAR kpar,
			const int icol,
			const unsigned imat) const = 0;

  virtual double avninter(const double E,
			  const pen_KPAR kpar,
			  const int icol,
			  const unsigned imat,
			  const bool calc_piecewise) const = 0;  

  virtual double getIF(const double forcer,
		       const pen_KPAR kpar,
		       const int icol,
		       const unsigned imat,
		       const bool calc_piecewise) const = 0;
  
  inline void getMatBaseArray(const abc_material* mats[constants::MAXMAT]) const {

    for(unsigned i = 0; i < constants::MAXMAT; i++){
      mats[i] = materials[i];
    }
  }

  inline double getMatEABS(const unsigned imat, const unsigned kpar){
    if(imat >= nMats){
      char error[300];
      sprintf(error,"getMatEABS: %d exceeds number of available materials (%d).",imat,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return materials[imat]->getEABS(kpar);
    }
  }
  inline double getEABS(const unsigned ibody, const unsigned kpar) const{
      unsigned index = kpar+ibody*constants::nParTypes;
      return maxEABS[index];
  }
  inline unsigned getDET(const unsigned ibody) const {
    return geometry->getDET(ibody);
  }
  inline double getDSMAX(const unsigned ibody) const {
    return geometry->getDSMAX(ibody);
  }
  
  int setGeometry(const wrapper_geometry* geoIn){
    if(geoIn == nullptr)
      return -1;

    //Clear previous geometry
    clearGeo();
    
    //Save pointer
    geometry = geoIn;

    //Get number of bodies in geometry
    geoBodies = geometry->getBodies();
    
    //Allocate memory for absorption energies of each body
    maxEABS = new double[(geoBodies+1)*constants::nParTypes];
    if(maxEABS == nullptr)
      return -2;
    
      
    return 0;
  }
  inline const wrapper_geometry* readGeometry() const {return geometry;}
  
  int updateEABS(){
    if(geometry == nullptr)
      return -1;

    //Check if the number of elements has been changed
    if(geoBodies != geometry->getBodies()){
      //Geometry has been changed, run set geometry again
      int err = setGeometry(geometry);
      if(err != 0)
	return -2;
    }

    //Fill absorption energies
    for(unsigned i = 0; i < geoBodies; i++){
      unsigned index0 = i*constants::nParTypes;

      //Get material index for this body
      int mat = geometry->getMat(i);
      if(mat < 0 || mat > (int)nMats){ //Check material and accept void
	//Index out of range
	return -3;
      }

      //Iterate over particle types
      for(unsigned j = 0; j < constants::nParTypes; j++){
	double Egeo = geometry->getEabs(i,j);
	if(mat < 1){ // Void
	  maxEABS[index0+j] = Egeo;
	}
	else{
	  double Emat = materials[mat-1]->getEABS(j);
	  if(Egeo > Emat){ maxEABS[index0+j] = Egeo;}
	  else{ maxEABS[index0+j] = Emat;}
	}
      }
    }
    

    return 0;
  }
  
  template<class derivedMat>
  int setMats(const unsigned M){

    //Check if materials has been already created
    if(matsSet)
      return -1;
    
    //Check index
    if(M > constants::MAXMAT || M < 1)
      return -2;

    //Clear materials (just in case)
    clearMats();
    
    //Create materials
    for(unsigned i = 0; i < M; i++){
      materials[i] = new derivedMat();
      if(materials[i] == nullptr){
	clearMats();
	return -3;
      }
    }
    nMats = M;
    matsSet = true;
    return 0;
  }

  unsigned getNMats() const {return nMats;}
  
  inline int setMatEABS(const unsigned M,
			const unsigned kpar,
			const double eabs){

    if(M >= constants::MAXMAT)
      return -1;

    if(materials[M] == nullptr)
      return -2;

    materials[M]->setEABS(kpar,eabs);
    return 0;
  }

  //Ensure that read and get material uses only convertible classes to baseMat
  template <class derivedMat>
  inline const typename std::enable_if<std::is_convertible<derivedMat,baseMat>::value, derivedMat&>::type
  readMaterial(const unsigned M) const{
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return *(static_cast<derivedMat*>(materials[M]));
    }
  }

  template <class derivedMat>
  inline const typename std::enable_if<std::is_convertible<derivedMat,baseMat>::value, derivedMat*>::type
  readMaterialPointer(const unsigned M) const{
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return static_cast<derivedMat*>(materials[M]);
    }
  }

  
  template <class derivedMat>
  inline typename std::enable_if<std::is_convertible<derivedMat,baseMat>::value, derivedMat&>::type
  getMaterial(const unsigned int M){
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return *(static_cast<derivedMat*>(materials[M]));
    }
  }

  inline const baseMat& readBaseMaterial(const unsigned M) const{
    
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
    
      return *materials[M];
    }
  }

  inline baseMat& getBaseMaterial(const unsigned M){
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return *materials[M];
    }
  }
  
  //Check material compatibility function
  template <class derivedMat>
  bool compatible() const{
    if(!matsSet)
      return false;

    derivedMat* pderived = nullptr;
    pderived = dynamic_cast<derivedMat*>(materials[0]);

    if(pderived == nullptr){
      return false;
    }
    return true;
  }
  
  virtual ~abc_context(){
    clear();
  };
};

//-------------------
// Energy grid
//-------------------

class abc_grid{

public:
  // ****  Energy grid and interpolation constants. The means of "Raw" calificative
  //       is that the variable has not been modifiqued by any transformation. For
  //       example, in a logarithmic scale, if the mean energy is 50eV, EL must be
  //       set to 50eV, not to log(50).
  
  // EMIN  : Minimum grid value
  // EL    : Raw lowest grid value (typically 0.99999*EMIN)
  // EU    : Raw last grid value
  // ET    : Raw Value of each bin (with no transformation)
  // DLEMP : Transformed ET (in a logarithmic scale, the "I" component must be log(ET[i]))
  // DLFC  : Inverse of distance between transofrmed bins
  // DLEMP1: Transformed EL (for example, in a logarithmic scale, log(EL))
  // 
  double EMIN, EL, EU, ET[constants::NEGP], DLEMP[constants::NEGP], DLEMP1, DLFC;

  bool initialized;

  abc_grid() : initialized(false){}
  
  virtual int init(double EMINu, double EMAXu) = 0;
  virtual void getInterval(const double E, int& KE, double& XEL, double& XE, double& XEK) const = 0;
  
  virtual ~abc_grid(){};
};

template<size_t dim>
class abc_genericGrid{
public:
  static const size_t size = dim;
  double EMIN, EL, EU, ET[size], DLEMP[size], DLEMP1, DLFC;

  bool initialized;

  abc_genericGrid() : initialized(false){}
  
  virtual int init(double EMINu, double EMAXu) = 0;
  virtual void getInterval(const double E, long int& KE,
			   double& XEL, double& XE, double& XEK) const = 0;  
  virtual ~abc_genericGrid(){};
};


//  Implementations
//-------------------

template<class stateType, class contextType, class materialType>
void abc_particle<stateType,contextType,materialType>::
vr_particleStack(const unsigned long long nhist,
		 pen_rand& random, const unsigned verbose){

  if(genericVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    genericVR->run_particleStack(nhist,kpar,KDET,
				 state,genericStates,
				 created,spaceAvailable,random);
      
    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Generic VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }
    stateType defaultState;
    for(unsigned i = 0; i < created; ++i){
      defaultState.copyBase(genericStates[i]);
      stack.store(defaultState);
    }
  }
  if(specificVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    specificVR->run_particleStack(nhist,kpar,KDET,
				  state,specificStates,
				  created,spaceAvailable,random);

    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Specific VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }

    for(unsigned i = 0; i < created; ++i){
      stack.store(specificStates[i]);
    }
  }
}

template<class stateType, class contextType, class materialType>
double abc_particle<stateType,contextType,materialType>::
vr_matChange(const unsigned long long nhist,
	     pen_rand& random, const unsigned verbose){

  double de = 0.0;
  if(genericVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    genericVR->run_matChange(nhist,kpar,MATL,
			     state,genericStates,
			     created,spaceAvailable,random);
      
    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Generic VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }
    stateType defaultState;
    for(unsigned i = 0; i < created; ++i){
      defaultState.copyBase(genericStates[i]);
      stack.store(defaultState);
      de += defaultState.E*defaultState.WGHT;
    }
  }
  if(specificVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    specificVR->run_matChange(nhist,kpar,MATL,
			      state,specificStates,
			      created,spaceAvailable,random);

    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Specific VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }

    for(unsigned i = 0; i < created; ++i){
      stack.store(specificStates[i]);
      de += specificStates[i].E*specificStates[i].WGHT;      
    }
  }
  return de;
}

template<class stateType, class contextType, class materialType>
double abc_particle<stateType,contextType,materialType>::
vr_interfCross(const unsigned long long nhist,
	       pen_rand& random, const unsigned verbose){

  double de = 0.0;
  if(genericVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    genericVR->run_interfCross(nhist,kpar,KDET,
			       state,genericStates,
			       created,spaceAvailable,random);

    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Generic VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }
    stateType defaultState;
    for(unsigned i = 0; i < created; ++i){
      defaultState.copyBase(genericStates[i]);
      stack.store(defaultState);
      de += defaultState.E*defaultState.WGHT;      
    }
  }
  if(specificVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    specificVR->run_interfCross(nhist,kpar,KDET,
				state,specificStates,
				created,spaceAvailable,random);

    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Specific VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }

    for(unsigned i = 0; i < created; ++i){
      stack.store(specificStates[i]);
      de += specificStates[i].E*specificStates[i].WGHT;      
    }
  }
  return de;
}

#endif
