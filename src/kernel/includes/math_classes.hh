 
//
//
//    Copyright (C) 2023-2024 Universitat de València - UV
//    Copyright (C) 2023-2024 Universitat Politècnica de València - UPV
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

#ifndef __PEN_MATH_CLASSES__
#define __PEN_MATH_CLASSES__

#include <cmath>
#include <numeric>
#include <functional>
#include <istream>
#include <string>
#include <limits>
#include <vector>
#include <array>
#include <sstream>

//--------------------------------
// Auxiliar structs and classes
//--------------------------------

template<class E, class T> struct container;

template<class T>
struct vector2D{
  T x,y;
    
  constexpr vector2D() {}
  constexpr vector2D(T xin, T yin) : x(xin), y(yin){}
    
  inline constexpr vector2D operator+(const vector2D& a) const{
    return vector2D(x + a.x, y + a.y);
  }
  inline vector2D& operator+=(const vector2D& a){
    x += a.x;
    y += a.y;
    return *this;
  }

  inline constexpr vector2D operator+(const T& a) const{
    return vector2D(x + a, y + a);
  }
  inline vector2D& operator+=(const T& a){
    x += a;
    y += a;
    return *this;
  }
  
  inline constexpr vector2D operator-(const vector2D& a) const{
    return vector2D(x - a.x, y - a.y);
  }
  inline vector2D& operator-=(const vector2D& a){
    x -= a.x;
    y -= a.y;
    return *this;
  }

  inline constexpr vector2D operator-(const T& a) const{
    return vector2D(x - a, y - a);
  }
  inline vector2D& operator-=(const T& a){
    x -= a;
    y -= a;
    return *this;
  }
  
  inline constexpr T operator*(const vector2D& a) const{
    return x*a.x + y*a.y;
  }

  inline constexpr vector2D operator*(const T& a) const{
    return vector2D(a*x, a*y);
  }

  inline constexpr vector2D operator/(const T& a) const{
    return vector2D(x/a, y/a);
  }
  inline vector2D& operator/=(const T& a){
    x /= a;
    y /= a;
    return *this;
  }  
    
  inline void add(const vector2D& a){
    x += a.x;
    y += a.y;
  }

  inline void add(const T& a){
    x += a;
    y += a;
  }
    
  inline void subs(const vector2D& a){
    x -= a.x;
    y -= a.y;
  }

  inline void subs(const T& a){
    x -= a;
    y -= a;
  }

  inline void pow(const unsigned a){
    x = pow(x,a);
    y = pow(y,a);
  }  
    
  inline constexpr T mod2() const{
    return x*x + y*y;
  }
    
  inline constexpr T mod() const{
    return sqrt(mod2());
  }

  inline constexpr T dist(const vector2D& a) const{
    return (a - *this).mod();
  }

  inline constexpr vector2D dir(const vector2D& a) const{
    vector2D normDir = (a - *this);
    normDir.normalize();
    return normDir;
  }
  
  inline void normalize(){
    T norm = mod();
    x /= norm;
    y /= norm;
  }

  inline std::string stringify() const{
    char str[60];
    sprintf(str,"(%.5E,%.5E)",
	    static_cast<double>(x),
	    static_cast<double>(y));
    return std::string(str);
  }
  
};

template<class T>
struct vector3D{
  T x,y,z;
    
  constexpr vector3D() {}
  constexpr vector3D(T xin, T yin, T zin) : x(xin), y(yin), z(zin){}
    
  inline constexpr vector3D operator+(const vector3D& a) const{
    return vector3D(x + a.x, y + a.y, z + a.z);
  }
  inline vector3D& operator+=(const vector3D& a){
    x += a.x;
    y += a.y;
    z += a.z;
    return *this;
  }

  inline constexpr vector3D operator+(const T& a) const{
    return vector3D(x + a, y + a, z + a);
  }
  inline vector3D& operator+=(const T& a){
    x += a;
    y += a;
    z += a;
    return *this;
  }  
  
  inline constexpr vector3D operator-(const vector3D& a) const{
    return vector3D(x - a.x, y - a.y, z - a.z);
  }
  inline vector3D& operator-=(const vector3D& a){
    x -= a.x;
    y -= a.y;
    z -= a.z;
    return *this;
  }  

  inline constexpr vector3D operator-(const T& a) const{
    return vector3D(x - a, y - a, z - a);
  }
  inline vector3D& operator-=(const T& a){
    x -= a;
    y -= a;
    z -= a;
    return *this;
  }
  
  inline constexpr T operator*(const vector3D& a) const{
    return x*a.x + y*a.y + z*a.z;
  }

  inline constexpr vector3D operator*(const T& a) const{
    return vector3D(a*x, a*y, a*z);
  }

  inline constexpr vector3D operator/(const T& a) const{
    return vector3D(x/a, y/a, z/a);
  }
  inline vector3D& operator/=(const T& a){
    x /= a;
    y /= a;
    z /= a;
    return *this;
  }
    
  inline constexpr vector3D operator^(const vector3D& a) const{
    return vector3D(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
  }
    
  inline void add(const vector3D& a){
    x += a.x;
    y += a.y;
    z += a.z;
  }

  inline void add(const T& a){
    x += a;
    y += a;
    z += a;
  }
    
  inline void subs(const vector3D& a){
    x -= a.x;
    y -= a.y;
    z -= a.z;
  }

  inline void subs(const T& a){
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
    
  inline constexpr T mod2() const{
    return x*x + y*y + z*z;
  }
    
  inline constexpr T mod() const{
    return sqrt(mod2());
  }

  inline constexpr T dist(const vector3D& a) const{
    return (a - *this).mod();
  }

  inline constexpr vector3D dir(const vector3D& a) const{
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

//Define a structures to save measures
namespace penred{
  namespace measurements{

      enum errors{
	SUCCESS = 0,
	DIMENSION_NOT_FOUND,
	DIMENSION_MISMATCH,
	DIMENSION_OUT_OF_RANGE,
	DIMENSION_REPEATED,
	EFFECTIVE_DIMENSIONS_MISMATCH,
	NUMBER_OF_BINS_NOT_FOUND,
	MIN_BIN_GREATER_THAN_MAX,
	BIN_OUT_OF_RANGE,
	INVALID_NUMBER_OF_BINS,
	LIMITS_NOT_FOUND,
	INVALID_LIMITS,
	SIGMA_NOT_FOUND,
	DATA_NOT_FOUND,
	CORRUPTED_DATA,
      };

    constexpr const char* errorToString( const unsigned i ){
      switch(i){
      case SUCCESS: return "Success";
      case DIMENSION_NOT_FOUND: return "Dimension not found";
      case DIMENSION_MISMATCH: return "Dimensions mismatch";
      case DIMENSION_OUT_OF_RANGE: return "Dimension is out of range";
      case DIMENSION_REPEATED: return "Dimension repeated";
      case EFFECTIVE_DIMENSIONS_MISMATCH: return "Effective dimensions mismatch";
      case NUMBER_OF_BINS_NOT_FOUND: return "Number of bins not found";
      case MIN_BIN_GREATER_THAN_MAX: return "Minimum bin is greater than maximum one";
      case BIN_OUT_OF_RANGE: return "Bin out of range";
      case INVALID_NUMBER_OF_BINS: return "Invalid number of bins";
      case LIMITS_NOT_FOUND: return "Limits not found";
      case INVALID_LIMITS: return "Invalid limits";
      case SIGMA_NOT_FOUND: return "Sigma not found";
      case DATA_NOT_FOUND: return "Data not found";
      case CORRUPTED_DATA: return "Corrupted data";
      default: return "Unknown error";
      };
    }

	
    constexpr size_t maxDims = 1000;

    typedef std::pair<double, double> limitsType;

    inline std::string triml(const std::string& strin){
      const std::string delimiters = " \n\r\t\f\v";
      size_t first = strin.find_first_not_of(delimiters);
      return (first == std::string::npos) ? "" : strin.substr(first);
    }

    inline std::string trimr(const std::string& strin){
      const std::string delimiters = " \n\r\t\f\v";
      size_t last = strin.find_last_not_of(delimiters);
      return (last == std::string::npos) ? "" : strin.substr(0,last+1);
    }

    inline std::string trim(const std::string& strin){
      return trimr(triml(strin));
    }

    template<size_t dim = 1>
    class multiDimension{

      static_assert(dim <= maxDims,
		    "penred::measurements: "
		    "Maximum number of dimensions exceeded");

    public:
      
    protected:

      // Variables
      
      std::array<unsigned long, dim> nBins;
      unsigned long totalBins;
      std::array<unsigned long, dim> binsPerIncrement;
      std::array<double, dim> binWidths;
      std::array<std::pair<double, double>, dim> limits;

      static constexpr const unsigned nHeaders = dim + 2;
      std::array<std::string, nHeaders> headers; 

      //Initialization functions
      
      template<size_t dimInit>
      inline std::enable_if_t< (dimInit < dim), int>
      initDims(const std::array<unsigned long, dimInit>& nBinsIn,
	       const std::array<std::pair<double, double>, dimInit>& limitsIn){

	//Init auxiliary arrays
	std::array<unsigned long, dim> nBinsAux;
	std::fill(nBinsAux.begin(), nBinsAux.end(), 1ul);

	std::array<std::pair<double, double>, dim> limitsAux;
	const std::pair<double, double> defaultLimit(-1.0e35, 1.0e35);
	std::fill(limitsAux.begin(), limitsAux.end(), defaultLimit);

	for(size_t i = 0; i < dimInit; ++i){
	  nBinsAux[i] = nBinsIn[i];
	}
	for(size_t i = 0; i < dimInit; ++i){
	  limitsAux[i] = limitsIn[i];
	}

	return initDims(nBinsAux, limitsAux);
      }

      template<size_t dimInit>
      inline std::enable_if_t< (dimInit > dim), int>
      initDims(const std::array<unsigned long, dimInit>& nBinsIn,
	       const std::array<std::pair<double, double>, dimInit>& limitsIn) = delete;

      template<size_t dimInit>
      std::enable_if_t< (dimInit == dim), int>
      initDims(const std::array<unsigned long, dimInit>& nBinsIn,
	       const std::array<std::pair<double, double>, dimInit>& limitsIn){

	//Calculate total number of bins
	totalBins = std::accumulate(nBinsIn.begin(),
				    nBinsIn.end(), 1,
				    std::multiplies<unsigned long>());
	if(totalBins == 0)
	  return errors::INVALID_NUMBER_OF_BINS;
    
	//Check limits
	for(size_t i = 0; i < dim; ++i){
	  if(limitsIn[i].first >= limitsIn[i].second){
	    return errors::INVALID_LIMITS;
	  }
	}

	//Save bins
	nBins = nBinsIn;

	//Calculate bins per increment in each dimension
	binsPerIncrement[0] = 1;
	for(size_t i = 1; i < dim; ++i){
	  binsPerIncrement[i] = binsPerIncrement[i-1]*nBins[i-1];
	}

	//Save limits
	limits = limitsIn;

	//Calculate bin widths
	for(size_t i = 0; i < dim; ++i){
	  binWidths[i] = (limits[i].second - limits[i].first)/static_cast<double>(nBins[i]);
	}

	return errors::SUCCESS;
      }

      template<size_t dimInit>
      inline int initDims(const multiDimension<dimInit>& c){
	return initDims(c.nBins, c.limits);
      }

      template<size_t binDims, size_t limitsDims>
      inline int initDims(const unsigned long(&nBinsIn)[binDims],
			  const std::pair<double, double>(&limitsIn)[limitsDims]){

	static_assert(binDims == limitsDims,
		      "Bins and limits dimensions mismatch.");

	//Create arrays
	std::array<unsigned long, binDims> auxBins;
	std::array<std::pair<double, double>, limitsDims> auxLimits;

	for(size_t i = 0; i < binDims; ++i){
	  auxBins[i] = nBinsIn[i];
	}
	for(size_t i = 0; i < limitsDims; ++i){
	  auxLimits[i] = limitsIn[i];
	}
	
	return initDims(auxBins,auxLimits);
	
      }

      inline int initDims(){
	std::array<unsigned long, 1> auxBins;
	std::array<std::pair<double, double>, 1> auxLimits;

	auxBins[0] = 1;
	auxLimits[0] = std::pair<double,double>(-1.0e35, 1.0e35);
	
	return initDims(auxBins, auxLimits);
      }
      
      //Copy functions
      inline void copyDimInfo(const multiDimension<dim>& o){
	*this = o;
      }
      
    public:

      //Variables
      std::string description;

      //Constructors
      inline multiDimension(){
	initHeaders();
	initDims();
      }
      
      //Loop functions
      void forEach(const std::function<void(const unsigned long,
		   const std::array<unsigned long, dim>&)>& f) const {

	//Iterates calling the function "f" for each bin, providing the global
	//bin index and the local index for each dimension
	
	//Create an array with local dimensions indexes
	std::array<unsigned long, dim> indexes;
	std::fill(indexes.begin(), indexes.end(), 0ul);
	
	for(unsigned long i = 0; i < totalBins; ++i){

	  //Call user function
	  f(i, indexes);

	  //Increase first dimension index
	  ++indexes[0];
	  //Check for other dimensions index increment
	  for(size_t idim = 0; idim < dim; ++idim){
	    if(indexes[idim] >= nBins[idim]){
	      
	      indexes[idim] = 0;
	      if(idim < dim-1){
		++indexes[idim+1];
	      }
	    }
	    else{
	      break;
	    }	    
	  }
	}
      }

      int forEach(const std::array<std::pair<unsigned long, unsigned long>, dim>& binLimits,
		  const std::function<void(const unsigned long,
		  const std::array<unsigned long, dim>&)>& f) const {

	//Iterates calling the function "f" for each bin in the specified limits,
	//providing the global bin index and the local index for each dimension

	//Check limits
	for(size_t idim = 0; idim < dim; ++idim){
	  if(binLimits[idim].first >= binLimits[idim].second){
	    return errors::MIN_BIN_GREATER_THAN_MAX;
	  }
	  if(binLimits[idim].second > nBins[idim])
	    return errors::BIN_OUT_OF_RANGE;	  
	}

	//Calculate the number of bins to skip in each dimension when the limit is reached
	std::array<unsigned long, dim> toSkip;
	for(size_t idim = 0; idim < dim; ++idim){
	  toSkip[idim] = (nBins[idim] - binLimits[idim].second) + binLimits[idim].first;
	  toSkip[idim] *= binsPerIncrement[idim];
	}

	//Store bin indexes in each dimension
	std::array<unsigned long, dim> indexes;

	//Init dimension bin to lower bin limits
	size_t globIndex = 0;
	for(size_t idim = 0; idim < dim; ++idim){
	  indexes[idim] = binLimits[idim].first;
	  globIndex += indexes[idim]*binsPerIncrement[idim];
	}

	//Iterate over bins in limits
	while(globIndex < totalBins){

	  //Call user function
	  f(globIndex, indexes);
	  
	  //Increase global index
	  ++globIndex;
	  
	  //Increase first dimension index
	  ++indexes[0];
	  
	  //Check for other dimensions index increment
	  for(size_t idim = 0; idim < dim; ++idim){
	    
	    if(indexes[idim] >= binLimits[idim].second){
	      indexes[idim] = binLimits[idim].first;

	      //Increase global index
	      globIndex += toSkip[idim];
	      
	      if(idim < dim-1){
		++indexes[idim+1];
	      }
	    }
	    else{
	      break;
	    }
	    
	  }
	}

	return errors::SUCCESS;
      }

      inline int forEach(const std::vector<std::array<unsigned long, 3>>& binLimits,
			 const std::function<void(const unsigned long,
			 const std::array<unsigned long, dim>&)>& f) const {

	std::array<std::pair<unsigned long, unsigned long>, dim> auxBinLimits;
	for(size_t i = 0; i < dim; ++i)
	  auxBinLimits[i] = std::pair<unsigned long, unsigned long>(0lu, nBins[i]);

	for(size_t i = 0; i < binLimits.size(); ++i){
	  unsigned long idim = binLimits[i][0];
	  if(idim >= dim){
	    return errors::DIMENSION_OUT_OF_RANGE;
	  }
	  auxBinLimits[idim] =
	    std::pair<unsigned long, unsigned long>(binLimits[i][1], binLimits[i][2]);
	}

	return forEach(auxBinLimits, f);
      }

      //Init functions
      inline void initHeaders(){
	for(size_t i = 0; i < dim; ++i){
	  headers[i] = std::string("Dimension ") + std::to_string(i);
	}
	
	headers[dim] = "Value";
	headers[dim+1] = "Sigma";	
      }

      //Headers functions
      inline const std::string& readDimHeader(const unsigned idim) const {
	if(idim > dim){
	  return "";
	}	
	return headers[idim];
      }
      inline const std::string& readValueHeader() const {
	return headers[dim];
      }
      inline const std::string& readSigmaHeader() const {
	return headers[dim+1];
      }

      inline int setDimHeader(const unsigned idim, const std::string& h){

	if(idim > dim){
	  return errors::DIMENSION_OUT_OF_RANGE;
	}
	
	headers[idim] = h;

	//Replace all '|'
	std::string::size_type pos = headers[idim].find('|');
	while(pos != std::string::npos){
	  headers[idim][pos] = '!';
	  pos = headers[idim].find('|', pos);
	}
	headers[idim] = trim(headers[idim]);

	return errors::SUCCESS;	
      }
      inline void setValueHeader(const std::string& h){
	headers[dim] = h;

	//Replace all '|'
	std::string::size_type pos = headers[dim].find('|');
	while(pos != std::string::npos){
	  headers[dim][pos] = '!';
	  pos = headers[dim].find('|', pos);
	}
	headers[dim] = trim(headers[dim]);
      }
      inline void setSigmaHeader(const std::string& h){
	headers[dim+1] = h;

	//Replace all '|'
	std::string::size_type pos = headers[dim+1].find('|');
	while(pos != std::string::npos){
	  headers[dim+1][pos] = '!';
	  pos = headers[dim+1].find('|', pos);
	}
	headers[dim+1] = trim(headers[dim+1]);
      }

      //Get functions
      const std::array<unsigned long, dim>& readDimBins() const { return nBins; }
      const std::array<std::pair<double, double>, dim>& readLimits() const {
	return limits;
      }
      inline unsigned long getNBins(const unsigned idim) const { return nBins[idim]; }
      inline unsigned long getNBins() const { return totalBins; }
      inline unsigned long effectiveDim() const {
	//Calculate effective dimensions
	unsigned long effectiveDims = dim;
	for(size_t i = 0; i < dim; ++i){
	  if(nBins[i] == 1){
	    --effectiveDims;
	  }
	}
	return effectiveDims;
      }

      //Index functions
      inline unsigned long getGlobalIndex(const std::array<unsigned long, dim>& dimIndexes)const{
	// Calculate the global index from local indexes in each dimension
	unsigned long globIndex = dimIndexes[0];
	for(size_t i = 1; i < dim; ++i){
	  globIndex += dimIndexes[i]*this->binsPerIncrement[i];
	}
	return globIndex;
      }
      
      //Stringify functions
      std::string stringifyInfo() const {

	std::stringstream buf;
	
	for(size_t i = 0; i < dim; ++i){
	  buf << "\n" << headers[i] << ":\n"
	      << "  + Number of bins:\n     "
	      << nBins[i] << "\n"
	      << "  + Limits:\n     [";
	  if(limits[i].first < -1.0e20)
	    buf << "-Inf, ";
	  else
	    buf << limits[i].first << ", ";

	  if(limits[i].second > 1.0e20)
	    buf << " Inf";
	  else
	    buf << limits[i].second;
	    
	  buf << ")" << std::endl;
	  if(nBins[i] > 1){
	    buf << "  + Bin width:\n     "
		<< binWidths[i] << std::endl;
	  }
	}
	return buf.str();
      }

      //Print functions
      void printDims(FILE* fout,
		     const bool printOnlyEffective = false) const{

	//Print description
	if(!description.empty()){
	  fprintf(fout, "#\n");

	  //Get next endLine in description string
	  std::string::size_type endLinePos = description.find('\n');
	  fprintf(fout, "# %s\n", description.substr(0,endLinePos).c_str());
	  while(endLinePos != std::string::npos){
	    std::string::size_type initPos = endLinePos;
	    endLinePos = description.find('\n',initPos+1);
	    fprintf(fout, "# %s\n", description.substr(initPos+1,endLinePos).c_str());
	  }
	}
	
	if(printOnlyEffective){
	  fprintf(fout, "# Dimensions:\n"
		  "# %lu\n"
		  "# Number of bins:\n"
		  "#",
		  effectiveDim());
	}
	else{
	  fprintf(fout, "# Dimensions:\n"
		  "# %lu\n"
		  "# Number of bins:\n"
		  "#",
		  static_cast<unsigned long>(dim));
	}
	
	for(size_t i = 0; i < dim; ++i){
	  if(printOnlyEffective && nBins[i] == 1){
	    continue;
	  }	  
	  fprintf(fout, " %lu", nBins[i]);
	}
	fprintf(fout, "\n"
		"# Limits [low,top): \n");
	for(size_t i = 0; i < dim; ++i){
	  if(printOnlyEffective && nBins[i] == 1){
	    continue;
	  }	  
	  fprintf(fout, "# %15.5E %15.5E  %s\n",
		  limits[i].first,
		  limits[i].second,
		  headers[i].c_str());
	}	
      }

      void printData(FILE* fout,
		     const std::function<
		     std::pair<double,double>(const unsigned long,
		     const std::array<unsigned long, dim>&)
		     > fvalue,
		     const unsigned nSigma,
		     const bool printCoordinates,
		     const bool printBinNumber,
		     const bool printOnlyEffective = false) const {

	//Print description and dimension information
	printDims(fout, printOnlyEffective);
	
	//Print sigma multiplier
	fprintf(fout, "# Printed sigmas:\n");
	fprintf(fout, "# %u\n", nSigma);

	//Calculate spaces required to print each dimension information
	constexpr unsigned long reqSpaceBins = 5;
	constexpr unsigned long reqSpaceCoor = 16;
	const unsigned long reqSpaceNums =
	  (printBinNumber ? reqSpaceBins : 0) + (printCoordinates ? reqSpaceCoor : 0);
	std::array<int, nHeaders> headersSpace;
	if(printBinNumber || printCoordinates){
	  for(size_t iHeader = 0; iHeader < dim; ++iHeader){
	    if(reqSpaceNums < headers[iHeader].length()){
	      headersSpace[iHeader] = headers[iHeader].length();
	    }else{
	      headersSpace[iHeader] = reqSpaceNums;
	    }
	  }
	}

	//Set header space for value and sigma
	headersSpace[dim]   = std::max(15u, static_cast<unsigned>(headers[dim].length()));
	headersSpace[dim+1] = std::max(13u, static_cast<unsigned>(headers[dim+1].length()));

	//Print dimensions headers
	fprintf(fout, "#");
	if(printBinNumber || printCoordinates){
	  for(unsigned long idim = 0; idim < dim; ++idim){
	    if(nBins[idim] > 1){
	      fprintf(fout, " %*s |",
		      headersSpace[idim], headers[idim].c_str());
	    }
	  }
	}
	//Print value and sigma headers
	for(unsigned long iHeader = dim; iHeader < nHeaders; ++iHeader){
	  fprintf(fout, " %*s |",
		  headersSpace[iHeader], headers[iHeader].c_str());	  
	}	  
	
	fprintf(fout, "\n");

	//Create an array to store local dim indexes
	std::array<unsigned long, dim> indexes;
	std::fill(indexes.begin(), indexes.end(), 0lu);

	for(unsigned long i = 0; i < totalBins; ++i){

	  fprintf(fout, " ");
	  //Print dimensions bin and coordinates information
	  for(size_t idim = 0; idim < dim; ++idim){
	    
	    if(nBins[idim] == 1){
	      //Skip dimensions with a single bin
	      continue;
	    }
  
	    //Print bin information
	    if(printBinNumber){
	      if(printCoordinates){
		fprintf(fout, " %4lu", indexes[idim]);
	      }
	      else{
		fprintf(fout, " %*lu  ", headersSpace[idim], indexes[idim]);
	      }
	    }

	    //Print position information
	    if(printCoordinates){
	      if(printBinNumber){
		fprintf(fout, " %*.5E  ",headersSpace[idim]-5,
			limits[idim].first + indexes[idim]*binWidths[idim]);
	      }
	      else{
		fprintf(fout, " %*.5E  ",headersSpace[idim],
			limits[idim].first + indexes[idim]*binWidths[idim]);
	      }
	    }
	  }

	  //Print value information and error
	  std::pair<double, double> val = fvalue(i, indexes);
	
	  fprintf(fout, " %*.5E   %*.2E\n",
		  headersSpace[dim], val.first,
		  headersSpace[dim+1], static_cast<double>(nSigma)*val.second);

	  //Increase first dimension index
	  ++indexes[0];
	  //Check for other dimensions index increment
	  for(size_t idim = 0; idim < dim; ++idim){
	    if(indexes[idim] >= nBins[idim]){
	      
	      if(nBins[idim] > 1){
		//Print a blank line to separate coordinate blocks
		fprintf(fout, "\n");
	      }
	      indexes[idim] = 0;
	      if(idim < dim-1){
		++indexes[idim+1];
	      }
	    }
	    else{
	      break;
	    }
	    
	  }
	}	
      }

      //Read function
      int loadData(std::vector<double>& data,
		   std::vector<double>& sigma,
		   std::istream& in) {

	//Clear data
	initDims();
	
	std::string line;

	// ** Dimensions
	unsigned long readDim;
    
	//Skip lines until "# Dimensions" is found
	bool found = false;
	bool textInDescription = false;
	while(std::getline(in, line)){
	  if(line.find("# Dimensions:") == 0){
	    found = true;
	    break;
	  }
	  line = trim(line);
	  if(line.size() > 0){
	    if(line[0] == '#'){
	      line.erase(0,1);
	      line = trim(line);
	    }
	  }

	  if(line.empty() && !textInDescription){
	    continue;
	  }

	  if(line.empty()){
	    description += '\n';
	  }
	  else{
	    if(textInDescription){
	      description += '\n' + line;
	    }
	    else{
	      description += line;
	      textInDescription = true;
	    }
	  }
	}

	if(!found){
	  return errors::DIMENSION_NOT_FOUND;
	}

	//Read dimension
	std::getline(in, line);
	if(sscanf(line.c_str(), " %*c %lu", &readDim) != 1){
	  return errors::DIMENSION_NOT_FOUND;
	}

	
	// ** Number of bins
	found = false;
	while(std::getline(in, line)){
	  if(line.find("# Number of bins:") == 0){
	    found = true;
	    break;
	  }
	}

	if(!found){
	  return errors::NUMBER_OF_BINS_NOT_FOUND;
	}

	//Remove comment character
	in.ignore(std::numeric_limits<std::streamsize>::max(), '#');
	
	//Read number of bins for each file dimension	
	std::vector<unsigned long> readNBins(readDim);
	for(size_t i = 0; i < readDim; ++i){
	  in >> readNBins[i];
	  if(!in){
	    return errors::NUMBER_OF_BINS_NOT_FOUND;
	  }
	  if(readNBins[i] == 0){
	    return errors::INVALID_NUMBER_OF_BINS;
	  }
	}

	//Count effective dimensions
	unsigned long readEffectiveDim = 0;
	for(size_t i = 0; i < readDim; ++i){
	  if(readNBins[i] > 1){
	    ++readEffectiveDim;
	  }
	}

	// * Check dimensions
	bool readEffective = false;
	if(readDim > dim){
	  if(readEffectiveDim <= dim)
	    readEffective = true;
	  else
	    return errors::DIMENSION_MISMATCH;
	}

	//Set number of bins
	if(readEffective){
	  size_t j = 0;
	  for(size_t i = 0; i < readDim; ++i){
	    if(readNBins[i] > 1){
	      nBins[j++] = readNBins[i];
	    }
	  }
	}
	else{
	  for(size_t i = 0; i < readDim; ++i){
	    nBins[i] = readNBins[i];
	  }
	}

	//Calculate total number of bins
	totalBins = std::accumulate(nBins.begin(),
				    nBins.end(), 1,
				    std::multiplies<unsigned long>());

	if(totalBins == 0)
	  return errors::INVALID_NUMBER_OF_BINS;    

	//Calculate bins per increment in each dimension
	binsPerIncrement[0] = 1;
	for(size_t i = 1; i < dim; ++i){
	  binsPerIncrement[i] = binsPerIncrement[i-1]*nBins[i-1];
	}
    
	//Skip line
	in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	// ** Limits
	found = false;
	while(std::getline(in, line)){
	  if(line.find("# Limits [low,top):") == 0){
	    found = true;
	    break;
	  }
	}

	if(!found){
	  return errors::LIMITS_NOT_FOUND;
	}

	//Read limits
	size_t j = 0;
	for(size_t i = 0; i < readDim; ++i){
	  if(readEffective && readNBins[i] == 1){
	    //Skip line
	    in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	  }
	  else{
	    //Read limits
	    std::getline(in, line);

	    //Create a buffer from this line
	    std::stringstream ssline;
	    ssline.str(line);

	    //Skip comment character
	    ssline.ignore(std::numeric_limits<std::streamsize>::max(), '#');

	    //Read low limit
	    ssline >> limits[j].first;
	    if(!ssline){
	      return errors::LIMITS_NOT_FOUND;
	    }

	    //Read top limit
	    ssline >> limits[j].second;
	    if(!ssline){
	      return errors::LIMITS_NOT_FOUND;
	    }

	    //Read dimension header
	    std::getline(ssline,headers[j]);

	    if(headers[j].length() > 2){
	      //Trim header right spaces
	      headers[j] = trim(headers[j]);
	    }
	    
	    if(limits[j].first >= limits[j].second){
	      return errors::INVALID_LIMITS;
	    }
	    ++j;
	  }
	}

	//Calculate bin widths
	for(size_t i = 0; i < dim; ++i){
	  binWidths[i] = (limits[i].second - limits[i].first)/static_cast<double>(nBins[i]);
	}    

	// ** Sigma
	found = false;
	while(std::getline(in, line)){
	  if(line.find("# Printed sigmas:") == 0){
	    found = true;
	    break;
	  }
	}

	if(!found){
	  return errors::SIGMA_NOT_FOUND;
	}

	//Read printed sigmas
	std::getline(in, line);
	int nSigmasAux;
	if(sscanf(line.c_str(), " %*c %d ", &nSigmasAux) != 1){
	  return errors::SIGMA_NOT_FOUND;
	}	

	double nSigmas = static_cast<double>(std::abs(nSigmasAux));

	// ** Headers

	//Read header line
	std::getline(in, line);

	//Parse value and sigma headers
	std::string::size_type lastBar = line.find_last_of('|');
	if(lastBar != std::string::npos){
	  std::string::size_type lastBar2 = line.find_last_of('|', lastBar-1);
	  if(lastBar2 != std::string::npos){
	    std::string::size_type lastBar3 = line.find_last_of('|', lastBar2-1);
	    if(lastBar3 != std::string::npos){
	      headers[dim] = line.substr(lastBar3+1, lastBar2-lastBar3-1);
	      headers[dim+1] = line.substr(lastBar2+1, lastBar-lastBar2-1);

	      //Trim header left spaces
	      headers[dim] = trim(headers[dim]);
	      headers[dim+1] = trim(headers[dim+1]);
	    }
	  }
	}
	
	// ** Data
	data.resize(totalBins);
	sigma.resize(totalBins);

	//Read the first non comment line (First data bin)
	found = false;
	while(std::getline(in, line)){
	  if(!line.empty()){
	    char firstChar;
	    if(sscanf(line.c_str(), " %c ", &firstChar) == 1){
	      if(firstChar != '#'){
		//First non empty/comment line found
		found = true;
		break;
	      }
	    }
	  }
	}

	if(!found){
	  return errors::DATA_NOT_FOUND;
	}

	//Count the number of numbers in the line
	std::stringstream ssline;
	ssline.str(line);

	size_t nNums = 0;
	while(true){
	  double aux;
	  ssline >> aux;
	  if(!ssline){
	    break;
	  }
	  ++nNums;
	}

	if(nNums < 2){
	  return errors::CORRUPTED_DATA;
	}
	
	size_t toSkip = nNums-2;	
	if(toSkip != 0 &&
	   toSkip != 2*readEffectiveDim &&
	   toSkip != readEffectiveDim){
	  return errors::CORRUPTED_DATA;
	}

	//Save first bin results
	ssline.str("");
	ssline.clear();
	ssline.str(line);
	//Skip coordinate and bin numbers
	for(size_t iskip = 0; iskip < toSkip; ++iskip){
	  double dummy;
	  ssline >> dummy;
	}

	//Read value and sigma
	ssline >> data[0];
	ssline >> sigma[0];
	sigma[0] /= nSigmas;

    
	for(unsigned long i = 1; i < totalBins; ++i){

	  //Skip coordinate and bin numbers
	  for(size_t iskip = 0; iskip < toSkip; ++iskip){
	    double dummy;
	    in >> dummy;
	  }

	  //Read value information and error
	  in >> data[i];
	  if(!in){
	    return errors::CORRUPTED_DATA;
	  }
	  in >> sigma[i];
	  if(!in){
	    return errors::CORRUPTED_DATA;
	  }
	  sigma[i] /= nSigmas;
	}

	return errors::SUCCESS;
      }
      
    };
    
    template<class type, size_t dim = 1>
    class results : public multiDimension<dim>{
      
    public:
      
      std::vector<type> data;
      std::vector<double> sigma;

      //Constructors
      inline results(){
	clear();
      }

      //Get functions
      const std::vector<type>& readData() const { return data; }
      const std::vector<double>& readSigam() const { return sigma; }

      template<size_t dimInit>
      int init(const std::array<unsigned long, dimInit>& nBinsIn,
	       const std::array<std::pair<double, double>, dimInit>& limitsIn){
	
	//Init dimensions
	this->initDims(nBinsIn, limitsIn);
    
	//Resize vectors
	data.resize(this->totalBins);
	std::fill(data.begin(), data.end(), static_cast<type>(0));

	sigma.resize(this->totalBins);
	std::fill(sigma.begin(), sigma.end(), 0.0);

	return errors::SUCCESS;    
      }

      template<size_t dimInit>
      inline int init(const results<type, dimInit>& c){
	return init(c.nBins, c.limits);
      }

      template<size_t binDims, size_t limitsDims>
      inline int init(const unsigned long(&nBinsIn)[binDims],
		      const std::pair<double, double>(&limitsIn)[limitsDims]){

	static_assert(binDims == limitsDims,
		      "Bins and limits dimensions mismatch.");

	//Create arrays
	std::array<unsigned long, binDims> auxBins;
	std::array<std::pair<double, double>, limitsDims> auxLimits;

	for(size_t i = 0; i < binDims; ++i){
	  auxBins[i] = nBinsIn[i];
	}
	for(size_t i = 0; i < limitsDims; ++i){
	  auxLimits[i] = limitsIn[i];
	}
	
	return init(auxBins,auxLimits);
	
      }

      inline int init(){
	std::array<unsigned long, 1> auxBins;
	std::array<std::pair<double, double>, 1> auxLimits;

	auxBins[0] = 1;
	auxLimits[0] = std::pair<double,double>(-1.0e35, 1.0e35);
	
	return init(auxBins, auxLimits);
      }
      
      //Clear function
      inline void clear(){
	init();	
	this->initHeaders();
      }


      //Profile functions
      template<size_t profDim>
      int profileByBins(const std::array<size_t, profDim>& profDimsIndex,
			const std::array<std::pair<unsigned long, unsigned long>, dim>& binLimits,
			results<type, profDim>& profile) const {

	//Create a profile at the dimension "ibin" using
	//the bin range [minBin, maxBin) in that dimension

	//Check profile dimension
	if(profDim >= dim){
	  return errors::DIMENSION_OUT_OF_RANGE;
	}
	for(size_t i = 0; i < profDim; ++i){
	  if(profDimsIndex[i] >= dim)
	    return errors::DIMENSION_OUT_OF_RANGE;
	}
	for(size_t i = 0; i < profDim-1; ++i){
	  for(size_t j = i+1; j < profDim; ++j){
	    if(profDimsIndex[i] == profDimsIndex[j])
	    return errors::DIMENSION_REPEATED;
	  }
	}

	//Calculate profile bins
	std::array<unsigned long, profDim> profileBins;
	for(size_t i = 0; i < profDim; ++i){
	  profileBins[i] =
	    binLimits[profDimsIndex[i]].second - binLimits[profDimsIndex[i]].first;
	}

	//Calculate profile limits
	std::array<std::pair<double,double>, profDim> profileLimits;
	for(size_t i = 0; i < profDim; ++i){
	  size_t profIndex = profDimsIndex[i];
	  profileLimits[i].first = this->limits[profIndex].first + 
	    static_cast<double>(binLimits[profIndex].first)*this->binWidths[profIndex];
	  profileLimits[i].second = this->limits[profIndex].first + 
	    static_cast<double>(binLimits[profIndex].second)*this->binWidths[profIndex];
	}
	
	//Init profile
	profile.init(profileBins,
		     profileLimits);

	//Set headers
	for(size_t i = 0; i < profDim; ++i){
	  profile.setDimHeader(i,this->headers[profDimsIndex[i]]);
	}
	profile.setValueHeader(this->headers[dim]);
	profile.setSigmaHeader(this->headers[dim+1]);

	int err = this->
	  forEach(binLimits,
		  [this, &profDimsIndex, &binLimits, &profile]
		  (const unsigned long globIndex,
		   const std::array<unsigned long, dim>& indexes) -> void{

		    //Calculate profile local indexes
		    std::array<unsigned long, profDim> profIndexes;
		    for(size_t i = 0; i < profDim; ++i){
		      profIndexes[i] =
			indexes[profDimsIndex[i]] - binLimits[profDimsIndex[i]].first;
		    }
		    
		    //Calculate profile global indexes
		    const unsigned long iprof = profile.getGlobalIndex(profIndexes);

		    //Add contribution
		    profile.data[iprof]  += data[globIndex];
		    profile.sigma[iprof] += sigma[globIndex]*sigma[globIndex];    
			    
		  });
	if(err != errors::SUCCESS)
	  return err;

	//Obtain final sigma for each bin
	for(unsigned long iprof = 0; iprof < profile.getNBins(); ++iprof){
	  profile.sigma[iprof] = sqrt(profile.sigma[iprof]);
	}
	
	return errors::SUCCESS;
      }

      template<size_t profDim>
      int profileByBins(const std::array<size_t, profDim>& profDimsIndex,
			const std::vector<std::array<unsigned long, 3>>& binLimits,
			results<type, profDim>& profile) const {

	std::array<std::pair<unsigned long, unsigned long>, dim> auxBinLimits;
	for(size_t i = 0; i < dim; ++i)
	  auxBinLimits[i] = std::pair<unsigned long, unsigned long>(0lu, this->nBins[i]);

	for(size_t i = 0; i < binLimits.size(); ++i){
	  unsigned long idim = binLimits[i][0];
	  if(idim >= dim){
	    return errors::DIMENSION_OUT_OF_RANGE;
	  }
	  auxBinLimits[idim] =
	    std::pair<unsigned long, unsigned long>(binLimits[i][1], binLimits[i][2]);
	}

	return profileByBins(profDimsIndex, auxBinLimits, profile);
      }

      inline int profile1D(const unsigned long profDim,
			   const std::array<std::pair<unsigned long, unsigned long>, dim>& binLimits,
			   results<type, 1>& profile) const{
	std::array<size_t, 1> aux;
	aux[0] = profDim;
	return profileByBins(aux, binLimits, profile);
      }

      inline int profile1D(const unsigned long profDim,
			   const std::vector<std::array<unsigned long, 3>>& binLimits,
			   results<type, 1>& profile) const{
	std::array<size_t, 1> aux;
	aux[0] = profDim;
	return profileByBins(aux, binLimits, profile);
      }
      
      //Print functions
      inline void print(FILE* fout,
			const unsigned nSigma,
			const bool printCoordinates,
			const bool printBinNumber,
			const bool printOnlyEffective = false) const {

	//Print description and dimension information
	this->
	  printData(fout,
		    [this](const unsigned long i, //Function to calculate value and sigma
			   const std::array<unsigned long, dim>&) -> std::pair<double,double> {
		      return {data[i], sigma[i]};
		    },
		    nSigma,
		    printCoordinates,
		    printBinNumber,
		    printOnlyEffective);
	
      }

      //Read functions
      inline int read(std::istream& in){

	//Clear previous data
	clear();

	//Read data
	return this->loadData(data, sigma, in);
      }
      
    };

    template<class type, size_t dim = 1>
    class measurement : public multiDimension<dim>{

      static_assert(dim <= maxDims,
		    "penred::measurements: "
		    "Maximum number of dimensions exceeded");
  
    private:

      std::vector<type> data;
      std::vector<type> data2;      
      std::vector<type> tmp;
      std::vector<unsigned long long> lastHist;
      
    public:
      
      static constexpr size_t dimensions = dim;

      //Constructors
      inline measurement(){
	clear();
      };
      
      //Get functions
      inline const std::vector<type>& readData() const { return data; }
      inline const std::vector<type>& readData2() const { return data2; }
      
      inline std::vector<type>& getData() { return data; }
      inline std::vector<type>& getData2() { return data2; }

      //Init functions
      template<size_t dimInit>
      int init(const std::array<unsigned long, dimInit>& nBinsIn,
	       const std::array<std::pair<double, double>, dimInit>& limitsIn){
	
	//Init dimensions
	this->initDims(nBinsIn, limitsIn);
    
	//Resize vectors
	data.resize(this->totalBins);
	std::fill(data.begin(), data.end(), static_cast<type>(0));

	data2.resize(this->totalBins);
	std::fill(data2.begin(), data2.end(), static_cast<type>(0));

	tmp.resize(this->totalBins);
	std::fill(tmp.begin(), tmp.end(), static_cast<type>(0));

	lastHist.resize(this->totalBins);
	std::fill(lastHist.begin(), lastHist.end(), 0ull);


	return errors::SUCCESS;    
      }

      template<size_t dimInit>
      inline int init(const measurement<type, dimInit>& c){
	return init(c.nBins, c.limits);
      }

      template<size_t binDims, size_t limitsDims>
      inline int init(const unsigned long(&nBinsIn)[binDims],
		      const std::pair<double, double>(&limitsIn)[limitsDims]){

	static_assert(binDims == limitsDims,
		      "Bins and limits dimensions mismatch.");

	//Create arrays
	std::array<unsigned long, binDims> auxBins;
	std::array<std::pair<double, double>, limitsDims> auxLimits;

	for(size_t i = 0; i < binDims; ++i){
	  auxBins[i] = nBinsIn[i];
	}
	for(size_t i = 0; i < limitsDims; ++i){
	  auxLimits[i] = limitsIn[i];
	}
	
	return init(auxBins,auxLimits);
	
      }

      inline int init(){
	std::array<unsigned long, 1> auxBins;
	std::array<std::pair<double, double>, 1> auxLimits;

	auxBins[0] = 1;
	auxLimits[0] = std::pair<double,double>(-1.0e35, 1.0e35);
	
	return init(auxBins, auxLimits);
      }

      //Clear function
      inline void clear(){
	init();	
	this->initHeaders();
      }

      //Measurement functions
      void add(const std::array<double, dim>& pos,
	       const type& value,
	       const unsigned long long hist){

	//Check limits
	for(size_t i = 0; i < dim; ++i){
	  if(pos[i] < this->limits[i].first || pos[i] >= this->limits[i].second){
	    return;
	  }
	}

	//Get indexes
	std::array<unsigned long, dim> index;
	for(size_t i = 0; i < dim; ++i){
	  index[i] = (pos[i] - this->limits[i].first)/this->binWidths[i];
	}


	//Get global index
	unsigned long globIndex = this->getGlobalIndex(index);
      
	if(hist > lastHist[globIndex]){

	  //Update counters
	  data[globIndex] += tmp[globIndex];
	  data2[globIndex] += tmp[globIndex]*tmp[globIndex];

	  //Restart tmp counter
	  tmp[globIndex] = value;

	  //Update last hist
	  lastHist[globIndex] = hist;
	}else{
	  tmp[globIndex] += value;
	}
      }

      int add(measurement<type,dim> toAdd){

	//Check number of bins
	for(size_t i = 0; i < dim; ++i){
	  if(this->nBins[i] != toAdd.nBins[i])
	    return -1;
	}
	for(size_t i = 0; i < this->totalBins; ++i){
	  data[i]  += toAdd.data[i];
	  data2[i] += toAdd.data2[i];
	}

	return 0;
      }
  
      void flush(){
	for(unsigned long i = 0; i < this->totalBins; ++i){
	  //Skip empty bins
	  if(lastHist[i] == 0){continue;}

	  //Update counters
	  data[i] += tmp[i];
	  data2[i] += tmp[i]*tmp[i];

	  //Reset tmp counter
	  tmp[i] = static_cast<type>(0);

	  //Reset last history to avoid recounting
	  lastHist[i] = 0;      
	}
      }

      void results(const unsigned long long nhists, results<type, dim>& res) const {

	res.init(this->nBins, this->limits);

	//Set headers
	for(size_t idim = 0; idim < dim; ++idim){
	  res.setDimHeader(idim, this->headers[idim]);
	}
	res.setValueHeader(this->headers[dim]);
	res.setSigmaHeader(this->headers[dim+1]);
	
	//Calculate normalization factor
	const double factor = 1.0/static_cast<double>(nhists);
    
	for(unsigned long i = 0; i < this->totalBins; ++i){

	  //Calculate resulting value and error
	  const double q = data[i]*factor;
	  double sigma = data2[i]*factor - q*q;
	  if(sigma > 0.0){
	    sigma = sqrt(sigma*factor);
	  }
	  else{
	    sigma = 0.0;
	  }

	  res.data[i] = q;
	  res.sigma[i] = sigma;
	}
      }

      //Print functions
      inline void print(FILE* fout,
			const unsigned long long nhists,
			const unsigned nSigma,
			const bool printCoordinates,
			const bool printBinNumber,
			const bool printOnlyEffective = false) const {

	//Calculate normalization factor
	const double factor = 1.0/static_cast<double>(nhists);

	//Print description and dimension information
	this->
	  printData(fout,
		    [this, factor] //Function to calculate value and sigma
		    (const unsigned long i, //Global index
		     const std::array<unsigned long, dim>&) -> std::pair<double,double> {

		      //Print value information and error
		      const double q = data[i]*factor;
		      double sigma = data2[i]*factor - q*q;
		      if(sigma > 0.0){
			sigma = sqrt(sigma*factor);
		      }
		      else{
			sigma = 0.0;
		      }

		      return {q, sigma};
		    },
		    nSigma,
		    printCoordinates,
		    printBinNumber,
		    printOnlyEffective);
	
      }

    };
    
  }; //namespace measurements
}; //namespace penred

#endif
