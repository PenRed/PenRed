
//
//
//    Copyright (C) 2020-2022 Universitat de València - UV
//    Copyright (C) 2020-2022 Universitat Politècnica de València - UPV
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

 
#ifndef __PEN_KERMA_TRACK_LENGTH_TALLY__
#define __PEN_KERMA_TRACK_LENGTH_TALLY__

#include "pen_constants.hh"
#include "pen_grids.hh"

namespace pen_tally_KTL{

  const double inf = 1.0e35;
  
  struct vect3d{
    double x,y,z;
    vect3d() : x(0.0),y(0.0),z(0.0){}          
    vect3d(const double a,const double b,const double c) : x(a),
							   y(b),
							   z(c){}      
  };
  struct vect3i{
    long int x,y,z;
    vect3i() : x(0),y(0),z(0){}          
    vect3i(const long int a,const long int b,const long int c) : x(a),
								 y(b),
								 z(c){}
  };

  inline double closestAlpha(const double alpha1, const double alpha2){
    //Take the closest alpha
    if(alpha1 >= 0.0 && alpha1 <= 1.0)
      if(alpha2 >= 0.0 && alpha2 <= 1.0)
	return std::min(alpha1,alpha2);
      else
	return alpha1;
    else if(alpha2 >= 0.0 && alpha2 <= 1.0)
      return alpha2;
    else
      return -1.0e35; //Both alpha outside the path    
  }
  inline double farestAlpha(const double alpha1, const double alpha2){
    //Take the closest alpha
    if(alpha1 >= 0.0 && alpha1 <= 1.0)
      if(alpha2 >= 0.0 && alpha2 <= 1.0)
	return std::max(alpha1,alpha2);
      else
	return alpha1;
    else if(alpha2 >= 0.0 && alpha2 <= 1.0)
      return alpha2;
    else
      return -1.0e35; //Both alpha outside the path    
  }
  
  inline void move(const double ds, const vect3d& dir,
		   vect3d& pos){
    pos.x += ds*dir.x;
    pos.y += ds*dir.y;
    pos.z += ds*dir.z;		        
  }

  inline bool moveIn(double pos,
		     const double dir,
		     double& ds2in,
		     double& ds2out,
		     const double max){

    //Check if a particle at position 'pos' and direction 'dir'
    //is in the geometry limits or will enter to it.  

    // This function asumes that geometry is in the interval [0,max)
    // Also, stores in 'ds' the distance until geometry is reached

    const double eps = 1.0e-8;

    if(std::signbit(pos) || pos >= max){
      // position component is out of geometry,
      // check if is moving to it
      if(dir == 0.0){
	//Particle is not moving on this axy, so never will reaches the mesh
	ds2in = ds2out = inf;	 
	return false;
      }

      if(std::signbit(dir) == std::signbit(pos)){
	//Particle is moving away from geometry
	ds2in = ds2out = inf;	 
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
    if(std::signbit(dir)){
      //Negative direction
      ds2out = pos/fabs(dir)-eps;
    }
    else{
      //Positive direction
      ds2out = (max-pos)/dir-eps;
    }
    return true; //Particle is in geometry limits
  }

  inline bool moveIn(vect3d& pos,
		     const vect3d& dir,
		     double& ds,
		     const vect3d& sizes){
    
    double ds2inx,ds2iny,ds2inz;
    double ds2outx,ds2outy,ds2outz;
    if(!moveIn(pos.x, dir.x, ds2inx, ds2outx, sizes.x) ||
       !moveIn(pos.y, dir.y, ds2iny, ds2outy, sizes.y) ||
       !moveIn(pos.z, dir.z, ds2inz, ds2outz, sizes.z)){
      return false;
    }

    double ds2allIn = std::max(std::max(ds2inx,ds2iny),ds2inz);
    if(ds2allIn > 0.0){
      double ds2someOut = std::min(std::min(ds2outx,ds2outy),ds2outz);
      if(ds2allIn >= ds2someOut){
	return false; //The particle dosn't reaches the mesh
      }
      else{
	move(ds2allIn,dir,pos);
      }
    }
    ds = ds2allIn;
    return true;
  }

  inline bool crossVox(double& ds,
		       const long unsigned nvoxAxis,
		       const long int dvox1D,
		       const long int dvox3D,
		       double& remaining,
		       long int& ivox,
		       long int& iaxis) {

    //
    // ds -> distance to next wall
    // imat -> initial material index
    // nvoxAxis -> number of voxels in the axis that will be crossed
    // dvox1D -> Index increment in crossed axis (-1 or +1 for backward or forward)
    // dvox3D -> Global index increment (+-1 for X, +- nx for Y and +- nx*ny for Z)
    // remaining -> Remaining distance to be traveled by the particle
    // ivox -> Actual global (3D) voxel index
    // iaxis -> Actual voxel index in crossed axis (1D) (ix,iy,iz)
  
    if(ds > 0.0E0){
      if(ds >= remaining){  // Maximum traveled distance reached
	ds = remaining;
	remaining = 0.0;
	return true;
      }
      // Update remaining traveling distance
      remaining -= ds;	
    }

    // Update axis voxel index
    iaxis += dvox1D;
    if(iaxis < 0 || iaxis >= (long int)nvoxAxis){ // Particle scapes from mesh	
      return true;
    }
    // Update global voxel index
    ivox += dvox3D;
    // The travel continues
    return false;
  }
  
  void cart2CylSph(const vect3d& pos,
		   vect3d& cyl,
		   vect3d& sph);

  inline void scoreTally(const double wght,
			 const double E,
			 const double muen,
			 const double l,
			 const long int ibin,
			 const unsigned long long nhist, 
			 double* score,
			 double* score2,
			 double* scoreTmp,
			 unsigned long long* lastHist){
    //Add travel to tallies
    double toScore = wght*E*muen*l;
    if(nhist > lastHist[ibin]){
      const double tmp = scoreTmp[ibin];
      score[ibin]  += tmp;
      score2[ibin] += tmp*tmp;
      scoreTmp[ibin] = toScore;
      lastHist[ibin] = nhist;
    }
    else{
      scoreTmp[ibin] += toScore;
    }
  }

}

class pen_tallyKermaTrackLength: public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_tallyKermaTrackLength,pen_particleState)

  public:
  static const size_t nbinmax = 1000;
  static const long int meshAxeMax = 1000;
  static const double PI05;
  
  private:

  pen_KPAR kparTrig;
  pen_tally_KTL::vect3d lastPos;

  pen_genericLogGrid<nbinmax> grid;
  double muen[constants::MAXMAT][nbinmax];
  bool activeMat[constants::MAXMAT];

  //Cartesian
  pen_tally_KTL::vect3i nbinsCart;
  pen_tally_KTL::vect3d minsCart;
  pen_tally_KTL::vect3d maxsCart;
  pen_tally_KTL::vect3d meshSizeCart;  
  pen_tally_KTL::vect3d dbinCart;

  double xPlanes[meshAxeMax+1];
  double yPlanes[meshAxeMax+1];
  double zPlanes[meshAxeMax+1];
  
  double* cartesian;
  double* cartesian2;
  double* cartesianTmp;
  unsigned long long* cartesianLastHist;
  double volumeCart;
  bool activeCart;
  
  //Cylindrical
  pen_tally_KTL::vect3i nbinsCyl;
  unsigned long nrphiCyl;
  double radCyl,radCylmin;
  double zminCyl,zmaxCyl, heightCyl;
  pen_tally_KTL::vect3d dbinCyl;

  double rPlanesCyl[meshAxeMax+1];
  double phiPlanesCyl[meshAxeMax+1];
  double zPlanesCyl[meshAxeMax+1];
  
  double* cylindrical;
  double* cylindrical2;
  double* cylindricalTmp;
  unsigned long long* cylindricalLastHist;
  double pvolumesCyl[meshAxeMax];
  bool activeCyl;
  
  //Spherical
  pen_tally_KTL::vect3i nbinsSph;
  long int halfnThetaCones; //Number of cones in one Z hemisphere
  long int planeThetaConeIndex; //Cone index whith z=0
  unsigned long nrthetaSph;
  double radSph,radSphmin;
  pen_tally_KTL::vect3d dbinSph;

  double rPlanesSph[meshAxeMax+1];
  double thetaPlanesSph[meshAxeMax+1];
  double phiPlanesSph[meshAxeMax+1];
  
  double* spherical;
  double* spherical2;
  double* sphericalTmp;
  unsigned long long* sphericalLastHist;
  double* pvolumesSph;
  bool activeSph;

public:
    
  pen_tallyKermaTrackLength() : pen_genericTally(USE_JUMP | USE_STEP),
				kparTrig(PEN_PHOTON),
				cartesian(nullptr),
				cartesian2(nullptr),
				cartesianTmp(nullptr),
				cartesianLastHist(nullptr),
				activeCart(false),
				cylindrical(nullptr),
				cylindrical2(nullptr),
				cylindricalTmp(nullptr),
				cylindricalLastHist(nullptr),
				activeCyl(false),
				spherical(nullptr),
				spherical2(nullptr),
				sphericalTmp(nullptr),
				sphericalLastHist(nullptr),
				pvolumesSph(nullptr),
				activeSph(false)
  {}
    
  void kermaTrackLengthCart(const unsigned long long nhist,
			    const double wght,
			    const double E,
			    const double muenVal,
			    const double dsef,
			    const pen_tally_KTL::vect3d p1,
			    const pen_tally_KTL::vect3d dir);

  void kermaTrackLengthCyl(const unsigned long long nhist,
			   const double wght,
			   const double E,
			   const double muenVal,
			   const double dsef,
			   const pen_tally_KTL::vect3d p1,
			   const pen_tally_KTL::vect3d dp,
			   const pen_tally_KTL::vect3d p1cyl,
			   const pen_tally_KTL::vect3d p2cyl);

  void kermaTrackLengthSph(const unsigned long long nhist,
			   const double wght,
			   const double E,
			   const double muenVal,
			   const double dsef,
			   const pen_tally_KTL::vect3d p1,
			   const pen_tally_KTL::vect3d dp,
			   const pen_tally_KTL::vect3d p1sph,
			   const pen_tally_KTL::vect3d p2sph);
  
  void tally_step(const unsigned long long nhist,
		  const pen_KPAR kpar,
		  const pen_particleState& state,
		  const tally_StepData& stepData);  

  void tally_jump(const unsigned long long /*nhist*/,
		const pen_KPAR /*kpar*/,
		const pen_particleState& state,
		const double ds);
  
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  void flush();
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_tallyKermaTrackLength& tally);

  inline const double* readCartesians() const { return cartesian; }
  inline const double* readCartesians2() const { return cartesian2; }
  inline bool enabledCart() const { return activeCart; }
  
  inline const double* readCyl() const { return cylindrical; }
  inline const double* readCyl2() const { return cylindrical2; }
  inline bool enabledCyl() const { return activeCyl; }
  
  inline const double* readSph() const { return spherical; }
  inline const double* readSph2() const { return spherical2; }
  inline bool enabledSph() const { return activeSph; }
  
  
  
  ~pen_tallyKermaTrackLength();
};

#endif
