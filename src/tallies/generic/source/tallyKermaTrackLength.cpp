
//
//
//    Copyright (C) 2020 Universitat de València - UV
//    Copyright (C) 2020 Universitat Politècnica de València - UPV
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


#include "tallyKermaTrackLength.hh"

#include <array>

inline void pen_tally_KTL::cart2CylSph(const vect3d& pos,
				       vect3d& cyl,
				       vect3d& sph){
  
  //Convert position to cylindrical and spherical coordinates
  vect3d pos2;
  pos2.x = pos.x*pos.x;
  pos2.y = pos.y*pos.y;
  pos2.z = pos.z*pos.z;
  double rho2 = pos2.x+pos2.y;
  double rho  = sqrt(rho2);
  double rad  = sqrt(rho2+pos2.z);
  double azim;
  double polar;
  if(fabs(pos.x) > 1.0e-10 || fabs(pos.y) > 1.0e-10){
    //Azimutal
    azim = atan2(pos.y,pos.x);
    if(std::signbit(azim))
      azim += 2.0*M_PI;
    //Polar
    polar = acos(pos.z/rad);
  }
  else{
    //Azimutal
    azim = 0.0;
    //Polar
    if(fabs(pos.z) < 1.0e-10)
      polar = 0.0;
    else
      polar = acos(pos.z/rad);
  }

  cyl.x = rho; cyl.y = azim ; cyl.z = pos.z;
  sph.x = rad; sph.y = polar; sph.z = azim ;
}

//Init static variables
const double pen_tallyKermaTrackLength::PI05 = M_PI/2.0;

inline void pen_tallyKermaTrackLength::kermaTrackLengthCart(const unsigned long long nhist,
							    const double wght,
							    const double E,
							    const double muenVal,
							    const double dsef,
							    const pen_tally_KTL::vect3d p1,
							    const pen_tally_KTL::vect3d dir){

  using namespace pen_tally_KTL;
  //Change the origin to mesh origin 
  vect3d in(p1.x-minsCart.x,p1.y-minsCart.y,p1.z-minsCart.z);
  //Try to move init position to the mesh
  double ds = 0.0;
  if(!moveIn(in,dir,ds,meshSizeCart)){
    // The particle doesn't reaches the mesh, leave
    return;
  }

  if(dsef <= ds){ //Particle has not reached the mesh
    return;
  }
  
  //Substract the distance traveled
  double toTravel = dsef-ds;

  //Get initial indexes
  vect3i index(in.x/dbinCart.x,in.y/dbinCart.y,in.z/dbinCart.z);
  
  long int ibin = nbinsCart.x*(index.z*nbinsCart.y + index.y) + index.x;
  vect3d dsVox;
  
  //Calculate the required distance to travel a single voxel on each direction
  vect3d idir;
  //X
  if(dir.x != 0.0E0){
    idir.x = 1.0/dir.x;
    dsVox.x = fabs(idir.x*dbinCart.x);
  }
  else{
    idir.x = dsVox.x = inf;
  }
  //Y
  if(dir.y != 0.0E0){
    idir.y = 1.0/dir.y;
    dsVox.y = fabs(idir.y*dbinCart.y);
  }
  else{
    idir.y = dsVox.y = inf;
  }
  //Z
  if(dir.z != 0.0E0){
    idir.z = 1.0/dir.z;
    dsVox.z = fabs(idir.z*dbinCart.z);
  }
  else{
    idir.z = dsVox.z = inf;
  }

  // Calculate the distance to the nex three voxel walls
  vect3d ds2vox;
  vect3i voxInc, voxIncGlob;
  long int nxy = nbinsCart.x*nbinsCart.y;
  //X
  if(std::signbit(idir.x)){
    // Moves backward on X
    ds2vox.x = (static_cast<double>(index.x)*dbinCart.x-in.x)*idir.x;
    voxInc.x = -1;
    voxIncGlob.x = -1;    
  }
  else{
    //Moves forward on X
    ds2vox.x = (static_cast<double>(index.x+1)*dbinCart.x-in.x)*idir.x;
    voxInc.x = +1;
    voxIncGlob.x = +1;    
  }
  //Y
  if(std::signbit(idir.y)){
    // Moves backward on Y
    ds2vox.y = (static_cast<double>(index.y)*dbinCart.y-in.y)*idir.y;
    voxInc.y = -1;    
    voxIncGlob.y = -nbinsCart.x;    
  }
  else{
    //Moves forward on Y
    ds2vox.y = (static_cast<double>(index.y+1)*dbinCart.y-in.y)*idir.y;
    voxInc.y = +1;
    voxIncGlob.y = +nbinsCart.x;    
  }
  //Z
  if(std::signbit(idir.z)){
    // Moves backward on Z
    ds2vox.z = (static_cast<double>(index.z)*dbinCart.z-in.z)*idir.z;
    voxInc.z = -1;    
    voxIncGlob.z = -nxy;    
  }
  else{
    //Moves forward on Z
    ds2vox.z = (static_cast<double>(index.z+1)*dbinCart.z-in.z)*idir.z;
    voxInc.z = +1;
    voxIncGlob.z = +nxy;
  }
  
  double remaining = toTravel;
  bool traveling = true;
  do{
    
    if(ds2vox.x < ds2vox.y && ds2vox.x < ds2vox.z){
      double l = std::max(ds2vox.x,0.0);
      const long int prevIbin = ibin;
      traveling = !crossVox(l,nbinsCart.x,voxInc.x,voxIncGlob.x,
			    remaining,ibin,index.x);
      scoreTally(wght,E,muenVal,l,prevIbin,nhist,
		 cartesian,cartesian2,cartesianTmp,cartesianLastHist);

      ds2vox.x  = dsVox.x;
      ds2vox.y -= l;
      ds2vox.z -= l;
    }
    else if(ds2vox.y < ds2vox.z){
      double l = std::max(ds2vox.y,0.0);
      const long int prevIbin = ibin;
      traveling = !crossVox(l,nbinsCart.y,voxInc.y,voxIncGlob.y,
			    remaining,ibin,index.y);
      scoreTally(wght,E,muenVal,l,prevIbin,nhist,
		 cartesian,cartesian2,cartesianTmp,cartesianLastHist);
      
      ds2vox.x -= l;
      ds2vox.y  = dsVox.y;
      ds2vox.z -= l;      
    }
    else{
      double l = std::max(ds2vox.z,0.0);
      const long int prevIbin = ibin;
      traveling = !crossVox(l,nbinsCart.z,voxInc.z,voxIncGlob.z,
			    remaining,ibin,index.z);
      scoreTally(wght,E,muenVal,l,prevIbin,nhist,
		 cartesian,cartesian2,cartesianTmp,cartesianLastHist);
      
      ds2vox.x -= l;
      ds2vox.y -= l;
      ds2vox.z  = dsVox.z;
    }
    
  }while(traveling);
  
}

void pen_tallyKermaTrackLength::kermaTrackLengthCyl(const unsigned long long nhist,
						    const double wght,
						    const double E,
						    const double muenVal,
						    const double dsef,
						    const pen_tally_KTL::vect3d p1,
						    const pen_tally_KTL::vect3d dp,
						    const pen_tally_KTL::vect3d p1cyl,
						    const pen_tally_KTL::vect3d p2cyl){

  //
  //   x = p1.x + alpha*dp.x
  //   y = p1.y + alpha*dp.y
  //   z = p1.z + alpha*dp.z
  //
  // rho = sqrt(x*x+y*y)
  // phi = atan2(y,x)
  //   z = z
  
  using namespace pen_tally_KTL;

  //Check if both points are inside the minimus radius
  if(p1cyl.x <= radCylmin && p2cyl.x <= radCylmin)
    return; //Nothing to score
  
  //Precalculate some values
  double dx2 = dp.x*dp.x;
  double dy2 = dp.y*dp.y;

  double drho2 = dx2 + dy2;
  double rho02 = p1.x*p1.x + p1.y*p1.y;

  double Brho  = 2.0*(p1.x*dp.x + p1.y*dp.y);
  double Brho2 = Brho*Brho;

  // Create auxiliary lambda functions to convert from coordinates
  // to alpha and from alpha to coordinates
  auto rho2alpha = [drho2,Brho,Brho2,rho02]
    (const double rho, double& alphaP, double& alphaM)
    -> void{
		     
		     if(drho2 < 1.0e-10){
		       alphaP = alphaM = -1.0e35;
		       return;
		     }
		     
		     const double Ax2 = 2.0*drho2;
		     const double Ax4 = 4.0*drho2;
		     double C = rho02 - rho*rho;
		     double sqrtArg = Brho2 - Ax4*C;
		     if(std::signbit(sqrtArg)){
		       alphaP = alphaM = -1.0e35;
		       return;
		     }
		     double sqrtRes = sqrt(sqrtArg);
		     alphaP = (-Brho+sqrtRes)/Ax2;
		     alphaM = (-Brho-sqrtRes)/Ax2;
		   };
  auto alpha2rho = [drho2,Brho,rho02](const double alpha) -> double{
		     const double alpha2 = alpha*alpha;
		     double rho2 = drho2*alpha2 + Brho*alpha + rho02;
		     return sqrt(rho2);
		   };
  auto phi2alpha = [dp,p1]
    (const double phi) -> double{
		    
		    double tanPhi = tan(phi);
		    double den = tanPhi*dp.x - dp.y;
		    if(fabs(den) < 1.0e-10)
		      return 1.0e35;
		    double num = p1.y-tanPhi*p1.x;
		    return num/den;
		  };
  auto alpha2phi = [dp,p1]
    (const double alpha) -> double{
		     double x = p1.x + alpha*dp.x;
		     double y = p1.y + alpha*dp.y;
		     if(x == 0.0 && y == 0.0)
		       return 0.0;
		     double phi = atan2(y,x);
		     if(std::signbit(phi))
		       phi += 2.0*M_PI;
		     return phi;
		   };
  auto z2alpha = [dp,p1]
    (const double z) -> double{
		   if(fabs(dp.z) < 1.0e-10)
		     return 1.0e35;
		   else
		     return (z-p1.z)/dp.z;
		 };
  auto alpha2z = [dp,p1]
    (const double alpha) -> double{
		   return p1.z + alpha*dp.z;
		 };
  
  //Check the minimum alpha value to enter the cylinder
  const double eps = 1.0e-8;
  //Radius
  double alphaRadIn  = 0.0;
  double alphaRadOut = 1.0;
  if(p1cyl.x >= radCyl){
    //The particle comes from outsied cylinder
    double alpha1,alpha2;
    rho2alpha(radCyl-eps,alpha1,alpha2);
    //If the particle enters the cylinder, both alpha values
    //must be positive
    if(std::signbit(alpha1))
      return; // Path never enters the cylinder
    if(alpha1 < alpha2){
      alphaRadIn  = alpha1;
      alphaRadOut = alpha2;
    }else{
      alphaRadIn  = alpha2;
      alphaRadOut = alpha1;
    }
    if(alphaRadIn > 1.0)
      return; // Path never enters the cylinder
  }
  else if(p2cyl.x >= radCyl){
    //The particle go from inside the cylinder to outside
    double alpha1,alpha2;
    rho2alpha(radCyl-eps,alpha1,alpha2);
    //One alpha value must be negative and the other positive
    if(std::signbit(alpha1))
      alphaRadOut = std::max(0.0,alpha2); //Avoid roinding errors on boundary
    else
      alphaRadOut = alpha1;
  }
  //Z
  double alphaZin  = 0.0;
  double alphaZout = 1.0;

  if(p1cyl.z < zminCyl){
    alphaZin = z2alpha(zminCyl+eps);
    if(alphaZin > 1.0)
      return; //Doesn't reaches the cylinder
  }
  else if(p1cyl.z >= zmaxCyl){
    alphaZin = z2alpha(zmaxCyl-eps);
    if(alphaZin > 1.0)
      return; //Doesn't reaches the cylinder
  }

  if(p2cyl.z < zminCyl)
    alphaZout = z2alpha(zminCyl+eps);
  else if(p2cyl.z >= zmaxCyl)
    alphaZout = z2alpha(zmaxCyl-eps);
  
  //Get the required path fraction to enter "p1cyl" to the cylinder
  double alphamin = std::max(alphaRadIn,alphaZin);

  //Get the required path fraction to enter "p2cyl" to the cylinder
  double alphamax = std::min(std::min(alphaRadOut,alphaZout),1.0);
  
  if(alphamin >= alphamax || fabs(alphamin-alphamax) < 1.0e-10){
    return; //Nothing to score
  }
  
  //Get the point at alphamin
  vect3d in;
  if(alphamin == 0.0){
    in.x = p1cyl.x; in.y = p1cyl.y; in.z = p1cyl.z;
  }
  else{
    in.x = alpha2rho(alphamin);
    in.y = alpha2phi(alphamin);
    in.z = alpha2z(alphamin);
  }

  //Get the point at alphamax
  vect3d out;
  if(alphamax == 1.0){
    out.x = p2cyl.x; out.y = p2cyl.y; out.z = p2cyl.z;
  }
  else{
    out.x = alpha2rho(alphamax);
    out.y = alpha2phi(alphamax);
    out.z = alpha2z(alphamax);
  }

  //Save current alpha
  double alpha = alphamin;
  
  //Get the corresponding bin index for each "in" coordinate
  vect3i ibin;
  if(radCylmin > 1.0e-9)
    if(in.x < radCylmin)
      ibin.x = 0;
    else
      ibin.x = 1+static_cast<long int>((in.x-radCylmin)/dbinCyl.x);
  else
    ibin.x = in.x/dbinCyl.x;
  ibin.y = in.y/dbinCyl.y;
  ibin.z = (in.z-zminCyl)/dbinCyl.z;

  //Get the corresponding bin index for each "out" coordinate
  vect3i obin;
  if(radCylmin > 1.0e-9)
    if(out.x < radCylmin)
      obin.x = 0;
    else
      obin.x = 1+(out.x-radCylmin)/dbinCyl.x;
  else
    obin.x = out.x/dbinCyl.x;
  obin.y = out.y/dbinCyl.y;
  obin.z = (out.z-zminCyl)/dbinCyl.z;
  
  long int globBin = nbinsCyl.x*(nbinsCyl.y*ibin.z + ibin.y) + ibin.x;

  //Create arrays to store all alpha cuts. Alpha cuts are specified using
  //the pair <alpha,index>, where index specify the cutted coordinate.
  std::array<std::pair<double,int>,3*(meshAxeMax+1)+1> cuts;

  //Add the final alpha value (alpha = 1.0)
  size_t ncuts = 0;
  //* Calculate rho cuts
  //**********************
  if(drho2 > 1.0e-10){ //Ensure that dx^2 + dy^2 is not zero
    //The particle crosses radial bins
    //Check if the path crosses to the previous bin
    if(ibin.x == 0){
      if(ibin.x != obin.x){
	//The path goes to outside bins. So, the same bin can be
	//crossed only one time. Compute the first cross.
	double alpha01Raw,alpha02Raw;
	rho2alpha(rPlanesCyl[0],alpha01Raw,alpha02Raw);
	double alpha01,alpha02;
	alpha01 = alpha01Raw-alphamin;
	alpha02 = alpha02Raw-alphamin;
	//One of the two crosses must be negative or both values could
	//be zero if the initial point is so close to bin limit 
	if(alpha01 < 0.0){
	  if(alpha02 < 0.0) //Avoid rounding errors
	    alpha01 = 0.0;
	  else
	    alpha01 = alpha02;
	}
	cuts[ncuts++] = std::pair<double,int>(alpha01Raw,1);
	//Compute the following crosses
	for(long int irho = 1; irho < obin.x; ++irho){
	  double alpha1Raw,alpha2Raw;
	  rho2alpha(rPlanesCyl[irho],alpha1Raw,alpha2Raw);
	  double alpha1,alpha2;
	  alpha1 = alpha1Raw-alphamin;
	  alpha2 = alpha2Raw-alphamin;	  
	  //One of the two crosses must be negative and the other positive
	  if(std::signbit(alpha1))
	    alpha1 = alpha2;
	  cuts[ncuts++] = std::pair<double,int>(alpha1Raw,1);
	}
      }
    }
    else{
      //Inner radius bins could be crossed by the path two times.
      //Check the first inner shell
      double alpha01Raw,alpha02Raw;
      rho2alpha(rPlanesCyl[ibin.x-1],alpha01Raw,alpha02Raw);
      double alpha01,alpha02;
      alpha01 = alpha01Raw-alphamin;
      alpha02 = alpha02Raw-alphamin;	        
      //If both alpha values are negative, this inner shell will not be
      //crossed. On the other hand, if both are positive, this bin will
      //be crossed. However, some of the two values could be positive 
      //and the other negative due rounding errors if the path initial
      //point is on the bin wall.
      bool innerShellCrossed = true;
      if(std::signbit(alpha01) != std::signbit(alpha02)){
	//Take the farest alpha value 
	double alphaAlpha = alpha01;
	if(fabs(alphaAlpha) < fabs(alpha02))
	  alphaAlpha = alpha02;

	if(std::signbit(alphaAlpha)){
	  //Dominant alpha is negative, inner shell is not crossed
	  innerShellCrossed = false;
	}
	else{
	  //Dominant alpha is positive, inner shell will be crossed
	  cuts[ncuts++] = std::pair<double,int>(alphamin,-1);
	  cuts[ncuts++] = std::pair<double,int>(alphaAlpha+alphamin,1);
	}
      }
      else if(std::signbit(alpha01)){
	//Both negative
	innerShellCrossed = false; //No inner shell crossed
      }
      else{
	//Both positive
	if(alpha01 < alpha02){
	  cuts[ncuts++] = std::pair<double,int>(alpha01Raw,-1); //Go inner
	  cuts[ncuts++] = std::pair<double,int>(alpha02Raw,1);  //Go outer
	}else{
	  cuts[ncuts++] = std::pair<double,int>(alpha02Raw,-1); //Go inner
	  cuts[ncuts++] = std::pair<double,int>(alpha01Raw,1);  //Go outer
	}
      }

      //Compute cuts
      if(innerShellCrossed){
	//Inner shells will be crossed, check all of them
	for(long int irho = ibin.x-2; irho >= 0; --irho){
	  double alpha1Raw,alpha2Raw;
	  rho2alpha(rPlanesCyl[irho],alpha1Raw,alpha2Raw);
	  double alpha1,alpha2;
	  alpha1 = alpha1Raw-alphamin;
	  alpha2 = alpha2Raw-alphamin;	  
	  
	  if(std::signbit(alpha1) || std::signbit(alpha2))
	    break;
	  if(alpha1 < alpha2){
	    cuts[ncuts++] = std::pair<double,int>(alpha1Raw,-1); //Go inner
	    cuts[ncuts++] = std::pair<double,int>(alpha2Raw,1);  //Go outer
	  }else{
	    cuts[ncuts++] = std::pair<double,int>(alpha2Raw,-1); //Go inner
	    cuts[ncuts++] = std::pair<double,int>(alpha1Raw,1);  //Go outer
	  }
	}
      }
      //Outside shells
      //Check first outside shell
      rho2alpha(rPlanesCyl[ibin.x],alpha01Raw,alpha02Raw);
      alpha01 = alpha01Raw-alphamin;
      alpha02 = alpha02Raw-alphamin;	        
      //One alpha must be positive and the other negative, since out shell
      //can be only crossed by one point. However, another time, due rounding
      //errors both could be positive or negative if the first point is on
      //the bin wall.
      if(std::signbit(alpha01) == std::signbit(alpha02)){ //Rounding error case
	if(std::signbit(alpha01)){ //Both negative, cross to next bin
	  if(ibin.x < nbinsCyl.x-1){
	    cuts[ncuts++] = std::pair<double,int>(alphamin,1);
	  }
	}
	else{ //Both positive, the cross to next bin is on greater alpha
	  if(ibin.x < nbinsCyl.x-1){
	    if(alpha01 < alpha02) alpha01 = alpha02;
	    cuts[ncuts++] = std::pair<double,int>(alpha01Raw,1);
	  }
	}
      }
      else{ //Normal case
	//Take the positive alpha cut
	if(std::signbit(alpha01))
	  cuts[ncuts++] = std::pair<double,int>(alpha02Raw,1);
	else
	  cuts[ncuts++] = std::pair<double,int>(alpha01Raw,1);
      }

      //Check other outer shells
      for(long int irho = ibin.x+1; irho < obin.x; ++irho){
	double alpha1Raw,alpha2Raw;
	rho2alpha(rPlanesCyl[irho],alpha1Raw,alpha2Raw);
	double alpha1,alpha2;
	alpha1 = alpha1Raw-alphamin;
	alpha2 = alpha2Raw-alphamin;	
	//One must be negative
	if(std::signbit(alpha1))
	  alpha1 = alpha2;
	cuts[ncuts++] = std::pair<double,int>(alpha1Raw,1);	  
      }
      
    }
  }

  //* Calculate phi cuts
  //**********************
  if(ibin.y != obin.y){
    //Obtain phi direction calculating the two possible first crosses taking
    //alphamin as origin point
    double alpha01Raw = phi2alpha(phiPlanesCyl[ibin.y]);
    double alpha01 = alpha01Raw-alphamin;
    long int nextBin = (ibin.y == nbinsCyl.y-1) ? 0 : ibin.y+1;
    double alpha02Raw = phi2alpha(phiPlanesCyl[nextBin]);
    double alpha02 = alpha02Raw-alphamin;
    //Both alphas can have different or equal sign.
    bool dirNeg = false;
    if(std::signbit(alpha01) == std::signbit(alpha02)){
      //Both have equal sign
      if(std::signbit(alpha01)){
	//Both are negative. As input and output points have different phi bin,
	//this must be caused due rounding errors.
	if(alpha01 < alpha02){
	  //alpha01 is farest from the init point. So, we are moving on
	  //positive bin direction
	  cuts[ncuts++] = std::pair<double,int>(alphamin,2);	  
	}
	else{
	  //alpha02 is farest from the init point. So, we are moving on
	  //negative bin direction
	  cuts[ncuts++] = std::pair<double,int>(alphamin,-2);
	  dirNeg = true;
	}	
      }
      else{
	if(alpha01 > alpha02){
	  //alpha01 is farest from the init point. So, we will cross first
	  //alpha02 in the positive direction
	  cuts[ncuts++] = std::pair<double,int>(alpha02Raw,2);	  
	}
	else{
	  //alpha02 is farest from the init point. So, we will cross first
	  //alpha01 in the negative direction
	  cuts[ncuts++] = std::pair<double,int>(alpha01Raw,-2);
	  dirNeg = true;
	}
      }
    }
    else{
      //Both signs are not equal. The path will cross the positive alpha
      if(std::signbit(alpha01))
	cuts[ncuts++] = std::pair<double,int>(alpha02Raw,2);
      else{
	cuts[ncuts++] = std::pair<double,int>(alpha01Raw,-2);
	dirNeg = true;
      }
    }
    
    //Obtain other cuts
    if(dirNeg){ //Negative direction
      if(obin.y < ibin.y){
	for(long int iphi = ibin.y-1; iphi > obin.y; --iphi){
	  cuts[ncuts++] = std::pair
	    <double,int>(phi2alpha(phiPlanesCyl[iphi]),-2);	  
	}
      }
      else{
	//Interval from ibin to 0
	for(long int iphi = ibin.y-1; iphi >= 0; --iphi)
	  cuts[ncuts++] = std::pair
	    <double,int>(phi2alpha(phiPlanesCyl[iphi]),-2);
	//Interval from 0 to obin
	for(long int iphi = nbinsCyl.y-1; iphi > obin.y; --iphi)
	  cuts[ncuts++] = std::pair
	    <double,int>(phi2alpha(phiPlanesCyl[iphi]),-2);	
      }
    }
    else{ //Positive direction
      if(obin.y != nextBin){
	if(obin.y > nextBin){
	  for(long int iphi = nextBin+1; iphi <= obin.y; ++iphi)
	    cuts[ncuts++] = std::pair
	      <double,int>(phi2alpha(phiPlanesCyl[iphi]),2);	  
	}
	else{
	  //Interval from ibin to nbins-1
	  for(long int iphi = nextBin+1; iphi < nbinsCyl.y; ++iphi)
	    cuts[ncuts++] = std::pair
	      <double,int>(phi2alpha(phiPlanesCyl[iphi]),2);	  
	  //Interval from 0 to obin
	  for(long int iphi = 0; iphi <= obin.y; ++iphi)
	    cuts[ncuts++] = std::pair
	      <double,int>(phi2alpha(phiPlanesCyl[iphi]),2);	  	
	}
      }
    }
    
  }
  
  //* Calculate Z cuts
  //**********************
  if(ibin.z != obin.z && fabs(dp.z) > 1.0e-10){

    //Obtain first cut
    if(ibin.z > obin.z){
      for(long int ipz = ibin.z; ipz > obin.z; --ipz){
	cuts[ncuts++] = std::pair<double,int>
	  ((zPlanesCyl[ipz]-p1.z)/dp.z,-3);
      }
    }
    else{
      for(long int ipz = ibin.z+1; ipz <= obin.z; ++ipz){
	cuts[ncuts++] = std::pair<double,int>
	  ((zPlanesCyl[ipz]-p1.z)/dp.z,3);
      }
    }
  }

  //Add final cut
  cuts[ncuts++] = std::pair<double,int>(alphamax,0);
  
  //Sort cuts
  std::sort(cuts.begin(),cuts.begin()+ncuts);

  //Travell all cuts
  double alphalimit = alphamax*0.999999;
  for(size_t icut = 0; icut < ncuts; ++icut){
    //Get next alpha cut
    std::pair<double,int> cut = cuts[icut];
    double toTravel = cut.first - alpha;
    if(cut.first >= alphalimit){
      //Reached finish point
      double l = (alphamax-alpha)*dsef;

      scoreTally(wght,E,muenVal,l,globBin,nhist,
		 cylindrical,cylindrical2,cylindricalTmp,cylindricalLastHist);
      return;
    }

    //Score traveled distance
    double l = toTravel*dsef;

    alpha = cut.first;
    
    if(l > 0.0)
      scoreTally(wght,E,muenVal,l,globBin,nhist,
		 cylindrical,cylindrical2,cylindricalTmp,cylindricalLastHist);
    //Update index
    if(cut.second == 1){ //External radius cross
      //Update bin index
      ibin.x += 1;
      globBin += 1;
      if(ibin.x >= nbinsCyl.x){
	long int prev = icut-1;
	prev = std::max(prev,0l);
	printf("kermaTrackLengthCyl: Unexpected termination: "
	       "last radial bin crossed "
	       "(alpha = %12.4E alpha min: %12.4E alpha max: %12.4E)\n"
	       " Traveled segment: %12.4E total distance %12.4E\n"
	       " This %s the last cut (Number of cuts: %lu)\n"
	       "(    Next cut: alpha = %12.4E flag = %d)\n"
	       "(Previous cut: alpha = %12.4E flag = %d)\n"
	       "p1(cart): (%12.4E,%12.4E,%12.4E)\n"
	       "dp(cart): (%12.4E,%12.4E,%12.4E)\n"
	       "p1(cyl) : (%12.4E,%12.4E,%12.4E)\n"
	       "p2(cyl) : (%12.4E,%12.4E,%12.4E)\n",
	       alpha,alphamin,alphamax,toTravel,dsef,
	       icut == ncuts-1 ? "is" : "is not",
	       static_cast<unsigned long>(ncuts),cuts[icut+1].first,
	       cuts[icut+1].second,cuts[prev].first,cuts[prev].second,
	       p1.x,p1.y,p1.z,dp.x,dp.y,dp.z,p1cyl.x,p1cyl.y,p1cyl.z,
	       p2cyl.x,p2cyl.y,p2cyl.z);
	return; //Escapes from the cylinder
      }
    }
    else if(cut.second == -1){ //Out radius cross
      //Update bin index
      ibin.x -= 1;
      globBin -= 1;
      if(ibin.x < 0){
	long int prev = icut-1;
	prev = std::max(prev,0l);
	printf("kermaTrackLengthCyl: Unexpected termination: "
	       "first radial bin 'crossed' "
	       "(alpha = %12.4E alpha min: %12.4E alpha max: %12.4E)\n"
	       " Traveled segment: %12.4E total distance %12.4E\n"
	       " This %s the last cut (Number of cuts: %lu)\n"
	       "(    Next cut: alpha = %12.4E flag = %d)\n"
	       "(Previous cut: alpha = %12.4E flag = %d)\n"
	       "p1(cart): (%12.4E,%12.4E,%12.4E)\n"
	       "dp(cart): (%12.4E,%12.4E,%12.4E)\n"
	       "p1(cyl) : (%12.4E,%12.4E,%12.4E)\n"
	       "p2(cyl) : (%12.4E,%12.4E,%12.4E)\n",
	       alpha,alphamin,alphamax,toTravel,dsef,
	       icut == ncuts-1 ? "is" : "is not",
	       static_cast<unsigned long>(ncuts),cuts[icut+1].first,
	       cuts[icut+1].second,cuts[prev].first,cuts[prev].second,
	       p1.x,p1.y,p1.z,dp.x,dp.y,dp.z,p1cyl.x,p1cyl.y,p1cyl.z,
	       p2cyl.x,p2cyl.y,p2cyl.z);
	return; //Escapes from the cylinder
      }
    }
    else if(cut.second == 2){ //Positive phi bin cross
      //Update bin index
      ibin.y += 1;
      globBin += nbinsCyl.x;
      if(ibin.y >= nbinsCyl.y){
	//Last bin crossed, go to first bin
	ibin.y = 0;
	globBin = nbinsCyl.x*(ibin.z*nbinsCyl.y + ibin.y) + ibin.x;
      }
    }
    else if(cut.second == -2){ //Negative phi bin cross
      //Update bin index
      ibin.y -= 1;
      globBin -= nbinsCyl.x;
      if(ibin.y < 0){
	//First bin crossed, go to last bin
	ibin.y = nbinsCyl.y-1;
	globBin = nbinsCyl.x*(ibin.z*nbinsCyl.y + ibin.y) + ibin.x;
      }
    }
    else if(cut.second == 3){ //Positive Z bin cross
      //Update bin index
      ibin.z += 1;
      globBin += nrphiCyl;
      if(ibin.z >= nbinsCyl.z){
	printf("kermaTrackLengthCyl: Unexpected termination: "
	       "last Z bin crossed "
	       "(alpha = %12.4E alpha max: %12.4E)\n",
	       alpha,alphamax);	
	return; //Escapes from the cylinder
      }
    }
    else if(cut.second == -3){ //Negative Z bin cross
      //Update bin index
      ibin.z -= 1;
      globBin -= nrphiCyl;
      if(ibin.z < 0){
	printf("kermaTrackLengthCyl: Unexpected termination: "
	       "first Z bin crossed "
	       "(alpha = %12.4E alpha max: %12.4E)\n",
	       alpha,alphamax);	
	return; //Escapes from the cylinder
      }
    }
  }
}

void pen_tallyKermaTrackLength::kermaTrackLengthSph(const unsigned long long nhist,
						    const double wght,
						    const double E,
						    const double muenVal,
						    const double dsef,
						    const pen_tally_KTL::vect3d p1,
						    const pen_tally_KTL::vect3d dp,
						    const pen_tally_KTL::vect3d p1sph,
						    const pen_tally_KTL::vect3d p2sph){

  //
  //   x = p1.x + alpha*dp.x
  //   y = p1.y + alpha*dp.y
  //   z = p1.z + alpha*dp.z
  //
  // rad   = sqrt(x*x+y*y+z*z)
  // phi   = atan2(y,x)
  // theta = atan2(sqrt(x*x+y*y),z)
  
  using namespace pen_tally_KTL;

  //Check if both points are inside the minimus radius
  if(p1sph.x <= radSphmin && p2sph.x <= radSphmin)
    return; //Nothing to score
  
  //Precalculate some values
  double dx2 = dp.x*dp.x;
  double dy2 = dp.y*dp.y;
  double dz2 = dp.z*dp.z;
  double p1z2 = p1.z*p1.z;
  double alphaZ0;
  if(fabs(dp.z) < 1.0e-10)
    alphaZ0 = -1.0e35;
  else
    alphaZ0 = -p1.z/dp.z;

  double drho2 = dx2 + dy2;
  double drad2 = drho2 + dz2;
  double rho02 = p1.x*p1.x + p1.y*p1.y;
  double rad02 = rho02 + p1z2;

  double Brho  = 2.0*(p1.x*dp.x + p1.y*dp.y);
  double Brho2Brad = 2.0*p1.z*dp.z;
  double Brad  = Brho + Brho2Brad;
  double Brad2 = Brad*Brad;

  // Create auxiliary lambda functions to convert from coordinates
  // to alpha and from alpha to coordinates
  auto rad2alpha = [drad2,Brad,Brad2,rad02]
    (const double rad, double& alphaP, double& alphaM)
    -> void{
		     
		     if(drad2 < 1.0e-10){
		       alphaP = alphaM = -1.0e35;
		       return;
		     }
		     
		     const double Ax2 = 2.0*drad2;
		     const double Ax4 = 4.0*drad2;
		     double C = rad02 - rad*rad;
		     double sqrtArg = Brad2 - Ax4*C;
		     if(std::signbit(sqrtArg)){
		       alphaP = alphaM = -1.0e35;
		       return;
		     }
		     double sqrtRes = sqrt(sqrtArg);
		     alphaP = (-Brad+sqrtRes)/Ax2;
		     alphaM = (-Brad-sqrtRes)/Ax2;
		   };
  auto alpha2rad = [drad2,Brad,rad02](const double alpha) -> double{
		     const double alpha2 = alpha*alpha;
		     double rad2 = drad2*alpha2 + Brad*alpha + rad02;
		     return sqrt(rad2);
		   };
  auto phi2alpha = [dp,p1]
    (const double phi) -> double{
		    
		    double tanPhi = tan(phi);
		    double den = tanPhi*dp.x - dp.y;
		    if(fabs(den) < 1.0e-10)
		      return 1.0e35;
		    double num = p1.y-tanPhi*p1.x;
		    return num/den;
		  };
  auto alpha2phi = [dp,p1]
    (const double alpha) -> double{
		     double x = p1.x + alpha*dp.x;
		     double y = p1.y + alpha*dp.y;
		     if(x == 0.0 && y == 0.0)
		       return 0.0;
		     double phi = atan2(y,x);
		     if(phi < 0.0)
		       phi += 2.0*M_PI;
		     return phi;
		   };
  auto theta2alpha = [drho2,dp,p1z2,alphaZ0,dz2,Brho,Brho2Brad,rho02]
    (const double theta, double& alphaP, double& alphaM)
    -> void{
		       if(fabs(theta-PI05) < 1.0e-8){
			 alphaP = alphaM = alphaZ0;
			 return;
		       }
		       
		       double tanTheta = tan(theta);
		       double tanTheta2 = tanTheta*tanTheta;
		       double A = drho2 - tanTheta2*dz2;
		       if(fabs(A) < 1.0e-10){
			 alphaP = alphaM = -1.0e35;
			 return;
		       }

		       double B = Brho - Brho2Brad*tanTheta2;
		       double C = rho02 - p1z2*tanTheta2;
		       double sqrtArg = B*B-4.0*A*C;
		       if(std::signbit(sqrtArg)){
			 alphaP = alphaM = -1.0e35;
			 return;
		       }
		       double sqrtRes = sqrt(sqrtArg);
		       double Ax2 = 2.0*A;
		       alphaP = (-B+sqrtRes)/Ax2;
		       alphaM = (-B-sqrtRes)/Ax2;
		     };
  auto alpha2theta = [dp,p1]
    (const double alpha) -> double{
		       double x = p1.x + alpha*dp.x;
		       double y = p1.y + alpha*dp.y;
		       double z = p1.z + alpha*dp.z;
		       
		       double num = sqrt(x*x + y*y);
		       if(num == 0.0 && z == 0.0)
			 return 0.0;
		       double theta = atan2(num,z);
		       if(std::signbit(theta)) //Just in case
			 return fabs(theta);
		       return theta;
		 };
    
  //Check the minimum alpha value to enter the sphere
  const double eps = 1.0e-8;
  //Radius
  double alphamin = 0.0;
  if(p1sph.x >= radSph){
    //Check if the particle will reach the sphere max radius
    double alpha1,alpha2;
    rad2alpha(radSph-eps,alpha1,alpha2);
    alphamin = closestAlpha(alpha1,alpha2);
    if(std::signbit(alphamin))
      return; // Path never enters the sphere    
  }

  //Check the maximum alpha value to enter the sphere
  //Radius
  double alphamax = 1.0;
  if(p2sph.x >= radSph){
    double alpha1, alpha2;
    rad2alpha(radSph-eps,alpha1,alpha2);
    alphamax = farestAlpha(alpha1,alpha2);
    if(std::signbit(alphamax)){
      double alpha1rhoM, alpha2rhoM;
      double alpha1rhoP, alpha2rhoP;
      rad2alpha(p1sph.x,alpha1rhoP,alpha1rhoM);
      rad2alpha(p2sph.x,alpha2rhoP,alpha2rhoM);
      
      printf("pen_tally_KTL::kermaTrackLengthSph: Unexpected error: "
	     "p1 reaches the sphere but p2 can't be moved to "
	     "the radial frontiner (r=%12.4E cm)\n"
	     "     Alpha min: %12.4E  Rad min: %12.4E (alphas = (%12.4E,%12.4E))\n"
	     " Alpha rad max: %12.4E  Rad max: %12.4E (alphas = (%12.4E,%12.4E))\n"
	     " Posible alpha values: %12.4E and %12.4E\n"
	     "      P1: (%12.4E,%12.4E,%12.4E)\n"
	     "      P2: (%12.4E,%12.4E,%12.4E)\n"
	     "      dP: (%12.4E,%12.4E,%12.4E)\n"
	     "   P1sph: (%12.4E,%12.4E,%12.4E)\n"
	     "   P2sph: (%12.4E,%12.4E,%12.4E)\n",
	     radSph,alphamin,p1sph.x,alpha1rhoP,alpha1rhoM,
	     alphamax,p2sph.x,alpha2rhoP,alpha2rhoM,
	     alpha1,alpha2,p1.x,p1.y,p1.z,p1.x+dp.x,p1.y+dp.y,p1.z+dp.z,
	     dp.x,dp.y,dp.z,p1sph.x,p1sph.y,p1sph.z,p2sph.x,p2sph.y,p2sph.z);
      return;
    }
  }
  
  if(alphamin >= alphamax || fabs(alphamin-alphamax) < 1.0e-10){
    return; //Nothing to score
  }
  
  //Get the point at alphamin
  vect3d in;
  if(alphamin == 0.0){
    in.x = p1sph.x; in.y = p1sph.y; in.z = p1sph.z;
  }
  else{
    in.x = alpha2rad(alphamin);
    in.y = alpha2theta(alphamin);
    in.z = alpha2phi(alphamin);
  }

  //Get the point at alphamax
  vect3d out;
  if(alphamax == 1.0){
    out.x = p2sph.x; out.y = p2sph.y; out.z = p2sph.z;
  }
  else{
    out.x = alpha2rad(alphamax);
    out.y = alpha2theta(alphamax);
    out.z = alpha2phi(alphamax);
  }

  //Save current alpha
  double alpha = alphamin;
  
  //Get the corresponding bin index for each "in" coordinate
  vect3i ibin;
  if(radSphmin > 1.0e-9)
    if(in.x < radSphmin)
      ibin.x = 0;
    else
      ibin.x = 1+static_cast<long int>((in.x-radSphmin)/dbinSph.x);
  else  
    ibin.x = in.x/dbinSph.x;
  ibin.y = in.y/dbinSph.y;
  ibin.z = in.z/dbinSph.z;

  //Get the corresponding bin index for each "out" coordinate
  vect3i obin;
  if(radSphmin > 1.0e-9)
    if(out.x < radSphmin)
      obin.x = 0;
    else
      obin.x = 1+(out.x-radSphmin)/dbinSph.x;
  else
    obin.x = out.x/dbinSph.x;
  obin.y = out.y/dbinSph.y;
  obin.z = out.z/dbinSph.z;

  long int globBin = nbinsSph.x*(nbinsSph.y*ibin.z + ibin.y) + ibin.x;

  //Create arrays to store all alpha cuts. Alpha cuts are specified using
  //the pair <alpha,index>, where index specify the cutted coordinate.
  std::array<std::pair<double,int>,3*(meshAxeMax+1)+1> cuts;

  //Add the final alpha value (alpha = 1.0)
  size_t ncuts = 0;
  //* Calculate rho cuts
  //**********************
  if(drad2 > 1.0e-10){ //Ensure that dx^2 + dy^2 + dz^2 is not zero
    //The particle crosses radial bins
    //Check if the path crosses to the previous bin
    if(ibin.x == 0){
      if(ibin.x != obin.x){
	//The path goes to outside bins. So, the same bin can be
	//crossed only one time. Compute the first cross.
	double alpha01Raw,alpha02Raw;
	rad2alpha(rPlanesSph[0],alpha01Raw,alpha02Raw);
	double alpha01,alpha02;
	alpha01 = alpha01Raw-alphamin;
	alpha02 = alpha02Raw-alphamin;
	//One of the two crosses must be negative or both values could
	//be zero if the initial point is so close to bin limit 
	if(alpha01 < 0.0){
	  if(alpha02 < 0.0) //Avoid rounding errors
	    alpha01 = 0.0;
	  else
	    alpha01 = alpha02;
	}
	cuts[ncuts++] = std::pair<double,int>(alpha01Raw,1);
	//Compute the following crosses
	for(long int irad = 1; irad < obin.x; ++irad){
	  double alpha1Raw,alpha2Raw;
	  rad2alpha(rPlanesSph[irad],alpha1Raw,alpha2Raw);
	  double alpha1,alpha2;
	  alpha1 = alpha1Raw-alphamin;
	  alpha2 = alpha2Raw-alphamin;	  
	  //One of the two crosses must be negative and the other positive
	  if(alpha1 < 0.0)
	    alpha1 = alpha2;
	  cuts[ncuts++] = std::pair<double,int>(alpha1Raw,1);
	}
      }
    }
    else{
      //Inner radius bins could be crossed by the path two times.
      //Check the first inner shell
      double alpha01Raw,alpha02Raw;
      rad2alpha(rPlanesSph[ibin.x-1],alpha01Raw,alpha02Raw);
      double alpha01,alpha02;
      alpha01 = alpha01Raw-alphamin;
      alpha02 = alpha02Raw-alphamin;	        
      //If both alpha values are negative, this inner shel will not be
      //crossed. On the other hand, if both are positive, this bin will
      //be crossed. However, some of the two values could be positive 
      //and the other negative due rounding errors if the path initial
      //point is on the bin wall.
      bool innerShellCrossed = true;
      if(std::signbit(alpha01) != std::signbit(alpha02)){
	//Take the farest alpha value 
	double alphaAlpha = alpha01;
	if(fabs(alphaAlpha) < fabs(alpha02))
	  alphaAlpha = alpha02;

	if(std::signbit(alphaAlpha)){
	  //Dominant alpha is negative, inner shell is not crossed
	  innerShellCrossed = false;
	}
	else{
	  //Dominant alpha is positive, inner shell will be crossed
	  cuts[ncuts++] = std::pair<double,int>(alphamin,-1);
	  cuts[ncuts++] = std::pair<double,int>(alphaAlpha+alphamin,1);
	}
      }
      else if(std::signbit(alpha01)){
	//Both negative
	innerShellCrossed = false; //No inner shell crossed
      }
      else{
	//Both positive
	if(alpha01 < alpha02){
	  cuts[ncuts++] = std::pair<double,int>(alpha01Raw,-1); //Go inner
	  cuts[ncuts++] = std::pair<double,int>(alpha02Raw,1);  //Go outer
	}else{
	  cuts[ncuts++] = std::pair<double,int>(alpha02Raw,-1); //Go inner
	  cuts[ncuts++] = std::pair<double,int>(alpha01Raw,1);  //Go outer
	}
      }

      //Compute cuts
      if(innerShellCrossed){
	//Inner shells will be crossed, check all of them
	for(long int irad = ibin.x-2; irad >= 0; --irad){
	  double alpha1Raw,alpha2Raw;
	  rad2alpha(rPlanesSph[irad],alpha1Raw,alpha2Raw);
	  double alpha1,alpha2;
	  alpha1 = alpha1Raw-alphamin;
	  alpha2 = alpha2Raw-alphamin;	  	  
	  if(std::signbit(alpha1) || std::signbit(alpha2))
	    break;
	  if(alpha1 < alpha2){
	    cuts[ncuts++] = std::pair<double,int>(alpha1Raw,-1); //Go inner
	    cuts[ncuts++] = std::pair<double,int>(alpha2Raw,1);  //Go outer
	  }else{
	    cuts[ncuts++] = std::pair<double,int>(alpha2Raw,-1); //Go inner
	    cuts[ncuts++] = std::pair<double,int>(alpha1Raw,1);  //Go outer
	  }
	}
      }
      //Outside shells
      //Check first outside shell
      rad2alpha(rPlanesSph[ibin.x],alpha01Raw,alpha02Raw);
      alpha01 = alpha01Raw-alphamin;
      alpha02 = alpha02Raw-alphamin;	        
      //One alpha must be positive and the other negative, since out shell
      //can be only crossed by one point. However, another time, due rounding
      //errors both could be positive or negative if the first point is on
      //the bin wall.
      if(std::signbit(alpha01) == std::signbit(alpha02)){ //Rounding error case
	if(std::signbit(alpha01)){ //Both negative, cross to next bin
	  if(ibin.x < nbinsSph.x-1){
	    cuts[ncuts++] = std::pair<double,int>(alphamin,1);
	  }
	}
	else{ //Both positive, the cross to next bin is on greater alpha
	  if(ibin.x < nbinsSph.x-1){
	    if(alpha01 < alpha02) alpha01 = alpha02;
	    cuts[ncuts++] = std::pair<double,int>(alpha01Raw,1);
	  }
	}
      }
      else{ //Normal case
	//Take the positive alpha cut
	if(std::signbit(alpha01))
	  cuts[ncuts++] = std::pair<double,int>(alpha02Raw,1);
	else
	  cuts[ncuts++] = std::pair<double,int>(alpha01Raw,1);
      }

      //Check other outer shells
      for(long int irad = ibin.x+1; irad < obin.x; ++irad){
	double alpha1Raw,alpha2Raw;
	rad2alpha(rPlanesSph[irad],alpha1Raw,alpha2Raw);
	double alpha1,alpha2;
	alpha1 = alpha1Raw-alphamin;
	alpha2 = alpha2Raw-alphamin;	
	//One must be negative
	if(std::signbit(alpha1))
	  alpha1 = alpha2;
	cuts[ncuts++] = std::pair<double,int>(alpha1Raw,1);	  
      }
      
    }
  }
  
  //* Calculate theta cuts
  //***********************
  if((ibin.y == 0 || ibin.y == nbinsSph.y-1) && (ibin.y == obin.y)){
    //Particle is in theta bin 0 or in the last bin, so no inner cones can be
    //crossed. In addition, input and output bins are equal, thus there are
    //no theta bins to cross.
  }
  else{

    //Get the cones that could be crossed
    long int initCone; //ibin.y cone index in the interval [0,halfnThetaCones]
    bool north; //Stores if the init point is on north semisphere
    if(ibin.y < halfnThetaCones){
      initCone = ibin.y;
      north = true;
    }
    else{
      initCone = (nbinsSph.y-1)-ibin.y;
      north = false;
    }

    long int finalCone; //obin.y cone index in the interval [0,halfnThetaCones]
    if(obin.y < halfnThetaCones){
      //Obin is on norh semisphere
      if(north){
	//Both, "in" and "out", on north semisphere
	finalCone = obin.y;
      }else{
	//"In" is on south and "out" on north
	finalCone = halfnThetaCones;
      }
    }
    else{
      //Obin is on south semisphere
      if(north){
	//"In" is on north and "out" on south
	finalCone = halfnThetaCones;
      } else{
	//Both, "in" and "out", on south semisphere	
	finalCone = (nbinsSph.y-1)-obin.y;
      }
    }
    
    //Create an array to store theta cuts
    std::array<double,meshAxeMax+1> thetaCuts;
    size_t nthetaCuts = 0;
    
    //Calculate alpha cuts for inner cones
    for(long int icone = initCone-1; icone >= 0; --icone){
      double alpha01Raw,alpha02Raw;
      theta2alpha(thetaPlanesSph[icone],alpha01Raw,alpha02Raw);
      if(alpha01Raw < alphamin && alpha02Raw < alphamin)
	break;
      if(alpha01Raw >= alphamin && alpha01Raw <= alphamax)
	thetaCuts[nthetaCuts++] = alpha01Raw;
      if(alpha02Raw >= alphamin && alpha02Raw <= alphamax)
	thetaCuts[nthetaCuts++] = alpha02Raw;
    }

    //Calculate alpha cuts for outer cones
    for(long int icone = initCone; icone < finalCone; ++icone){
      double alpha01Raw,alpha02Raw;
      theta2alpha(thetaPlanesSph[icone],alpha01Raw,alpha02Raw);
      if(alpha01Raw >= alphamin && alpha01Raw <= alphamax)
	thetaCuts[nthetaCuts++] = alpha01Raw;
      if(icone != planeThetaConeIndex){ //Avoid duplicate cut at z=0 plane
	if(alpha02Raw >= alphamin && alpha02Raw <= alphamax)
	  thetaCuts[nthetaCuts++] = alpha02Raw;
      }
    }
    

    //Sort theta cuts
    std::sort(thetaCuts.begin(),thetaCuts.begin()+nthetaCuts);
    thetaCuts[nthetaCuts++] = alphamax;
    
    double actualAlpha = thetaCuts[0];
    for(size_t icut = 1; icut < nthetaCuts; ++icut){
      double cut = thetaCuts[icut];
      double midAlpha = actualAlpha + (cut - actualAlpha)*0.5;
      double theta = alpha2theta(midAlpha);
      //Add a constant to theta bin to interpret that value as an
      //absolute bin index and not an increment.
      long int thetaBin = theta/dbinSph.y+10; 
      cuts[ncuts++] = std::pair<double,int>(actualAlpha,thetaBin);
      actualAlpha = cut; //Go to next alpha
    }
  }

  //* Calculate phi cuts
  //**********************
  if(ibin.z != obin.z){
    //Obtain phi direction calculating the two possible first crosses taking
    //alphamin as origin point
    double alpha01Raw = phi2alpha(phiPlanesSph[ibin.z]);
    double alpha01 = alpha01Raw-alphamin;
    long int nextBin = (ibin.z == nbinsSph.z-1) ? 0 : ibin.z+1;
    double alpha02Raw = phi2alpha(phiPlanesSph[nextBin]);
    double alpha02 = alpha02Raw-alphamin;
    //Both alphas can have different or equal sign.
    bool dirNeg = false;
    if(std::signbit(alpha01) == std::signbit(alpha02)){
      //Both have equal sign
      if(std::signbit(alpha01)){
	//Both are negative. As input and output points have different phi bin,
	//this must be caused due rounding errors.
	if(alpha01 < alpha02){
	  //alpha01 is farest from the init point. So, we are moving on
	  //positive bin direction
	  cuts[ncuts++] = std::pair<double,int>(alphamin,3);	  
	}
	else{
	  //alpha02 is farest from the init point. So, we are moving on
	  //negative bin direction
	  cuts[ncuts++] = std::pair<double,int>(alphamin,-3);
	  dirNeg = true;
	}	
      }
      else{
	if(alpha01 > alpha02){
	  //alpha01 is farest from the init point. So, we will cross first
	  //alpha02 in the positive direction
	  cuts[ncuts++] = std::pair<double,int>(alpha02Raw,3);	  
	}
	else{
	  //alpha02 is farest from the init point. So, we will cross first
	  //alpha01 in the negative direction
	  cuts[ncuts++] = std::pair<double,int>(alpha01Raw,-3);
	  dirNeg = true;
	}
      }
    }
    else{
      //Both signs are not equal. The path will cross the positive alpha
      if(std::signbit(alpha01))
	cuts[ncuts++] = std::pair<double,int>(alpha02Raw,3);
      else{
	cuts[ncuts++] = std::pair<double,int>(alpha01Raw,-3);
	dirNeg = true;
      }
    }
    
    //Obtain other cuts
    if(dirNeg){ //Negative direction
      if(obin.z < ibin.z){
	for(long int iphi = ibin.z-1; iphi > obin.z; --iphi){
	  cuts[ncuts++] = std::pair
	    <double,int>(phi2alpha(phiPlanesSph[iphi]),-3);	  
	}
      }
      else{
	//Interval from ibin to 0
	for(long int iphi = ibin.z-1; iphi >= 0; --iphi)
	  cuts[ncuts++] = std::pair
	    <double,int>(phi2alpha(phiPlanesSph[iphi]),-3);
	//Interval from 0 to obin
	for(long int iphi = nbinsSph.z-1; iphi > obin.z; --iphi)
	  cuts[ncuts++] = std::pair
	    <double,int>(phi2alpha(phiPlanesSph[iphi]),-3);	
      }
    }
    else{ //Positive direction
      if(obin.z != nextBin){
	if(obin.z > nextBin){
	  for(long int iphi = nextBin+1; iphi <= obin.z; ++iphi)
	    cuts[ncuts++] = std::pair
	      <double,int>(phi2alpha(phiPlanesSph[iphi]),3);	  
	}
	else{
	  //Interval from ibin to nbins-1
	  for(long int iphi = nextBin+1; iphi < nbinsSph.z; ++iphi)
	    cuts[ncuts++] = std::pair
	      <double,int>(phi2alpha(phiPlanesSph[iphi]),3);	  
	  //Interval from 0 to obin
	  for(long int iphi = 0; iphi <= obin.z; ++iphi)
	    cuts[ncuts++] = std::pair
	      <double,int>(phi2alpha(phiPlanesSph[iphi]),3);	  	
	}
      }
    }
    
  }

  //Add final cut
  cuts[ncuts++] = std::pair<double,int>(alphamax,0);
  
  //Sort cuts
  std::sort(cuts.begin(),cuts.begin()+ncuts);
  
  //Travell all cuts
  double alphalimit = alphamax*0.999999;
  for(size_t icut = 0; icut < ncuts; ++icut){
    //Get next alpha cut
    std::pair<double,int> cut = cuts[icut];
    double toTravel = cut.first - alpha;
    if(cut.first >= alphalimit){
      //Reached finish point
      double l = (alphamax-alpha)*dsef;
      
      scoreTally(wght,E,muenVal,l,globBin,nhist,
		 spherical,spherical2,sphericalTmp,sphericalLastHist);
      return;
    }

    //Score traveled distance
    double l = toTravel*dsef;

    //Update alpha
    alpha = cut.first;
    
    if(l > 0.0)
      scoreTally(wght,E,muenVal,l,globBin,nhist,
		 spherical,spherical2,sphericalTmp,sphericalLastHist);
    //Update index
    if(cut.second == 1){ //External radius cross
      //Update bin index
      ibin.x += 1;
      globBin += 1;
      if(ibin.x >= nbinsSph.x){
	long int prev = icut-1;
	prev = std::max(prev,0l);
	printf("kermaTrackLengthSph: Unexpected termination: "
	       "last radial bin crossed "
	       "(alpha = %12.4E alpha min: %12.4E alpha max: %12.4E)\n"
	       " Traveled segment: %12.4E total distance %12.4E\n"
	       " This %s the last cut (Number of cuts: %lu)\n"
	       "(    Next cut: alpha = %12.4E flag = %d)\n"
	       "(Previous cut: alpha = %12.4E flag = %d)\n"
	       "p1(cart): (%12.4E,%12.4E,%12.4E)\n"
	       "dp(cart): (%12.4E,%12.4E,%12.4E)\n"
	       "p1(sph) : (%12.4E,%12.4E,%12.4E)\n"
	       "p2(sph) : (%12.4E,%12.4E,%12.4E)\n",
	       alpha,alphamin,alphamax,toTravel,dsef,
	       icut == ncuts-1 ? "is" : "is not",
	       static_cast<unsigned long>(ncuts),cuts[icut+1].first,
	       cuts[icut+1].second,cuts[prev].first,cuts[prev].second,
	       p1.x,p1.y,p1.z,dp.x,dp.y,dp.z,p1sph.x,p1sph.y,p1sph.z,
	       p2sph.x,p2sph.y,p2sph.z);
	return; //Escapes from the sphere
      }
    }
    else if(cut.second == -1){ //Out radius cross
      //Update bin index
      ibin.x -= 1;
      globBin -= 1;
      if(ibin.x < 0){
	long int prev = icut-1;
	prev = std::max(prev,0l);
	printf("kermaTrackLengthSph: Unexpected termination: "
	       "first radial bin 'crossed' "
	       "(alpha = %12.4E alpha min: %12.4E alpha max: %12.4E)\n"
	       " Traveled segment: %12.4E total distance %12.4E\n"
	       " This %s the last cut (Number of cuts: %lu)\n"
	       "(    Next cut: alpha = %12.4E flag = %d)\n"
	       "(Previous cut: alpha = %12.4E flag = %d)\n"
	       "p1(cart): (%12.4E,%12.4E,%12.4E)\n"
	       "dp(cart): (%12.4E,%12.4E,%12.4E)\n"
	       "p1(sph) : (%12.4E,%12.4E,%12.4E)\n"
	       "p2(sph) : (%12.4E,%12.4E,%12.4E)\n",
	       alpha,alphamin,alphamax,toTravel,dsef,
	       icut == ncuts-1 ? "is" : "is not",
	       static_cast<unsigned long>(ncuts),cuts[icut+1].first,
	       cuts[icut+1].second,cuts[prev].first,cuts[prev].second,
	       p1.x,p1.y,p1.z,dp.x,dp.y,dp.z,p1sph.x,p1sph.y,p1sph.z,
	       p2sph.x,p2sph.y,p2sph.z);
	return; //Escapes from the sphere
      }
    }
    else if(cut.second == 2){ //Positive theta bin cross
      //Update bin index
      ibin.y += 1;
      globBin += nbinsSph.x;
      if(ibin.y >= nbinsSph.y){
	printf("kermaTrackLengthSph: Unexpected termination: "
	       "last theta bin crossed "
	       "(alpha = %12.4E alpha max: %12.4E)\n",
	       alpha,alphamax);	
	return; //Escapes from the sphere	
      }
    }
    else if(cut.second == -2){ //Negative theta bin cross
      //Update bin index
      ibin.y -= 1;
      globBin -= nbinsSph.x;
      if(ibin.y < 0){	
	printf("kermaTrackLengthSph: Unexpected termination: "
	       "first theta bin crossed "
	       "(alpha = %12.4E alpha max: %12.4E)\n",
	       alpha,alphamax);	
	return; //Escapes from the sphere	
      }
    }
    else if(cut.second == 3){ //Positive phi bin cross
      //Update bin index
      ibin.z += 1;
      globBin += nrthetaSph;
      if(ibin.z >= nbinsSph.z){
	//Last bin crossed, go to first phi bin
	ibin.z = 0;
	globBin = nbinsSph.x*(ibin.z*nbinsSph.y + ibin.y) + ibin.x;
      }
    }
    else if(cut.second == -3){ //Negative phi bin cross
      //Update bin index
      ibin.z -= 1;
      globBin -= nrthetaSph;
      if(ibin.z < 0){
	//First bin crossed, go to last phi bin
	ibin.z = nbinsSph.z-1;
	globBin = nbinsSph.x*(ibin.z*nbinsSph.y + ibin.y) + ibin.x;
      }
    }
    else if(cut.second >= 10){ //Theta bin cross with explicit bin index
      //Update theta bin
      ibin.y = cut.second-10;
      //Update global bin
      globBin = nbinsSph.x*(ibin.z*nbinsSph.y + ibin.y) + ibin.x;      
    }
  }
}

void pen_tallyKermaTrackLength::flush(){

  unsigned long nbins;

  if(activeCart){
    nbins = nbinsCart.x*nbinsCart.y*nbinsCart.z;
    for(unsigned long i = 0; i < nbins; ++i){
      if(cartesianLastHist[i] == 0){ continue;}  // Skip void counters
      double tmp = cartesianTmp[i];
      cartesian[i] += tmp;
      cartesian2[i] += tmp*tmp;
      cartesianTmp[i] = 0.0;
      cartesianLastHist[i] = 0;
    }
  }

  if(activeCyl){
    nbins = nbinsCyl.x*nbinsCyl.y*nbinsCyl.z;
    for(unsigned long i = 0; i < nbins; ++i){
      if(cylindricalLastHist[i] == 0){ continue;}  // Skip void counters
      double tmp = cylindricalTmp[i];
      cylindrical[i] += tmp;
      cylindrical2[i] += tmp*tmp;
      cylindricalTmp[i] = 0.0;
      cylindricalLastHist[i] = 0;
    }
  }
  
  if(activeSph){
    nbins = nbinsSph.x*nbinsSph.y*nbinsSph.z;
    for(unsigned long i = 0; i < nbins; ++i){
      if(sphericalLastHist[i] == 0){ continue;}  // Skip void counters
      double tmp = sphericalTmp[i];
      spherical[i] += tmp;
      spherical2[i] += tmp*tmp;
      sphericalTmp[i] = 0.0;
      sphericalLastHist[i] = 0;
    }
  }  
}


void pen_tallyKermaTrackLength::tally_jump(const unsigned long long /*nhist*/,
					   const pen_KPAR kpar,
					   const pen_particleState& state,
					   const double /*ds*/){
  //Check kpar
  if(kpar == kparTrig){
    //Save last position
    lastPos.x = state.X;
    lastPos.y = state.Y;
    lastPos.z = state.Z;
  }
}

void pen_tallyKermaTrackLength::tally_step(const unsigned long long nhist,
					   const pen_KPAR kpar,
					   const pen_particleState& state,
					   const tally_StepData& stepData){

  using namespace pen_tally_KTL;

  if(kpar != kparTrig)
    return;
  if(!activeMat[stepData.originMAT])
    return;
  //Check energy to avoid count energies in the range (emin-ebin,emin)
  if(state.E < grid.EL || state.E >= grid.EU){
    return;
  }

  //Get energy interval
  long int KE;
  double XEL, XE, XEK;
  grid.getInterval(state.E,KE,XEL,XE,XEK);

  double logmuenKe = muen[stepData.originMAT][KE];
  double logmuenKeNext = muen[stepData.originMAT][KE+1];  
  double emuen = exp(logmuenKe + (logmuenKeNext - logmuenKe)*XEK);

  vect3d dp;
  dp.x = state.U*stepData.dsef;
  dp.y = state.V*stepData.dsef;
  dp.z = state.W*stepData.dsef;

  vect3d endPos = lastPos;
  endPos.x += dp.x;
  endPos.y += dp.y;
  endPos.z += dp.z;

  //Score paths
  if(activeCart){
    vect3d dir(state.U,state.V,state.W);
    kermaTrackLengthCart(nhist,state.WGHT,state.E,emuen,stepData.dsef,
			 lastPos,dir);
  }
  
  if(activeCyl || activeSph){

    vect3d pos1Cyl, pos2Cyl;
    vect3d pos1Sph, pos2Sph;
    cart2CylSph(lastPos,pos1Cyl,pos1Sph);
    cart2CylSph(endPos ,pos2Cyl,pos2Sph);
      
    if(activeCyl){
      kermaTrackLengthCyl(nhist,state.WGHT,state.E,emuen,stepData.dsef,
			  lastPos,dp,pos1Cyl,pos2Cyl);
    }
  
    if(activeSph){
      kermaTrackLengthSph(nhist,state.WGHT,state.E,emuen,stepData.dsef,
			  lastPos,dp,pos1Sph,pos2Sph);
    }
  }
}
  
int pen_tallyKermaTrackLength::configure(const wrapper_geometry& /*geometry*/,
				      const abc_material* const /*materials*/[constants::MAXMAT],
				      const pen_parserSection& config,
				      const unsigned verbose){
  int err;

 //Energies
 //****************
 
  // Minimum energy
  //***************    
  double emin, emax;
  err = config.read("emin", emin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	     "'emin' in configuration. Double expected\n");
    }
    return -1;
  }

  // Maximum energy
  //***************    
    
  err = config.read("emax", emax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	     "'emax' in configuration. Double expected\n");
    }
    return -2;
  }

  if(emin >= emax){
    if(verbose > 0){
      printf("pen_tallyKermaTrackLength:configure: Error: 'emin' must "
	     "be greater than 'emax'.\n"
	     "           emin: %12.4E\n"
	     "           emax: %12.4E\n",
	     emin,emax);
    }
    return -3;
  }
    
  if(verbose > 1){
    printf("Energy limits [Emin,Emax] (eV) and no. bins:\n");
    printf(" %12.5E %12.5E %lu \n\n",
	   emin,emax,static_cast<unsigned long>(nbinmax));
  }

  // Set energy grid
  grid.init(emin,emax);
  
  // Absorption coefficient data filenames
  //***************************************
  std::vector<double> muenData;
  std::vector<double> EData;
  std::vector<double> A;
  std::vector<double> B;
  std::vector<double> C;
  std::vector<double> D;
  EData.reserve(nbinmax);
  muenData.reserve(nbinmax);
  A.reserve(nbinmax);  B.reserve(nbinmax);
  C.reserve(nbinmax);  D.reserve(nbinmax);

  bool someMat = false;
  activeMat[0] = false;
  for(unsigned imat = 1; imat < constants::MAXMAT; ++imat){

    activeMat[imat] = false;
    
    std::string key("dataFiles/");
    key += std::to_string(imat);
    std::string filename;
    err = config.read(key, filename);
    if(err == INTDATA_SUCCESS){

      FILE* fin = nullptr;
      fin = ::fopen(filename.c_str(),"r");
      if(fin == nullptr){
	if(verbose > 0){
	  printf("pen_tallyKermaTrackLength:configure: Error: Unable to open "
		 "data file '%s'\n",filename.c_str());
	}
	return -8;
      }

      char line[1000];
      unsigned long actualLine = 0;
      unsigned long nread = 0;
      while(pen_getLine(fin,1000,line,nread) == 0){
	actualLine += nread;
	double E,muenAux;
	if(sscanf(line,"%lf %lf",&E,&muenAux) != 2){
	  if(verbose > 0){
	    printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
		   "Energy and muen from of bins at data file '%s'\n"
		   "       line: %s\n"
		   "line number: %lu",filename.c_str(),line,actualLine);
	  }
	  fclose(fin);
	  return -9;
	}
	//Convert to eV
	//E *= 1000.0;
	if(exp(EData.back()) >= E || E <= 0.0 || muenAux <= 0.0){
	  if(verbose > 0){
	    printf("pen_tallyKermaTrackLength:configure: Error: Energy bins "
		   "must be in increasing order and both vaules positive\n"
		   "   line: %s\n",line);
	  }
	  fclose(fin);
	  return -9;	  
	}
	EData.push_back(log(E));
	muenData.push_back(log(muenAux));
      }
      fclose(fin);

      size_t auxNbins = EData.size();
      
      if(verbose > 1){
	printf(" Coeficient bins for material %u: %lu\n",
	       imat,static_cast<unsigned long>(auxNbins));
      }

      //Check if there are bins to ajust
      if(auxNbins <= 10){
	if(verbose > 0)
	  printf("pen_tallyKermaTrackLength:configure: Error: Insuficient "
		 "bins at data file '%s'\n",filename.c_str());
	return -10;	
      }

      //Check if the specified energy rank is in the provided data 
      if(grid.EL < exp(EData[0]) || grid.EU > exp(EData.back())){
	if(verbose > 0)
	  printf("pen_tallyKermaTrackLength:configure: Error: energy rank "
		 "not included in the provided data file '%s'\n",
		 filename.c_str());
	return -11;		
      }
      
      A.resize(auxNbins); B.resize(auxNbins);
      C.resize(auxNbins); D.resize(auxNbins);
      SPLINE(EData.data(),muenData.data(),A.data(),B.data(),
	     C.data(),D.data(),0.0,0.0,auxNbins);

      size_t nextSplin = 0;
      FILE* fmuen = nullptr;
      if(verbose > 2){
	printf("#Muen data of '%s':\n",filename.c_str());
	printf("# %12s  %12s\n","Energy (eV)","mu_en(cm^2/g)\n");
	std::string filenameOut("processed-muen-");
	filenameOut += filename;
	fmuen = fopen(filenameOut.c_str(),"w");
	if(fmuen != nullptr)
	  fprintf(fmuen,"# %11s     %12s\n","E(eV)","mu_en g/cm^2");
      }
      double* pmuen = muen[imat];
      for(size_t j = 0; j < nbinmax; ++j){
	double nextE = grid.DLEMP[j];
	while((EData[nextSplin] > nextE || EData[nextSplin+1] < nextE) &&
	      nextSplin < auxNbins-1){
	  ++nextSplin;
	}
	pmuen[j] = A[nextSplin] +
	  nextE*(B[nextSplin] + nextE*(C[nextSplin] + nextE*D[nextSplin]));
	if(pmuen[j] <= 0.0){
	  double dmu = muenData[nextSplin+1]-muenData[nextSplin];
	  double dE = EData[nextSplin+1]-EData[nextSplin];
	  double m = dmu/dE;	  
	  pmuen[j] = muenData[nextSplin]+ m*(nextE-EData[nextSplin]);
	}
	if(fmuen != nullptr)
	  fprintf(fmuen," %12.4E  %12.4E\n",grid.ET[j],exp(pmuen[j]));
      }
      if(fmuen != nullptr)
	fclose(fmuen);      
      activeMat[imat] = true;
      someMat = true;
    }

  }

  if(!someMat){
    if(verbose > 0)
      printf("pen_tallyKermaTrackLength:configure: Error: No material "
	     "data file provided\n");
    return -12;		      
  }

  // Get meshes information
  //*************************

  //Cartesian
  if(config.isSection(std::string("cartesian"))){
    activeCart = true;
    int auxInt;
    //Number of bins
    err = config.read("cartesian/nx", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/nx' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsCart.x = auxInt;
    
    err = config.read("cartesian/ny", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/ny' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsCart.y = auxInt;
    
    err = config.read("cartesian/nz", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/nz' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsCart.z = auxInt;
    
    if(verbose > 1){
      printf("Cartesian number of bins (nx,ny,nz):"
	     " %ld %ld %ld\n",nbinsCart.x,nbinsCart.y,nbinsCart.z);
    }
    if(nbinsCart.x <= 0 || nbinsCart.y <= 0 || nbinsCart.z <= 0  ||
       nbinsCart.x > meshAxeMax || nbinsCart.y > meshAxeMax ||
       nbinsCart.z > meshAxeMax){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: number of "
	       "cartesian bins must be greater than zero and lesser "
	       "than %ld\n",meshAxeMax);
      }
      return -13;      
    }

    //Limits
    err = config.read("cartesian/xmin", minsCart.x);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/xmin' in configuration. Double expected\n");
      }
      return -14;
    }
    err = config.read("cartesian/xmax", maxsCart.x);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/xmax' in configuration. Double expected\n");
      }
      return -14;
    }
    err = config.read("cartesian/ymin", minsCart.y);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/ymin' in configuration. Double expected\n");
      }
      return -14;
    }
    err = config.read("cartesian/ymax", maxsCart.y);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/ymax' in configuration. Double expected\n");
      }
      return -14;
    }
    err = config.read("cartesian/zmin", minsCart.z);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/zmin' in configuration. Double expected\n");
      }
      return -14;
    }
    err = config.read("cartesian/zmax", maxsCart.z);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cartesian/zmax' in configuration. Double expected\n");
      }
      return -14;
    }

    if(verbose > 1){
      printf("Cartesian ranges (cm):\n"
	     "   x: %14.4E - %14.4E\n"
	     "   y: %14.4E - %14.4E\n"
	     "   z: %14.4E - %14.4E\n",
	     minsCart.x,maxsCart.x,minsCart.y,maxsCart.y,minsCart.z,maxsCart.z);
    }

    if(minsCart.x >= maxsCart.x ||
       minsCart.y >= maxsCart.y ||
       minsCart.z >= maxsCart.z){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Cartesian mesh "
	       "ranges must be positive\n");
      }
      return -14;      
    }

    meshSizeCart.x = maxsCart.x-minsCart.x;
    meshSizeCart.y = maxsCart.y-minsCart.y;
    meshSizeCart.z = maxsCart.z-minsCart.z;

    dbinCart.x = (meshSizeCart.x)/static_cast<double>(nbinsCart.x);
    dbinCart.y = (meshSizeCart.y)/static_cast<double>(nbinsCart.y);
    dbinCart.z = (meshSizeCart.z)/static_cast<double>(nbinsCart.z);

    if(verbose > 1){
      printf("Cartesian mesh sizes (cm):\n"
	     "   x: %14.4E\n"
	     "   y: %14.4E\n"
	     "   z: %14.4E\n",
	     dbinCart.x,dbinCart.y,dbinCart.z);
    }
    
    for(long int i = 0; i < nbinsCart.x+1; ++i)      
      xPlanes[i] = minsCart.x + static_cast<double>(i)*dbinCart.x;
    for(long int i = 0; i < nbinsCart.y+1; ++i)      
      yPlanes[i] = minsCart.y + static_cast<double>(i)*dbinCart.y;
    for(long int i = 0; i < nbinsCart.z+1; ++i)      
      zPlanes[i] = minsCart.z + static_cast<double>(i)*dbinCart.z;
  }

  //Cylindrical
  if(config.isSection(std::string("cylindrical"))){
    activeCyl = true;
    //Number of bins
    int auxInt;
    err = config.read("cylindrical/nr", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cylindrical/nr' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsCyl.x = auxInt;
    
    err = config.read("cylindrical/nphi", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cylindrical/ntheta' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsCyl.y = auxInt;
      
    err = config.read("cylindrical/nz", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cylindrical/nz' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsCyl.z = auxInt;

    if(verbose > 1){
      printf("Cylindrical number of bins (nr,ntheta,nz):"
	     " %ld %ld %ld\n",nbinsCyl.x,nbinsCyl.y,nbinsCyl.z);
    }
    if(nbinsCyl.x <= 0 || nbinsCyl.y <= 0 || nbinsCyl.z <= 0 ||
       nbinsCyl.x > meshAxeMax || nbinsCyl.y > meshAxeMax ||
       nbinsCyl.z > meshAxeMax){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: number of "
	       "cylindrical bins must be greater than zero and lesser "
	       "than %ld\n",meshAxeMax);
      }
      return -13;      
    }

    //Limits
    err = config.read("cylindrical/rmax", radCyl);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cylindrical/rmax' in configuration. Double expected\n");
      }
      return -14;
    }

    err = config.read("cylindrical/rmin", radCylmin);
    if(err != INTDATA_SUCCESS){
      radCylmin = 0.0;
    }
    long int effrbins = nbinsCyl.x;
    if(radCylmin < 1.0e-8)
      radCylmin = 0.0;
    else
      ++nbinsCyl.x; //Use an extra fake bin for the interval [0.0,radmin)

    
    err = config.read("cylindrical/zmin", zminCyl);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cylindrical/zmin' in configuration. Double expected\n");
      }
      return -14;
    }
    err = config.read("cylindrical/zmax", zmaxCyl);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'cylindrical/zmax' in configuration. Double expected\n");
      }
      return -14;
    }

    if(verbose > 1){
      printf("Cylindrical ranges:\n"
	     "    r: %14.4E - %14.4E (cm)\n"
	     "  phi: %14f - %14f (DEG)\n"
	     "    z: %14.4E - %14.4E (cm)\n",
	     radCylmin,radCyl,0.0,360.0,zminCyl,zmaxCyl);
    }

    if(radCyl < 1.0e-8 || radCylmin >= radCyl ||
       zminCyl >= zmaxCyl){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Cylindrical mesh "
	       "ranges must be positive\n");
      }
      return -14;      
    }

    heightCyl = zmaxCyl - zminCyl;
    
    dbinCyl.x = (radCyl-radCylmin)/static_cast<double>(effrbins);
    dbinCyl.y = 2.0*M_PI/static_cast<double>(nbinsCyl.y);
    dbinCyl.z = heightCyl/static_cast<double>(nbinsCyl.z);

    if(verbose > 1){
      printf("Cylindrical mesh bin sizes:\n"
	     "   r: %14.4E (cm)\n"
	     " phi: %14.4E (RAD)\n"
	     "   z: %14.4E (cm)\n",
	     dbinCyl.x,dbinCyl.y,dbinCyl.z);
    }

    if(radCylmin > 1.0e-9){
      rPlanesCyl[0] = radCylmin;
    }
    else{
      rPlanesCyl[0] = dbinCyl.x;
    }
    double rhomin = rPlanesCyl[0];
    for(long int i = 1; i < nbinsCyl.x; ++i)
      rPlanesCyl[i] = rhomin+static_cast<double>(i)*dbinCyl.x;
    for(long int i = 0; i < nbinsCyl.y+1; ++i)
      phiPlanesCyl[i] = static_cast<double>(i)*dbinCyl.y;
    for(long int i = 0; i < nbinsCyl.z+1; ++i)
      zPlanesCyl[i] = zminCyl + static_cast<double>(i)*dbinCyl.z;
    
    nrphiCyl =
      static_cast<long int>(nbinsCyl.x)*static_cast<long int>(nbinsCyl.y);
  }

  //Spherical
  if(config.isSection(std::string("spherical"))){
    activeSph = true;
    //Number of bins
    int auxInt;
    err = config.read("spherical/nr", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'spherical/nr' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsSph.x = auxInt;
    
    err = config.read("spherical/ntheta", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'spherical/ntheta' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsSph.y = auxInt;
    
    err = config.read("spherical/nphi", auxInt);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'spherical/nphi' in configuration. Integrer expected\n");
      }
      return -13;
    }
    nbinsSph.z = auxInt;

    if(verbose > 1){
      printf("Spherical number of bins (nr,ntheta,nphi):"
	     " %ld %ld %ld\n",nbinsSph.x,nbinsSph.y,nbinsSph.z);
    }
    if(nbinsSph.x <= 0 || nbinsSph.y <= 0 || nbinsSph.z <= 0 ||
       nbinsSph.x > meshAxeMax || nbinsSph.y > meshAxeMax ||
       nbinsSph.z > meshAxeMax){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: number of "
	       "spherical bins must be greater than zero and lesser "
	       "than %ld\n",meshAxeMax);
      }
      return -13;      
    }

    //Limits
    err = config.read("spherical/rmax", radSph);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Unable to read "
	       "'spherical/rmax' in configuration. Double expected\n");
      }
      return -14;
    }
    err = config.read("spherical/rmin", radSphmin);
    if(err != INTDATA_SUCCESS){
      radSphmin = 0.0;
    }

    long int effrbins = nbinsSph.x;
    if(radSphmin < 1.0e-8)
      radSphmin = 0.0;
    else
      ++nbinsSph.x; //Use an extra fake bin for the interval [0.0,radmin)

    if(verbose > 1){
      printf("Spherical ranges:\n"
	     "    r: %14.4E - %14.4E (cm)\n"
	     "theta: %14f - %14f (DEG)\n"
	     "  phi: %14f - %14f (DEG)\n",
	     radSphmin,radSph,0.0,180.0,0.0,360.0);
    }

    if(radSph < 1.0e-8 || radSphmin >= radSph){
      if(verbose > 0){
	printf("pen_tallyKermaTrackLength:configure: Error: Spherical mesh "
	       "radius must be greater than zero and range positive.\n");
      }
      return -14;      
    }

    dbinSph.x = (radSph-radSphmin)/static_cast<double>(effrbins);
    dbinSph.y = M_PI/static_cast<double>(nbinsSph.y);
    dbinSph.z = 2.0*M_PI/static_cast<double>(nbinsSph.z);

    if(verbose > 1){
      printf("Spherical bin sizes:\n"
	     "    r: %14.4E (cm)\n"
	     "theta: %14.4E (RAD)\n"
	     "  phi: %14.4E (RAD)\n",
	     dbinSph.x,dbinSph.y,dbinSph.z);
    }

    if(radSphmin > 1.0e-9){
      rPlanesSph[0] = radSphmin;
    }
    else{
      rPlanesSph[0] = dbinSph.x;
    }
    double radmin = rPlanesSph[0];
    for(long int i = 1; i < nbinsSph.x; ++i)
      rPlanesSph[i] = radmin+static_cast<double>(i)*dbinSph.x;
    for(long int i = 0; i < nbinsSph.y; ++i)      
      thetaPlanesSph[i] = static_cast<double>(i+1)*dbinSph.y;
    for(long int i = 0; i < nbinsSph.z+1; ++i)
      phiPlanesSph[i] = static_cast<double>(i)*dbinSph.z;

    if(nbinsSph.y % 2){
      halfnThetaCones = nbinsSph.y/2;
      planeThetaConeIndex = -1;
    }
    else{
      halfnThetaCones = nbinsSph.y/2;
      planeThetaConeIndex = nbinsSph.y/2-1;
    }
	
    nrthetaSph =
      static_cast<long int>(nbinsSph.x)*static_cast<long int>(nbinsSph.y);    
  }

  //Check if some meash has been selected
  if(!activeCart && !activeCyl && !activeSph){
    if(verbose > 0)
      printf("pen_tallyKermaTrackLength:configure: Error: No mesh "
	     "selected\n");    
    return 15; 
  }

  //Allocate required memory
  const double B2MB = 1.0/(1024.0*1024.0);
  if(activeCart){
    unsigned long long nxy = nbinsCart.x*nbinsCart.y;
    unsigned long long nbins =
      static_cast<unsigned long long>(nbinsCart.z)*nxy;
    if(verbose > 1)
      printf("Number of cartesian bins: %llu (%6.2f MB)\n",
	     nbins,static_cast<double>(4*nbins*sizeof(double))*B2MB);
    cartesian         = static_cast<double*>(malloc(nbins*sizeof(double)));
    cartesian2        = static_cast<double*>(malloc(nbins*sizeof(double)));
    cartesianTmp      = static_cast<double*>(malloc(nbins*sizeof(double)));
    cartesianLastHist = static_cast<unsigned long long*>(malloc(nbins*sizeof(unsigned long long)));

    for(unsigned long int i = 0; i < nbins; ++i) cartesian[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) cartesian2[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) cartesianTmp[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) cartesianLastHist[i] = 0;
    
    //Precalcuate element volume
    volumeCart = dbinCart.x*dbinCart.y*dbinCart.z;

    if(verbose > 1)
      printf("Cartesian element volume (cm^3): %12.4E\n",volumeCart);
    
    dump.toDump(cartesian,nbins);
    dump.toDump(cartesian2,nbins);
  }

  if(activeCyl){
    unsigned long long nbins =
      static_cast<unsigned long long>(nbinsCyl.z)*nrphiCyl;
    if(verbose > 1)
      printf("Number of cylindrical bins: %llu (%6.2f MB)\n",
	     nbins,
	     static_cast<double>((4*nbins+meshAxeMax)*sizeof(double))*B2MB);
    cylindrical         = static_cast<double*>(malloc(nbins*sizeof(double)));
    cylindrical2        = static_cast<double*>(malloc(nbins*sizeof(double)));
    cylindricalTmp      = static_cast<double*>(malloc(nbins*sizeof(double)));
    cylindricalLastHist = static_cast<unsigned long long*>(malloc(nbins*sizeof(unsigned long long)));

    for(unsigned long int i = 0; i < nbins; ++i) cylindrical[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) cylindrical2[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) cylindricalTmp[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) cylindricalLastHist[i] = 0;
    
    //Precalculate elements volumes
    double dyz = dbinCyl.y*dbinCyl.z;
    for(long int i = 0; i < nbinsCyl.x; ++i){
      double rho1;
      if(i == 0)
	rho1 = 0.0;
      else
	rho1 = rPlanesCyl[i-1];
      double rho2 = rPlanesCyl[i];

      pvolumesCyl[i] = (rho2*rho2 - rho1*rho1)/2.0 * dyz;
    }
    
    dump.toDump(cylindrical,nbins);
    dump.toDump(cylindrical2,nbins);
  }  

  if(activeSph){
    unsigned long long nbins =
      static_cast<unsigned long long>(nbinsSph.z)*nrthetaSph;
    if(verbose > 1)
      printf("Number of spherical bins: %llu (%6.2f MB)\n",
	     nbins,
	     static_cast<double>((4*nbins+nrthetaSph)*sizeof(double))*B2MB);
    spherical         = static_cast<double*>(malloc(nbins*sizeof(double)));
    spherical2        = static_cast<double*>(malloc(nbins*sizeof(double)));
    sphericalTmp      = static_cast<double*>(malloc(nbins*sizeof(double)));
    sphericalLastHist = static_cast<unsigned long long*>(malloc(nbins*sizeof(unsigned long long)));

    for(unsigned long int i = 0; i < nbins; ++i) spherical[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) spherical2[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) sphericalTmp[i] = 0.0;
    for(unsigned long int i = 0; i < nbins; ++i) sphericalLastHist[i] = 0;  
    
    pvolumesSph  = static_cast<double*>(malloc(nrthetaSph*sizeof(double)));
    
    //Precalculate elements volumes
    long long int t = 0;
    for(long int j = 0; j < nbinsSph.y; ++j){
      double theta1 = static_cast<double>(j)*dbinSph.y;
      double theta2 = theta1 + dbinSph.y;    
      for(long int i = 0; i < nbinsSph.x; ++i){
	double rad1;
	if(i == 0)
	  rad1 = 0.0;
	else
	  rad1 = rPlanesSph[i-1];
	double rad2 = rPlanesSph[i];
	pvolumesSph[t++] =
	  dbinSph.z*(cos(theta1)-cos(theta2))*(pow(rad2,3)-pow(rad1,3))/3.0;
      }
    }
    
    dump.toDump(spherical,nbins);
    dump.toDump(spherical2,nbins);
  }  
  
  return 0;
}

pen_tallyKermaTrackLength::~pen_tallyKermaTrackLength(){
  if(activeCart){
    free(cartesian);
    free(cartesian2);
    free(cartesianTmp);
    free(cartesianLastHist);
  }
  if(activeCyl){
    free(cylindrical);
    free(cylindrical2);
    free(cylindricalTmp);
    free(cylindricalLastHist);
  }
  if(activeSph){
    free(spherical);free(spherical2);free(sphericalTmp);
    free(sphericalLastHist);free(pvolumesSph);
  }
}

void pen_tallyKermaTrackLength::saveData(const unsigned long long nhist) const{

  const char* filenameCart = "kermaTrackLength-cartesian.dat";
  const char* filenameCyl = "kermaTrackLength-cylindrical.dat";
  const char* filenameSph = "kermaTrackLength-spherical.dat";
    
  double invn = 1.0/static_cast<double>(nhist);

  // Cartesian mesh
  if(activeCart){
    FILE* fout = nullptr;
    fout = fopen(filenameCart,"w");
    if(fout == nullptr){
      printf(" *********************************************\n");
      printf(" pen_tallyKermaTrackLength:saveData:ERROR: cannot "
	     "open output data '%s'\n",filenameCart);
      printf(" *********************************************\n");
      return;      
    }

    fprintf(fout,"#------------------------------------------------------------\n");
    fprintf(fout,"# PenRed: Cartesian kerma track length\n");
    fprintf(fout,"#\n");
    fprintf(fout,"# Its units are eV/(g*history).\n");
    fprintf(fout,"#\n");
    fprintf(fout,"# Number of bins (x,y,z):\n");
    fprintf(fout,"#  %ld %ld %ld\n",nbinsCart.x,nbinsCart.y,nbinsCart.z);
    fprintf(fout,"# Bin sizes (dx,dy,dz):\n");
    fprintf(fout,"#  %14.4E %14.4E %14.4E\n",dbinCart.x,dbinCart.y,dbinCart.z);
    fprintf(fout,"#\n");
    fprintf(fout,"#   x   |   y   |   z   |    value    |   2*sigma  \n");

    unsigned long t = 0;
    for(long int k = 0; k < nbinsCart.z; ++k){
      for(long int j = 0; j < nbinsCart.y; ++j){
	for(long int i = 0; i < nbinsCart.x; ++i){
	  double q = cartesian[t]*invn;
	  double sigma = (cartesian2[t]*invn - q*q)*invn;
	  sigma = sqrt((sigma > 0.0 ? sigma : 0.0));

	  q /= volumeCart;
	  sigma /= volumeCart;
	  
	  fprintf(fout," %6ld  %6ld  %6ld  %12.5E %12.5E\n",
		  i,j,k,q,2.0*sigma);
	  ++t;
	}
	fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
    }
    fclose(fout);
  }
  
  // Cylindrical mesh
  if(activeCyl){
    FILE* fout = nullptr;
    fout = fopen(filenameCyl,"w");
    if(fout == nullptr){
      printf(" *********************************************\n");
      printf(" pen_tallyKermaTrackLength:saveData:ERROR: cannot "
	     "open output data '%s'\n",filenameCyl);
      printf(" *********************************************\n");
      return;      
    }
    
    long int modI = 0;
    if(radCylmin > 1.0e-9){
      modI = -1;
    }
    
    fprintf(fout,"#------------------------------------------------------------\n");
    fprintf(fout,"# PenRed: Cylindrical kerma track length\n");
    fprintf(fout,"#\n");
    fprintf(fout,"# Its units are eV/(g*history).\n");
    fprintf(fout,"#\n");
    fprintf(fout,"# Number of bins (r,phi,z):\n");
    fprintf(fout,"#  %ld %ld %ld\n",nbinsCyl.x+modI,nbinsCyl.y,nbinsCyl.z);
    fprintf(fout,"# Bin sizes (dr,dphi,dz):\n");
    fprintf(fout,"#  %14.4E %14.4E %14.4E\n",dbinCyl.x,dbinCyl.y,dbinCyl.z);
    fprintf(fout,"# Radius range (rmin,rmax):\n");
    fprintf(fout,"#  %14.4E %14.4E\n",radCylmin,radCyl);
    fprintf(fout,"#\n");
    fprintf(fout,"#  rho  |   rlow(cm)   |  phi  |  philow(DEG) "
	    "|   z   |   zlow(cm)   |    Dose     |   2*sigma  \n");

    unsigned long t = 0;
    for(long int k = 0; k < nbinsCyl.z; ++k){
      for(long int j = 0; j < nbinsCyl.y; ++j){
	for(long int i = 0; i < nbinsCyl.x; ++i){
	  double rlow;
	  if(i == 0){
	    if(radCylmin > 1.0e-9){
	      ++t;
	      continue; //Skip [0.0,rmin) bin
	    }
	    rlow = 0.0;
	  }else{
	    rlow = rPlanesCyl[i-1];
	  }	  
	  double q = cylindrical[t]*invn;
	  double sigma = (cylindrical2[t]*invn - q*q)*invn;
	  sigma = sqrt((sigma > 0.0 ? sigma : 0.0));

	  double volume = pvolumesCyl[i];
	  q /= volume;
	  sigma /= volume;
	  
	  fprintf(fout," %6ld  %12.4E  %6ld   %12.4E  %6ld   %12.4E  "
		  "%12.5E %12.5E\n",
		  i+modI,rlow,j,phiPlanesCyl[j],k,zPlanesCyl[k],
		  q,2.0*sigma);
	  ++t;
	}
	if(nbinsCyl.y > 1)fprintf(fout,"\n");
      }
      if(nbinsCyl.z > 1)fprintf(fout,"\n");
    }
    fclose(fout);
  }

  // Spherical mesh
  if(activeSph){
    FILE* fout = nullptr;
    fout = fopen(filenameSph,"w");
    if(fout == nullptr){
      printf(" *********************************************\n");
      printf(" pen_tallyKermaTrackLength:saveData:ERROR: cannot "
	     "open output data '%s'\n",filenameSph);
      printf(" *********************************************\n");
      return;      
    }

    long int modI = 0;
    if(radSphmin > 1.0e-9){
      modI = -1;
    }    
    
    fprintf(fout,"#------------------------------------------------------------\n");
    fprintf(fout,"# PenRed: Spherical kerma track length\n");
    fprintf(fout,"#\n");
    fprintf(fout,"# Its units are eV/(g*history).\n");
    fprintf(fout,"#\n");
    fprintf(fout,"# Number of bins (r,theta,phi):\n");
    fprintf(fout,"#  %ld %ld %ld\n",nbinsSph.x+modI,nbinsSph.y,nbinsSph.z);
    fprintf(fout,"# Bin sizes (dr,dtheta,dphi):\n");
    fprintf(fout,"#  %14.4E %14.4E %14.4E\n",dbinSph.x,dbinSph.y,dbinSph.z);
    fprintf(fout,"#\n");
    fprintf(fout,"#  rho  |   rlow(cm)   | theta | thetalow(DEG) "
	    "|  phi  |  philow(DEG) |    Dose     |   2*sigma  \n");

    for(long int i = 0; i < nbinsSph.x; ++i){
      double rlow;
      if(i == 0){
	if(radSphmin > 1.0e-9)
	  continue; //Skip [0.0,rmin) bin
	rlow = 0.0;
      }else{
	rlow = rPlanesSph[i-1];
      }
      for(long int j = 0; j < nbinsSph.y; ++j){
	long int itheta = j*nbinsSph.x;
	for(long int k = 0; k < nbinsSph.z; ++k){
	  long int ibin = nbinsSph.x*(k*nbinsSph.y + j) + i;
	  double q = spherical[ibin]*invn;
	  double sigma = (spherical2[ibin]*invn - q*q)*invn;
	  sigma = sqrt((sigma > 0.0 ? sigma : 0.0));
	  
	  double volume = pvolumesSph[itheta+i];
	  q /= volume;
	  sigma /= volume;
      
	  fprintf(fout," %6ld  %12.4E  %6ld   %12.4E   %6ld   %12.4E  "
		  "%12.5E %12.5E\n",
		  i+modI,rlow,j,thetaPlanesSph[j],k,phiPlanesSph[k],
		  q,2.0*sigma);
	}
	fprintf(fout,"\n");
      }
      fprintf(fout,"\n");
    }
    fclose(fout);
  }
  
}


int pen_tallyKermaTrackLength::sumTally(const pen_tallyKermaTrackLength& tally){

  if(activeCart != tally.activeCart ||
     activeCyl  != tally.activeCyl  ||
     activeSph  != tally.activeSph )
    return -1;

  if(activeCart){
    if(nbinsCart.x != tally.nbinsCart.x ||
       nbinsCart.y != tally.nbinsCart.y ||
       nbinsCart.z != tally.nbinsCart.z )
      return -2;
  }

  if(activeCyl){
    if(nbinsCyl.x != tally.nbinsCyl.x ||
       nbinsCyl.y != tally.nbinsCyl.y ||
       nbinsCyl.z != tally.nbinsCyl.z )
      return -2;
  }

  if(activeSph){
    if(nbinsSph.x != tally.nbinsSph.x ||
       nbinsSph.y != tally.nbinsSph.y ||
       nbinsSph.z != tally.nbinsSph.z )
      return -2;
  }

  //Cartesian
  if(activeCart){
    unsigned long int nbins = nbinsCart.x*nbinsCart.y*nbinsCart.z;
    for(unsigned long i = 0; i < nbins; ++i)
      cartesian[i] += tally.cartesian[i];
    for(unsigned long i = 0; i < nbins; ++i)
      cartesian2[i] += tally.cartesian2[i];    
  }

  //Cylindrical
  if(activeCyl){
    unsigned long int nbins = nbinsCyl.x*nbinsCyl.y*nbinsCyl.z;
    for(unsigned long i = 0; i < nbins; ++i)
      cylindrical[i] += tally.cylindrical[i];
    for(unsigned long i = 0; i < nbins; ++i)
      cylindrical2[i] += tally.cylindrical2[i];    
  }

  //Spherical
  if(activeSph){
    unsigned long int nbins = nbinsSph.x*nbinsSph.y*nbinsSph.z;
    for(unsigned long i = 0; i < nbins; ++i)
      spherical[i] += tally.spherical[i];
    for(unsigned long i = 0; i < nbins; ++i)
      spherical2[i] += tally.spherical2[i];    
  }
    
  return 0;
  
}

REGISTER_COMMON_TALLY(pen_tallyKermaTrackLength, KERMA_TRACK_LENGTH)
