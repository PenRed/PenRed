
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


#include "cone_directionSampling.hh"

const double cone_directionSampling::deg2rad = constants::PI/180.0;

void cone_directionSampling::directionSampling(double dir[3], pen_rand& random) const{

  const double TWOPI = 2.0*constants::PI;
  
  double UT,VT,WT;
  double DF;
  double SUV;
  double UF,VF,WF;

  // Define a direction relative to the z-axis
  WT  = CAPER + (1.0-CAPER)*random.rand();
  DF  = TWOPI*random.rand();
  SUV = sqrt(1.0-WT*WT);
  UT  = SUV*cos(DF);
  VT  = SUV*sin(DF);
  // Rotate to the beam axis direction
  UF  = CPCT*UT-SPHI*VT+CPST*WT;
  VF  = SPCT*UT+CPHI*VT+SPST*WT;
  WF  =-STHE*UT+CTHE*WT;
  // Ensure normalisation
  double DXY  = UF*UF+VF*VF;
  double DXYZ = DXY+WF*WF;
  if(fabs(DXYZ-1.0) > 1.0e-14){
    double FNORM = 1.0/sqrt(DXYZ);
    dir[0] = FNORM*UF;
    dir[1] = FNORM*VF;
    dir[2] = FNORM*WF;
  }
  else{
    dir[0] = UF;
    dir[1] = VF;
    dir[2] = WF;    
  }
}

int cone_directionSampling::configure(const pen_parserSection& config, const unsigned verbose){

  int err;
  double theta,phi,alpha;
  //Store cosines (u,v,w)
  err = config.read("theta",theta);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("coneDirection:configure:unable to read 'theta' in configuration. Real number expected\n");
    }
    return -1;
  }
  err = config.read("phi",phi);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("coneDirection:configure:unable to read 'phi' in configuration. Real number expected\n");
    }
    return -1;
  }
  err = config.read("alpha",alpha);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("coneDirection:configure:unable to read 'alpha' in configuration. Real number expected\n");
    }
    return -1;
  }

  if(verbose > 1){
    printf("Theta: %12.4E DEG\n",theta);
    printf("Phi  : %12.4E DEG\n",phi);
    printf("Alpha: %12.4E DEG\n",alpha);
  }

  theta *= deg2rad;
  phi   *= deg2rad;
  alpha *= deg2rad;
  
  CPCT  = cos(phi)*cos(theta);
  CPST  = cos(phi)*sin(theta);
  SPCT  = sin(phi)*cos(theta);
  SPST  = sin(phi)*sin(theta);
  SPHI  = sin(phi);
  CPHI  = cos(phi);
  STHE  = sin(theta);
  CTHE  = cos(theta);
  CAPER = cos(alpha);
  
  return 0;
  
}

REGISTER_SAMPLER(cone_directionSampling, CONE)
