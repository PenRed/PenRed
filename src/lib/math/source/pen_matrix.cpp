
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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


#include "pen_matrix.hh"

void matmul3D(const double A[9], double B[3]){

  double C[3];
  
  for(int j = 0; j < 3; j++){
    C[j] = 0.0;
    for(int i = 0; i < 3; i++){
      C[j] += A[3*j+i]*B[i];
    }
  }
  for(int i = 0; i < 3; i++)
    B[i] = C[i];  
  
}

void matmul3D(const float A[9], float B[3]){

  float C[3];
  
  for(int j = 0; j < 3; j++){
    C[j] = 0.0;
    for(int i = 0; i < 3; i++){
      C[j] += A[3*j+i]*B[i];
    }
  }
  for(int i = 0; i < 3; i++)
    B[i] = C[i];  
  
}

void createRotationZYZ(const double omega, const double theta, const double phi, double rotation[9])
{
  //*******************************************************************
  //*    Computes the rotation matrix determined by the three Euler   *
  //*    angles as defined in the PENELOPE manual.                    *
  //*                                                                 *
  //*    Input:                                                       *
  //*      omega -> rotation angle around the z axis (rad)            *
  //*      theta -> rotation angle around the y axis (rad)            *
  //*      phi   -> rotation angle around the z axis (rad)            *
  //*    Output:                                                      *
  //*      rot -> 3x3 matrix of the 3D rotation                       *
  //*******************************************************************

  double somega,comega,stheta,ctheta,sphi,cphi;
  
  somega = sin(omega);
  comega = cos(omega);
  stheta = sin(theta);
  ctheta = cos(theta);
  sphi   = sin(phi);
  cphi   = cos(phi);
    
  rotation[0] = cphi*ctheta*comega-sphi*somega;
  rotation[1] = -cphi*ctheta*somega-sphi*comega;
  rotation[2] = cphi*stheta;
    
  rotation[3] = sphi*ctheta*comega+cphi*somega;
  rotation[4] = -sphi*ctheta*somega+cphi*comega;
  rotation[5] = sphi*stheta;
    
  rotation[6] = -stheta*comega;
  rotation[7] = stheta*somega;
  rotation[8] = ctheta;
}

void rollAlign(const double u,
	       const double v,
	       const double w,
	       const double omega,
	       double rotation[9]){
  //*******************************************************************
  //*    Computes the required rotation matrix to align the vector    *
  //*    (0,0,1) to the direction (u,v,w) with a initial Z rotation   *
  //*    specified by "omega".                                        *
  //*                                                                 *
  //*    Input:                                                       *
  //*      u -> direction x component                                 *
  //*      v -> direction y component                                 *
  //*      w -> direction z component                                 *
  //*      omega -> initial rotation in the Z axis                    *
  //*    Output:                                                      *
  //*      rot -> 3x3 matrix of the 3D rotation                       *
  //*******************************************************************

  //To solve:
  //
  // cphi * stheta = u
  // sphi * stheta = v
  // ctheta = w
  
  double somega,comega,stheta,ctheta,sphi,cphi;
  
  ctheta = w;
  
  if(fabs(ctheta) < 0.99999999999){
    double phi = atan2(v,u);
    cphi = cos(phi);
    sphi = sin(phi);
    
    stheta = fabs(cphi) > fabs(sphi) ? u/cphi : v/sphi;
    //sphi = v/stheta;
  }
  else{
    stheta = 0.0;
    ctheta = 1.0;
    cphi = 1.0;
    sphi = 0.0;
  }
  
  somega = sin(omega);
  comega = cos(omega);
    
  rotation[0] = cphi*ctheta*comega-sphi*somega;
  rotation[1] = -cphi*ctheta*somega-sphi*comega;
  rotation[2] = cphi*stheta;
    
  rotation[3] = sphi*ctheta*comega+cphi*somega;
  rotation[4] = -sphi*ctheta*somega+cphi*comega;
  rotation[5] = sphi*stheta;
    
  rotation[6] = -stheta*comega;
  rotation[7] = stheta*somega;
  rotation[8] = ctheta;
}

void rollAlignf(const float u,
		const float v,
		const float w,
		const float omega,
		float rotation[9]){
  //*******************************************************************
  //*    Computes the required rotation matrix to align the vector    *
  //*    (0,0,1) to the direction (u,v,w) with a initial Z rotation   *
  //*    specified by "omega".                                        *
  //*                                                                 *
  //*    Input:                                                       *
  //*      u -> direction x component                                 *
  //*      v -> direction y component                                 *
  //*      w -> direction z component                                 *
  //*      omega -> initial rotation in the Z axis                    *
  //*    Output:                                                      *
  //*      rot -> 3x3 matrix of the 3D rotation                       *
  //*******************************************************************

  //To solve:
  //
  // cphi * stheta = u
  // sphi * stheta = v
  // ctheta = w
  
  float somega,comega,stheta,ctheta,sphi,cphi;
  
  ctheta = w;
  
  if(fabs(ctheta) < 0.99999999999){
    float phi = atan2(v,u);
    cphi = cos(phi);
    sphi = sin(phi);
    
    stheta = fabs(cphi) > fabs(sphi) ? u/cphi : v/sphi;
    //sphi = v/stheta;
  }
  else{
    stheta = 0.0;
    ctheta = 1.0;
    cphi = 1.0;
    sphi = 0.0;
  }
  
  somega = sin(omega);
  comega = cos(omega);
    
  rotation[0] = cphi*ctheta*comega-sphi*somega;
  rotation[1] = -cphi*ctheta*somega-sphi*comega;
  rotation[2] = cphi*stheta;
    
  rotation[3] = sphi*ctheta*comega+cphi*somega;
  rotation[4] = -sphi*ctheta*somega+cphi*comega;
  rotation[5] = sphi*stheta;
    
  rotation[6] = -stheta*comega;
  rotation[7] = stheta*somega;
  rotation[8] = ctheta;
}
