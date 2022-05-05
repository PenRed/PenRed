
//
//
//    Copyright (C) 2019-2022 Universitat de València - UV
//    Copyright (C) 2019-2022 Universitat Politècnica de València - UPV
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


#include "sphereSection_directionSampling.hh"

const double sphereSection_directionSampling::deg2rad = constants::PI/180.0;

void sphereSection_directionSampling::directionSampling(double dir[3], pen_rand& random) const{

  double polarSin, phi;

  //Sample on vector (0.0,0.0,1.0)
  dir[2] = polarcos0+dpolarcos*random.rand();
  polarSin = sqrt(1.0-dir[2]*dir[2]);
  phi = phi0 + dphi*random.rand();
  dir[1] = polarSin*sin(phi);
  dir[0] = polarSin*cos(phi);

  //Rotate direction
  matmul3D(rotation,dir);
}

int sphereSection_directionSampling::configure(const pen_parserSection& config, const unsigned verbose){

  int err;
  //Store cosines (u,v,w)
  double dir[3];
  err = config.read("u",dir[0]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("sphereSection:configure:unable to read 'u' in "
	     "configuration. Double expected\n");
    }
    return -1;
  }
  err = config.read("v",dir[1]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("sphereSection:configure:unable to read 'v' in "
	     "configuration. Double expected\n");
    }
    return -1;
  }
  err = config.read("w",dir[2]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("sphereSection:configure:unable to read 'w' in "
	     "configuration. Double expected\n");
    }
    return -1;
  }
  
  //Normalize director
  double norm = sqrt(pow(dir[0],2) + pow(dir[1],2) + pow(dir[2],2));
  if(norm < 1.0e-16){
    if(verbose > 0){
      printf("sphereSection:configure: Invalid direction.\n");
    }
    return -2;
  }
  dir[0] /= norm;
  dir[1] /= norm;
  dir[2] /= norm;

  //Create rotation matrix to convert (0.0,0.0,1.0) to dir
  if(fabs(dir[0]) > 0.0 || fabs(dir[1]) > 0.0){
    createRotationZYZ(0.0,acos(dir[2]),atan2(dir[1],dir[0]),rotation);
  } else {
    createRotationZYZ(0.0,acos(dir[2]),0.0,rotation);    
  }

  //Get polar angle interval
  double theta0,theta1;

  err = config.read("theta0",theta0);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("sphereSection:configure:unable to read 'theta0' in "
	     "configuration. Double expected\n");
    }
    return -3;
  }
  err = config.read("theta1",theta1);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("sphereSection:configure:unable to read 'theta1' "
	     "in configuration. Double expected\n");
    }
    return -3;
  }

  theta0 *= deg2rad;
  theta1 *= deg2rad;

  
  //Store azimutal angle interval
  err = config.read("phi0",phi0);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("sphereSection:configure:unable to read 'phi0' in "
	     "configuration. Double expected\n");
    }
    return -4;
  }
  err = config.read("dphi",dphi);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("sphereSection:configure:unable to read 'dphi' "
	     "in configuration. Double expected\n");
    }
    return -4;
  }
  
  phi0 *= deg2rad;
  dphi *= deg2rad;

  if(phi0 < 0.0 || phi0 >= 2.0*constants::PI || dphi < 0.0 || dphi > 2.0*constants::PI){
    if(verbose > 0){
      printf("sphereSection:configure: Invalid azimutal interval.\n");
    }
    return -5;
  }

  //Calculate polar initial and interval cosinus 
  polarcos0 = cos(theta0);
  dpolarcos = cos(theta1)-polarcos0;

  if(verbose > 1){
    printf("Polar interval (rad)   : %12.4E - %12.4E\n", theta0, theta1);
    printf("Azimutal angle (rad)   : %12.4E\n", phi0);
    printf("Azimutal overture (rad): %12.4E\n", dphi);
  }
    
  return 0;
  
}

REGISTER_SAMPLER(sphereSection_directionSampling, SOLID_ANGLE)
