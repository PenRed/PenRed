
//
//
//    Copyright (C) 2026 Vicent Giménez Alventosa
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
//    
//


#include "circle_spatialSampling.hh"

void circle_spatialSampling::geoSampling(double pos[3], pen_rand& random) const{
  
  constexpr double pi = 3.141592653589793;
  constexpr double pi2 = 2.0*pi;
  
  double r = beamRad * sqrt(random.rand());
  double theta = random.rand() * pi2; 
  
  pos[0] = r * cos(theta);
  pos[1] = r * sin(theta);;
  pos[2] = 0.0;
  
}

int circle_spatialSampling::configure(const pen_parserSection& config, const unsigned verbose){

  int err;

  err = config.read("radius",beamRad);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("circleSpatial:configure:unable to read 'radius' in configuration. Double expected\n");
    }
    return -1;
  }

  beamRad = std::max(beamRad, 0.0);

  err = config.read("position/x",translation[0]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){      
      printf("circleSpatial:configure:unable to read 'position/x' in configuration. Double expected\n");
    }
    return -2;
  }

  err = config.read("position/y",translation[1]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("circleSpatial:configure:unable to read 'position/y' in configuration. Double expected\n");    
    }
    return -2;
  }

  err = config.read("position/z",translation[2]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("circleSpatial:configure:unable to read 'position/z' in configuration. Double expected\n");    
    }
    return -2;
  }

  //Read rotation
  bool toRotate = false;
  double omega, theta, phi;
  err = config.read("euler/omega",omega);
  if(err != INTDATA_SUCCESS){
    omega = 0.0;
  }else toRotate = true;

  err = config.read("euler/theta",theta);
  if(err != INTDATA_SUCCESS){
    theta = 0.0;
  }else toRotate = true;

  err = config.read("euler/phi",phi);
  if(err != INTDATA_SUCCESS){
    phi = 0.0;
  }else toRotate = true;

  //Create the rotation
  if(toRotate)
    setRotationZYZ(omega,theta,phi);

  if(verbose > 1){
    printf("Circle center (x,y,z) :\n %12.4E %12.4E %12.4E\n",
           translation[0],translation[1],translation[2]);
    printf("Circle radius         :\n %12.4E %12.4E %12.4E\n",beamRad);    
    printf("Rotation angles(Z,Y,Z):\n %12.4E %12.4E %12.4E\n",
           omega,theta,phi);
  }
  
  return 0;
}

REGISTER_SAMPLER(circle_spatialSampling,CIRCLE)
