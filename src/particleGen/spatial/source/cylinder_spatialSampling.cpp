
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es               (Sandra Oliver Gil)
//    
//


#include "cylinder_spatialSampling.hh"

const double cylinder_spatialSampling::PI = 3.141592653589793;
const double cylinder_spatialSampling::PI2 = 2.0*cylinder_spatialSampling::PI;


void cylinder_spatialSampling::geoSampling(double pos[3], pen_rand& random) const{

  //Sample theta
  double theta = random.rand() * PI2;
  
  //Select cylinder shell region to be sampled
  double p = random.rand();
  if(p < p1){ //Top cover

    //Sample radious and height
    double r = rmax*std::sqrt(random.rand());
    pos[0] = r*std::cos(theta);
    pos[1] = r*std::sin(theta);
    pos[2] = dzmin05 + random.rand()*coverHeight;
    
  }else if(p < p2){ //Lateral shell
  
    double r = std::sqrt( rmin2 + dr2 * random.rand() );
  
    pos[0] = r*std::cos(theta);
    pos[1] = r*std::sin(theta);
    pos[2] = dzmin*random.rand()-dzmin05;
    
  }else{ //Bottom cover

    //Sample radious and height
    double r = rmax*std::sqrt(random.rand());
    pos[0] = r*std::cos(theta);
    pos[1] = r*std::sin(theta);
    pos[2] = -dzmin05 - random.rand()*coverHeight;
    
  }
  
}

int cylinder_spatialSampling::configure(const pen_parserSection& config, const unsigned verbose){

  int err;

  err = config.read("size/rmin",rmin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("cylinderSpatial:configure:unable to read 'size/rmin' in configuration. Double expected\n");
    }
    return -1;
  }

  err = config.read("size/rmax",rmax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("cylinderSpatial:configure:unable to read 'size/rmax' in configuration. Double expected\n");
    }
    return -1;
  }

  err = config.read("size/dzmax",dzmax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("cylinderSpatial:configure:unable to read 'size/dzmax' in configuration. Double expected\n");    
    }
    return -1;
  }  
  
  err = config.read("size/dzmin",dzmin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("cylinderSpatial:configure:unable to read 'size/dzmin' in configuration. Double expected\n");    
    }
    return -1;
  }  
  
  if(rmin < 0.0 || rmax < 0.0 || rmin > rmax ||
     dzmax < 0.0 || dzmin < 0.0 || dzmax < dzmin){
    if(verbose > 0){
      printf("cylinderSpatial:configure: Invalid size values. rmin, rmax and dz"
	     " must be postive with rmin <= rmax and dzmin <= dzmax\n");    
    }
    return -2;
  }

  //Precalculate values
  dzmin05 = dzmin/2.0;
  dzmax05 = dzmax/2.0;
  coverHeight = dzmax05-dzmin05;
  rmin2 = rmin*rmin;
  dr2 = rmax*rmax-rmin2;

  err = config.read("position/x",translation[0]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){      
      printf("cylinderSpatial:configure:unable to read 'position/x' in configuration. Double expected\n");
    }
    return -2;
  }

  err = config.read("position/y",translation[1]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("cylinderSpatial:configure:unable to read 'position/y' in configuration. Double expected\n");    
    }
    return -2;
  }

  err = config.read("position/z",translation[2]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("cylinderSpatial:configure:unable to read 'position/z' in configuration. Double expected\n");    
    }
    return -2;
  }    
  
  //Calculate the sampling probabilities in each region of the cylinder shell
  //
  // 1- Top cylindrical cover region
  // 2- Cylinder lateral shell region
  // 3- Bottom cylindrical cover region
  //
  //                 ___________
  //                |_____1_____|
  //                | |       | |
  //                | |       | |
  //                | |       | |
  //                |2|       |2|
  //                | |       | |
  //                | |       | |
  //                |_|_______|_|
  //                |_____3_____|
  //
  //
  double vCover = PI*rmax*rmax*coverHeight;
  double vLateral = PI*dzmin*(rmax*rmax - rmin*rmin);
  double vTotal = 2.0*vCover + vLateral;

  p1 = vCover/vTotal;
  p2 = p1 + vLateral/vTotal;
  p3 = 1.0;

  if(verbose > 1){
    printf("Cylinder center (x,y,z):\n %12.4E %12.4E %12.4E\n",translation[0],translation[1],translation[2]);
    printf("Cylinder radious (rmin,rmax):\n %12.4E %12.4E\n",rmin,rmax);    
    printf("Cylinder heights (dzmin,dzmax):\n %12.4E %12.4E\n",dzmin,dzmax);
    printf("Cylinder sampling cumulative probabilities:\n");
    printf("      -Top cover    : %E\n", p1);
    printf("      -Lateral shell: %E\n", p2);
    printf("      -Bottom cover : %E\n", p3);
  }
  
  return 0;
}

REGISTER_SAMPLER(cylinder_spatialSampling,CYLINDER)
