
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


#include "box_spatialSampling.hh"

void box_spatialSampling::geoSampling(double pos[3], pen_rand& random) const{

  pos[0] = dx*random.rand()-dx05;
  pos[1] = dy*random.rand()-dy05;
  pos[2] = dz*random.rand()-dz05;
  
}

int box_spatialSampling::configure(const pen_parserSection& config, const unsigned verbose){

  int err;

  err = config.read("size/dx",dx);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("boxSpatial:configure:unable to read 'size/dx' in configuration. Double expected\n");
    }
    return -1;
  }

  err = config.read("size/dy",dy);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("boxSpatial:configure:unable to read 'size/dy' in configuration. Double expected\n");
    }
    return -1;
  }

  err = config.read("size/dz",dz);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("boxSpatial:configure:unable to read 'size/dz' in configuration. Double expected\n");    
    }
    return -1;
  }  
  
  if(dx < 0.0 || dy < 0.0 || dz < 0.0){
    return -2;
  }

  dx05 = dx*0.5;
  dy05 = dy*0.5;
  dz05 = dz*0.5;

  err = config.read("position/x",translation[0]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){      
      printf("boxSpatial:configure:unable to read 'position/x' in configuration. Double expected\n");
    }
    return -2;
  }

  err = config.read("position/y",translation[1]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("boxSpatial:configure:unable to read 'position/y' in configuration. Double expected\n");    
    }
    return -2;
  }

  err = config.read("position/z",translation[2]);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("boxSpatial:configure:unable to read 'position/z' in configuration. Double expected\n");    
    }
    return -2;
  }    

  if(verbose > 1){
    printf("Box center (x,y,z):\n %12.4E %12.4E %12.4E\n",translation[0],translation[1],translation[2]);
    printf("Box size (dx,dy,dz):\n %12.4E %12.4E %12.4E\n",dx,dy,dz);    
  }
  
  return 0;
}

REGISTER_SAMPLER(box_spatialSampling,BOX)
