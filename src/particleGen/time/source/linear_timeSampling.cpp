
//
//
//    Copyright (C) 2026 Universitat Politècnica de València - UPV
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


#include "linear_timeSampling.hh"

void linear_timeSampling::timeSampling(double& time, pen_rand& random) const{
  time = tmin+random.rand()*dt;
}

int linear_timeSampling::configure(const pen_parserSection& config, const unsigned verbose){

  int err;

  //Get time window
  err = config.read("tmin",tmin);
  if(err != INTDATA_SUCCESS){
    tmin = 0.0;
  }
  
  err = config.read("tmax",tmax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("linearTime:configure:unable to read 'tmax' in configuration. Double expected\n");
    }
    return -1;    
  }

  if(tmin < 0.0 || tmax < 0.0 || tmax < tmin){
    if(verbose > 0){
      printf("linearTime:configure: Negative time intervals are not allowed.\n");
    }
    return -2;
  }

  //Calculate range
  dt = tmax - tmin;

  if(verbose > 1){
    printf("tmin (s) : %12.4E\n",tmin);
    printf("tmax (s) : %12.4E\n",tmax);
  }  
  return 0;
}

REGISTER_SAMPLER(linear_timeSampling, LINEAR)
