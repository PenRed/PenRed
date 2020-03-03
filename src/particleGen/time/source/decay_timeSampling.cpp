
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


#include "decay_timeSampling.hh"

const double decay_timeSampling::LOG2 = log(2);

void decay_timeSampling::timeSampling(double& time, pen_rand& random) const{

  double rand = rand0+random.rand()*drand;
  time = -tau*log(1.0-rand);
}

int decay_timeSampling::configure(const pen_parserSection& config, const unsigned verbose){

  int err;
  //Store activity and half life
  double halfLife;
  err = config.read("halfLife",halfLife);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("decayTime:configure:unable to read 'halfLife' in configuration. Double expected\n");
    }
    return -1;    
  }

  if(halfLife <= 0.0)
    return -2;

  //Get time window
  double time0;
  double time1;
  err = config.read("time/time0",time0);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){    
      printf("decayTime:configure:unable to read 'time/time0' in configuration. Double expected\n");
    }
    return -3;    
  }
  
  err = config.read("time/time1",time1);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("decayTime:configure:unable to read 'time/time1' in configuration. Double expected\n");
    }
    return -3;    
  }

  if(time0 < 0.0 || time1 < 0.0 || time1 < time0){
    if(verbose > 0){
      printf("decayTime:configure: Negative time intervals are not allowed.\n");
    }
    return -4;
  }

  tau = halfLife/LOG2;

  //Calculate range of randoms
  rand0 = 1.0-exp(-time0/tau);
  drand = 1.0-exp(-time1/tau);

  if(drand >= 1.0)
    drand = 1.0-1.0e-16;

  drand -= rand0;

  if(verbose > 1){
    printf("T1/2 (s): %12.4E\n",halfLife);
    printf("tau (s) : %12.4E\n",tau);
  }
  
  return 0;
}

REGISTER_SAMPLER(decay_timeSampling, DECAY)
