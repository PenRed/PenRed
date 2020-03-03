
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


#include "monoenergetic.hh"


void monoenergetic::energySampling(double& energy, pen_rand& /*random*/) const{
  // Set energy
  energy = E;
}

int monoenergetic::configure(double& Emax, const pen_parserSection& config, const unsigned verbose){

  int err;

  err = config.read("energy",E);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("monoenergetic:configure:unable to read 'energy' in configuration. Real number expected.\n");
    }
    return -1;
  }

  if(verbose > 1)
    printf("Energy: %12.4E\n",E);
  
  Emax = E;
  
  return 0;
}

REGISTER_SAMPLER(monoenergetic, MONOENERGETIC)
