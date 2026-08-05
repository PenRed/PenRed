
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


#include "measure_energySampling.hh"


void measure_energySampling::energySampling(double& energy,
						 pen_rand& random) const{

  std::array<double, 1> sampled = sampler.samplePositions(random);
  energy = sampled[0];
}

int measure_energySampling::configure(double& Emax,
                                      const pen_parserSection& config,
                                      const unsigned verbose){

  int err;
  std::string filename;
  err = config.read("filename",filename);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("measure_energySampling:configure:Error: Unable to read "
             "'filename' in configuration. String expected\n");
    }
    return -1;
  }

  //Read measure file
  penred::measurements::results<double, 1> distrib;
  std::ifstream fin(filename, std::ifstream::in);
  if(!fin){
    if(verbose > 0){
      printf("measure_energySampling:configure:Error: Unable to "
             "open distribution file '%s'",
             filename.c_str());
    }
    return -1;
  }
  
  err = distrib.read(fin);
  if(err != 0){
    if(verbose > 0){
      printf("measure_energySampling:configure:Error: Unable to "
             "read distribution file '%s'.\n  Error code: %d\n",
             reader.filename.c_str(), err);
    }
    return -2;
  }

  
  //Get maximum energy
  Emax = distrib.readLimits()[0]->second;

  if(verbose > 1){
    printf("%s\n",distrib.stringifyInfo().c_str());
  }

  //Init sampling
  err = sampler.init(distrib.readData(),
		     distrib.readDimBins(),
		     distrib.readLimits());

  if(err != 0){
    if(verbose > 0){
      printf("measure_energySampling:configure:Error: Unable to "
	     "init energy distribution with walker's aliasing.\n  Error code: %d\n",
	     err);
    }
    return -3;
  }
  
  return 0;
}

REGISTER_SAMPLER(measure_energySampling, MEASURE_ENERGY)
