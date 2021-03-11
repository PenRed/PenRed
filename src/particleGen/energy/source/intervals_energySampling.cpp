
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


#include "intervals_energySampling.hh"


void intervals_energySampling::energySampling(double& energy, pen_rand& random) const{

  // Get random number
  double rand = random.rand();

  // Get interval
  unsigned interval = seeki(cummulative,rand,nIntervals);

  // Calculate sampled energy
  energy = energies[interval]+(rand-cummulative[interval])*dE[interval];
}

int intervals_energySampling::configure(double& Emax, const pen_parserSection& config, const unsigned verbose){

  int err;
  int nintervals;
  err = config.read("nintervals",nintervals);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("energy_intervals:configure:unable to read 'nintervals' in configuration. Integer expected\n");
    }
    return -1;
  }

  nIntervals = unsigned(nintervals);
    
  if(nIntervals > maxIntervals){
    if(verbose > 0){
      printf("intervals_energySampling: Invalid number of intervals %d, the maximum is %d.",nIntervals,maxIntervals);
    }
    return -5;
  }
  
  //Get configuration arrays
  pen_parserArray configElow;
  pen_parserArray configEtop;
  pen_parserArray configCummul;
  err = config.read("lowE",configElow);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("energy_intervals:configure:unable to read 'lowE' in configuration. Array expected.\n");
    }
    return -2;
  }
  err = config.read("topE",configEtop);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("energy_intervals:configure:unable to read 'topE' in configuration. Array expected.\n");
    }
    return -2;
  }
  
  err = config.read("probabilities",configCummul);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("energy_intervals:configure:unable to read 'cummulative' in configuration. Array expected\n");
    }
    return -2;
  }

  // Store data
  cummulative[0] = 0.0; //Set lower cummulative prob to 0
  double Etop[maxIntervals];
  for(unsigned i = 0; i < nIntervals; i++){
    int err1,err2;
    err1 = configElow[i].read(energies[i]);
    err2 = configEtop[i].read(Etop[i]);
    if(err1 != INTDATA_SUCCESS || err2 != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("energy_intervals:configure:unable to read energy interval in configuration at position %d. Doubles expected\n",i);
      }
      return -3;
    }
    err = configCummul[i].read(cummulative[i+1]);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("energy_intervals:configure:unable to read 'cummulative' in configuration at position %d. Double expected\n",i);
      }
      return -3;
    }

    if(energies[i] < 0.0 || Etop[i] < 0.0 || energies[i] > Etop[i] || cummulative[i+1] <= 0.0){
      if(verbose > 0){
	printf("energy_intervals:configure: Invalid energy and probability interval at position %d\n",i);
      }
      return i;
    }
    
    cummulative[i+1] += cummulative[i];
  }
  
  //Check probability
  if(cummulative[nIntervals] <= 0.0){
    if(verbose > 0){
      printf("energy_intervals:configure: Null probability\n");
    }
    return -4;
  }
  
  //Normalize prob
  for(unsigned j = 1; j <= nIntervals; j++){
    cummulative[j] /= cummulative[nIntervals];
  }

  //Prepare array for sampling
  for(unsigned j = 0; j < nIntervals; j++){
    dE[j] = (Etop[j]-energies[j])/(cummulative[j+1]-cummulative[j]);
  }

  //Get maximum possible value of sampled energy
  Emax = Etop[0];
  for(unsigned j = 1; j < nIntervals; j++)
    if(Emax < Etop[j])
      Emax = Etop[j];
  
  //Add an extra interval for the seeki function because we
  //add a 0.0 value at the beginning of the cummulative array
  nIntervals++;

  if(verbose > 1){
    printf("      Elow          Etop       LowCummul     topCummul\n");
    for(unsigned j = 0; j < nIntervals-1; j++)
      printf(" %12.4E  %12.4E  %12.4E  %12.4E\n",energies[j],Etop[j],cummulative[j],cummulative[j+1]);

    printf("Maximum possible energy: %12.4E eV\n",Emax);
  }
  
  return 0;
}

REGISTER_SAMPLER(intervals_energySampling, INTERVALS)
