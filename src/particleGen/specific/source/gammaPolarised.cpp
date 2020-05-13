
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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

#include "gammaPolarised.hh" 

void gammaPolarised_specificSampler::sample(pen_state_gPol& state,
					    pen_KPAR& /*genKpar*/,
					    unsigned long long& /*dhist*/,
					    pen_rand& /*random*/){

  state.IPOL = IPOL0;
  state.SP1  = SP10;
  state.SP2  = SP20;
  state.SP3  = SP30;
}

int gammaPolarised_specificSampler::configure(double& /*Emax*/,
					      const abc_spatialSampler* /*pSpatial*/,
					      const abc_directionSampler* /*pDirection*/,
					      const abc_energySampler* /*pEnergy*/,
					      const abc_timeSampler* /*pTime*/,
					      const pen_parserSection& config,
					      const unsigned verbose){

  int err = 0;

  // Read Stokes parameters
  //*************************
  pen_parserArray stokesparam;
  err = config.read("stokes", stokesparam);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("gammaPolarised:configure: Error: Unable to read 'stokes' in configuration. Array with 3 elements expected.\n");
    }
    return -1;
  }

  err = stokesparam[0].read(SP10);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("gammaPolarised:configure: Error: Unable to read 'SP1' in array at position 0. Double expected\n");
    }
    return -2;
  }

  err = stokesparam[1].read(SP20);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("gammaPolarised:configure: Error: Unable to read 'SP2' in array at position 1. Double expected\n");
    }
    return -3;
  }

  err = stokesparam[2].read(SP30);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("gammaPolarised:configure: Error: Unable to read 'SP3' in array at position 2. Double expected\n");
    }
    return -4;
  }

  if(fabs(SP10) > 1.0E-16 || fabs(SP20) > 1.0E-16 || fabs(SP30) > 1.0E-16){
    IPOL0 = 1;
  }
  else{
    IPOL0 = 0;
  }

  if(verbose > 1){
    printf("Polarised primary photons. Stokes Parameters:\n");
    printf("      P1 = %13.6E (linear polarisation at 45 deg azimuth)\n",SP10);
    printf("      P2 = %13.6E (circular polarisation)\n",SP20);
    printf("      P3 = %13.6E (linear polarisation at zero azimuth)\n",SP30);
    printf("      Polarisation %s\n", IPOL0 == 1 ? "enabled" : "disabled");
      
  }
  
  return 0;
}

REGISTER_SPECIFIC_SAMPLER(gammaPolarised_specificSampler,pen_state_gPol, GAMMA_POL)
