
//
//
//    Copyright (C) 2020 Universitat de València - UV
//    Copyright (C) 2020 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

#include "VRsplitting.hh"

int pen_VRradialSplitting::configure(const pen_parserSection& config,
			       const wrapper_geometry& /*geometry*/,
			       const unsigned verbose){

  if(verbose > 1)
    printf("*** Radial splitting configuration ***\n\n");

  // Get weight window
  //*******************
  if(config.read("minWght",minWght) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: Missing minimum "
	     "weight ('minWght'). Double expected.\n");
    }
    return -1;
  }
  if(config.read("maxWght",maxWght) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: Missing maximum "
	     "weight ('maxWght'). Double expected.\n");
    }
    return -2;
  }

  if(minWght < 0.0)
    minWght = 0.0;

  if(maxWght <= minWght){
    if(verbose > 0)
      printf("pen_VRradialSplitting: configure: Error: Maximum weight must be "
	     "greater than minimum.\n"
	     "          minimum weight: %.5E\n"
	     "          maximum weight: %.5E\n",
	     minWght,maxWght);
    return -1;
  }

  if(verbose > 1)
    printf("\n Weight window: [%.5E,%.5E) \n",minWght,maxWght);

  // Get center coordinates
  //************************
  if(config.read("x",x0) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: Missing 'X' "
	     "coordinate. Double expected.\n");
    }
    return -3;
  }
  if(config.read("y",y0) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: Missing 'Y' "
	     "coordinate. Double expected.\n");
    }
    return -4;
  }
  if(config.read("z",z0) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: Missing 'Z' "
	     "coordinate. Double expected.\n");
    }
    return -5;
  }

  if(verbose > 1)
    printf("\n Splitting sphere center: (%.5E,%.5E,%.5E) \n",x0,y0,z0);
  
  // Get minimum and maximum radius
  //********************************
  if(config.read("rmin",rmin) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: Missing 'rmin' "
	     "parameter. Double expected.\n");
    }
    return -6;
  }
  if(config.read("rmax",rmax) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: Missing 'rmax' "
	     "parameter. Double expected.\n");
    }
    return -7;
  }

  if(rmin <= 0.0){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: 'rmin' must be "
	     "greater than 0.\n"
	     "        provided 'rmin' value: %.5E\n",rmin);
    }
    return -8;
  }
  
  if(rmin >= rmax){
    if(verbose > 0){
      printf("pen_VRradialSplitting: configure: Error: 'rmax' must be "
	     "greater than 'rmin'.\n"
	     "        provided 'rmin' value: %.5E\n"
	     "        provided 'rmax' value: %.5E\n",rmin,rmax);
    }
    return -9;
  }

  rmin2 = rmin*rmin;
  rmax2 = rmax*rmax;

  if(verbose > 1)
    printf("\n Radial interval: [%.5E,%.5E) \n",rmin,rmax);

  //Check if the user has selected linear multiplier
  if(config.read("linear",linear) != INTDATA_SUCCESS){
    linear = false;
  }

  if(verbose > 1){
    if(linear)
      printf("\n Linear splitting factor scale selected\n");
    else
      printf("\n Quadratic splitting factor scale selected\n");
  }
  
  return 0;
}


void pen_VRradialSplitting::vr_particleStack(const unsigned long long /*nhist*/,
					     const pen_KPAR /*kpar*/,
					     const unsigned /*kdet*/,
					     pen_particleState& state,
					     std::vector<pen_particleState>& stack,
					     unsigned& created,
					     const unsigned available,
					     pen_rand& /*random*/) const{

  //Check weight window
  if(state.WGHT < minWght || state.WGHT >= maxWght)
    return;

  //Calculate distance to the center
  double dx = state.X - x0;
  double dy = state.Y - y0;
  double dz = state.Z - z0;

  double r2 = dx*dx + dy*dy + dz*dz;

  //Check if is out of the splitting zone
  if(r2 < rmin2 || r2 >= rmax2)
    return;

  unsigned toSplit;
  if(linear){
    toSplit = sqrt(r2)/rmin;
  }
  else
    toSplit = r2/rmin2;

  if(toSplit <= 1)
    return;
  
  //Check the available space at the stack
  unsigned freeSpace = available - created;
  unsigned nsplit = std::min(toSplit,freeSpace);
  
  //Reduce the weight according to splitting factor
  state.WGHT /= (double) nsplit;

  
  //Clone the state
  for(unsigned isplit = 1; isplit < nsplit; ++isplit){
    stack[created++] = state;
  }
  
  
}

REGISTER_VR(pen_VRradialSplitting, pen_particleState, RADIAL_SPLITTING)
