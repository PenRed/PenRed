
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

int pen_VRsplitting::configure(const pen_parserSection& config,
			       const wrapper_geometry& geometry,
			       const unsigned verbose){

  if(verbose > 1)
    printf("*** Splitting configuration ***\n\n");

  int err;

  //Get materials used by the current geometry
  bool usedMat[constants::MAXMAT+1];
  geometry.usedMat(usedMat);

  //Get weight window
  if(config.read("minWght",minWght) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRsplitting: configure: Error: Missing minimum "
	     "weight ('minWght'). Double expected.\n");
    }
    return -1;
  }
  if(config.read("maxWght",maxWght) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRsplitting: configure: Error: Missing maximum "
	     "weight ('maxWght'). Double expected.\n");
    }
    return -2;
  }  

  if(minWght < 0.0)
    minWght = 0.0;

  if(maxWght <= minWght){
    if(verbose > 0)
      printf("pen_VRsplitting: configure: Error: Maximum weight must be "
	     "greater than minimum.\n"
	     "          minimum weight: %.5E\n"
	     "          maximum weight: %.5E\n",
	     minWght,maxWght);
    return -3;
  }

  if(verbose > 1)
    printf("\n Weight window: [%.5E,%.5E) \n",minWght,maxWght);
  
  // Materials
  //************

  std::vector<std::string> mats;
  config.ls("materials",mats);

  if(mats.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Material splitting:\n\n");
      printf("  Mat  | splitting\n");      
    }
    
    for(unsigned imname = 0; imname < mats.size(); imname++){

      //Get ibody name section
      pen_parserSection matSection;
      std::string matSecKey = std::string("materials/") + mats[imname];
      if(config.readSubsection(matSecKey,matSection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_VRsplitting: configure: '%s' is not a "
		 "section, skip this material.\n",matSecKey.c_str());
	}
	continue;
      }

      // Material index
      //***********************

      //Read index
      int imat;
      err = matSection.read("mat-index",imat);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_VRsplitting: configure: Error: Material index "
		 "not specified for material '%s'. Integer expected\n",
		 mats[imname].c_str());
	}
	return -4;
      }

      //Check if material index is in range
      if(imat < 1 || imat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_VRsplitting: configure: Error: Specified index "
		 "(%d) for material '%s' out of range.\n",
		 imat,mats[imname].c_str());
	}
	return -5;
      }

      //Check if material is used at current geometry
      if(!usedMat[imat]){
	if(verbose > 0){
	  printf("pen_VRsplitting: configure: Error: Specified index "
		 "(%d) for material '%s' is not used at current geometry.\n",
		 imat,mats[imname].c_str());
	}
	return -6;
      }

      //Get splitting factor
      int splitting;
      err = matSection.read("splitting",splitting);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("pen_VRsplitting: configure: Error: Unable to read "
		 "field 'splitting' for x-ray splitting on "
		 "material '%s'. Integer expected.\n",
		 mats[imname].c_str());
	return -7;
      }

      //Check splitting factor
      if(splitting < 1){
	if(verbose > 0)
	  printf("pen_VRsplitting: configure: Error: Invalid "
		 "splitting factor (%d).\n",splitting);
	return -8;
      }

      //Set splitting factor for bodies with this material index
      for(unsigned ibody = 0; ibody < geometry.getBodies(); ibody++){
      
	if(geometry.getMat(ibody) != (unsigned)imat) continue;

	if(ibody >= pen_geoconst::NB){
	  if(verbose > 0){
	    printf("pen_VRsplitting: configure: Error: Maximum body "
		   "index reached (%u)\n",pen_geoconst::NB);
	  }
	  return -9;
	}
	
	ISPL[ibody] = (unsigned)splitting;
	LSPL[ibody] = true;
      }

      //Print configuration
      if(verbose > 1){
	printf(" %5d   %5d\n", imat,splitting);
      }
      
    }
  }
  else if(verbose > 1){
    printf("No material with splitting enabled.\n");
  }
  
  // Bodies
  //************

  std::vector<std::string> bodies;
  config.ls("bodies",bodies);

  if(bodies.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Body splitting:\n\n");
      printf(" Body  | splitting\n");      
    }
    
    for(unsigned ibname = 0; ibname < bodies.size(); ibname++){

      //Get ibody name section
      pen_parserSection bodySection;
      std::string bodySecKey = std::string("bodies/") + bodies[ibname];
      if(config.readSubsection(bodySecKey,bodySection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_VRsplitting: configuration: '%s' is not a "
		 "section, skip this body.\n",bodySecKey.c_str());
	}
	continue;
      }

      // Body index
      //***********************

      //Check if the specified body exists
      unsigned ibody = geometry.getIBody(bodies[ibname].c_str());
      if(ibody >= geometry.getBodies()){
	if(verbose > 0){
	  printf("pen_VRsplitting: configure: Error: Body '%s' doesn't "
		 "exists in the loaded geometry.\n",
		 bodies[ibname].c_str());
	}
	return -28;
      }
      else if(ibody >= pen_geoconst::NB){
	if(verbose > 0){
	  printf("pen_VRsplitting: configure: Error: Maximum body "
		 "index is (%u)\n",pen_geoconst::NB);
	  printf("                     specified index: %u\n",ibody);
	}
	return -28;
      }

      // Splitting factor
      //***********************

      //Get splitting factor
      int splitting;
      err = bodySection.read("splitting",splitting);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("pen_VRsplitting: configure: Error: Unable to "
		 "read field 'splitting' for splitting on "
		 "body '%s'. Integer expected.\n",
		 bodies[ibname].c_str());
	return -29;
      }

      //Check splitting factor
      if(splitting < 1){
	if(verbose > 0)
	  printf("pen_VRsplitting: configure: Error: Invalid "
		 "splitting factor (%d).\n",splitting);
	return -30;
      }

      //Set splitting factor for specified body
      ISPL[ibody] = (unsigned)splitting;	  
      LSPL[ibody] = true;	  

      //Print configuration
      if(verbose > 1){
	printf(" %5d   %5d\n", ibody,splitting);
      }
      
    }
  }
  else if(verbose > 1){
    printf("No bodies with splitting enabled.\n");
  }

  if(verbose > 1){
    printf("\n\nFinal splitting:\n\n");
    printf(" Body  | splitting\n");
    for(unsigned ibody = 0; ibody < pen_geoconst::NB; ibody++){      
      
      if(LSPL[ibody])
	printf(" %5u   %5u\n", ibody,ISPL[ibody]);
    }
    printf("\n\n");
  }

  return 0;
}


void pen_VRsplitting::vr_interfCross(const unsigned long long /*nhist*/,
				     const pen_KPAR kpar,
				     const unsigned /*kdet*/,
				     pen_particleState& state,
				     std::array<pen_particleState,constants::NMS>& stack,
				     unsigned& created,
				     const unsigned available,
				     pen_rand& /*random*/) const{

  //Check if splitting is enabled in this body
  if(LSPL[state.IBODY]){

    if(state.WGHT < minWght || state.WGHT >= maxWght)
      return;
    
    //Check the available space at the stack
    unsigned freeSpace = available - created;
    unsigned nsplit = std::min(ISPL[state.IBODY],freeSpace);

    if(nsplit <= 1)
        return;
    
    //Reduce the weight according to splitting factor
    state.WGHT /= (double) nsplit;

    //Clone the state
    for(unsigned isplit = 1; isplit < nsplit; ++isplit){
      stack[created++] = state;
    }
  }
  
}

REGISTER_VR(pen_VRsplitting,pen_particleState, SPLITTING)
