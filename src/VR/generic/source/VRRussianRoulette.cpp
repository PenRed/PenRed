
//
//
//    Copyright (C) 2020-2021 Universitat de València - UV
//    Copyright (C) 2020-2021 Universitat Politècnica de València - UPV
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

#include "VRRussianRoulette.hh"

int pen_VRRussianRoulette::configure(const pen_parserSection& config,
				     const wrapper_geometry& geometry,
				     const unsigned verbose){

  if(verbose > 1)
    printf("*** Russian roulette configuration ***\n\n");

  int err;

  //Get materials used by the current geometry
  bool usedMat[constants::MAXMAT+1];
  geometry.usedMat(usedMat);

  //Get weight window
  if(config.read("minWght",minWght) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRRussianRoulette: configure: Error: Missing minimum "
	     "weight ('minWght'). Double expected.\n");
    }
    return -1;
  }
  if(config.read("maxWght",maxWght) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_VRRussianRoulette: configure: Error: Missing maximum "
	     "weight ('maxWght'). Double expected.\n");
    }
    return -2;
  }  

  if(minWght < 0.0)
    minWght = 0.0;

  if(maxWght <= minWght){
    if(verbose > 0)
      printf("pen_VRRussianRoulette: configure: Error: Maximum weight must be "
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
      printf("\n\n **** Material russian roulette:\n\n");
      printf("  Mat  | survival prob\n");      
    }
    
    for(unsigned imname = 0; imname < mats.size(); imname++){

      //Get ibody name section
      pen_parserSection matSection;
      std::string matSecKey = std::string("materials/") + mats[imname];
      if(config.readSubsection(matSecKey,matSection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_VRRussianRoulette: configure: '%s' is not a "
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
	  printf("pen_VRRussianRoulette: configure: Error: Material index "
		 "not specified for material '%s'. Integer expected\n",
		 mats[imname].c_str());
	}
	return -4;
      }

      //Check if material index is in range
      if(imat < 1 || imat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_VRRussianRoulette: configure: Error: Specified index "
		 "(%d) for material '%s' out of range.\n",
		 imat,mats[imname].c_str());
	}
	return -5;
      }

      //Check if material is used at current geometry
      if(!usedMat[imat]){
	if(verbose > 0){
	  printf("pen_VRRussianRoulette: configure: Error: Specified index "
		 "(%d) for material '%s' is not used at current geometry.\n",
		 imat,mats[imname].c_str());
	}
	return -6;
      }

      //Get survival probability factor
      double survivalProb;
      err = matSection.read("prob",survivalProb);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("pen_VRRussianRoulette: configure: Error: Unable to read "
		 "field 'prob' for russian roulette on "
		 "material '%s'. Double expected.\n",
		 mats[imname].c_str());
	return -7;
      }

      //Check survival probability
      if(survivalProb < 1.0e-5 || survivalProb > 1.0){
	if(verbose > 0)
	  printf("pen_VRRussianRoulette: configure: Error: Invalid "
		 "survival probability (%.5E).\n",survivalProb);
	return -8;
      }

      //Set survival probability for bodies with this material index
      for(unsigned ibody = 0; ibody < geometry.getBodies(); ibody++){
      
	if(geometry.getMat(ibody) != (unsigned)imat) continue;

	if(ibody >= pen_geoconst::NB){
	  if(verbose > 0){
	    printf("pen_VRRussianRoulette: configure: Error: Maximum body "
		   "index reached (%u)\n",pen_geoconst::NB);
	  }
	  return -9;
	}
	
	DRR[ibody] = survivalProb;
	LRR[ibody] = true;
      }

      //Print configuration
      if(verbose > 1){
	printf(" %5d   %.5E\n", imat,survivalProb);
      }
      
    }
  }
  else if(verbose > 1){
    printf("No material with russian roulette enabled.\n");
  }
  
  // Bodies
  //************

  std::vector<std::string> bodies;
  config.ls("bodies",bodies);

  if(bodies.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Body roussian roulette:\n\n");
      printf(" Body  | survival prob\n");      
    }
    
    for(unsigned ibname = 0; ibname < bodies.size(); ibname++){

      //Get ibody name section
      pen_parserSection bodySection;
      std::string bodySecKey = std::string("bodies/") + bodies[ibname];
      if(config.readSubsection(bodySecKey,bodySection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_VRRussianRoulette: configuration: '%s' is not a "
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
	  printf("pen_VRRussianRoulette: configure: Error: Body '%s' doesn't "
		 "exists in the loaded geometry.\n",
		 bodies[ibname].c_str());
	}
	return -28;
      }
      else if(ibody >= pen_geoconst::NB){
	if(verbose > 0){
	  printf("pen_VRRussianRoulette: configure: Error: Maximum body "
		 "index is (%u)\n",pen_geoconst::NB);
	  printf("                     specified index: %u\n",ibody);
	}
	return -28;
      }

      // Survival probability
      //***********************

      //Get survival probability
      double survivalProb;
      err = bodySection.read("prob",survivalProb);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("pen_VRRussianRoulette: configure: Error: Unable to "
		 "read field 'prob' for russian roulette on "
		 "body '%s'. double expected.\n",
		 bodies[ibname].c_str());
	return -29;
      }

      //Check survival prob
      if(survivalProb < 1.0e-5 || survivalProb > 1.0){
	if(verbose > 0)
	  printf("pen_VRRussianRoulette: configure: Error: Invalid "
		 "survival probability (%.5E).\n",survivalProb);
	return -30;
      }

      //Set survival probability for specified body
      DRR[ibody] = survivalProb;	  
      LRR[ibody] = true;	  

      //Print configuration
      if(verbose > 1){
	printf(" %5d   %.5E\n", ibody,survivalProb);
      }
      
    }
  }
  else if(verbose > 1){
    printf("No bodies with russian roulette enabled.\n");
  }

  if(verbose > 1){
    printf("\n\nFinal russian roulette:\n\n");
    printf(" Body  | survival prob\n");
    for(unsigned ibody = 0; ibody < pen_geoconst::NB; ibody++){      
      
      if(LRR[ibody])
	printf(" %5u   %.5E\n", ibody,DRR[ibody]);
    }
    printf("\n\n");
  }

  return 0;
}


void pen_VRRussianRoulette::vr_interfCross(const unsigned long long /*nhist*/,
					   const pen_KPAR /*kpar*/,
					   const unsigned /*kdet*/,
					   pen_particleState& state,
					   std::array<pen_particleState,constants::NMS>& /*stack*/,
					   unsigned& /*created*/,
					   const unsigned /*available*/,
					   pen_rand& random) const{

  //Check if russian roulette is enabled in this body
  if(LRR[state.IBODY]){

    if(state.WGHT < minWght || state.WGHT >= maxWght)
      return;

    if(random.rand() > DRR[state.IBODY]){
      //Kill the particle
      state.E = 0.0;
    }
    else{
      //The particle survive
      state.WGHT /= DRR[state.IBODY];
    }
  }
}

REGISTER_VR(pen_VRRussianRoulette,pen_particleState, RUSSIAN_ROULETTE)
