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
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
// 

#include "VRxraySplitting.hh"

int pen_VRxraysplitting::configure(const pen_parserSection& config,
				   const wrapper_geometry& geometry,
				   const unsigned verbose){

  if(verbose > 1)
    printf("*** X-Ray splitting configuration ***\n\n");

  int err;

  //Get materials used by the current geometry
  bool usedMat[constants::MAXMAT+1];
  geometry.usedMat(usedMat);
  
  // Materials
  //************

  std::vector<std::string> xRayMat;
  config.ls("materials",xRayMat);

  if(xRayMat.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Material x-ray splitting:\n\n");
      printf("  Mat  | x-ray splitting\n");      
    }
    
    for(unsigned imname = 0; imname < xRayMat.size(); imname++){

      //Get ibody name section
      pen_parserSection matSection;
      std::string matSecKey = std::string("materials/") + xRayMat[imname];
      if(config.readSubsection(matSecKey,matSection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_VRxraysplitting: configure: '%s' is not a "
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
	  printf("pen_VRxraysplitting: configure: Error: Material index "
		 "not specified for material '%s'. Integer expected\n",
		 xRayMat[imname].c_str());
	}
	return -23;
      }

      //Check if material index is in range
      if(imat < 1 || imat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_VRxraysplitting: configure: Error: Specified index "
		 "(%d) for material '%s' out of range.\n",
		 imat,xRayMat[imname].c_str());
	}
	return -24;
      }

      //Check if material is used at current geometry
      if(!usedMat[imat]){
	if(verbose > 0){
	  printf("pen_VRxraysplitting: configure: Error: Specified index "
		 "(%d) for material '%s' is not used at current geometry.\n",
		 imat,xRayMat[imname].c_str());
	}
	return -25;
      }

      //Get splitting factor
      int splitting;
      err = matSection.read("splitting",splitting);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("pen_VRxraysplitting: configure: Error: Unable to read "
		 "field 'splitting' for x-ray splitting on "
		 "material '%s'. Integer expected.\n",
		 xRayMat[imname].c_str());
	return -26;
      }

      //Check splitting factor
      if(splitting < 1){
	if(verbose > 0)
	  printf("pen_VRxraysplitting: configure: Error: Invalid x-ray "
		 "splitting factor (%d).\n",splitting);
	return -27;
      }

      //Set splitting factor for bodies with this material index
      for(unsigned ibody = 0; ibody < geometry.getBodies(); ibody++){
      
	if(geometry.getMat(ibody) != (unsigned)imat) continue;

	if(ibody >= pen_geoconst::NB){
	  if(verbose > 0){
	    printf("pen_VRxraysplitting: configure: Error: Maximum body "
		   "index for IF reached (%u)\n",pen_geoconst::NB);
	  }
	  return -27;
	}
	
	IXRSPL[ibody] = (unsigned)splitting;
	LXRSPL[ibody] = true;
      }

      //Print configuration
      if(verbose > 1){
	printf(" %5d   %5d\n", imat,splitting);
      }
      
    }
  }
  else if(verbose > 1){
    printf("No material with x-ray splitting enabled.\n");
  }
  
  // Bodies
  //************

  std::vector<std::string> xRayBodies;
  config.ls("bodies",xRayBodies);

  if(xRayBodies.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Body x-ray splitting:\n\n");
      printf(" Body  | x-ray splitting\n");      
    }
    
    for(unsigned ibname = 0; ibname < xRayBodies.size(); ibname++){

      //Get ibody name section
      pen_parserSection bodySection;
      std::string bodySecKey = std::string("bodies/") + xRayBodies[ibname];
      if(config.readSubsection(bodySecKey,bodySection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_VRxraysplitting: configuration: '%s' is not a "
		 "section, skip this body.\n",bodySecKey.c_str());
	}
	continue;
      }

      // Body index
      //***********************

      //Check if the specified body exists
      unsigned ibody = geometry.getIBody(xRayBodies[ibname].c_str());
      if(ibody >= geometry.getBodies()){
	if(verbose > 0){
	  printf("pen_VRxraysplitting: configure: Error: Body '%s' doesn't "
		 "exists in the loaded geometry.\n",
		 xRayBodies[ibname].c_str());
	}
	return -28;
      }
      else if(ibody >= pen_geoconst::NB){
	if(verbose > 0){
	  printf("pen_VRxraysplitting: configure: Error: Maximum body "
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
	  printf("pen_VRxraysplitting: configure: Error: Unable to "
		 "read field 'splitting' for x-ray splitting on "
		 "body '%s'. Integer expected.\n",
		 xRayBodies[ibname].c_str());
	return -29;
      }

      //Check splitting factor
      if(splitting < 1){
	if(verbose > 0)
	  printf("pen_VRxraysplitting: configure: Error: Invalid x-ray "
		 "splitting factor (%d).\n",splitting);
	return -30;
      }

      //Set splitting factor for specified body
      IXRSPL[ibody] = (unsigned)splitting;	  
      LXRSPL[ibody] = true;	  

      //Print configuration
      if(verbose > 1){
	printf(" %5d   %5d\n", ibody,splitting);
      }
      
    }
  }
  else if(verbose > 1){
    printf("No bodies with x-ray splitting enabled.\n");
  }

  if(verbose > 1){
    printf("\n\nFinal x-ray splitting:\n\n");
    printf(" Body  | x-ray splitting\n");
    for(unsigned ibody = 0; ibody < pen_geoconst::NB; ibody++){      
      
      if(LXRSPL[ibody])
	printf(" %5u   %5u\n", ibody,IXRSPL[ibody]);
    }
    printf("\n\n");
  }

  return 0;
}

void pen_VRxraysplitting::vr_particleStack(const unsigned long long /*nhist*/,
					   const pen_KPAR /*kpar*/,
					   const unsigned /*kdet*/,
					   pen_state_gPol& state,
					   std::vector<pen_state_gPol>& stack,
					   unsigned& created,
					   const unsigned available,
					   pen_rand& random) const{
  
  if(LXRSPL[state.IBODY] && state.ILB[3] > 0){
    //Is a characteristic x-ray in a body with x-ray splitting enabled
    if(state.ILB[0] == 2 && state.ILB[2] < 9){
      // Unsplitted 2nd generation photon
      unsigned freeSpace = available - created;
      unsigned nsplit = std::min(IXRSPL[state.IBODY],freeSpace);
      state.WGHT /= (double) nsplit;
      state.ILB[2] = 9; //Labels split x rays
	    
      //Store 'nsplit' states
      pen_state_gPol stateSplit;
      stateSplit = state;

      for(unsigned isplit = 1; isplit < nsplit; ++isplit){
	      
	stateSplit.W = -1.0 + 2.0 * random.rand();
	double SDTS = sqrt(1.0 - stateSplit.W*stateSplit.W);
	double DF = constants::TWOPI*random.rand();
	stateSplit.U = cos(DF)*SDTS;
	stateSplit.V = sin(DF)*SDTS;
	stack[created++] = stateSplit;
      }
    }
  }
  
}

REGISTER_VR(pen_VRxraysplitting, pen_state_gPol, XRAY_SPLITTING)
