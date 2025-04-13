
//
//
//    Copyright (C) 2020-2023 Universitat de València - UV
//    Copyright (C) 2020-2023 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#include "CTsource.hh"

void ct_specificSampler::skip(const unsigned long long dhists){
  if(!genericSource)
    psf.skip(dhists);
}

int ct_specificSampler::configure(double& Emax,
				  const pen_parserSection& config,
				  const unsigned nthreads,
				  const unsigned verbose){
  
  //First, initialize the phase space file sampler
  psf.setThread(getThread());
  
  //Check if the CT source will use generic samplers or a psf
  if(spatial() != nullptr || direction() != nullptr || energy() != nullptr || time() != nullptr){
    //Generic samplers will be used
    
    genericSource = true;
    if(spatial() == nullptr || direction() == nullptr || energy() == nullptr){
      
      if(verbose > 0){
	printf("ctSource:configure: Error: Spatial, direction and energy sources are "
	       "mandatory to use generic sources within the CT source. Specified generic sources:\n"
	       " Spatial   : %s\n"
	       " Direction : %s\n"
	       " Energy    : %s\n"
	       " Time      : %s\n",
	       spatial() == nullptr ? "Null" : spatial()->readID(),
	       direction() == nullptr ? "Null" : direction()->readID(),
	       energy() == nullptr ? "Null" : energy()->readID(),
	       time() == nullptr ? "Null" : time()->readID());
      }
      return 1;
    }
    
    //All mandatory samplers have been specified, read the number of particles per history
    
    // Number of particles per history
    //**********************************
    int partPerHistAux;
    int err = config.read("nSecondaries", partPerHistAux);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("ctSource:configure: Error: Unable to read 'nSecondaries' "
	       "in configuration. Integer expected\n");
      }
      return -2;
    }

    if(partPerHistAux < 1){
      if(verbose > 0){
	printf("ctSource:configure: Error: The number of secondary "
	       "particles per history must be greater than zero\n");
      }
      return -3;
    }
    
    partPerHist = static_cast<unsigned int>(partPerHistAux);
    
    if(verbose > 1){
      printf("\nNumber of secondary particles per history: %u\n",partPerHist);
    }
    
  }
  else{
    
    //No generic sampler specified, a psf source will be used
    //genericSource = false;
    
    if(verbose > 1)
      printf("CT with phase space file source enabled.\n");
    
    //Get psf subsection
    pen_parserSection psfSection;
    int err = config.readSubsection("psf",psfSection);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0)
	printf("ctSource:configure: Error: Unable to read 'psf' configuration section.\n"
	       "               error code: %d\n",err);
      return err;
    }
    
    err = psf.configure(Emax,psfSection,nthreads,verbose);
    if(err != 0){
      if(verbose > 0)
	printf("ctSource:configure: Error: Unable to configure psf.\n"
	       "               error code: %d\n",err);
      return err;
    }

    if(verbose > 1)
      printf("Phase space file source configured.\n");
    
  }

  //Continue with common CT configuration for both source types

  double  phi0, phif, dphi;
  int err = config.read("rad", r);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("ctSource:configure: Error: Unable to read 'rad' in "
	     "configuration. Double expected\n");
    }
    return -1;
  }
  
  err = config.read("phi-ini", phi0);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("CTsource:configure: Error: Unable to read 'phi-ini' "
	     "in configuration. Double expected\n");
    }
    return -2;

  }

  
  err = config.read("phi-end", phif);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("CTsource:configure: Error: Unable to read 'phi-end' "
	     "in configuration. Double expected\n");
    }
    return -3;

  }

  if(phi0 >= phif)
    {
      if(verbose > 0){
	printf("CTsource:configure: Error: 'phi-ini' "
	       "value must be lower than 'phi-end' value.\n");
      }
      return -4;
    }

  int auxNProj;
  err = config.read("nProjections", auxNProj);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("CTsource:configure: Error: Unable to read 'angularStep' "
	     "in configuration. Double expected\n");
    }
    return -5;

  }

  if(auxNProj < 1){
    if(verbose > 0){
      printf("CTsource:configure: Error: Invalid number of projections (%d). "
	     "One or more projections is required\n", auxNProj);
    }
    return -5;    
  }

  nphi = static_cast<unsigned long>(auxNProj);

  //Convert to rad
  phif *= M_PI/180.0;
  phi0 *= M_PI/180.0;
  dphi = (phif-phi0)/static_cast<double>(nphi);
  
  if(verbose > 1){
    printf("Number of total angular positions (projections):, \n");
    printf(" %lu \n\n",nphi);
  }
    
  double ctOrigin[3];
  int errx = config.read("CTx0",ctOrigin[0]);
  int erry = config.read("CTy0",ctOrigin[1]);
  int errz = config.read("CTz0",ctOrigin[2]);

  if(errx != INTDATA_SUCCESS ||
     erry != INTDATA_SUCCESS ||
     errz != INTDATA_SUCCESS){
    if(verbose > 0)
      printf("ctsource:configure:Error: Unable to read 'CTx0,CTy0,CTz0'. "
	     " Doubles expected.\n");
    return -8;
  }

  //After the rotation, we need to move the psf to
  //the CT center:
  // 3- Translate psf to CT center (ctOrigin) 
  origin2CT.x = ctOrigin[0];
  origin2CT.y = ctOrigin[1];
  origin2CT.z = ctOrigin[2];
  
  
  if(verbose > 1){
    printf("CT origin:\n"
	   "(%12.5E, %12.5E, %12.5E) cm \n\n",
	   origin2CT.x ,origin2CT.y, origin2CT.z);
  }
  
  // Time window
  //************
  err = config.read("tmin", tmin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("ctSource:configure: Error: Unable to read 'tmin' "
	     "in configuration. Double expected\n");
    }
    return -9;
  }

  
  err = config.read("dt", dt);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("ctSource:configure: Error: Unable to read 'dt' "
	     "in configuration. Double expected\n");
    }
    return -10;
  }

  
  //pre-calculate all required rotation cosinus and sinus
  rotations.resize(nphi);
  for(size_t i = 0; i < nphi; ++i){
    rotations[i].c = cos(phi0+static_cast<double>(i)*dphi);
    rotations[i].s = sin(phi0+static_cast<double>(i)*dphi);
  }
  
  //Allocate initial vector memory
  histStates.reserve(10000);
  
  //Important!! Initialize these variables or the sample function
  //will not work properly
  savedDhist = 0;
  actualCTpos = 0;
  //If a generic sampler is used, force the first sampled particle to be a new history
  actualState = genericSource ? partPerHist : 0;
  
  return 0;  
}

void ct_specificSampler::sample(pen_particleState& state,
				pen_KPAR& genKpar,
				unsigned long long& dhist,
				pen_rand& random){
  
  //Check if the source is using a generic or a psf sampler
  if(genericSource){
      
    //Generic source      
    sampleGeneric(state, random);
    state.ILB[0] = 2; //Flag particles as secondaries
      
    //Check if the number of particles sampled in this history reaches the limit
    if(actualState >= partPerHist){
      //Increase history number
      dhist = 1;
      //Increment CT position
      ++actualCTpos;
      if(actualCTpos >= nphi){
	actualCTpos = 0;
      }
      //Reset number of particles sampled in this history
      actualState = 0;
    }
    else{
      dhist = 0;
    }
    //Increase number of particles sampled in this history
    ++actualState;
      
  }
  else{
      
    //PSF sources. In this sampling method, each history will be simulated on each CT position
      
    //Check if no history increment has been saved, i.e. the next history must be read
    if(savedDhist == 0){
      //Read all particles belonging to the next history
      for(;;){
	unsigned long long localdhist = 0;
	pen_particleState localstate;
	pen_KPAR localgenKpar;
	psf.sample(localstate,localgenKpar,localdhist,random);
	      
	//Check if the history or the psf ends
	if(localdhist != 0 || localgenKpar == ALWAYS_AT_END){

	  //Check if the buffer is empty
	  if(histStates.size() == 0){
	    //The buffer is empty because the first psf particle
	    //skips one or more histories. Save that new history
	    histStates.push_back(std::pair<pen_particleState,pen_KPAR>
				 (localstate,localgenKpar));
	    continue; //Continue to read the whole history
	    
	  }
	  
	  //Save the first state of the next history, because has been already read
	  nextHistFirstState.first  = localstate;
	  nextHistFirstState.second = localgenKpar;
	  //Save the history increment to be applied on history finish
	  savedDhist = localdhist;
	  if(savedDhist == 0)
	    savedDhist = 1;
	  break;
	}
	else{
	  //Save particle state in the buffer and read the next state
	  histStates.push_back(
			       std::pair<pen_particleState,pen_KPAR>
			       (localstate,localgenKpar));
	}
      }
    }
      
    //Get the next particle state from the buffer
    if(actualState < histStates.size()){
      //Get next state
      state   = histStates[actualState].first;
      genKpar = histStates[actualState].second;
      dhist = 0;
      ++actualState;
    }
    else{
      //All states in the buffer have been simulated on this
      //CT position. Therefore, increase the CT position
      ++actualCTpos;
      dhist = savedDhist;
      //Check if we are on the last CT position
      if(actualCTpos < nphi){
	//Change projection and reset buffer read
	state   = histStates[0].first;
	genKpar = histStates[0].second;
	actualState = 1;
      }
      else{
	//Finished buffer and projections, get next hist first state
	if(nextHistFirstState.second == ALWAYS_AT_END){
	  //The psf end has been reached, finish the sampling
	  genKpar = ALWAYS_AT_END;
	  return;
	}
	    
	//The next history must be read, return the first history state and read the
	//whole history on the next sample call
	state   = nextHistFirstState.first;
	genKpar = nextHistFirstState.second;

	//Clear previous buffer
	histStates.clear();
	//Add next hist first state
	histStates.push_back(nextHistFirstState);
	    
	//Reset sampling state values
	actualState = 1;        
	savedDhist = 0;
	actualCTpos = 0;
      }
    }
      
  }
    
  //Once the particle state has been read/sampled, rotate and translate it to fit the CT movement
    
  //Get the CT position
  const unsigned long CTpos = actualCTpos;

  //Add the corresponding time
  state.PAGE += tmin+CTpos*dt;
  state.LAGE = true;
      
  //Apply the corresponding translation to R x+=r
  state.X += r;
  
  //Apply the corresponding rotation according to CT position
  rotations[CTpos].rotate(state);

  //Finally, move the particle to the CT
  origin2CT.translate(state);

}

REGISTER_SPECIFIC_SAMPLER(ct_specificSampler,pen_particleState, CT)
