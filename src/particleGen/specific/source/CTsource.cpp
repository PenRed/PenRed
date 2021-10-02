
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#include "CTsource.hh"

void ct_specificSampler::skip(const unsigned long long dhists){
  psf.skip(dhists);
}

int ct_specificSampler::configure(double& Emax,
				  const abc_spatialSampler* /*pSpatial*/,
				  const abc_directionSampler* /*pDirection*/,
				  const abc_energySampler* /*pEnergy*/,
				  const abc_timeSampler* /*pTime*/,
				  const pen_parserSection& config,
				  const unsigned verbose){

  //First, initialize the phase space file sampler
  psf.setThread(getThread());
  int err = psf.configure(Emax,nullptr,nullptr,nullptr,nullptr,config,verbose);
  if(err != 0){
    if(verbose > 0)
      printf("ctSource:configure: Error: Unable to configure psf.\n"
	     "               error code: %d\n",err);
    return err;
  }

  if(verbose > 1)
    printf("Phase space file source configured.\n");

  double  phi0, phif, dphi;
  err = config.read("rad", r);
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

  err = config.read("angularStep", dphi);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("CTsource:configure: Error: Unable to read 'angularStep' "
	     "in configuration. Double expected\n");
    }
    return -5;

  }
  

  //Convert to rad
  phif *= M_PI/180.0;
  phi0 *= M_PI/180.0;
  dphi *= M_PI/180.0;
    
  nphi = (phif-phi0)/dphi+1;
  
  if(verbose > 1){
    printf("Number of total angular positions (projections):, \n");
    printf(" %lu \n\n",nphi);
  }

  double errx, erry, errz;
  double psfOrigin[3];
  errx = config.read("psfx0",psfOrigin[0]);
  erry = config.read("psfy0",psfOrigin[1]);
  errz = config.read("psfz0",psfOrigin[2]);
  

  if(errx != INTDATA_SUCCESS ||
     erry != INTDATA_SUCCESS ||
     errz != INTDATA_SUCCESS){
    if(verbose > 0)
      printf("ctsource:configure:Error: Unable to read 'psfx0,psfy0,psfz0'. "
	     " Doubles expected.\n");
    return -6;
  }
  
  double errDirx, errDiry, errDirz;
  double psfOriginDir[3];
  errDirx = config.read("psfDirx0",psfOriginDir[0]);
  errDiry = config.read("psfDiry0",psfOriginDir[1]);
  errDirz = config.read("psfDirz0",psfOriginDir[2]);
  
  if(errDirx != INTDATA_SUCCESS ||
     errDiry != INTDATA_SUCCESS ||
     errDirz != INTDATA_SUCCESS){
    if(verbose > 0)
      printf("ctsource:configure:Error: Unable to read 'psfDirx0,psfDiry0,psfDirz0'. "
	     " Doubles expected.\n");
    return -7;
  }
  
  //Normalize direction of the PSF
  double psfOriginDirModule = sqrt(pow(psfOriginDir[0],2) + pow(psfOriginDir[1],2) + pow(psfOriginDir[2],2));
  
  double psfOriginDirNorm[3];
  psfOriginDirNorm[0] = psfOriginDir[0]/psfOriginDirModule;
  psfOriginDirNorm[1] = psfOriginDir[1]/psfOriginDirModule;
  psfOriginDirNorm[2] = psfOriginDir[2]/psfOriginDirModule;
  
  //thetaX=acos(-psfOriginDirNorm[0]/psfOriginDirModule);
  
  //Obtain polar and azimutal angle of the original PSF direction
  if(fabs(psfOriginDirNorm[0]) > 1.0e-10 || fabs(psfOriginDirNorm[1]) > 1.0e-10)
  {
      //Azimutal
      azim = atan2(psfOriginDirNorm[1],psfOriginDirNorm[0]);
      if(std::signbit(azim))
        azim += 2.0*M_PI;
      
      //Polar
      polar=acos(psfOriginDirNorm[2]);
  }
  else
  {
      //Azimutal
      azim = 0.0;
      //Polar
      if(fabs(psfOriginDirNorm[2]) < 1.0e-10)
          polar = 0.0;
      else
      {
          polar = acos(psfOriginDirNorm[2]);
      }
      
  }
  
  //printf("%14.5E, %14.5E degrees \n\n",polar*180.0/M_PI, azim*180.0/M_PI);
  
  //Angles to rotate the vector of the PSF
  polarRot = M_PI*0.5 - polar; //If polar > PI/2, the rotation will be negative
  azimRot = - azim; //If azim > PI, the rotation will be negative

  createRotationZYZ(azimRot,polarRot,M_PI,particleRot);
  
  double ctOrigin[3];
  errx = config.read("CTx0",ctOrigin[0]);
  erry = config.read("CTy0",ctOrigin[1]);
  errz = config.read("CTz0",ctOrigin[2]);

  if(errx != INTDATA_SUCCESS ||
     erry != INTDATA_SUCCESS ||
     errz != INTDATA_SUCCESS){
    if(verbose > 0)
      printf("ctsource:configure:Error: Unable to read 'CTx0,CTy0,CTz0'. "
	     " Doubles expected.\n");
    return -8;
  }  

  //Calculate the total translation
  // 1- Move particle of psf to origin (using psf Origin information)
  part2psfOrigin.x = -psfOrigin[0];
  part2psfOrigin.y = -psfOrigin[1];
  part2psfOrigin.z = -psfOrigin[2];
  
 if(verbose > 1){
     printf("Translation values of the PSF particle to origin\n");
     printf("(%12.5E, %12.5E, %12.5E) \n\n",part2psfOrigin.x,part2psfOrigin.y, part2psfOrigin.z);
 }
 
  if(verbose > 1){
     printf("Normalized original direction of the PSF\n");
     printf("(%12.5E, %12.5E, %12.5E) \n\n",psfOriginDirNorm[0],psfOriginDirNorm[1],psfOriginDirNorm[2]);
     printf("Rotation needed for -x PSF direction, polar and azimutal angles\n");
     printf("%14.5E, %14.5E degrees \n\n",polarRot*180.0/M_PI, azimRot*180.0/M_PI);
 }  

  //After the rotation, we need to move the psf to
  //the CT center:
  // 3- Translate psf to CT center (ctOrigin) 
  origin2CT.x = ctOrigin[0];
  origin2CT.y = ctOrigin[1];
  origin2CT.z = ctOrigin[2];
  
  
    if(verbose > 1){
      printf("Translation values of the PSF to the CT origin\n");
      printf("(%12.5E, %12.5E, %12.5E) \n\n",origin2CT.x ,origin2CT.y, origin2CT.z);
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
  actualState = 0;
  
  return 0;  
}

void ct_specificSampler::sample(pen_particleState& state,
				 pen_KPAR& genKpar,
				 unsigned long long& dhist,
				 pen_rand& random){

    if(savedDhist == 0){
        //Read first history
        for(;;){
            unsigned long long localdhist = 0;
            pen_particleState localstate;
            pen_KPAR localgenKpar;
            psf.sample(localstate,localgenKpar,localdhist,random);
            if(localdhist != 0 || localgenKpar == ALWAYS_AT_END){
                nextHistFirstState.first  = localstate;
                nextHistFirstState.second = localgenKpar;
                savedDhist = localdhist;
                if(savedDhist == 0)
                    savedDhist = 1;
                break;
            }
            else{
                histStates.push_back(
                    std::pair<pen_particleState,pen_KPAR>
                    (localstate,localgenKpar));
            }
        }
    }
  //First, get a particle from buffer
  if(actualState < histStates.size()){
      //Get next state
      state   = histStates[actualState].first;
      genKpar = histStates[actualState].second;
      dhist = 0;
      ++actualState;
  }
  else{
      //Increment CT position
      ++actualCTpos;
      dhist = savedDhist;
      if(actualCTpos < nphi){
          //Change projection and reset buffer read
          state   = histStates[0].first;
          genKpar = histStates[0].second;
          actualState = 1;
      }
      else{
        //Finished buffer and projections, get next hist first state
        if(nextHistFirstState.second == ALWAYS_AT_END){
            genKpar = ALWAYS_AT_END;
            return;
        }
        state   = nextHistFirstState.first;
        genKpar = nextHistFirstState.second;

        //Clear previous buffer
        histStates.clear();
        //Add next hist first state
        histStates.push_back(nextHistFirstState);
        
        actualState = 1;        
        savedDhist = 0;
        actualCTpos = 0;
      }
  }
  
  //Get the CT position
  const unsigned long CTpos = actualCTpos;

  //Add the corresponding time
  state.PAGE += tmin+CTpos*dt;
  state.LAGE = true;
  
  //printf(out,"#CTpos = %lu\n", CTpos);
  
  //printf(out,"#Original state:   state.X   state.Y   state.Z    state.U    state.V    state.W\n");
  //printf(out,"                  %14.5E    %14.5E    %14.5E     %14.5E    %14.5E    %14.5E\n", state.X, state.Y, state.Z,  state.U, state.V, state.W);
  
  
  //Apply the corresponding translations and rotations to the particle
  //First, move it to the psf center 
  part2psfOrigin.translate(state);
  

  
  //Apply the corresponding rotation in y axis due the angle between x direction
  //of the psf and their module
  double pos[3] = {state.X,state.Y,state.Z};
  double dir[3] = {state.U,state.V,state.W};
  matmul3D(particleRot,pos);
  matmul3D(particleRot,dir);

  state.X = pos[0];
  state.Y = pos[1];
  state.Z = pos[2];

  state.U = dir[0];
  state.V = dir[1];
  state.W = dir[2];
  
  //Apply the corresponding translation to R x+=r
  state.X += r;
  
  
  //Apply the corresponding rotation according to CT position
  rotations[CTpos].rotate(state);

  
  
  //Finally, move the particle to the CT
  origin2CT.translate(state);
}

REGISTER_SPECIFIC_SAMPLER(ct_specificSampler,pen_particleState, CT)
