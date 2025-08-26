
//
//
//    Copyright (C) 2020-2021 Universitat de València - UV
//    Copyright (C) 2020-2021 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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

#include "tallyCTsinogram.hh"

void pen_CTsinogram::flush(){
    
  //*******************************//
  //*      Particle counters      *//
  //*******************************//
       
  for(unsigned long int i = 0; i < sinoDim; i++) {    
    if(lasthist[i] < 0.5){continue;}
    //Reset last visited to avoid recounting in next report
    lasthist[i] = 0.0;
    //Spectrum particle counters
    if(sinotemp[i] == 0.0){continue;}  // Skip void counters
    sino[i] += sinotemp[i];
    sino2[i] += sinotemp[i]*sinotemp[i];
    //Reset counters
    sinotemp[i] = 0.0;
  }
  
  for(unsigned long int i = 0; i < sinoDim; i++) {    
    if(lasthistNorm[i] < 0.5){continue;}
    //Reset last visited to avoid recounting in next report
    lasthistNorm[i] = 0.0;
    //Spectrum particle counters
    if(sinotempNorm[i] == 0.0){continue;}  // Skip void counters
    sinoNorm[i] += sinotempNorm[i];
    sino2Norm[i] += sinotempNorm[i]*sinotempNorm[i];
    //Reset counters
    sinotempNorm[i] = 0.0;
  }  
  
  for(auto& tally : edepMat)
      tally.flush();
}


bool pen_CTsinogram::detIndexes(const pen_particleState& state,
                                const CTsinogram_point finalPoint, 
                                const double partPage,
                                unsigned long& iphi, 
                                unsigned long& ipix){
  
  double phiMin, phiMax; 
  unsigned long phiBin;    

  //Check if the particle is in the time window
  if(partPage < tmin || partPage >= tmax)
  return false;

  
  
  //printf("Inside time intervale: %14.5E\n", partPage);
  
  //Obtain the phi position of the detector
  phiBin = (partPage-tmin)/dt; 

  //Minimum and maximum angle of the detector. Angular detector window
  phiMin = phi0 - phi*0.5 + dphi*phiBin;
  phiMax = phi + phiMin;
  
  //printf("phiMin = %14.5E, phiMax = %14.5E\n", phiMin, phiMax);
  
  //Values of angles must be in the interval [0,360] degrees
  //No negative values are possible
  if(phiMin < 0.0)
      phiMin = 2.0*M_PI+phiMin;
  if(phiMax < 0.0)
      phiMax = 2.0*M_PI+phiMax;
  
  //No values larger than 360 deg are possible
  if(phiMin == 2.0*M_PI)
      phiMin=0.0;
  if(phiMax == 2.0*M_PI)
      phiMax=0.0;
  if(phiMin > 2.0*M_PI){
      double extra = phiMin/(2.0*M_PI)-1.0;
      phiMin = 2.0*M_PI*extra;
  }
  if(phiMax > 2.0*M_PI){
      double extra = phiMax/(2.0*M_PI)-1.0;
      phiMax = 2.0*M_PI*extra;
  }
  
  /*if(phiMin >= 2.0*M_PI){
      unsigned extra = phiMin/(2.0*M_PI);
      phiMin -= 2.0*M_PI*extra;
  }
  if(phiMax >= 2.0*M_PI){
      unsigned extra = phiMax/(2.0*M_PI);
      phiMax -= 2.0*M_PI*extra;
  }
  */
  //Elements of the quadratic equation to solve y = ax² + bx +c
  double a,b,c; 
  double inv2a, l, x, y, z;
  
  x = finalPoint.X - xOrigin;
  y = finalPoint.Y - yOrigin;
  z = finalPoint.Z - zOrigin;
  a = pow(state.U,2) + pow(state.V,2);

  if(a<1.0E-10)
  return false;
  
  //printf(" x = %14.5E\n y = %14.5E\n z = %14.5E\n", x, y, z);  

  b = state.U*x + state.V*y;
  b *= 2.0;
  c = pow(x,2) + pow(y,2) - ri2;
  inv2a = 0.5/a;

  double ac4 = 4.0*a*c;
  double b2 = pow(b,2);
  
  if(ac4 > b2)
  return false;
  
  //printf("ac4 < b2: %14.5E < %14.5E \n", ac4, b2);  
  
  //Only positive l value of the solution of the equation has physical sense
  
  double l1 = (-b + sqrt(b2-ac4))*inv2a;
  double l2 = (-b - sqrt(b2-ac4))*inv2a;

  if(std::signbit(l1)){
      if(std::signbit(l2))
        return false;
      else
        l = l2;
  }else if(std::signbit(l2)){
      l = l1;
  }
  else{
      l = std::max(l1,l2);
  }
  
  
  //Coordenates where the particle crosses the detector
  double Xfin = x + state.U*l;
  double Yfin = y + state.V*l;
  double Zfin = z + state.W*l;
  //Angle between the particle direction 
  //and z axis to determine the sinogram bin
  double phiPart;
  
  //Check if escaped particle crosses the detector and in which pixel occurs
 
  //printf("Xfin = %14.5E\n", Xfin);
  //printf("Yfin = %14.5E\n", Yfin);
  //printf("Zfin = %14.5E\n", Zfin);
  if(Zfin >= -semiDet && Zfin <= semiDet) //Ring detector intersected
  {
      //printf("Ring intersected\n");
      
      //Check if Xfin and Yfin are zero. In that case, phiPart==0.0
      if(fabs(Xfin)<1.0E-10 && fabs(Yfin)<1.0E-10)
      {
          phiPart=0.0;
      }
      else
      {
          phiPart=atan2(Yfin,Xfin);
          //Check if the angular value is negative
          if(std::signbit(phiPart))
          {
              phiPart+=2.0*M_PI;
          }
      }
      
    
      //printf("phiPart = %14.5E rad  %14.5E graus\n", phiPart, phiPart*180.0/M_PI);
      
     
     
     unsigned long binPixel;
     if(phiMin < phiMax && (phiPart >= phiMin && phiPart < phiMax)){
            binPixel = (phiPart - phiMin)/phipx;
     }
     else if(phiMin > phiMax && (phiPart >= phiMin || phiPart < phiMax))
     {
         if(phiPart >= phiMin){
            binPixel = (phiPart - phiMin)/phipx;
         }
         else{
            binPixel = (phiPart + 2.0*M_PI - phiMin)/phipx; 
         }
     }
     else{
        //Exit from function
        return false;
     }
     
     
    iphi = phiBin;
    ipix  = binPixel;
     
    return true; //The particle crosses the detector
 }
    return false;
}

void pen_CTsinogram::countDetector(const unsigned long long nhist, 
                                   const pen_particleState& state){
    unsigned long iphi, ipix;
    if(detIndexes(state,lastPos,lastPage,iphi,ipix)){
        unsigned long sinoBin = ipix+iphi*npixs;
        
        if(nhist > lasthist[sinoBin]){
            sino[sinoBin] += sinotemp[sinoBin];
            sino2[sinoBin] += sinotemp[sinoBin]*sinotemp[sinoBin];
            sinotemp[sinoBin] = state.WGHT;
            lasthist[sinoBin] = nhist + 0.5;
        }
        else{
            sinotemp[sinoBin] += state.WGHT;
        }        
    }
    
}

void pen_CTsinogram::countSource(const unsigned long long nhist, 
                                 const pen_particleState& state){
    unsigned long iphi, ipix;
    CTsinogram_point origin;
    origin.X = state.X;
    origin.Y = state.Y;
    origin.Z = state.Z;
    if(detIndexes(state,origin,state.PAGE,iphi,ipix)){
        unsigned long sinoBin = ipix+iphi*npixs;
        
        if(nhist > lasthistNorm[sinoBin]){
            sinoNorm[sinoBin] += sinotempNorm[sinoBin];
            sino2Norm[sinoBin] += sinotempNorm[sinoBin]*sinotempNorm[sinoBin];
            sinotempNorm[sinoBin] = state.WGHT;
            lasthistNorm[sinoBin] = nhist + 0.5;
        }
        else{
            sinotempNorm[sinoBin] += state.WGHT;
        }        
    }
    
}

void pen_CTsinogram::tally_localEdep(const unsigned long long nhist,
				  const pen_KPAR kpar,
				  const pen_particleState& state,
				  const double dE){

  if(iproj >= 0)
    edepMat[iproj].tally_localEdep(nhist,kpar,state,dE);
}

void pen_CTsinogram::tally_beginPart(const unsigned long long nhist,
				     const unsigned kdet,
				     const pen_KPAR kpar,
				     const pen_particleState& state){
  inGeo = true;
  
  iproj = getProjection(state);
  if(iproj >= 0)
    edepMat[iproj].tally_beginPart(nhist,kdet,kpar,state);  
}

void pen_CTsinogram::tally_sampledPart(const unsigned long long nhist,
				       const unsigned long long dhist,
				       const unsigned kdet,
				       const pen_KPAR kpar,
				       const pen_particleState& state){
    
  iproj = getProjection(state);
  if(iproj >= 0)
    edepMat[iproj].tally_sampledPart(nhist,dhist,kdet,kpar,state);
}

void pen_CTsinogram::tally_step(const unsigned long long nhist,
			     const pen_KPAR kpar,
			     const pen_particleState& state,
			     const tally_StepData& stepData){
  if(iproj >= 0)
    edepMat[iproj].tally_step(nhist,kpar,state,stepData);    
}
void pen_CTsinogram::tally_move2geo(const unsigned long long nhist,
				 const unsigned kdet,
				 const pen_KPAR kpar,
				 const pen_particleState& state,
				 const double dsef,
				 const double dstot){
  iproj = getProjection(state);
  if(iproj >= 0){
    edepMat[iproj].tally_move2geo(nhist,kdet,kpar,state,dsef,dstot);
    if(state.E > emin && state.E <= emax && kpar == part){
        countSource(nhist, state);
        moved2geo = true;
    }
  }
}

void pen_CTsinogram::tally_endHist(const unsigned long long nhist){
  if(iproj >= 0)
    edepMat[iproj].tally_endHist(nhist);    
}

void pen_CTsinogram::tally_jump(const unsigned long long /*nhist*/,
				const pen_KPAR kpar,
				const pen_particleState& state,
				const double /*ds*/){

  if(kpar == part){
    //Save last particle relevant state
    lastPos.X = state.X;
    lastPos.Y = state.Y;
    lastPos.Z = state.Z;
    lastPage  = state.PAGE;
  }
  
}

void pen_CTsinogram::tally_knock(const unsigned long long /*nhist*/,
		 const pen_KPAR /*kpar*/,
		 const pen_particleState& /*state*/,
		 const int /*icol*/){
    knocked = true;
}

void pen_CTsinogram::tally_endPart(const unsigned long long nhist,
				   const pen_KPAR kpar,
				   const pen_particleState& state){

  if(!scatter){ 
    //Check if the particle has been scattered
    if(knocked || !moved2geo){
        knocked = false;
        moved2geo = false;
        inGeo = false;
        lastPos.X = lastPos.Y = lastPos.Z = lastPage = 1.0e35;
        //Particle has been scattered, so we don't add it to the detector counter
        //clear the boolean variables and inicializate last positions, and exit
        return;
    }
    //Clear boolean variables, check has been set before
    knocked = false;
    moved2geo = false;
  }
    
  //Check if the particle reaches the geometry
    
  if(!inGeo)
  return;
  
  //Particle reaches de geometry system. Clear the inGeo flag
  //for next particle and check particle type
  inGeo = false;
 
  if(kpar != part){
    lastPos.X = lastPos.Y = lastPos.Z = lastPage = 1.0e35;
    return;
  }

  //This particle is of the desired type,
  //check if it has escaped from geometry
  if(state.MAT == 0){
    //Particle escaped from the geometry.
    //Check if is in the energy window
    if(state.E > emin && state.E <= emax)
      {
        //Count the particle
        countDetector(nhist,state);
      }
  }

  
}



int pen_CTsinogram::configure(const wrapper_geometry& geometry,
			      const abc_material* const materials[constants::MAXMAT],
			      const pen_parserSection& config,
			      const unsigned verbose
			      ){
  int err;
  
  //Init last position outside the CT cylinder
  lastPos.X = 1.0e35;
  lastPos.Y = 1.0e35;
  lastPos.Z = 1.0e35;
  lastPage  = 1.0e35;
  inGeo = false;
 
  //Energies
  //****************//
   
  err = config.read("emin", emin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'emin' in configuration. Double expected\n");
    }
    return -1;
  }


  // Maximum energy
  //***************//    
    
  err = config.read("emax", emax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read "
	     "'emax' in configuration. Double expected\n");
    }
    return -2;
  }
    
  // Number of pixels of the detector
  //*********************************//    
  err = config.read("npixels", npixs);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read "
	     "'npixels' in configuration. Integrer expected\n");
    }
    return -3;
  }
    
    
  if(npixs > npxmax){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Invalid entry, "
	     "npixels must be lower than %d.\n", npxmax);
    }
    return -4;
  }
    
  if(npixs < 1){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Invalid entry, npixels must be 1 at least.\n");
    }
    return -5;
  }
    
  if(verbose > 1){
    printf("Sinogram energy limits [Emin,Emax] (eV)");
    printf(" %12.5E %12.5E \n\n",emin,emax);
  }


  

  err = config.read("pixel-depth", dz);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'pixel-depth' in configuration. Double expected\n");
    }
    return -6;
  }
  
  if(dz <= 0.0){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Invalid entry, pixel-depth must be greater than zero.\n");
    }
    return -7;
  }
  
  //Calcualte once the semi depth of the pixel
  semiDet = dz*0.5; 
    
  if(verbose > 1){
    printf("Pixel information height dimensions and no. pixels at the detector:\n");
    printf(" %12.5E %d \n\n",dz,npixs);
  }


  // Detector radius
  //****************//    
    
  err = config.read("r-inner", ri);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'r-inner' in configuration. Double expected\n");
    }
    return -8;
  }
  
   if(ri <= 0.0)
    {
      if(verbose > 0){
	printf("pen_CTsinogram:configure: Error: ri value must be greater than zero.\n");
      }
      return -9;
    }
    
  //Calculate once the constant value of ri²
  ri2 = pow(ri,2); 

  /*
  err = config.read("r-outer", ro);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'r-outer' in configuration. Double expected\n");
    }
    return -10;

  }
    
  if(ri >= ro)
    {
      if(verbose > 0){
	printf("pen_CTsinogram:configure: Error: ri value must be lower than ro value.\n");
      }
      return -11;
    }
  */

  err = config.read("xOrigin", xOrigin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'xOrigin' in configuration. Double expected\n");
    }
    return -12;

  }

  err = config.read("yOrigin", yOrigin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'yOrigin' in configuration. Double expected\n");
    }
    return -13;

  }

  err = config.read("zOrigin", zOrigin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'zOrigin' in configuration. Double expected\n");
    }
    return -14;

  }

  if(verbose > 1){
    printf("Detector information inner radious:\n");
    printf(" %12.5E \n\n",ri);
  } 

  if(verbose > 1){
    printf("Detector origin (x, y, z):\n");
    printf(" %12.5E %12.5E %12.5E \n\n",xOrigin,yOrigin,zOrigin);
  }


  // Angular information
  //********************//

  err = config.read("phi-ini", phi0);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'phi-ini' in configuration. Double expected\n");
    }
    return -15;

  }

  
  err = config.read("phi-end", phif);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'phi-end' in configuration. Double expected\n");
    }
    return -16;

  }

  if(phi0 >= phif)
    {
      if(verbose > 0){
	printf("pen_CTsinogram:configure: Error: 'phi-ini' value must be lower than 'phi-end' value.\n");
      }
      return -17;
    }

  int auxProj;
  err = config.read("nProjections", auxProj);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'nProjections' "
	     "in configuration. Integer expected\n");
    }
    return -18;
    
  }

  if(auxProj <= 0){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Invalid number of projections %d. "
	     "One or more projections is required\n", auxProj);
    }
    return -18;
  }
  
  //Number of total phi positions/projections 
  // depending of the angula step
  nphi = static_cast<unsigned long>(auxProj);

  //Calculate angular step
  dphi = (phif-phi0)/static_cast<double>(nphi);

  if(verbose > 1){
    printf("Angular information angle initial, angle final and angular step for each projection:\n");
    printf(" %12.5E %12.5E %12.5E \n\n",phi0,phif,dphi);
  }

  //Detector size, as an angle
  err = config.read("aperture", phi);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'aperture' in configuration. Double expected\n");
    }
    return -19;

  }
  
  //Angular values change it to  radians
  dphi *= M_PI/180.0;
  phi0 *= M_PI/180.0;
  phif *= M_PI/180.0;
  phi  *= M_PI/180.0;

  // Some precalculate variables
  //******************************//
  
  //Value of each angular size of a pixel of the detector
  phipx=phi/double(npixs);
  //Inner arc of the detector
  si = phipx*ri*double(npixs);

  //Dimension of the sinogram, number of possible angles projections 
  //by the number of pixels of the detector that is the same in each projection
  sinoDim = nphi*npixs;



  

  
  sino     = (double*) malloc(sizeof(double)*sinoDim);
  sino2    = (double*) malloc(sizeof(double)*sinoDim);
  sinotemp = (double*) malloc(sizeof(double)*sinoDim);
  lasthist = (double*) malloc(sizeof(double)*sinoDim);

  sinoNorm     = (double*) malloc(sizeof(double)*sinoDim);
  sino2Norm    = (double*) malloc(sizeof(double)*sinoDim);
  sinotempNorm = (double*) malloc(sizeof(double)*sinoDim);
  lasthistNorm = (double*) malloc(sizeof(double)*sinoDim);
  
  if(verbose > 1){
    printf("Number of total angular positions (projections):, \n");
    printf(" %lu \n\n",nphi);
  }
  
  //Create a material deposition tally for each projection
  edepMat.resize(nphi);
  //Create a dummy configure section for edep tallies
  bool usedMats[constants::MAXMAT+1];
  geometry.usedMat(usedMats);
  nmat = 0;
  for(unsigned i = 0; i < constants::MAXMAT+1; ++i){
      if(usedMats[i]){
          nmat = i;
      }
  }
  if(nmat == 0){
      if(verbose > 0)
      printf("pen_CTsinogram:configure: Error: No materials used at "
	     " specified geometry.\n");
    return -20;          
  }
  else if(verbose > 1){
    printf("Number of geometry used materials: \n");
    printf(" %d \n\n",nmat);      
  }
  
  pen_parserSection configMat;
  configMat.set("nmat",nmat);
  for(auto& tally : edepMat){
      //Set tally name
      tally.setName(readName());
      if(tally.configure(geometry,materials,configMat,1) != 0){
        if(verbose > 0)
        printf("pen_CTsinogram:configure: Error at material deposition "
            "configuration.\n");
        return -21;                    
      }
  }
    
  //Particle that will be counted at the sinogram
  //*********************************************//

  std::string partName;
  err = config.read("particle", partName);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read particle name "
	     " ('particle') in configuration. String expected\n");
    }
    return -22;
  }

  part = static_cast<pen_KPAR>(particleID(partName.c_str()));
  if(!isKpar(part)){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unkown particle name %s\n",
	     partName.c_str());
    }
    return -23;    
  }

  // Time window
  //************//
  err = config.read("tmin", tmin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'tmin' in configuration. Double expected\n");
    }
    return -24;
  }

  
    err = config.read("timeInterval", dt);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'timeInterval' in configuration. Double expected\n");
    }
    return -25;
  }

  err = config.read("scattered", scatter);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_CTsinogram:configure: Error: Unable to read 'scattered' in configuration. Boolean expected\n");
    }
    return -26;
  }

  
  
  //Calculate maximum measure time
  tmax = tmin + dt*static_cast<double>(nphi);
  
  if(verbose > 1){
    printf("Time information,    tmin,     tmax,      dt:, \n");
    printf("                %12.5E  %12.5E  %12.5E \n\n",tmin,tmax,dt);
  }
    
 
  //Init arrays
  for(size_t j = 0; j < sinoDim; j++){        
    //Sinogram counters
    sino[j] = 0.0;
  } 
  for(size_t j = 0; j < sinoDim; j++){        
    //Sinogram counters
    sino2[j] = 0.0;
  } 
  for(size_t j = 0; j < sinoDim; j++){        
    //Sinogram counters
    sinotemp[j] = 0.0;
  } 
  for(size_t j = 0; j < sinoDim; j++){        
    //Sinogram counters
    lasthist[j] = 0.0;
  } 
  
  for(size_t j = 0; j < sinoDim; j++){        
    //Sinogram counters
    sinoNorm[j] = 0.0;
  } 
  for(size_t j = 0; j < sinoDim; j++){        
    //Sinogram counters
    sino2Norm[j] = 0.0;
  } 
  for(size_t j = 0; j < sinoDim; j++){        
    //Sinogram counters
    sinotempNorm[j] = 0.0;
  } 
  for(size_t j = 0; j < sinoDim; j++){        
    //Sinogram counters
    lasthistNorm[j] = 0.0;
  } 
  
    
  //Register data to dump    
  dump.toDump(sino,sinoDim);
  dump.toDump(sino2,sinoDim);
  dump.toDump(sinotemp,sinoDim);
  dump.toDump(lasthist,sinoDim);

  dump.toDump(sinoNorm,sinoDim);
  dump.toDump(sino2Norm,sinoDim);
  dump.toDump(sinotempNorm,sinoDim);
  dump.toDump(lasthistNorm,sinoDim);
    
  
  /*
  pen_particleState stateOrigin;
  stateOrigin.PAGE = 60.0;
  
  FILE* fout = fopen("testCTtally.dat","w");
  fprintf(fout, "# phiMin     phiMax      phiCentral      phiBin      state.PAGE\n");
  for(size_t i = 0; i < nphi; ++i){
      
    //Obtain the phi position of the detector
    long int phiBin = (stateOrigin.PAGE-tmin)/dt; 

    //Minimum and maximum angle of the detector. Angular detector window
    double phiMin = phi0 - phi*0.5 + dphi*phiBin;
    double phiMax = phi + phiMin;
    
    
    //Values of angles must be in the interval [0,360] degrees
    //No negative values are possible
    if(phiMin < 0.0)
        phiMin = 2.0*M_PI+phiMin;
    if(phiMax < 0.0)
        phiMax = 2.0*M_PI+phiMax;
    
    //No values larger than 360 deg are possible
    if(phiMin == 2.0*M_PI)
        phiMin=0.0;
    if(phiMax == 2.0*M_PI)
        phiMax=0.0;
    if(phiMin > 2.0*M_PI){
        double extra = phiMin/(2.0*M_PI)-1.0;
        phiMin = 2.0*M_PI*extra;
    }
    if(phiMax > 2.0*M_PI){
        double extra = phiMax/(2.0*M_PI)-1.0;
        phiMax = 2.0*M_PI*extra;
    }
    
    phiBin=int(phiBin);
    
        fprintf(fout," %12.4E   %12.4E    %12.4E    %12.4E    %12.4E",
            phiMin*180.0/M_PI,phiMax*180.0/M_PI,(phiMin+phiMax)*0.5*180.0/M_PI,phiBin,stateOrigin.PAGE);
        
        
    fprintf(fout,"\n");
    fprintf(fout,"\n");
    stateOrigin.PAGE +=dt;
  }
  fclose(fout);
  */
  return 0;

}

void pen_CTsinogram::tally_endSim(const unsigned long long /*nhist*/){
  
  flush();
}

void pen_CTsinogram::saveData(const unsigned long long nhist) const{

  const char* filenameSINOGRAM = "CTsinogram.dat";
  
  
  double numberHist = (static_cast<double>(nhist)/static_cast<double>(nphi));
  double invn = 1.0/numberHist;
  

     
  //*************//
  //  Sinogram   //
  //************//

  FILE* outSINO = NULL;
  outSINO = fopen(filenameSINOGRAM, "w");
  if(outSINO == NULL){
    
    printf(" *********************************************\n");
    printf(" pen_CTsinogram:saveData:ERROR: cannot open output "
	   "data sinogram file\n");
    printf(" *********************************************\n");
    return;
  }
     
  fprintf(outSINO, "#------------------------------------------------------------\n");
  fprintf(outSINO, "# PenRed: CT sinogram report\n");
  fprintf(outSINO, "# Units are counts/history for fluence and"
                    " no dimensions for mu*ds\n");
  fprintf(outSINO, "#\n");
        
    
  fprintf(outSINO, "# Energy       min    ,     max\n");
  fprintf(outSINO, "#         %11.5e     %10.5e \n", emin, emax);

  
  fprintf(outSINO,"#\n");
    
  fprintf(outSINO, "#     nBin    nProj    nPix    sourceFlu  : +-2sigma detectedFlu  : +-2sigma mu*ds(average)  : +-2sigma : \n");
  
    
    
    
    //Limits of bins of each projection ex: 0 to 672, 673 to 1345 etc
    unsigned long int initProjectionBin = 0; 
    unsigned long int finalProjectionBin = npixs;
    for(unsigned long i=0; i < nphi; i++)
    {
        // Initialize pixel number
        unsigned long pix = 0;
        for(unsigned long int j = initProjectionBin; j < finalProjectionBin ; j++){
            double q, q2,sigma;
            double qNorm, q2Norm,sigmaNorm;
            double sigmaMux;
            
            q = sino[j]*invn;
            q2 = sino2[j]*invn;
            sigma = (q2-(q*q))*invn;
            if(sigma > 0.0){ sigma = sqrt(sigma);}
            else{sigma = 0.0;}
            
            qNorm = sinoNorm[j]*invn;
            q2Norm = sino2Norm[j]*invn;
            sigmaNorm = (q2Norm-(qNorm*qNorm))*invn;
            if(sigmaNorm > 0.0){ sigmaNorm = sqrt(sigmaNorm);}
            else{sigmaNorm = 0.0;}
            
            //Estimate mu*x as ln(qNorm/q)
            //
            //  N=N0*e^(-mu*x) -> ln(N/N0) = -mu*x -> ln(N0/N) = mu*x
            //
            double mux = 0.0E0;
            if(q > 0.0E0){
                mux = qNorm/q;
                if(mux > 0.0){
                    mux = log(mux);
                } else{
                    mux = 0.0;
                }   
                sigmaMux = sqrt(pow(sigmaNorm/qNorm,2) + pow(sigma/q,2));
            }
            else{
                sigmaMux = 0.0;
            }
            
           //double sigmaMux = sqrt(pow(sigmaNorm/qNorm,2) + pow(sigma/q,2));

            fprintf(outSINO,"   %6lu  %6lu  %6lu  %12.5E    %12.5E  %12.5E    %12.5E  %12.5E    %12.5E\n",j, i, pix, qNorm, 2.0*sigmaNorm, q, 2.0*sigma, mux, 2.0*sigmaMux);
            pix++;
        }
        //Update limits of bins of each projection
        initProjectionBin = finalProjectionBin;
        finalProjectionBin += npixs;
        fprintf(outSINO,"\n");      
    }


        
  fprintf(outSINO,"\n");
    
    
  fclose(outSINO);    
  
  //Save material energy deposition
  FILE* foutEdep = nullptr;
  foutEdep = fopen("CTedep.dat","w");
  if(foutEdep == nullptr){    
    printf(" *********************************************\n");
    printf(" pen_CTsinogram:saveData:ERROR: cannot open output "
	   "data energy deposition file\n");
    printf(" *********************************************\n");
    return;
  }

  fprintf(foutEdep, "#------------------------------------------------------------\n");
  fprintf(foutEdep, "# PenRed: CT energy deposition report\n");
  fprintf(foutEdep, "# Units are eV/history\n");
  fprintf(foutEdep, "#\n");

  fprintf(foutEdep,"#  Projection  ");
  for(int i = 1; i <= nmat; ++i)
      fprintf(foutEdep,"|       Material   %3d       ",i);
  fprintf(foutEdep, "\n");
  
  fprintf(foutEdep,"#              ");
  for(int i = 1; i <= nmat; ++i)
    fprintf(foutEdep,"|   Energy        +-2sigma   ");
  fprintf(foutEdep, "\n");
  
  unsigned long projIndex = 0;
  double edep, twoSigma;
  double totalEdep[constants::MAXMAT];
  double total2Sigma[constants::MAXMAT];
  
  for(int i=0; i<nmat; ++i)
  {
      totalEdep[i] = 0.0;
      total2Sigma[i] = 0.0;
  }
  
  for(auto& tally : edepMat){
      tally.saveData(numberHist);
      
      FILE* fin = nullptr;
      fin = fopen("materialEnergyDeposition.dat", "r");
      if(fin == nullptr){    
            printf(" *********************************************\n");
            printf(" pen_CTsinogram:saveData:ERROR: cannot open input "
            "data energy deposition file\n");
            printf(" *********************************************\n");
            return;
      }
      fprintf(foutEdep,"      %3lu      ", projIndex); 
      
      for(int i=0; i<5; ++i)
      {
          fscanf(fin, "%*[^\n]\n");
      }
      
      for(int i=0; i<nmat; ++i)
      {
         fscanf(fin,"    %*d      %le      %le",&edep,&twoSigma);
         
         fprintf(foutEdep, "  %12.5E    %8.1E  ", edep, twoSigma); 
         totalEdep[i]    += edep;
         total2Sigma[i]  += twoSigma;
          
      }
      fclose(fin); 
      
      fprintf(foutEdep, "\n");
      
      ++projIndex;
  }
  
  fprintf(foutEdep, "# Total Edep \n"); 
  fprintf(foutEdep,"               ");
  for(int i = 0; i < nmat; ++i)
    fprintf(foutEdep,"  %12.5E    %8.1E  ",totalEdep[i], total2Sigma[i]);
  
  fprintf(foutEdep, "\n");
  
  
  fclose(foutEdep);
  
}
  




int pen_CTsinogram::sumTally(const pen_CTsinogram& tally){

    
  if(sinoDim != tally.sinoDim)
    return -1;

  if(edepMat.size() != tally.edepMat.size())
    return -2;    
   
  //Sinogram
  //*********// 
  for(unsigned int j = 0; j < sinoDim; ++j){
    sino[j] += tally.sino[j];
  }
  for(unsigned int j = 0; j < sinoDim; ++j){
    sino2[j] += tally.sino2[j];
  }

  //Input fluence
  //***************// 
  for(unsigned int j = 0; j < sinoDim; ++j){
    sinoNorm[j] += tally.sinoNorm[j];
  }
  for(unsigned int j = 0; j < sinoDim; ++j){
    sino2Norm[j] += tally.sino2Norm[j];
  }
  
  
  //Material energy deposition
  for(size_t i = 0; i < edepMat.size(); ++i){
    edepMat[i].sumTally(tally.edepMat[i]);
  }
  
  return 0;
}

REGISTER_COMMON_TALLY(pen_CTsinogram)
