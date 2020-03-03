
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


#include "tallyImpactDetector.hh"

void pen_ImpactDetector::flush(){
    
  //*******************************
  //* Particle dependent counters *
  //*******************************
       
  for(unsigned int k = 0; k < constants::nParTypes; k++)
    {
      if(fln){
        //* Fluence particle counters    
        //*******************************
                  
        for(int i = 0; i < nbin; i++) {          
            
	  if(flutmp[k][i] == 0.0){ continue;}  // Skip void counters
	  flu[k][i]    += flutmp[k][i];
	  flu2[k][i]   += flutmp[k][i]*flutmp[k][i];
	  flutmp[k][i] = 0.0;
        }
      }
       
      if(spc){ 
        //* Spectrum particle counters    
        //*******************************
        
        for(int i = 0; i < nbin; i++) {          
	  //Spectrum particle counters
	  if(spectrumTmp[k][i] == 0.0){continue;}  // Skip void counters
	  spectrum[k][i] += spectrumTmp[k][i];
	  spectrum2[k][i]   += spectrumTmp[k][i]*spectrumTmp[k][i];
	  spectrumTmp[k][i] = 0.0;
        }
      }
    }

  //*******************
  //* Global counters *
  //*******************

  if(fln){
    //* Fluence
    //***********
        
    for(int i = 0; i < nbin; i++){          
      //Fluence particle counters      
      if(flutmptotal[i] == 0.0){continue;}  // Skip void counters
      flutotal[i] += flutmptotal[i];
      flu2total[i] += flutmptotal[i]*flutmptotal[i];
      flutmptotal[i] = 0.0;	
    }
  }

  if(spc){
    //* Spectrum
    //***********
    
    for(int i = 0; i < nbin; i++){
      //Spectrum particle counters
      if(spectrumTotalTmp[i] == 0.0){continue;} // Skip void counters
      spectrumTotal[i] += spectrumTotalTmp[i];
      spectrum2Total[i] += spectrumTotalTmp[i]*spectrumTotalTmp[i];
      spectrumTotalTmp[i] = 0.0;
    }
  }

  if(ageActive){
    //* Age
    //***********
    
    for(int i = 0; i < nbinAge; i++){
      //Age counters
      if(agetmp[i] == 0.0){continue;}  //Skip void counters
      age[i] += agetmp[i];
      age2[i] += agetmp[i]*agetmp[i];
      agetmp[i] = 0.0;
    }
  }

  if(eDepActive){
    //* Deposited energy
    //*******************
    
    //Deposited energy at detector
    edepCounter += edepCounterTmp;
    edep2Counter += edepCounterTmp*edepCounterTmp;
        
    if(edepCounterTmp >1.0e-5)
      {
	int binEdep;
	if(!isLinScaleDep)
	  {
	    binEdep = (log(edepCounterTmp) - emin)*iebin;
	  }
	else
	  {
	    binEdep = (edepCounterTmp - emin)*iebin;
	  }
            
	//Check if we are inside of energy bin range
	if(binEdep >= 0 && binEdep < nbin)
	  {
	    detected[binEdep] += 1.0;
	  }
      }  
    
    edepCounterTmp = 0.0;
    
        
  }
}


//Fluence distribution of particles within the
// detector. Discrete collisions only.
void pen_ImpactDetector::discreteTrackL(const pen_KPAR kpar,
					const double ds,
					const double energy){

  //Check energy to avoid count energies in the range (emin-ebin,emin)
  if(energy < configEmin)
    return;
  
  int binhigh;
  //Energy grid in no Linear Scale
  if(!isLinScaleFlu)
    {
      //Lowest E bin that scores
      binhigh = int((log(energy) - emin)*iebin);
    }
  else
    {
      //Lowest E bin that scores
      binhigh = int((energy - emin)*iebin);
    }
        
  if(binhigh >= 0 && binhigh < nbin)
    {   
      flutmptotal[binhigh] += ds;
      flutmp[kpar][binhigh] += ds;            
    }
  return;  
}


//Fluence distribution of particles within the
//detector. Continuous slowing down
void pen_ImpactDetector::continuousTrackL(const pen_KPAR kpar,
					  const pen_particleState& state,
					  const double ds,
					  const double de,
					  const double energy){

  //Check energy to avoid count energies in the range (emin-ebin,emin)
  if(energy < configEmin)
    return;
  
  int binhigh, binlow;
    
  double eContinuous = energy;
  if(eContinuous > emax)
    {
      eContinuous = emax;
    }
  double eFinal = energy - de;
  if(eFinal < emin)
    {
      eFinal = emin;
    }
  //Out of energy range
  if(eFinal > eContinuous){return;}
    
  //Energy grid in no Linear Scale
  if(!isLinScaleFlu)
    {
      binhigh = int((log(eContinuous) - emin)*iebin);
      binlow = int((log(eFinal) - emin)*iebin);
    }
  else
    {
      binhigh = int((eContinuous - emin)*iebin);
      binlow = int((eFinal - emin)*iebin);
    }
    
  double fact = ds/de;
    
  for(int i = binlow; i <= binhigh; i++)
    {
      double EA = emin + double(i)*ebin;
      EA = EA < eFinal ? eFinal : EA;

      double EB = emin + double(i+1)*ebin;
      EB = EB > eContinuous ? eContinuous : EB;
        
      double trackBin = (EB-EA)*fact;
      flutmptotal[i] += state.WGHT*trackBin;	  	
      flutmp[kpar][i] += state.WGHT*trackBin;	
    }
}




//Energy spectrum of entering particles
void pen_ImpactDetector::countSpectrum(const pen_KPAR kpar,
					   const pen_particleState& state){
  
  int binE;
  //Calculate energy bin
  if(!isLinScaleSpc)
    {
      binE = (log(state.E) - emin)*iebin;
    }
  else
    {
      binE = (state.E - emin)*iebin;
    }
    
  //Check if we are inside of energy bin range
  if(binE >= 0 && binE < nbin)
    {
      //Total counters
      spectrumTotalTmp[binE] += state.WGHT;
      
      //Particle kpar counters
      spectrumTmp[kpar][binE] += state.WGHT;
    }
}



//Age distribution
void pen_ImpactDetector::ageSpectrum(const pen_particleState& state){

   
  //Check energy to avoid count energies in the range (emin-ebin,emin)
  if(state.E < configEmin)
    return; 
    
  //Check age to avoid count ages in the range (ageMin-ageBinW,ageBinW)
  if(state.PAGE < configAgeMin)
    return;
  
  int binAge;
  //Calculate energy bin
  if(!isLinScaleAge)
    {
      if(state.PAGE <= 0.0)
	return;
      binAge = (log(state.PAGE) - ageMin)*iageBinW;
    }
  else
    {
      binAge = (state.PAGE - ageMin)*iageBinW;
    }

  if(binAge >= 0 && binAge < nbinAge)
    {
      agetmp[binAge] += state.WGHT; 
    }
}


void pen_ImpactDetector::tally_step(const double /*nhist*/,
				    const pen_KPAR kpar,
				    const pen_particleState& state,
				    const tally_StepData& stepData){
                  

  //Check if the flight is inside detection material
  if(inside){

    edepCounterTmp += stepData.softDE*state.WGHT;

    if(fln){
      if(stepData.softDE > 1.0e-16)
	{
	  // Store deposited energy
	  double eini = state.E+stepData.softDE;
	  continuousTrackL(kpar,state,stepData.dsef,stepData.softDE,eini);
	}
      else
	{
	  discreteTrackL(kpar,stepData.dsef*state.WGHT,state.E);
	}
    }
  }
}


void pen_ImpactDetector::tally_localEdep(const double /*nhist*/,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state,
		       const double dE){
    
    if(!inside){return;}

    // Store deposited energy
    edepCounterTmp += dE*state.WGHT;
    
}

void pen_ImpactDetector::tally_interfCross(const double /*nhist*/,
					       const unsigned kdet,
					       const pen_KPAR kpar,
					       const pen_particleState& state){
        
  //Check if, after crossing interface, we are in the detector now
  if(idet == kdet){
    if(!inside){
      //Check if is in energy window
      if(state.E >= configEmin && state.E < configEmax){
	if(spc){
	  countSpectrum(kpar, state);
	}
	if(ageActive){
	  ageSpectrum(state);
	}
      }
      inside = true;
    }
  }
  else{
    inside = false;
  }
}


void pen_ImpactDetector::tally_move2geo(const double /*nhist*/,
					    const unsigned kdet,
					    const pen_KPAR kpar,
					    const pen_particleState& state,
					    const double /*dsef*/,
					    const double /*dstot*/){
    
  //Check if, after move the particle to the geometry, we are in the detector now
  if(idet == kdet){
    
    //Check if is in energy window
    if(state.E >= configEmin && state.E < configEmax){
      if(spc){
	countSpectrum(kpar, state);
      }
      if(ageActive){
	ageSpectrum(state);
      }
    }
    // Store deposited energy
    edepCounterTmp += state.E*state.WGHT;
  
    inside = true;
  }
  else{
    inside = false;
  }
}


void pen_ImpactDetector::tally_beginPart(const double /*nhist*/,
					     const unsigned kdet,
					     const pen_KPAR /*kpar*/,
					     const pen_particleState& state){
        
  if(idet == kdet)   //Inside det body
    {
        inside = true;   //Set flag
       
        // Store deposited energy
        edepCounterTmp -= state.E*state.WGHT;
       
    }
  else                 //Outside det material
    {
      inside = false;
    }
}




void pen_ImpactDetector::tally_endHist(const double /*nhist*/){
    
    flush();
}





int pen_ImpactDetector::configure(const wrapper_geometry& /*geometry*/,
				      const abc_material* const /*materials*/[constants::MAXMAT],
				      const pen_parserSection& config,
				      const unsigned verbose
				      ){
  int err;
  

  //Activate fluence
  //***************
  err = config.read("fluence", fln);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'fluence' activated. Results of fluence are disabled. \n\n");
    }
    fln = false;
  }
  
  
  //Activate spectrum
  //***************
  err = config.read("spectrum", spc);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'spectrum' activated. Results of energy spectrum are disabled. \n\n");
    }
    spc = false;
  }
  
  
  //Activate age
  //***************
   err = config.read("age", ageActive);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'age' active. Results of age distribution are disabled\n\n");
    }
    ageActive = false;
  }
  
  
  
  //Activate energy deposition detector
  //***************************************
   err = config.read("energy-dep", eDepActive);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'energy-dep' active. Results of energy deposition distribution are disabled\n\n");
    }
    eDepActive = false;
  }
    
  
 // Detector
 //*************** 
  int auxkDet;
  err = config.read("detector", auxkDet);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_ImpactDetector:configure: Error: Unable to read 'detector' in configuration. Integrer expected\n");
    }
    return -1;
  }

  if(auxkDet <= 0){
    if(verbose > 0){
      printf("pen_ImpactDetector:configure: Error: Index of detector must be greater than 0.\n");
    }
    return -2;
  }
    
  idet = (unsigned)auxkDet;
  if(verbose > 1){
    printf("Detector:\n");
    printf(" %u \n\n",idet);
  }
 
 
 
 //Energies
 //****************
 
 if(spc || fln || eDepActive){
  // Minimum energy
  //***************    
  err = config.read("emin", emin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_ImpactDetector:configure: Error: Unable to read 'emin' in configuration. Double expected\n");
    }
    return -4;
  }

  // Maximum energy
  //***************    
    
  err = config.read("emax", emax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_ImpactDetector:configure: Error: Unable to read 'emax' in configuration. Double expected\n");
    }
    return -5;
  }
    
  // Number of bins
  //***************    
    
  err = config.read("nbin-energy", nbin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_ImpactDetector:configure: Error: Unable to read 'nbin' in configuration. Integrer expected\n");
    }
    return -6;
  }
    
    
  if(nbin > nbinmax){
    if(verbose > 0){
      printf("pen_ImpactDetector:configure: Error: Invalid entry, nbin must be lower than %d.\n", nbinmax);
    }
    return -7;
  }
    
  if(nbin < 1){
    if(verbose > 0){
      printf("pen_pen_ImpactDetector:configure: Error: Invalid entry, nbin must be 1 at least.\n");
    }
    return -8;
  }
    
  if(verbose > 1){
    printf("Spectrum limits [Emin,Emax] (eV) and no. bins:\n");
    printf(" %12.5E %12.5E %d \n\n",emin,emax,nbin);
  }

  //Save original emin
  configEmin = emin;
  configEmax = emax;
 }
    
 // Age interval
 //***************  
 if(ageActive)
   {
     err = config.read("age-min", ageMin);
     if(err != INTDATA_SUCCESS){
       if(verbose > 0){
	 printf("pen_ImpactDetector:configure: Error: Unable to read 'age-min' in configuration. Double expected\n");
       }
       return -10;
     }
    
    
     err = config.read("age-max", ageMax);
     if(err != INTDATA_SUCCESS){
       if(verbose > 0){
	 printf("pen_ImpactDetector:configure: Error: Unable to read 'age-max' in configuration. Double expected\n");
       }
       return -11;
     }
 
 
     if(ageMin < 1.0e-20)
       {
	 ageMin = 1.0e-20;
       }
    
     if(ageMin < 0.0 || ageMax < 0.0)
       {
	 if(verbose > 0){
	   printf("pen_ImpactDetector:configure: Error: Age values must be greater than zero.\n");
	 }
	 return -12;
       }

    
     if(ageMin >= ageMax)
       {
	 if(verbose > 0){
	   printf("pen_ImpactDetector:configure: Error: age-min value must be lower than age-max value.\n");
	 }
	 return -13;
       }
        
     if(ageMax < ageMin + 1.0e-19)
       {
	 if(verbose > 0){
	   printf("pen_ImpactDetector:configure: Error: Incorrect age limits. Age window must be at least 1.0e-19.\n");
	 }
	 return -14;
       }
    
     // Number of age bins
     //***************    
    
     err = config.read("nbin-age", nbinAge);
     if(err != INTDATA_SUCCESS){
       if(verbose > 0){
	 printf("pen_ImpactDetector:configure: Error: Unable to read 'nbin-age' in configuration. Integrer expected\n");
       }
       return -15;
     }
    
    
     if(nbinAge > nbinmax){
       if(verbose > 0){
	 printf("pen_ImpactDetector:configure: Error: Invalid entry, nbin- age must be lower than %d.\n", nbinmax);
       }
       return -16;
     }
    
     if(nbinAge < 1){
       if(verbose > 0){
	 printf("pen_pen_ImpactDetector:configure: Error: Invalid entry, nbin-age must be 1 at least.\n");
       }
       return -17;
     }
     if(verbose > 1){
       printf("Age distribution is enable. Result file with age distribution data will be shown.\n");
            
       printf("Particle age window [ageMin,ageMax] (seconds) and no. bins:\n");
       printf(" %12.5E %12.5E %d\n",ageMin, ageMax, nbinAge);
  
     }

     //Save configuration ageMin
     configAgeMin = ageMin;
   }
        

  // Scale type
  //***************    
   
  //Fluence
  if(fln){
  err = config.read("linearScale-fln", isLinScaleFlu);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'linearScale-fln' specified, linear scale will be set.\n");
    }
    isLinScaleFlu = true;
  }
    
    
  if(isLinScaleFlu){
    if(verbose > 1){
      printf("Linear scale for fluence. \n\n");
    }  
  }
  else{
    if(verbose > 1){
      printf("Logarithmic scale for fluence. \n\n");
    }
    //Default values due to the Logarithmic scale
    if(emin == 0.0){
      emin = 50.0;
    }

    // Store values in log scale
    emin = log(emin); 
    emax = log(emax);
  }
}
  
  
  
  //Spectrum
  if(spc){
    err = config.read("linearScale-spc", isLinScaleSpc);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'linearScale-spc' specified, linear scale will be set.\n");
    }
    isLinScaleSpc = true;
  }
    
    
  if(isLinScaleSpc){
    if(verbose > 1){
      printf("Linear scale for spectrum. \n\n");
    }  
  }
  else{
    if(verbose > 1){
      printf("Logarithmic scale for spectrum. \n\n");
    }
    //Default values due to the Logarithmic scale
    if(emin == 0.0){
      emin = 50.0;
    }

    // Store values in log scale
    emin = log(emin); 
    emax = log(emax);
  }
  
}

  //Energy deposition
  if(eDepActive){
    err = config.read("linearScale-edep", isLinScaleDep);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'linearScale-edep' specified, linear scale will be set.\n");
    }
    isLinScaleDep = true;
  }
    
    
  if(isLinScaleDep){
    if(verbose > 1){
      printf("Linear scale for spectrum. \n\n");
    }  
  }
  else{
    if(verbose > 1){
      printf("Logarithmic scale for spectrum. \n\n");
    }
    //Default values due to the Logarithmic scale
    if(emin == 0.0){
      emin = 50.0;
    }

    // Store values in log scale
    emin = log(emin); 
    emax = log(emax);
  }
  
}

  //Age
  if(ageActive){
    err = config.read("linearScale-age", isLinScaleAge);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'linearScale-age' specified, linear scale will be set.\n");
    }
    isLinScaleAge = true;
  }
    
    
  if(isLinScaleAge){
    if(verbose > 1){
      printf("Linear scale for age. \n\n");
    }  
  }
  else{
    if(verbose > 1){
      printf("Logarithmic scale for age. \n\n");
    }
    
    if(ageMin == 0.0){
      ageMin = 1.0e-20;
    }

    // Store values in log scale
    ageMin = log(ageMin);
    ageMax = log(ageMax);
      
  }
}
    
  //For energy bin width
  ebin = (emax-emin)/(double)nbin;
  iebin = 1.0/ebin;
    
  // Set energy grid
  for(int i = 0; i<nbin; i++){
    egrid[i] = emin+double(i)*ebin;
    ebingrd[i] = ebin;
  }
    
  //For age bin width
  ageBinW = (ageMax - ageMin)/(double)nbinAge;
  iageBinW = 1.0/ageBinW;
 
  // Clear counters:
  for(unsigned int i = 0; i < constants::nParTypes; i++){
    for(int j = 0; j < nbinmax; j++){
            
      //Fluence counters
      flu[i][j] = 0.0;
      flu2[i][j] = 0.0;
      flutmp[i][j] = 0.0;
            
      //Spectrum counters
      spectrum[i][j] = 0.0;
      spectrum2[i][j] = 0.0;
      spectrumTmp[i][j] = 0.0;
    }
  }
    
  for(int j = 0; j < nbinmax; j++){
        
    //Fluence counters
    flutmptotal[j] = 0.0;
    flutotal[j] = 0.0;
    flu2total[j] = 0.0;
        
    //Spectrum counters
    spectrumTotalTmp[j] = 0.0;
    spectrumTotal[j] = 0.0;
    spectrum2Total[j] = 0.0;
    
    //Energy deposited detector counters
    detected[j] = 0.0;
        
    //Age counters
    agetmp[j] = 0.0; 
    age[j] = 0.0;
    age2[j] = 0.0;
  } 
  
  
  //Energy deposited detector counters
  edepCounterTmp = 0.0;
  edepCounter = 0.0;
  edep2Counter = 0.0;
    
  //Register data to dump
   for(unsigned k = 0; k < constants::nParTypes; k++)
    {
        dump.toDump(flu[k],nbin);
        dump.toDump(flu2[k],nbin);
        dump.toDump(spectrum[k],nbin);
        dump.toDump(spectrum2[k],nbin);
    }
    
  dump.toDump(flutotal,nbin);
  dump.toDump(flu2total,nbin);
  dump.toDump(spectrumTotal,nbin);
  dump.toDump(spectrum2Total,nbin);
  
  dump.toDump(detected,nbin);
  dump.toDump(&edepCounter,1);
  dump.toDump(&edep2Counter,1);
  
  dump.toDump(age,nbinAge);
  dump.toDump(age2,nbinAge);
    
  return 0;

}

void pen_ImpactDetector::saveData(const double nhist) const{

  char filenameFLU[400];
  char filenameSPC[400];
  char filenameEDEP[400];
  char filenameAGE[400];
  
  sprintf(filenameFLU, "fluenceTrackLength-%u.dat",idet);
  sprintf(filenameSPC, "spectrum-impdet-%u.dat",idet);
  sprintf(filenameEDEP,"energyDeposition-impdet-%u.dat",idet);
  sprintf(filenameAGE, "age-impdet-%u.dat",idet);
  
  
  double fact, emiddle, invn;
  double qFLU, sigmaFLU,
    qSPC, sigmaSPC,
    qEDEP, sigmaEDEP,
    qAGE, sigmaAGE;
    
  invn = 1.0/(double)nhist;
    
if(fln){
  FILE* outFLU = NULL;
        
  //Fluence
  outFLU = fopen(filenameFLU, "w");
  if(outFLU == NULL){
    
    printf(" *********************************************\n");
    printf(" pen_ImpactDetector:saveData:ERROR: cannot open output data fluence file\n");
    printf(" *********************************************\n");
    return;
  }
    
   
  fprintf(outFLU,"#------------------------------------------------------------\n");
  fprintf(outFLU,"# PenRed: Integrated fluence spectrum\n");
  fprintf(outFLU,"#\n");
  fprintf(outFLU,"# The spectral fluence is integrated over the detector volume. Its units are cm/eV.\n");
  fprintf(outFLU,"#\n");
  fprintf(outFLU,"# Detector:\n");
  fprintf(outFLU,"#  %u\n",idet);
  fprintf(outFLU,"# Number of energy bins:\n");
  fprintf(outFLU,"#  %d\n", nbin);
  
  
  if (isLinScaleFlu){
    fprintf(outFLU,"# Linear energy scale, bin width (eV):\n");
    fprintf(outFLU,"#  %12.5E\n",ebin);
  }
  else{
    fprintf(outFLU,"# Log energy scale, E_{i+1}/E_i:\n");
    fprintf(outFLU,"#  %12.5E\n",exp(ebin));
  }
      
  fprintf(outFLU,"#\n");
  fprintf(outFLU,"#\n");
  fprintf(outFLU,"#\n");

  //Print headers
  fprintf(outFLU,"#%33.33s|%6.6s%15.15s%6.6s|",
	  " "," ", "Total fluence", " ");
  

  for(unsigned ipart = 0; ipart < constants::nParTypes; ipart++){
    fprintf(outFLU,"%6.6s%15.15s%6.6s|"," ", particleName(ipart), " ");
  }
  fprintf(outFLU,"\n");
  
  
  fprintf(outFLU,"#%33.33s|"," ");
  for(unsigned ipart = 0; ipart <= constants::nParTypes; ipart++){
      fprintf(outFLU, "%6.6s%15.15s%6.6s|"," ", " ", " ");
  }
  fprintf(outFLU, "\n");

  fprintf(outFLU,"# nbin :  Elow(eV)  : Emiddle(eV) | ");
  for(unsigned ipart = 0; ipart <= constants::nParTypes; ipart++){
    fprintf(outFLU,"fluence(cm/eV) : +-2sigma | ");
  }
  fprintf(outFLU,"\n");


  //Print data

  //*************
  // Fluence
  //*************
  
  for(int i = 0; i < nbin; i++){

    //Total Fluence
    qFLU = flutotal[i]*invn;
    sigmaFLU = (flu2total[i]*invn - qFLU*qFLU)*invn;
    sigmaFLU = sqrt((sigmaFLU > 0.0 ? sigmaFLU : 0.0));
    fact = 1.0/ebingrd[i];
    qFLU = qFLU*fact;
    sigmaFLU = sigmaFLU*fact;
        
    // Middle energy of the bin
    emiddle = egrid[i] + ebin*0.5;
          
    fprintf(outFLU," %5d  %12.5E %12.5E    %12.5E    %8.1E   ",i,egrid[i],emiddle,qFLU,2.0*sigmaFLU);

    //Particle fluence
    for(unsigned ipart = 0; ipart < constants::nParTypes; ipart++){    
      qFLU = flu[ipart][i]*invn;
      sigmaFLU = (flu2[ipart][i]*invn - qFLU*qFLU)*invn;
      sigmaFLU = sqrt((sigmaFLU > 0.0 ? sigmaFLU : 0.0));
      fact = 1.0/ebingrd[i];
      qFLU = qFLU*fact;
      sigmaFLU = sigmaFLU*fact;
        
      // Middle energy of the bin
      emiddle = egrid[i] + ebin*0.5;
          
      fprintf(outFLU," %12.5E    %8.1E   ",qFLU,2.0*sigmaFLU);
    }
    fprintf(outFLU,"\n");    
  }
    
  fclose(outFLU);
}
     
  //*************
  // Spectrum
  //*************
  if(spc){
    FILE* outSPC = NULL;
    outSPC = fopen(filenameSPC, "w");
    if(outSPC == NULL){
    
      printf(" *********************************************\n");
      printf(" pen_ImpactDetector:saveData:ERROR: cannot open output data spectrum file\n");
      printf(" *********************************************\n");
      return;
    }
     
    fprintf(outSPC, "#------------------------------------------------------------\n");
    fprintf(outSPC, "# PenRed: Energy spectrum report\n");
    fprintf(outSPC, "# detector number: %2d\n",idet);
    fprintf(outSPC, "# Units are 1/(eV*history)\n");
    fprintf(outSPC, "#\n");
        
    
    if(isLinScaleSpc){
      fprintf(outSPC, "# Linear scale\n");
      fprintf(outSPC, "# Energy     min    ,    max    , number of bins\n");
      fprintf(outSPC, "#        %11.5e %10.5e  %5d\n", emin, emax, nbin);
    }
    else{
      fprintf(outSPC, "# Logarithmic scale\n");
      fprintf(outSPC, "# Energy     min    ,    max    , number of bins\n");
      fprintf(outSPC, "#        %11.5e %10.5e  %5d\n", exp(emin), exp(emax), nbin);
    }
    
  
    fprintf(outSPC,"#\n");
    
    fprintf(outSPC, "#%13.13s|%5.5s%15.15s%5.5s|",
	    " "," ", "Total spectrum", " ");
   
    for(unsigned ipart = 0; ipart < constants::nParTypes; ipart++){
      fprintf(outSPC, "%5.5s%15.15s%5.5s|"," ", particleName(ipart), " ");
    }
    fprintf(outSPC, "\n");
  
    fprintf(outSPC,"#%13.13s|"," ");
    for(unsigned ipart = 0; ipart <= constants::nParTypes; ipart++){
      fprintf(outSPC, "%5.5s%15.15s%5.5s|"," ", " ", " ");
    }
    fprintf(outSPC, "\n");
  

    fprintf(outSPC, "# Energy (eV) |"); 
    for(unsigned ipart = 0; ipart <= constants::nParTypes; ipart++){
      fprintf(outSPC, "     prob    : +-2sigma  |"); 
    }
    fprintf(outSPC, "\n");
    
    
    for(int j = 0; j < nbin; j++)
      {
	double energy;
	double localiebin = iebin;
	if(isLinScaleSpc){
	  energy = emin + ((double)j+0.5)*ebin;
	}
	else{
	  double xLow = exp(emin +(double)j*ebin);
	  double xUp = exp(emin + (double)(j+1)*ebin);
	  energy = 0.5*(xUp + xLow);
	  localiebin = 1.0/(xUp - xLow);
	}
        
	fprintf(outSPC," %12.5E ",energy);

	//Total spectrum
	qSPC = spectrumTotal[j]*invn;
	double q2SPC = spectrum2Total[j]*invn;
	sigmaSPC = (q2SPC-qSPC*qSPC)*invn;
	if (sigmaSPC > 0.0){sigmaSPC = sqrt(sigmaSPC);}
	else{sigmaSPC = 0.0;}
    
	fprintf(outSPC," %12.5E   %8.3E ",
		qSPC*localiebin,2.0*sigmaSPC*localiebin);

	//Particle spectrum
	for(unsigned ipart = 0; ipart < constants::nParTypes; ipart++){
	
	  qSPC = spectrum[ipart][j]*invn;
	  q2SPC = spectrum2[ipart][j]*invn;
	  sigmaSPC = (q2SPC-qSPC*qSPC)*invn;
	  if (sigmaSPC > 0.0){sigmaSPC = sqrt(sigmaSPC);}
	  else{sigmaSPC = 0.0;}
    
	  fprintf(outSPC," %12.5E   %8.3E ",
		  qSPC*localiebin,2.0*sigmaSPC*localiebin);
	}
	
	fprintf(outSPC,"\n");
      }
    fclose(outSPC);
  }
    

//*******************
// Energy deposition
//*******************
  if(eDepActive){
      FILE* outEDEP = NULL;
      outEDEP = fopen(filenameEDEP, "w");
      
      if(outEDEP == NULL){
        printf(" *********************************************\n");
        printf(" pen_ImpactDetector:saveData:ERROR: cannot open output data energy deposition file\n");
        printf(" *********************************************\n");
        return;
      }
     
    fprintf(outEDEP, "#------------------------------------------------------------\n");
    fprintf(outEDEP, "# PenRed: Deposited energy spectrum\n");
    fprintf(outEDEP, "# detector number: %2d\n",idet);
    
    double eDepMean = edepCounter*invn;
    double sigmaeDepMean = fabs(edep2Counter - edepCounter*edepCounter)*invn;
    
    fprintf(outEDEP, "# Mean energy deposited: %11.5e  +- %8.3e\n", eDepMean, 2.0*invn*sqrt(sigmaeDepMean));
    fprintf(outEDEP, "# Units are 1/(eV*particle)\n");
    fprintf(outEDEP, "#\n");
    
    if(isLinScaleDep){
      fprintf(outEDEP, "# Linear scale\n");
      fprintf(outEDEP, "# Energy     min    ,    max    , number of bins\n");
      fprintf(outEDEP, "#       %11.5e  %10.5e    %5d\n", emin, emax, nbin);
    }
    else{
      fprintf(outEDEP, "# Logarithmic scale\n");
      fprintf(outEDEP, "# Energy     min    ,    max    , number of bins\n");
      fprintf(outEDEP, "#       %11.5e  %10.5e    %5d\n", exp(emin), exp(emax), nbin);
    }
    
    fprintf(outEDEP, "#\n");
    
    fprintf(outEDEP, "# deposited E (eV) :     prob     : +-2sigma\n");
    
    
    for(int j = 0; j < nbin; j++)
    {
        double energy;
	double localiebin = iebin;
        if(isLinScaleDep){
            energy = emin + ((double)j+0.5)*ebin;
        }
        else
        {
            double xLow = exp(emin +(double)j*ebin);
            double xUp = exp(emin + (double)(j+1)*ebin);
            energy = 0.5*(xUp + xLow);
            localiebin = 1.0/(xUp - xLow);
        }
        
        fprintf(outEDEP," %12.5E ",energy);
        
        //Deposited E
        qEDEP = detected[j]*invn;
        sigmaEDEP = fabs(qEDEP*(1.0 - qEDEP))*invn;
        if (sigmaEDEP > 0.0){sigmaEDEP = sqrt(sigmaEDEP);}
        else{sigmaEDEP = 0.0;}
    
        fprintf(outEDEP," %12.5E   %8.3E ",
		qEDEP*localiebin,2.0*sigmaEDEP*localiebin);
        
        fprintf(outEDEP,"\n");
    }
    
    fclose(outEDEP);    
  }
  
//*************
// Age
//*************
  if(ageActive){
    FILE* outAGE = NULL;
    outAGE = fopen(filenameAGE, "w");
    if(outAGE == NULL){
    
      printf(" *********************************************\n");
      printf(" pen_ImpactDetector:saveData:ERROR: cannot open output data age file\n");
      printf(" *********************************************\n");
      return;
    }
        
    fprintf(outAGE, "#------------------------------------------------------------\n");
    fprintf(outAGE, "# PenRed: Age distribution report\n");
    fprintf(outAGE, "# detector number: %2d\n",idet);
    fprintf(outAGE, "#\n");
    
   
    if(isLinScaleAge){
      fprintf(outAGE, "# Linear scale\n");
      fprintf(outAGE, "# Age        min    ,    max    , number of bins\n");
      fprintf(outAGE, "#        %11.5e %10.5e  %5d\n",ageMin,ageMax,nbinAge);
    }
    else
      {
	fprintf(outAGE, "# Logarithmic scale\n");
	fprintf(outAGE, "# Age        min    ,    max    , number of bins\n");
	fprintf(outAGE, "#        %11.5e %10.5e  %5d\n",exp(ageMin),exp(ageMax),nbinAge);	
      }
         

         
    fprintf(outAGE, "# Units are 1/(s*history)\n");      
    fprintf(outAGE, "#\n");
    fprintf(outAGE, "#    age (s)   :     prob     : +-2sigma\n");
        
    double sec = 0.0;
        
    for(int j = 0; j < nbinAge; j++)
      {
	double iagebin = iageBinW;
	if(isLinScaleAge)
	  {
	    sec = ageMin + ((double)j+0.5)*ageBinW;
	  }
	else
	  {
	    double xLow = exp(ageMin + (double)j*ageBinW);
	    double xUp = exp(ageMin + (double)(j+1)*ageBinW);
	    sec = 0.5*(xUp + xLow);
	    iagebin = 1.0/(xUp - xLow); 
	  }
            
	qAGE = age[j]*invn;
	double q2AGE = age2[j]*invn;
	sigmaAGE = (q2AGE - qAGE*qAGE)*invn;
	if(sigmaAGE > 0.0){sigmaAGE = sqrt(sigmaAGE);}
	else{sigmaAGE = 0.0;}
	fprintf(outAGE," %14.5E %14.5E %8.1E\n",
		sec,qAGE*iagebin,2.0*sigmaAGE*iagebin);
      }

    fclose(outAGE);
  }
}





int pen_ImpactDetector::sumTally(const pen_ImpactDetector& tally){

    
  if(nbin != tally.nbin)
    return -1;
  
  if(nbinAge != tally.nbinAge)
    return -2;
    
   
  //Spectrum
  //*********
  for(unsigned int i = 0; i < constants::nParTypes; ++i)
    {
      for(int j = 0; j < nbin; ++j){
	spectrum[i][j] += tally.spectrum[i][j];
      }
    }
    
  for(unsigned int i = 0; i < constants::nParTypes; ++i)
    {
      for(int j = 0; j < nbin; ++j){
	spectrum2[i][j] += tally.spectrum2[i][j];
      }
    }
    
    
  for(int j = 0; j < nbin; ++j){
    spectrumTotal[j] += tally.spectrumTotal[j];
  }
    
    
  for(int j = 0; j < nbin; ++j){
    spectrum2Total[j] += tally.spectrum2Total[j];
  }
    
    
    
  //Fluence
  //************
  for(unsigned int i = 0; i < constants::nParTypes; ++i)
    {
      for(int j = 0; j < nbin; ++j){
	flu[i][j] += tally.flu[i][j];
      }
    }
    
  for(unsigned int i = 0; i < constants::nParTypes; ++i)
    {
      for(int j = 0; j < nbin; ++j){
	flu2[i][j] += tally.flu2[i][j];
      }
    }
    
    
  for(int j = 0; j < nbin; ++j){
    flutotal[j] += tally.flutotal[j];
  }
  
    
  for(int j = 0; j < nbin; ++j){
    flu2total[j] += tally.flu2total[j];
  }

  for(int i = 0; i < nbin; ++i){
    detected[i] += tally.detected[i];
  }

  edepCounter += tally.edepCounter;
  edep2Counter += tally.edep2Counter;
    
  //Age
  //***********
  for(int i = 0; i < nbinAge; ++i){
    age[i] += tally.age[i];
  }
    
  for(int i = 0; i < nbinAge; ++i){
    age2[i] += tally.age2[i];
  }
    
  return 0;
  
}

REGISTER_COMMON_TALLY(pen_ImpactDetector, IMPACT_DET)
