
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

 
#include "tallyEmergingParticlesDistribution.hh"


void pen_EmergingPartDistrib::scapedParticle(const unsigned long long nhist,
					    const pen_KPAR kpar,
					    const pen_particleState& state){

  //Calculate the energy bin
  const double TWOPI = 2.0*constants::PI;
  int particleEbin;
    
  if(energyLogscale)
    {
      double auxE = log(state.E) - emin;
      //Out of range
      if(auxE < 0.0){particleEbin = -1;}  
      else{particleEbin = auxE*iebin;}
    }
  else
    {
      double auxE = state.E - emin;
      //Out of range
      if(auxE < 0.0){particleEbin = -1;}  
      else{particleEbin = auxE*iebin;}
        
    }
    
  double theta, phi;
    
  theta = acos(state.W);
  if(fabs(state.V) > 1.0e-16 || fabs(state.U) > 1.0e-16)
    {
      phi = atan2(state.V,state.U);
    }
  else
    {
      phi = 0.0;
    }
    
  if(phi < 0.0)
    {
      phi = TWOPI + phi;
    }
    
    
  int particleThetaBin,particlePhiBin;
    
  if(angularLogscale)
    {
      if(theta < 1.0e-12)
        {
	  theta = 1.0e-12;
        }
        
      particleThetaBin = int((log(theta)-tmin)/tbin);
    }
  else
    {
      particleThetaBin = int(theta/tbin);
    }

  particlePhiBin = int(phi/pbin);
    
  //Increment counters
    
  //Energy spectrum upbound and downbound
    
  //Check if the energy of the particle is in energy range
  if(state.E >= configEmin && particleEbin >= 0 && particleEbin < nBinsE) 
  {
      //Downbound case
      if(state.W <= 0.0)  
      {
        // Transfer partial to totals when a new history visits
        if(nhist > enlastDown[kpar][particleEbin])  
        {
	      //Energy-down
	      eHistDown[kpar][particleEbin]  += eHistTmpDown[kpar][particleEbin];
	      eHist2Down[kpar][particleEbin] += eHistTmpDown[kpar][particleEbin]*eHistTmpDown[kpar][particleEbin];
	      eHistTmpDown[kpar][particleEbin] = state.WGHT;
	      // Add 1/2 to avoid roundoff errors
          enlastDown[kpar][particleEbin] = nhist;    
            }
	  else
            {
	      eHistTmpDown[kpar][particleEbin] += state.WGHT;
            }
        }
      else   
        {
            //Upbound case
            // Transfer partial to totals when a new history visits
            if(nhist > enlastUp[kpar][particleEbin])  
            {
                //Energy-up
                eHistUp[kpar][particleEbin]  += eHistTmpUp[kpar][particleEbin];
                eHist2Up[kpar][particleEbin] += eHistTmpUp[kpar][particleEbin]*eHistTmpUp[kpar][particleEbin];
                eHistTmpUp[kpar][particleEbin] = state.WGHT;
                // Add 1/2 to avoid roundoff errors
                enlastUp[kpar][particleEbin] = nhist;    
            }
	  else
            {      
	      eHistTmpUp[kpar][particleEbin] += state.WGHT;
            }
        }
    }
    
  if(particleThetaBin >= nBinsTheta){particleThetaBin = nBinsTheta-1;}
  if(particlePhiBin >= nBinsPhi){particlePhiBin = nBinsPhi-1;}
  
  int angularBin = particleThetaBin*nBinsPhi + particlePhiBin; 
  if(nhist > angnlast[kpar][angularBin]) 
    {
        // Transfer partial to totals when a new history visits
        angHist[kpar][angularBin]  += angHistTmp[kpar][angularBin];
        angHist2[kpar][angularBin] += angHistTmp[kpar][angularBin]*angHistTmp[kpar][angularBin];
        angHistTmp[kpar][angularBin] = state.WGHT;
        // Add 1/2 to avoid roundoff errors
        angnlast[kpar][angularBin] = nhist;   
    }
  else
    {        
      angHistTmp[kpar][angularBin] += state.WGHT;
    }
}


void pen_EmergingPartDistrib::tally_move2geo(const unsigned long long nhist,
					    const unsigned /*kdet*/,
					    const pen_KPAR kpar,
					    const pen_particleState& state,
					    const double /*dsef*/,
					    const double /*dstot*/){
     
  if(state.MAT == 0){   
    //Particle scape the geometry
    //The particle doesn't reached the material system    
    scapedParticle(nhist,kpar,state);
  }
    
}
 
 
void pen_EmergingPartDistrib::tally_matChange(const unsigned long long nhist,
					     const pen_KPAR kpar,
					     const pen_particleState& state,
					     const unsigned /*prevMat*/){

  if(state.MAT == 0){  
    //Particle scape the geometry
    //The particle scaped from the material system    
    scapedParticle(nhist,kpar,state);    
  }  
}


void pen_EmergingPartDistrib::flush(){
    
  // Dump temp counters and obtain max score:
  for(unsigned k = 0; k < constants::nParTypes; k++)
    {
      for(int i = 0; i < nBinsE; i++)
        {
	  //Up counters
	  //*************
            
	  // Skip empty bins
	  if (enlastUp[k][i] == 0){ continue;}
	  // Transfer temp counter
	  eHistUp[k][i]  += eHistTmpUp[k][i];
	  eHist2Up[k][i] += eHistTmpUp[k][i]*eHistTmpUp[k][i];
	  // Reset counter
	  eHistTmpUp[k][i] = 0.0;
	  // Reset last visited to avoid recounting in next report
	  enlastUp[k][i] = 0;
	}

      for(int i = 0; i < nBinsE; i++)
	{	
	  //Down counters
	  //*************
            
	  // Skip empty bins
	  if (enlastDown[k][i] == 0){ continue;}
	  // Transfer temp counter
	  eHistDown[k][i]  += eHistTmpDown[k][i];
	  eHist2Down[k][i] += eHistTmpDown[k][i]*eHistTmpDown[k][i];
	  // Reset counter
	  eHistTmpDown[k][i] = 0.0;
	  // Reset last visited to avoid recounting in next report
	  enlastDown[k][i] = 0;
	}
        
      for(int j = 0; j < nAngBins; j++)
        {
	  //Angular counters
	  //*************
            
	  // Skip empty bins
	  if (angnlast[k][j] == 0){ continue;}
	  // Transfer temp counter
	  angHist[k][j]  += angHistTmp[k][j];
	  angHist2[k][j] += angHistTmp[k][j]*angHistTmp[k][j];
	  // Reset counter
	  angHistTmp[k][j] = 0.0;
	  // Reset last visited to avoid recounting in next report
	  angnlast[k][j] = 0;
        }
    }
    
}


void pen_EmergingPartDistrib::tally_endSim(const unsigned long long /*nhist*/){        
  flush();
}

int pen_EmergingPartDistrib::configure(const wrapper_geometry& /*geometry*/,
				      const abc_material* const /*materials*/[constants::MAXMAT],
				      const pen_parserSection& config,
				      const unsigned verbose){
                   
  int err;

    
  err = config.read("emin", emin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("EmergingPartDistrib:configure: Error: unable to read 'emin' in configuration. Double expected\n");
    }
    return -1;
  }
    
  err = config.read("emax", emax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("EmergingPartDistrib:configure: Error: unable to read 'emax' in configuration. Double expected\n");
    }
    return -2;
  }

  
  if(emin >= emax)
    {
      if(verbose > 0){
	printf("EmergingPartDistrib:configure: Error: minimum energy (%.5E eV) can't be greater than maximum energy (%.5E eV).\n", emin, emax);
      }
      return -3;
    }
    
    
  if(verbose > 1){
    printf("Energy interval:\n");
    printf(" %.5E - %.5E\n",emin,emax);
  }

  //Save configuration emin
  configEmin = emin;
    
  err = config.read("nBinsE", nBinsE);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("EmergingPartDistrib:configure: Error: unable to read 'nBinsE' in configuration. Integrer expected\n");
    }
    return -4;
  }
    
  if(nBinsE < 1)
    {
      if(verbose > 0){
	printf("EmergingPartDistrib:configure: Error: The number of energy bins can't be 0.\n");
      }
      return -5;   
    }
    
    
    
  if(verbose > 1)
    {
      printf("Energy bins:\n");
      printf(" %d\n",nBinsE);
    }
    
  if(nBinsE >= nbinmax)
    {
      if(verbose > 0){
        printf("EmergingPartDistrib:configure: Error: The maximum number of energy bins is %d.\n",nbinmax);
      }
      return -6;
    }
    
    
  //Check if user select log scale
    
  err = config.read("energyLogscale", energyLogscale);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'energyLogscale' specified, linear scale will be set.\n\n");
    }
    energyLogscale = false;
  }
    
  if(!energyLogscale)
    {
      if(verbose > 1)
        {
	  printf("Linear energy scale\n");
        }
    }
  else
    {
      if(verbose > 1)
        {
	  printf("Logarithmic energy scale\n");
        }

      //Store values in log scale
      emin = log(emin); 
      emax = log(emax); 
    }
    
  ebin = (emax - emin)/(double)nBinsE;
  iebin = 1.0/ebin;
    
    
  //Angular bins
    
  err = config.read("nBinsTheta", nBinsTheta);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("EmergingPartDistrib:configure: Error: unable to read 'nBinsTheta' in configuration. Integrer expected\n");
    }
    return -7;
  }
  

  if(nBinsTheta < 1)
    {
      if(verbose > 0)
        {
	  printf("EmergingPartDistrib:configure: Error: The number of polar bins must be greater than 0.\n");
        }
      return -8;
    }
    
  err = config.read("nBinsPhi", nBinsPhi);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("EmergingPartDistrib:configure: Error: unable to read 'nBinsPhi' in configuration. Integrer expected.\n");
    }
    return -9;
  }
    

  if(nBinsPhi <= 0)
    {
      if(verbose > 0)
        {
	  printf("EmergingPartDistrib:configure: Error:  The number of azimuthal must be greater than 0.\n");
        }
      return -10;
    }
    
  nAngBins = nBinsTheta*nBinsPhi;
  
    
  if(verbose > 1)
    {
      printf("Angular bins (Theta, Phi, Total):\n");
      printf(" %5d %5d %5d\n",nBinsTheta,nBinsPhi,nAngBins);
    }
           
    
  //Check if user select log scale
    
  err = config.read("angularLogscale", angularLogscale);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No 'angularLogscale' specified, linear scale will be set.\n\n");
    }
    angularLogscale = false;
  }

  
    if(!angularLogscale)
    {
        if(verbose > 1)
        {
            printf("Linear angular scale\n");
        }
        tmin = 0.0;
        tmax = constants::PI;
    }
    else
    {
        if(verbose > 1)
        {
            printf("Logarithmic angular scale\n");
        }
        //Store values in log scale
        tmin = log(1.0e-2);  
        tmax = log(constants::PI);

    }
    

    
  tbin = (tmax-tmin)/(double)nBinsTheta;
  pbin = (pmax-pmin)/(double)nBinsPhi;
    
  itbin = 1.0/tbin;
  ipbin = 1.0/pbin;
    
    
  //Init arrays
    
  for(unsigned int i = 0; i < constants::nParTypes; i++)
    {
      for(int j=0;j<nbinmax;j++)
        {
	  //Energy up
	  eHistUp[i][j]    = 0.0;
	  eHistTmpUp[i][j] = 0.0;
	  eHist2Up[i][j]   = 0.0;
	  enlastUp[i][j]   = 0;

	  //Energy down
	  eHistDown[i][j]    = 0.0;
	  eHistTmpDown[i][j] = 0.0;
	  eHist2Down[i][j]   = 0.0;	  
	  enlastDown[i][j]   = 0;
        }

      angHist[i] = (double*) calloc(nAngBins,sizeof(double));
      angHistTmp[i] = (double*) calloc(nAngBins,sizeof(double));
      angHist2[i] = (double*) calloc(nAngBins,sizeof(double));
      angnlast[i] = (unsigned long long*)
	calloc(nAngBins,sizeof(unsigned long long));

      double* p = angHist[i];
      for(int j = 0; j < nAngBins; ++j)
	p[j] = 0.0;
      p = angHistTmp[i];
      for(int j = 0; j < nAngBins; ++j)
	p[j] = 0.0;
      p = angHist2[i];
      for(int j = 0; j < nAngBins; ++j)
	p[j] = 0.0;
      unsigned long long* p2 = angnlast[i];
      for(int j = 0; j < nAngBins; ++j)
	p2[j] = 0ull;
      
    }
    
    
    //Register data to dump
    for(unsigned k = 0; k < constants::nParTypes; k++)
    {
        dump.toDump(eHistUp[k],nBinsE);
        dump.toDump(eHist2Up[k],nBinsE);
        dump.toDump(eHistDown[k],nBinsE);
        dump.toDump(eHist2Down[k],nBinsE);
        dump.toDump(angHist[k],nAngBins);
        dump.toDump(angHist2[k],nAngBins);
        dump.toDump(&nBinsE,1);  
        dump.toDump(&nAngBins,1);
    }
    
    return 0;
  
}
               


void pen_EmergingPartDistrib::saveData(const unsigned long long nhist) const{

  FILE* outUp = 0;
  FILE* outDown = 0;
  FILE* outAng = 0;
  FILE* outPolar = 0;
  int i;
  double invn;
  const double rad2deg = 180.0/constants::PI;

  invn = 1.0/static_cast<double>(nhist);
     
  outUp   = fopen("emerging-upbound.dat", "w");
  outDown = fopen("emerging-downbound.dat", "w");
  outAng = fopen("emerging-angle.dat", "w");
  outPolar = fopen("emerging-polar-angle.dat", "w");
     
  if(outUp == NULL || outDown == NULL || outAng == NULL || outPolar == NULL)
    {
      printf(" *********************************************\n");
      printf(" EmergingPartDistrib:saveData:ERROR: cannot open output data file\n");
      printf(" *********************************************\n");
      return;
    }
     
  fprintf(outUp, "#------------------------------------------------------------\n");
  fprintf(outDown, "#------------------------------------------------------------\n");
  fprintf(outAng, "#------------------------------------------------------------\n");
  fprintf(outPolar, "#------------------------------------------------------------\n");

  fprintf(outUp, "# PenRed: Emerging upbound particles distribution\n");
  fprintf(outDown, "# PenRed: Emerging downbound particles distribution\n");
  fprintf(outAng, "# PenRed: Emerging angular particles distribution\n");
  fprintf(outPolar, "# PenRed: Emerging polar particles distribution\n");
     
     
  //Energy report
     
  //Up
  fprintf(outUp, "# Units are particles/(eV history)\n");
  fprintf(outUp, "# Energy window (eV), energy bins:\n");
  if(energyLogscale)
    {
      fprintf(outUp, "# %12.5E %12.5E %5d\n",exp(emin),exp(emax),nBinsE);
    }
  else
    {
      fprintf(outUp, "# %12.5E %12.5E %5d\n",emin,emax,nBinsE);
    }
     
  fprintf(outUp, "#\n");
  fprintf(outUp, "# Energy (eV) :     e-     : +-2sigma :   photons  : +-2sigma :     e+     : +-2sigma :\n");  
     
  //Down
  fprintf(outDown, "# Units are particles/(eV history)\n");
  fprintf(outDown, "# Energy minimum, maximum, energy bins:\n");
  if(energyLogscale)
    {
      fprintf(outDown, "# %12.5E %12.5E %5d\n",exp(emin),exp(emax),nBinsE);
    }
  else
    {
      fprintf(outDown, "# %12.5E %12.5E %5d\n",emin,emax,nBinsE);
    }      
  fprintf(outDown, "#\n");
  fprintf(outDown, "# Energy (eV) :     e-     : +-2sigma :   photons  : +-2sigma :     e+     : +-2sigma :\n"); 
     
  double idbin = iebin;
  for(i = 0; i < nBinsE; i++)
    {
      double binE;
      if(energyLogscale)
	{
	  double binElow = exp(emin + double(i)*ebin);
	  double binEtop = exp(emin + double(i+1)*ebin);
	  binE = 0.5*(binEtop + binElow);
	  idbin = 1.0/(binEtop - binElow);
	}
      else
	{
	  binE = emin + (double(i) + 0.5)*ebin;
	}
         
      fprintf(outUp  ," %13.5E",binE);
      fprintf(outDown," %13.5E",binE);
         
      for(unsigned int j = 0; j < constants::nParTypes; j++)
	{
	  double qUp  = eHistUp[j][i]*invn;
	  double q2Up = eHist2Up[j][i]*invn;
	  double sigmaUp = (q2Up-qUp*qUp)*invn;
	  if(sigmaUp > 0.0){ sigmaUp = sqrt(sigmaUp);}
	  else{sigmaUp = 0.0;}
             
	  fprintf(outUp," %12.5E %10.3E", qUp*idbin, 2.0*sigmaUp*idbin);
             
	  double qDown  = eHistDown[j][i]*invn;
	  double q2Down = eHist2Down[j][i]*invn;
	  double sigmaDown = (q2Down-qDown*qDown)*invn;
	  if(sigmaDown > 0.0){ sigmaDown = sqrt(sigmaDown);}
	  else{sigmaDown = 0.0;}
             
	  fprintf(outDown," %12.5E %10.3E", qDown*idbin,2.0*sigmaDown*idbin);
	}
         
      fprintf(outUp  ,"\n");
      fprintf(outDown,"\n");
    }
     
  fclose(outUp);
  fclose(outDown);
     
     
     
     
  //Angular report
     
    
  fprintf(outAng, "# Units are particles/(sr history)\n");
  fprintf(outAng, "# Polar angle minimum, maximum, polar bins:\n");
  if(angularLogscale)
    {
      fprintf(outAng, "# %12.5E %12.5E %5d\n",exp(tmin)*rad2deg,exp(tmax)*rad2deg,nBinsTheta);
    }
  else
    {
      fprintf(outAng, "# %12.5E %12.5E %5d\n",tmin*rad2deg,tmax*rad2deg,nBinsTheta);
    }
     
  fprintf(outAng, "# Azimuthal angle minimum, maximum, azimuthal bins:\n");
  fprintf(outAng, "# %12.5E %12.5E %5d\n",pmin*rad2deg,pmax*rad2deg,nBinsPhi);
  fprintf(outAng, "#\n");
  fprintf(outAng, "# theta (deg) :  phi (deg) :     e-     : +-2sigma :   photons  : +-2sigma :     e+     : +-2sigma :\n");
  
     
     
  fprintf(outPolar, "# Units are particles/(sr history)\n");
  fprintf(outPolar, "#\n");
  fprintf(outPolar, "# Polar angle minimum, maximum, polar bins:\n");
  if(angularLogscale)
    {
      fprintf(outPolar, "# %12.5E %12.5E %5d\n",exp(tmin)*rad2deg,exp(tmax)*rad2deg,nBinsTheta);
    }
  else
    {
      fprintf(outPolar, "# %12.5E %12.5E %5d\n",tmin*rad2deg,tmax*rad2deg,nBinsTheta);
    }            
  fprintf(outPolar, "#\n");
  fprintf(outPolar, "# theta (deg) :     e-     : +-2sigma :   photons  : +-2sigma :     e+     : +-2sigma :\n"); 
     
     
  for(i = 0; i < nBinsTheta; i++)
    {
      double binTlow = tmin+double(i)*tbin;
      double binTtop = tmin+double(i+1)*tbin;
      if(angularLogscale)
	{
	  binTlow = exp(binTlow);
	  binTtop = exp(binTtop);
	}
        
      double binT = 0.5*(binTtop + binTlow);
        
      fprintf(outPolar, " %13.5E", binT*rad2deg);
        
      double qPolar[constants::nParTypes], q2Polar[constants::nParTypes]; 
        
      for(unsigned int j = 0; j < constants::nParTypes; j++)
	{
	  qPolar[j] = 0.0;
	  q2Polar[j] = 0.0;
	}
        
      double dsang = (cos(binTlow) - cos(binTtop))*pbin;
      double idsang = 1.0/dsang;
        
      for(int k = 0; k < nBinsPhi; k++)
	{
	  double binP = (k + 0.5)*pbin;
	  fprintf(outAng," %13.5E %12.5E", binT*rad2deg, binP*rad2deg);
        
	  for(unsigned int j = 0; j < constants::nParTypes; j++)
	    {		  
	      double qAng  = angHist[j][i*nBinsPhi+k]*invn;
	      double q2Ang = angHist2[j][i*nBinsPhi+k]*invn;

	      qPolar[j]  += qAng;
	      q2Polar[j] += q2Ang;
	      

	      double sigmaAng = (q2Ang-qAng*qAng)*invn;
	      if(sigmaAng > 0.0){ sigmaAng = sqrt(sigmaAng);}
	      else{sigmaAng = 0.0;}
	  
	      fprintf(outAng," %12.5E %10.3E", qAng*idsang, 2.0*sigmaAng*idsang);
	    }
             
	  fprintf(outAng,"\n");
	}
    
      //Let blank line between theta bins
      fprintf(outAng,"\n"); 
      dsang = (cos(binTlow)-cos(binTtop))*pmax;
      idsang = 1.0/dsang;
    
      for(unsigned int j = 0; j < constants::nParTypes; j++)
	{
	  double sigmaPolar = (q2Polar[j]-qPolar[j]*qPolar[j])*invn;
	  if(sigmaPolar > 0.0){ sigmaPolar = sqrt(sigmaPolar);}
	  else{sigmaPolar = 0.0;}
	  
	  fprintf(outPolar," %12.5E %10.3E", qPolar[j]*idsang,2.0*sigmaPolar*idsang);
	}
    
      fprintf(outPolar,"\n");
    }
    
  fclose(outAng);
  fclose(outPolar); 
}


int pen_EmergingPartDistrib::sumTally(const pen_EmergingPartDistrib& tally){

    if(nBinsE != tally.nBinsE)
        return -1;
    
    if(nAngBins != tally.nAngBins)
        return -2;
    
    
    
    //Upbound
    //********
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        for(int i = 0; i < nBinsE; ++i){
            eHistUp[k][i] += tally.eHistUp[k][i];
        }
    }
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        for(int i = 0; i < nBinsE; ++i){
            eHist2Up[k][i] += tally.eHist2Up[k][i];
        }
    }
    
    //Downbound
    //********
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        for(int i = 0; i < nBinsE; ++i){
            eHistDown[k][i] += tally.eHistDown[k][i];
        }
    }
    
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        for(int i = 0; i < nBinsE; ++i){
            eHist2Down[k][i] += tally.eHist2Down[k][i];
        }
    }    
    
    //Angular
    //*********
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        for(int i = 0; i < nAngBins; ++i){
            angHist[k][i] += tally.angHist[k][i];
        }
    }
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        for(int i = 0; i < nAngBins; ++i){
            angHist2[k][i] += tally.angHist2[k][i];
        }
    }
    
    return 0;
}


REGISTER_COMMON_TALLY(pen_EmergingPartDistrib, EMERGING_PART_DISTRIB)
 
        
    
    
    
    
    
    
    
    
    
      


      
      
      
      
      
      
      
      
      
      
      
