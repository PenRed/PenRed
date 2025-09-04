
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//

 
#include "tallyAngularDetector.hh"



void pen_AngularDet::scapedParticle(const pen_KPAR kpar,
                                    const pen_particleState& state){
    //Calculate the energy bin
    const double TWOPI = 2.0*constants::PI;
    int particleEbin;
    
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
        phi += TWOPI;
    }

    
    if(phi1 < 0.0 && phi > constants::PI)
    {
        phi -= TWOPI;
    }
    
    if(phi >= phi1 && phi <= phi2 && theta >= theta1 && theta <= theta2)
    {
        if(state.E >= emin && state.E < emax)
        {
            double auxE;
              if(!isLinScale)
              {
                auxE = log(state.E) - emin;
              }
              else
              {
                auxE = state.E - emin;
              }
            particleEbin = auxE*iebin;              
    
        
            //Chek if the energy of the particle is in the energy range
            if(particleEbin >= 0 && particleEbin < nBinsE)
            {
                angDetSpcTmp[kpar][particleEbin] += state.WGHT;
            }
        }
    }
}

void pen_AngularDet::flush(){
    
    for(unsigned int k = 0; k < constants::nParTypes; k++)
    {
        for(int i = 0; i < nBinsE; i++)
        {
            //Spectrum energy angular detector counters
            if(angDetSpcTmp[k][i] == 0.0){continue;} // Skip void counters
            angDetSpc[k][i]  += angDetSpcTmp[k][i];
            angDetSpc2[k][i] += angDetSpcTmp[k][i]*angDetSpcTmp[k][i];
            angDetSpcTmp[k][i] = 0.0;
        }
    }
}

void pen_AngularDet::tally_move2geo(const unsigned long long /*nhist*/,
					    const unsigned /*kdet*/,
					    const pen_KPAR kpar,
					    const pen_particleState& state,
					    const double /*dsef*/,
					    const double /*dstot*/){
     
  if(state.MAT == 0){   
    //Particle scape the geometry
    //The particle scaped from the material system    
    scapedParticle(kpar,state);
  }
    
}

void pen_AngularDet::tally_matChange(const unsigned long long /*nhist*/,
                         const pen_KPAR kpar,
					     const pen_particleState& state,
					     const unsigned /*prevMat*/){

        if(state.MAT == 0){
            //Particle scape the geometry
            //The particle scaped from the material system
            scapedParticle(kpar,state);
        }
}


void pen_AngularDet::tally_endHist(const unsigned long long /*nhist*/){
    flush();
}


int pen_AngularDet::configure(const wrapper_geometry& /*geometry*/,
				      const abc_material* const /*materials*/[constants::MAXMAT],
				      const pen_parserSection& config,
				      const unsigned verbose){
    
    int err;
    
    
    //Minimum energy
    //**************
    err = config.read("emin", emin);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
        printf("pen_AngularDet:configure: Error: unable to read 'emin' in configuration. Double expected\n");
        }
        return -1;
    }
    
    
    //Maximum energy
    //**************
    err = config.read("emax", emax);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
        printf("pen_AngularDet:configure: Error: unable to read 'emax' in configuration. Double expected\n");
        }
        return -2;
    }
    
    
    if(emin >= emax)
    {
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: minimum energy (%.5E eV) can't be greater than maximum energy (%.5E eV).\n", emin, emax);
        }
        return -3;
    }
    
    if(verbose > 1){
        printf("Energy interval:\n");
        printf(" %.5E - %.5E\n",emin,emax);
    }
    
    
    
    //Number of bins
    //**************
    err = config.read("nBinsE", nBinsE);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: unable to read 'nBinsE' in configuration. Integrer expected\n");
        }
        return -4;
  }
    
    
    if(nBinsE < 1)
    {
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: The number of energy bins can't be 0.\n");
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
            printf("pen_AngularDet:configure: Error: The maximum number of energy bins is %d.\n",nbinmax);
        }
        return -6;
    }
    
    
    
    //Minimum angles
    //**************
    
    const double degree2rad = constants::PI/180.0;
    
    err = config.read("theta1", theta1);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: unable to read 'theta1' in configuration. Double expected\n");
        }
        return -7;
    }
    
    
    err = config.read("phi1", phi1);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: unable to read 'phi1' in configuration. Double expected\n");
        }
        return -8;
    }
    
    
    
    //Maximum angles
    //**************
    err = config.read("theta2", theta2);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: unable to read 'theta2' in configuration. Double expected\n");
        }
        return -9;
    }
    
    
    err = config.read("phi2", phi2);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: unable to read 'phi2' in configuration. Double expected\n");
        }
        return -10;
    }
    

    
    
    //Theta restrictions
    //*******************
    if(theta1 >= theta2)
    {
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: Invalid values for theta1 = %11.4E and theta2 = %11.4E. Theta2 must be greater than theta1.\n", theta1, theta2);
        }
        return -11;
    }
    
    
    //Theta interval must be (0,180)
    if(theta1 < 0.0 || theta1 > 180.0)
    {
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: Invalid value for theta1 = %11.4E. Theta1 must be in the interval (0,180).\n", theta1);
        }
        return -12;
    }
    
    
    
    if(theta2 < 0.0 || theta2 > 180.0)
    {
        if(verbose > 0){
            printf("pen_AngularDet:configure: Error: Invalid value for theta2 = %11.4E. Theta2 must be in the interval (0,180).\n", theta2);
        }
        return -13;
    }
    
    
    
    //Phi restrictions
    //*******************    
    if(phi1 >= phi2)
    {
         if(verbose > 0){
            printf("pen_AngularDet:configure: Error: Invalid values for phi1 = %11.4E and phi2 = %11.4E. Phi2 must be greater than phi1.\n", phi1, phi2);
        }
        return -14;
    }
    
    
    //Phi interval must be (0,360) or (-180,180)
    if(phi1 < 0.0)
    {
        if(phi1 < -180.0)
        {
            if(verbose > 0){
            printf("pen_AngularDet:configure: Error: Invalid value for phi1 = %11.4E. Phi1 must be larger than -180.\n", phi1);
            }
            return -15;
        }
        if(phi2 > 180.0)
        {
            if(verbose > 0){
            printf("pen_AngularDet:configure: Error: Invalid value for phi2 = %11.4E. Phi2 must be lower than 180.\n", phi2);
            }   
            return -16;
        }
    }
    else
    {
        if(phi1 > 360.0)
        {
            if(verbose > 0){
            printf("pen_AngularDet:configure: Error: Invalid value for phi1 = %11.4E. Phi1 must be lower than 360.\n", phi1);
            }
            return -17;
        }
        if(phi2 > 360.0)
        {
            if(verbose > 0){
            printf("pen_AngularDet:configure: Error: Invalid value for phi2 = %11.4E. Phi2 must be lower than 360.\n", phi2);
            }
            return -18;
        }
        
    }
    
    //Convert angles in radians to operate properly
    theta1 = theta1*degree2rad;
    theta2 = theta2*degree2rad;
    phi1 = phi1*degree2rad;
    phi2 = phi2*degree2rad;
    
    
    //Scale type
    //**************
    err = config.read("linearScale", isLinScale);
    if(err != INTDATA_SUCCESS){
        if(verbose > 1){
            printf("No 'linearScale' specified, linear scale will be set.\n\n");
        }
        isLinScale = true;
    }
    
    if(isLinScale){
        if(verbose > 1){
            printf("Linear scale\n");
        }
    }
    else{
        if(verbose > 1){
            printf("Logarithmic scale\n");
        }
        //Default values due to the Logarithmic scale
        if(emin == 0.0){
            emin = 50.0;
        }
        
        //Store values in log scale
        emin = log(emin);
        emax = log(emax);
    }
    
    ebin = (emax - emin)/(double)nBinsE;
    iebin = 1.0/ebin;
    
    
    
  // Detector
  //***************
  int auxkDet;
  err = config.read("detector", auxkDet);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_AngularDet:configure: Error: Unable to read 'detector' in configuration. Integrer expected\n");
    }
    return -19;
  }

  if(auxkDet <= 0){
    if(verbose > 0){
      printf("pen_AngularDet:configure: Error: Index of detector must be greater than 0.\n");
    }
    return -20;
  }
    
  idet = (unsigned)auxkDet;
  if(verbose > 1){
    printf("Detector:\n");
    printf(" %u\n",idet);
  }
    
    
    
    
    
    
    //Init arrays
    
    for(unsigned int i = 0; i < constants::nParTypes; i++)
    {
        for(int j = 0; j < nbinmax; j++)
        {
            angDetSpc[i][j]    = 0.0;
            angDetSpcTmp[i][j] = 0.0;
            angDetSpc2[i][j]   = 0.0;
        }
    }
    
    //Register data to dump
    for(unsigned k = 0; k < constants::nParTypes; k++)
    {
        dump.toDump(angDetSpc[k],nBinsE);
        dump.toDump(angDetSpc2[k],nBinsE);
        dump.toDump(&nBinsE,1); 
    }
    
    
    return 0;
}


void pen_AngularDet::saveData(const unsigned long long nhist) const{

  FILE* out = NULL;

  const double rad2degree = 180.0/constants::PI;
  
  char filename[400];
  
  //Prepare output file
  double invn = 1.0/static_cast<double>(nhist);
  
  sprintf(filename, "spc-angdet-%u.dat",idet);
  out = fopen(filename, "w");
  
  if(out == NULL)
  {
      printf(" *********************************************\n");
      printf(" AngularDet:saveData:ERROR: cannot open output data file\n");
      printf(" *********************************************\n");
      return;
  }
  
  
  fprintf(out, "#-------------------------------------------------------\n");
  fprintf(out, "# PenRed: Angular detector energy spectrum\n");
  
  //Angular intervals
  fprintf(out,"# Angular intervals : theta1 = %13.6E,    theta2 = %13.6E\n", theta1*rad2degree, theta2*rad2degree);
  fprintf(out,"#                        phi1 = %13.6E,    phi2 = %13.6E\n", phi1*rad2degree, phi2*rad2degree);
  fprintf(out,"# Energy window (eV), energy bins:\n");
  if(isLinScale)
  {
      fprintf(out,"# %12.5E %12.5E %5d\n", emin, emax, nBinsE); 
  }
  else
  {
      fprintf(out,"# %12.5E %12.5E %5d\n", exp(emin), exp(emax), nBinsE);
  }
  
  fprintf(out,"# Energy spectra of emerging particles (1/(eV*sr*particle)).\n");
  fprintf(out,"#\n");
    
   fprintf(out,"#             |");
   for(unsigned ipart = 0; ipart < constants::nParTypes; ipart++){
      fprintf(out, "%5.5s%15.15s%5.5s|"," ", particleName(ipart), " ");
   }
   fprintf(out,"\n");
    
   fprintf(out,"#             |");
   for(unsigned ipart = 0; ipart < constants::nParTypes; ipart++){
   fprintf(out,"%5.5s%15.15s%5.5s|"," "," ", " ");
   }
   fprintf(out,"\n");
   

   fprintf(out, "# Energy (eV) |");
   for(unsigned ipart = 0; ipart < constants::nParTypes; ipart++){
        fprintf(out, "     prob    : +-2sigma  |"); 
    }
   fprintf(out, "\n");
    
   
   double solidAngle = (phi2 - phi1)*(cos(theta1) - cos(theta2));
   
    double deBin = ebin;
   for(int j = 0; j < nBinsE; j++)
   {
       	double energy;
        if(isLinScale){
            energy = emin + ((double)j+0.5)*ebin;
        }
        else{
            double xLow = exp(emin +(double)j*ebin);
            double xUp = exp(emin + (double)(j+1)*ebin);
            energy = 0.5*(xUp + xLow);
            deBin = (xUp - xLow);
        }
        
   
        fprintf(out, " %12.5E ",energy);
        
        double factor = 1.0/(deBin*solidAngle);
   
        //Particle spectrum
        double q, q2, sigma;
        for(unsigned ipart = 0; ipart < constants::nParTypes; ipart++)
        {
            q = angDetSpc[ipart][j]*invn;
            q2 = angDetSpc2[ipart][j]*invn;
            sigma = (q2 - q*q)*invn;
            if(sigma > 0.0){sigma = sqrt(sigma);}
            else{sigma = 0.0;}
            
            fprintf(out," %12.5E   %8.3E ", q*factor, 2.0*sigma*factor);
        }
        fprintf(out, "\n");
   }
   fclose(out);
}


int pen_AngularDet::sumTally(const pen_AngularDet& tally){
    
    
    if(nBinsE != tally.nBinsE)
        return -1;    
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        for(int i = 0; i < nBinsE; ++i){
            angDetSpc[k][i] += tally.angDetSpc[k][i];
        }
    }
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        for(int i = 0; i < nBinsE; ++i){
            angDetSpc2[k][i] += tally.angDetSpc2[k][i];
        }
    }

  return 0;
    
    
}

penred::measurements::results<double, 1> pen_AngularDet::generateResults(const pen_KPAR kpar,
									 const unsigned long long nhists){

  if(!isLinScale){
    penred::measurements::results<double, 1> results;
    results.description = "Error: Results can only be generated for linear scale";
    return results;
  }
  
  const double rad2deg = 180.0/constants::PI;
  const double solidAngle = (phi2 - phi1)*(cos(theta1) - cos(theta2));
  const double deBin = ebin;
  const double factor = 1.0/(deBin*solidAngle);
  double invn = 1.0/static_cast<double>(nhists);
  const unsigned ebins = static_cast<unsigned>(nBinsE);
	  
  std::string description("PenRed: Angular detector energy spectrum for ");
  description += particleName(kpar);
  description += "\n";
  description += "  - Angular intervals :  theta1 = " + std::to_string(theta1*rad2deg);
  description += ",  theta2 = " + std::to_string(theta2*rad2deg) + "\n";
  description += "                         phi1 = " + std::to_string(phi1*rad2deg);
  description += ",  phi2 = " + std::to_string(phi2*rad2deg) + "\n";
  description += " Histories simulated: " + std::to_string(nhists) + "\n";
  description += " Energy spectra of emerging particles (1/(eV*sr*particle)).\n";
  
  //Create results
  penred::measurements::results<double, 1> results;
  results.initFromLists({ebins},
			{penred::measurements::limitsType(emin, emax)});
	  
  results.description = description;
  
  results.setDimHeader(0, "Energy (eV)");
  results.setValueHeader("Particles (1/hist)");

  for(unsigned j = 0; j < ebins; ++j)
    {
      double q = angDetSpc[kpar][j]*invn;
      double q2 = angDetSpc2[kpar][j]*invn;
      double sigma = (q2 - q*q)*invn;
      if(sigma > 0.0){sigma = sqrt(sigma);}
      else{sigma = 0.0;}

      results.data[j] = q*factor;
      results.sigma[j] = sigma*factor;
    }
	  
  return results;    
}


REGISTER_COMMON_TALLY(pen_AngularDet)






























