
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

 
#include "tallySphericalDoseDistrib.hh"

void pen_SphericalDoseDistrib::updateEdepCounters(const double dE,
						  const unsigned long long nhist,
						  const double X,
						  const double Y,
						  const double Z,
						  const double WGHT){
    int i;
    double r2;
    
    r2 = X*X + Y*Y + Z*Z;
    //Not in ROI
    if (r2 < r2min || r2 >= r2max){return;}
    //Bin index where score will be added
    i = (int)((sqrt(r2) - rmin)*idr);                        
    
    //Transfer partial to totals when a new history visits
    if (nhist > nlast[i])    
    {
            edep[i]  += edptmp[i];
            edep2[i] += edptmp[i]*edptmp[i];
            edptmp[i] = dE*WGHT;
            // Add 1/2 to avoid roundoff errors
            nlast[i]  = nhist;             
    }
    else
    {
        //Same history as before, keep adding to temp counter
        edptmp[i] += dE*WGHT;       
    }
}

 void pen_SphericalDoseDistrib::flush(){
       // Dump temp counters and obtain max score:
    for(int i = 0; i < nbin; ++i)
    {
      // Skip empty bins
      if (nlast[i] < 0.5){ continue;}      
      // Transfer temp counter
      edep[i]  += edptmp[i];   
      edep2[i] += edptmp[i]*edptmp[i];
      // Reset counter
      edptmp[i]= 0.0;                  
      // Reset last visited to avoid recounting in next report
      nlast[i] = 0.0;                  
    }
 }
 
 
 

void pen_SphericalDoseDistrib::tally_localEdep(const unsigned long long nhist,
					       const pen_KPAR /*kpar*/,
					       const pen_particleState& state,
					       const double dE){
    
   if(dE == 0.0){return;}  
   //Nothing to deposit.
   //Energy deposited at material
   updateEdepCounters(dE, nhist, state.X, state.Y, state.Z, state.WGHT);
}

void pen_SphericalDoseDistrib::tally_beginPart(const unsigned long long nhist,
					       const unsigned /*kdet*/,
					       const pen_KPAR /*kpar*/,
					       const pen_particleState& state){
 
    //Extract energy from material to create new particle
    updateEdepCounters(-state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
}

void pen_SphericalDoseDistrib::tally_beginHist(const unsigned long long nhist,
					       const unsigned /*kdet*/,
					       const pen_KPAR /*kpar*/,
					       const pen_particleState& state){

    if(state.MAT > 0){
        //Particle created at non void volume. Add particle energy to compensate
        //substracted one when beginPart will be called.
        updateEdepCounters(state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
    }
}

void pen_SphericalDoseDistrib::tally_step(const unsigned long long nhist,
					  const pen_KPAR /*kpar*/,
					  const pen_particleState& state,
					  const tally_StepData& stepData){

   if(stepData.softDE == 0.0){return;}  //Nothing to deposit.
   //Energy deposited
   updateEdepCounters(stepData.softDE, nhist,
		      stepData.softX, stepData.softY, stepData.softZ,
		      state.WGHT);  
}


void pen_SphericalDoseDistrib::tally_move2geo(const unsigned long long nhist,
					      const unsigned /*kdet*/,
					      const pen_KPAR /*kpar*/,
					      const pen_particleState& state,
					      const double /*dsef*/,
					      const double /*dstot*/){

    //Primary particle has been created at void volume. Check if
    //after step call particle reached non void volume.
    if(state.MAT > 0){
        //Non void volume reached. Add particle energy to compensate
        //substracted one when beginPart will be called.
        updateEdepCounters(state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
    }
}



int pen_SphericalDoseDistrib::configure(const wrapper_geometry& geometry,
					const abc_material* const materials[constants::MAXMAT],
					const pen_parserSection& config,
                    const unsigned verbose
){

    double rmax;
    // A bit more than one
    const double oneplus = 1.0+1.0E-12;  
    // A bit less than one
    const double oneminus= 1.0-1.0E-12;  
    int i, err; 
    
    err = config.read("rmin", rmin);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_SphericalDoseDistrib:configure:unable to read 'rmin' in configuration. Double expected\n");
        }
        return -1;
    }
    
    err = config.read("rmax", rmax);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_SphericalDoseDistrib:configure:unable to read 'rmax' in configuration. Double expected\n");
        }
        return -2;
    }
    err = config.read("nbin", nbin);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_SphericalDoseDistrib:configure:unable to read 'nbin' in configuration. Integrer expected\n");
        }
        return -3;
    }
    

    
    if(rmin < 0){
        if(verbose > 0){
            printf("pen_SphericalDoseDistrib:configure:Invalid entry, rmin must be greater than 0.\n");
        }
        return -4;
    }
    
     if(rmin >= rmax){
        if(verbose > 0){
            printf("pen_SphericalDoseDistrib:configure:Invalid entry, rmin must be lower than rmax.\n");
        }
        return -5;
    }
    
     if(nbin < 1){
        if(verbose > 0){
            printf("pen_SphericalDoseDistrib:configure:Invalid entry, nbin must be at least 1.\n");
        }
        return -6;
    }
    if(nbin > nbinmax){
        if(verbose > 0){
            printf("pen_SphericalDoseDistrib:configure:Too many bins. Max no.bins is %d\n", nbinmax);
        }
        return -7;
    }

    
    err = config.read("print-xyz", prtxyz);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_SphericalDoseDistrib:configure:unable to read 'print-xyz' in configuration. Bool expected\n");
        }
        return -8;
    }
    
    
    if(verbose > 1){
        if(prtxyz)
        {
            printf("pen_SphericalDoseDistrib:configure: print-xyz =  true. Output file print rBinIndex,  rLow(cm) and  rAve(cm)\n");
        }
        else{
            printf("pen_SphericalDoseDistrib:configure: print-xyz = false. Output file doesn't print rBinIndex,  rLow(cm) and  rAve(cm). It occurs if print-xyz = true  \n");
        }
    }
    
    dr = (rmax - rmin)/(double)nbin;
    idr = 1.0/dr;
    // (1+-eps) tolerances avoid out of range bin indexes in tally routine
    r2min = rmin*rmin*oneplus;  
    r2max = rmax*rmax*oneminus;
    
    
     // Init arrays:
    for(i = 0; i < nbinmax; i++)
    {
      edptmp[i] = 0.0;
      edep[i]   = 0.0;
      edep2[i]  = 0.0;
      nlast[i]  = 0.0;
      idens[i]  = 0.0;
    }
    
    pen_particleState state;
    
    // For locate() to operate properly
    state.U = 0.0;  
    state.V = 0.0;
    state.W = 1.0;
    state.X = 0.0;
    state.Y = 0.0;
    
    
    for(i = 0; i < nbin; i++)
    {
         //This is to locate a point and find its material
        state.Z = rmin + dr*(i+0.5); 
        geometry.locate(state);
        if(state.MAT > 0)
        {
            //Local mass density
            double localdens = materials[state.MAT-1]->readDens(); 
            if (localdens > 0.0){idens[i] = 1.0/localdens;}
        }
    }

   if(verbose > 1){
    printf("Number of radial bins: %d\n",nbin);
    printf("Minimum r value (cm): %12.5E\n",rmin);
    printf("Bin width (cm): %12.5E\n",dr);
  }  
    
  //Register data to dump
  dump.toDump(nlast,nbin);
  dump.toDump(idens,nbin);
  dump.toDump(edep,nbin);
  dump.toDump(edep2,nbin);
  dump.toDump(&nbin,1);  
    
    
    
  return 0;   
    
    
}

void pen_SphericalDoseDistrib::tally_endSim(const unsigned long long /*nhist*/){
        
  flush();
}

void pen_SphericalDoseDistrib::saveData(const unsigned long long nhist) const{
   
  FILE*out = 0;
  int i;
  double q,sigma,fact,r,rave,rnext;
  double invdr3,invn;
  const double threefour = 3.0/4.0;
  const double pi = constants::PI;
  const double inv43pi =3.0/(4.0*pi);

    
  //Prepare output files:
  out = fopen("sphericalDoseDistrib.dat", "w");
  if (out == 0){
    printf("\n*********************************************\n");
    printf("SphericalDoseDistrib:saveData:ERROR: cannot open output data file");
    printf("\n*********************************************\n");
    return;
  }
  
  // Write header:
  fprintf(out,"#------------------------------------------------------------\n");
  fprintf(out,"# PenRed: Spherical dose distribution\n");
  fprintf(out,"# Dose units are eV/g per history\n");
  fprintf(out,"#\n");
  fprintf(out,"# No. of radial (r) bins:\n");
  fprintf(out,"#  %d\n",nbin);
  fprintf(out,"# Min r and bin width (cm):\n");
  fprintf(out,"#  %12.5E %12.5E\n", rmin, dr);
  if (prtxyz == 1){
    fprintf(out,"#\n");
    fprintf(out,"# Coordinate print enabled. Two radial values for each bin will be printed, the low and an average radius. the average radial value is a weighted proportionally to r^2.\n");
  }
    
  fprintf(out,"#\n");
  fprintf(out,"# ");
  if (prtxyz == 1){
    fprintf(out,"rBinIndex :  rLow(cm)  :  rAve(cm)  : ");
  }
  fprintf(out,"dose(eV/g) :+-2sigma: edep(eV/cm) :+-2sigma\n");
                
  // Write data
  invn = 1.0/static_cast<double>(nhist);  
    


  for(i = 0; i < nbin; i++)
    {
      r = rmin+dr*(i);
      rnext = r+dr;
      invdr3 = 1.0/(rnext*rnext*rnext-r*r*r);
      // This is 1/Delta_mass
      fact = idens[i]*inv43pi*invdr3;             
      if (prtxyz == 1){
	// Average with weight(r)~r^2
	rave = threefour*(pow(rnext,4)-pow(r,4))*invdr3; 
	fprintf(out," %6d      %12.5E %12.5E",i,r,rave);
      }
      q = edep[i]*invn;
      sigma = (edep2[i]*invn-(q*q))*invn;
      if(sigma > 0.0)
        {
	  sigma = sqrt(sigma);
        }
      else
        {
	  sigma = 0.0;
        }
        
      double sigma2 = 2.0*sigma;
      fprintf(out," %12.5E %8.1E %12.5E  %8.1E\n",q*fact,sigma2*fact,q*idr,sigma2*idr);
      

    }
        
  fclose(out); 
}

int pen_SphericalDoseDistrib::sumTally(const pen_SphericalDoseDistrib& tally){
    
    if(nbin != tally.nbin)
        return -1;
    
    for(int i = 0; i < nbin; ++i){
        edep[i] += tally.edep[i];
    }
    for(int i = 0; i < nbin; ++i){
        edep2[i] += tally.edep2[i];
    }

  return 0;
    
    
}


REGISTER_COMMON_TALLY(pen_SphericalDoseDistrib, SPHERICAL_DOSE_DISTRIB)
