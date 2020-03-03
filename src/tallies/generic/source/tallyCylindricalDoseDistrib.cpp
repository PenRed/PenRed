
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

 
#include "tallyCylindricalDoseDistrib.hh"

void pen_CylindricalDoseDistrib::updateEdepCounters(const double dE,
						    const double nhist,
						    const double X,
						    const double Y,
						    const double Z,
						    const double WGHT){

  //Avoid count range (zmin-dz,zmin)
  if(Z < zmin) return;
  
    int i,k;
    double r2;
    int bin = 0;
    
    k = (Z - zmin)*idz;
    if(k < 0 || k >= nz){return;}

    r2 = X*X + Y*Y;
    if (r2 < r2min || r2 >= r2max){return;}

    // Bin index where score will be added
    i = (int)((sqrt(r2) - rmin)*idr);
    // Maps i,k into a single index
    bin = i + k*nr;
    
    //Transfer partial to totals when a new history visits
    if (nhist > nlast[bin])
    {
            edep[bin]  += edptmp[bin];
            edep2[bin] += edptmp[bin]*edptmp[bin];
            edptmp[bin] = dE*WGHT;
            // Add 1/2 to avoid roundoff errors
            nlast[bin]  = nhist+0.5;
    }
    else
    {
        //Same history as before, keep adding to temp counter
        edptmp[bin] += dE*WGHT;
    }
}

 void pen_CylindricalDoseDistrib::flush(){
    // Dump temp counters and obtain max score:
    for(int k = 0; k < nz; k++)
    {
        for(int i = 0; i < nr; i++)
        {
            //Maps i,k into a single index
            int bin = i + k*nr;  
            // Skip empty bins
            if (nlast[bin] < 0.5){ continue;}
            // Transfer temp counter
            edep[bin]  += edptmp[bin];   
            edep2[bin] += edptmp[bin]*edptmp[bin];
            // Reset counter
            edptmp[bin]= 0.0;                  
            // Reset last visited to avoid recounting in next report
            nlast[bin] = 0.0;                  
        }
    }
 }
 
 
 

void pen_CylindricalDoseDistrib::tally_localEdep(const double nhist,
						 const pen_KPAR /*kpar*/,
						 const pen_particleState& state,
						 const double dE){
    
  //Nothing to deposit.
  if(dE == 0.0){return;}  
  //Energy deposited at material
  updateEdepCounters(dE, nhist, state.X, state.Y, state.Z, state.WGHT);
}

void pen_CylindricalDoseDistrib::tally_beginPart(const double nhist,
						 const unsigned /*kdet*/,
						 const pen_KPAR /*kpar*/,
						 const pen_particleState& state){
 
    //Extract energy from material to create new particle
    updateEdepCounters(-state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
}

void pen_CylindricalDoseDistrib::tally_beginHist(const double nhist,
						 const unsigned /*kdet*/,
						 const pen_KPAR /*kpar*/,
						 const pen_particleState& state){

    if(state.MAT > 0){
        //Particle created at non void volume. Add particle energy to compensate
        //substracted one when beginPart will be called.
        updateEdepCounters(state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
    }
}

void pen_CylindricalDoseDistrib::tally_step(const double nhist,
					    const pen_KPAR /*kpar*/,
					    const pen_particleState& state,
					    const tally_StepData& stepData){

  if(stepData.softDE > 0.0){  
    //Energy deposited
    updateEdepCounters(stepData.softDE, nhist,
		       stepData.softX, stepData.softY, stepData.softZ,
		       state.WGHT);
  }
}


void pen_CylindricalDoseDistrib::tally_move2geo(const double nhist,
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



int pen_CylindricalDoseDistrib::configure(const wrapper_geometry& geometry,
					const abc_material* const materials[constants::MAXMAT],
					const pen_parserSection& config,
                    const unsigned verbose
){
    int k;
    double rmax, zmax;
    // A bit more than one
    const double oneplus = 1.0+1.0E-12;  
    // A bit less than one
    const double oneminus= 1.0-1.0E-12;  
    int i, err; 
    
    err = config.read("rmin", rmin);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'rmin' in configuration. Double expected\n");
        }
        return -1;
    }
    
    err = config.read("rmax", rmax);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'rmax' in configuration. Double expected\n");
        }
        return -2;
    }
    err = config.read("nbinsr", nr);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'nr', number of r bins, in configuration. Integrer expected\n");
        }
        return -3;
    }
    

    
    if(rmin < 0){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:Invalid entry, rmin must be greater than 0.\n");
        }
        return -4;
    }
    
     if(rmin >= rmax){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:Invalid entry, rmin must be lower than rmax.\n");
        }
        return -5;
    }
    
     if(nr < 1){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:Invalid entry, nr must be at least 1.\n");
        }
        return -6;
    }
    
    
    err = config.read("print-xyz", prtxyz);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'print-xyz' in configuration. Bool expected\n");
        }
        return -7;
    }
    
    
    if(verbose > 1){
        if(prtxyz)
        {
            printf("pen_CylindricalDoseDistrib:configure: print-xyz =  true. Output file print rBinIndex,  rLow(cm) and  rAve(cm)\n");
        }
        else{
            printf("pen_CylindricalDoseDistrib:configure: print-xyz = false. Output file doesn't print rBinIndex,  rLow(cm) and  rAve(cm). It occurs if print-xyz = true  \n");
        }
    }
    
    
    dr = (rmax - rmin)/nr;
    idr = 1.0/dr;
    // (1+-eps) tolerances avoid out of range bin indexes in tally routine
    r2min = rmin*rmin*oneplus;  
    r2max = rmax*rmax*oneminus;
    
        
    err = config.read("zmin", zmin);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'zmin' in configuration. Double expected\n");
        }
        return -8;
    }
    
    err = config.read("zmax", zmax);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'zmax' in configuration. Double expected\n");
        }
        return -9;
    }
    err = config.read("nbinsz", nz);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'nz', number of z bins, in configuration. Integrer expected\n");
        }
        return -10;
    }
    

    
    if(zmin > zmax){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:Invalid entry, zmin must be lower than zmax.\n");
        }
        return -11;
    }
    
    if(nz < 0){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:Invalid entry, nz must be at least 0.\n");
        }
        return -12;
    }
    
    if(nz > 0 && zmin == zmax){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:Invalid entry, when nz > 0, zmin must be lower than zmax.\n");
        }
        return -13;
    }
    
    if (nz == 0)
    {
        dz = zmax - zmin;
        idz = 0.0;           
        nz = 1;
    }
    else
    {
        dz = (zmax - zmin)/nz;
        idz = 1.0/dz;
    }
    
    rnbin = nr*nz; 
   
    if(rnbin > nbinmax)
    {
        if(verbose > 0){
        printf("pen_CylindricalDoseDistrib:configure:Too many bins. Max no.bins is %d\n", nbinmax);
        }
        return -14;
    }
    
 
    
     // Init arrays:
    for(int unsigned j = 0; j < nbinmax; j++)
    {
      edptmp[j] = 0.0;
      edep[j]   = 0.0;
      edep2[j]  = 0.0;
      nlast[j]  = 0.0;
      idens[j]  = 0.0;
    }
    
    pen_particleState state;
    
    // For locate() to operate properly
    state.U = 0.0;  
    state.V = 0.0;
    state.W = 1.0;
    state.X = 0.0;
    state.Y = 0.0;
    
    for(k = 0; k < nz; k++)
    {
        state.Z = zmin + dz*(k + 0.5);
        for(i = 0; i < nr; i++)
        {
            //This is to locate a point and find its material
            state.X = rmin + dr*(i+0.5);  
            int bin = i + k*nr;
            geometry.locate(state);
            if(state.MAT > 0)
            {
                //Local mass density
                double localdens = materials[state.MAT-1]->readDens(); 
                if (localdens > 0.0){idens[bin] = 1.0/localdens;}
            }
        }
    }

   if(verbose > 1){
    printf("Number of radial bins: %d\n",nr);
    printf("Minimum r value (cm): %12.5E\n",rmin);
    printf("Bin width radial (cm): %12.5E\n",dr);
    printf("Number of depth bins: %d\n",nz);
    printf("Minimum z value (cm): %12.5E\n",zmin);
    printf("Bin width depth (cm): %12.5E\n",dz);
  }  
    
  //Register data to dump
  dump.toDump(idens,rnbin);
  dump.toDump(edep,rnbin);
  dump.toDump(edep2,rnbin);
  dump.toDump(&rnbin,1);    
    
  return 0;   
    
    
}

void pen_CylindricalDoseDistrib::tally_endSim(const double /*nhist*/){
        
  flush();
}

void pen_CylindricalDoseDistrib::saveData(const double nhist) const{
   
  char buffer[81];
  FILE*out = 0;
  int i,k,nzef;
  double q,sigma,fact,r,z,rave,zmiddle;
  double invdr2,invn,factidz;
  const double twothird=2.0/3.0;
  const double invpi=1.0/constants::PI;

  int bin = 0;
  //Prepare z factors and filename
  if(idz == 0.0)   // No z bins
    {
      nzef = 0;
      factidz = 1.0;
      strcpy(buffer,"radialDoseDistrib.dat");
    }
  else
    {
      nzef = nz;
      factidz = idz;
      strcpy(buffer,"cylindricalDoseDistrib.dat");
    }
     
  //Prepare output file:
  out = fopen(buffer, "w");
  if (out == NULL){
    printf("\n*********************************************\n");
    printf("CylindricalDoseDistrib:saveData:ERROR: cannot open output data file");
    printf("\n*********************************************\n");
    fclose(out); //Just in case
    return;
  }
  
  // Write header:
  fprintf(out,"#-------------------------------------------------------\n");
  fprintf(out,"# PenRed: %s dose distribution report\n",idz == 0 ? "Radial" : "Cylindrical");
  if(nzef != 0)
    { 
      fprintf(out,"# Dose units are eV/g per history\n");
    }
  else
    {
      fprintf(out,"# Dose units are eV.cm/g per history\n");
    }
  fprintf(out,"#\n");
  fprintf(out,"# Radial (r) info -> Number of bins, minimum value, bin width (cm):\n");
  fprintf(out,"# %4d %12.5E %12.5E\n",nr,rmin,dr);
  fprintf(out,"# Depth (z) info -> Number of bins, minimum value, bin width (cm):\n");
  fprintf(out,"# %4d %12.5E %12.5E\n",nz,zmin,dz);
    
  if (prtxyz == 1){
    fprintf(out,"#\n");
    fprintf(out,"# Enabled coordinate print. The low and average value will be printed for each coordinate of each bin.\n");
    fprintf(out,"# For the z coordinate the average is the middle point of the bin.\n");
    fprintf(out,"# For the r coordinate the average is weighted  with a weight proportional to the radius r.\n");
  }
    
  fprintf(out,"#\n");
  fprintf(out,"# ");
  if (prtxyz == 1){
    fprintf(out,"rBinIndex :  rLow(cm)  :  rAve(cm)  : ");
    if (nzef != 0)
      {
	fprintf(out,"zBinIndex : zLow(cm) : zMiddle(cm) : ");
      }
  }
  fprintf(out,"   dose   : +-2sigma\n");
                
  // Write data:    
  invn = 1.0/nhist;  
    
  for(k = 0; k < nz; k++)
    {
      z = zmin + dz*k;
      zmiddle = z + dz*0.5;
      // Since z is not written, give at least a summary
      if(nzef != 0 && prtxyz != 1)    
        {
	  fprintf(out,"# zBinIndex=%-3d zMiddle(cm)=%12.5E\n",k,zmiddle);
        }
      for( i = 0; i < nr; i++)
        {
	  r = rmin + dr*i;
	  // Note: dr2 = (r+dr)^2-r^2
	  invdr2 = 1.0/(dr*(dr+2.0*r));   
	  if (prtxyz == 1)
            {
	      // Average with weight(r)~r
	      rave = twothird*((r+dr)*(r+dr)*(r+dr)-r*r*r)*invdr2; 
	      fprintf(out,"   %5d    %12.5E %12.5E",i,r,rave);
            }
	  if (nzef != 0 && prtxyz == 1)
            {
	      fprintf(out,"   %5d    %12.5E %12.5E",k,z,zmiddle);
            }
	  bin = i+k*nr;
	  fact = factidz*idens[bin]*invpi*invdr2;
	  q = edep[bin]*invn;
	  sigma = sqrt((edep2[bin]*invn-q*q)*invn > 0.0 ? (edep2[bin]*invn-q*q)*invn : 0.0)*fact;
	  q = q*fact;
	  fprintf(out," %12.5E  %8.1E\n",q,2.0*sigma);

        }

      // Separate data blocks
      if (nr > 1){ fprintf(out,"  \n");}                  

    }
        
  fclose(out); 
}

int pen_CylindricalDoseDistrib::sumTally(const pen_CylindricalDoseDistrib& tally){
    
    
    if(rnbin != tally.rnbin)
        return -1;
    
    for(int i = 0; i < rnbin; ++i){
        edep[i] += tally.edep[i];
    }
    for(int i = 0; i < rnbin; ++i){
        edep2[i] += tally.edep2[i];
    }

  return 0;
    
    
}


REGISTER_COMMON_TALLY(pen_CylindricalDoseDistrib, CYLINDRICAL_DOSE_DISTRIB)
