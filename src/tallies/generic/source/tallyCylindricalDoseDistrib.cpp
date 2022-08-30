
//
//
//    Copyright (C) 2019-2022 Universitat de València - UV
//    Copyright (C) 2019-2022 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

 
#include "tallyCylindricalDoseDistrib.hh"

void pen_CylindricalDoseDistrib::updateEdepCounters(const double dE,
						    const unsigned long long nhist,
						    const double X,
						    const double Y,
						    const double Z,
						    const double WGHT){

  //Avoid count range (zmin-dz,zmin)
  if(Z < zmin) return;
  
    long int i,k;
    double r2;
    
    k = (Z - zmin)*idz;
    if(k < 0 || k >= nz){return;}

    r2 = X*X + Y*Y;
    if (r2 < r2min || r2 >= r2max){return;}
    
    // Bin index where score will be added
    i = static_cast<long int>((sqrt(r2) - rmin)*idr);

    long int phiBin = 0;
    if(nphi > 1){
      if(fabs(X) > 1.0e-10 || fabs(Y) > 1.0e-10){
	double phi = atan2(Y,X);
	if(std::signbit(phi))
	  phi += 2.0*M_PI;
	phiBin = phi*idphi;
	if(phiBin >= nphi) --phiBin;
      }
    }
    
    // Maps i,k into a single index
    long int bin = nr*(k*nphi + phiBin) + i;
    
    //Transfer partial to totals when a new history visits
    if (nhist > nlast[bin]){
      edep[bin]  += edptmp[bin];
      edep2[bin] += edptmp[bin]*edptmp[bin];
      edptmp[bin] = dE*WGHT;
      // Add 1/2 to avoid roundoff errors
      nlast[bin]  = nhist;
    }
    else{
      //Same history as before, keep adding to temp counter
      edptmp[bin] += dE*WGHT;
    }
}

 void pen_CylindricalDoseDistrib::flush(){
    // Dump temp counters and obtain max score:
   for(long int bin = 0; bin < nbins; ++bin)
     {
       // Skip empty bins
       if (nlast[bin] == 0){ continue;}
       // Transfer temp counter
       edep[bin]  += edptmp[bin];   
       edep2[bin] += edptmp[bin]*edptmp[bin];
       // Reset counter
       edptmp[bin]= 0.0;                  
       // Reset last visited to avoid recounting in next report
       nlast[bin] = 0;                  
     }    
 }
 
 
 

void pen_CylindricalDoseDistrib::tally_localEdep(const unsigned long long nhist,
						 const pen_KPAR /*kpar*/,
						 const pen_particleState& state,
						 const double dE){
    
  //Nothing to deposit.
  if(dE == 0.0){return;}  
  //Energy deposited at material
  updateEdepCounters(dE, nhist, state.X, state.Y, state.Z, state.WGHT);
}

void pen_CylindricalDoseDistrib::tally_beginPart(const unsigned long long nhist,
						 const unsigned /*kdet*/,
						 const pen_KPAR /*kpar*/,
						 const pen_particleState& state){
 
    //Extract energy from material to create new particle
    updateEdepCounters(-state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
}

void pen_CylindricalDoseDistrib::tally_sampledPart(const unsigned long long nhist,
						   const unsigned long long /*dhist*/,
						   const unsigned /*kdet*/,
						   const pen_KPAR /*kpar*/,
						   const pen_particleState& state){

    if(state.MAT > 0){
        //Particle created at non void volume. Add particle energy to compensate
        //substracted one when beginPart will be called.
        updateEdepCounters(state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
    }
}

void pen_CylindricalDoseDistrib::tally_step(const unsigned long long nhist,
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


void pen_CylindricalDoseDistrib::tally_move2geo(const unsigned long long nhist,
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
    double rmax, zmax;
    // A bit more than one
    const double oneplus = 1.0+1.0E-12;  
    // A bit less than one
    const double oneminus= 1.0-1.0E-12;  
    int err; 
    
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
    int auxInt;
    err = config.read("nbinsr", auxInt);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'nr', number of r bins, in configuration. Integrer expected\n");
        }
        return -3;
    }
    nr = auxInt;
    

    
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
    err = config.read("nbinsz", auxInt);
    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_CylindricalDoseDistrib:configure:unable to read 'nz', number of z bins, in configuration. Integrer expected\n");
        }
        return -10;
    }
    nz = auxInt; 

    
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

    err = config.read("nbinsPhi", auxInt);
    if(err != INTDATA_SUCCESS){
      nphi = 1;
    }
    else{
      if(auxInt <= 0 || auxInt >= nbinmax){
        if(verbose > 0){
	  printf("pen_CylindricalDoseDistrib:configure: Number of phi bins "
		 "must be greater than zero and lesser than %ld\n"
		 "    Specified phi bins: %d\n", nbinmax,auxInt);
        }
        return -14;
      }
      nphi = auxInt;
    }
    dphi = 2.0*M_PI/static_cast<double>(nphi);
    idphi = 1.0/dphi;
    
    nbins = nr*nphi*nz; 
   
    if(nbins > nbinmax)
    {
        if(verbose > 0){
        printf("pen_CylindricalDoseDistrib:configure:Too many bins. "
	       "Max no.bins is %ld\n", nbinmax);
        }
        return -14;
    }
    
 
    
     // Init arrays:
    for(unsigned long j = 0; j < nbinmax; j++)
    {
      edptmp[j] = 0.0;
      edep[j]   = 0.0;
      edep2[j]  = 0.0;
      nlast[j]  = 0;
      imass[j]  = 0.0;
    }
    
    pen_particleState state;
    
    // For locate() to operate properly
    state.U = 0.0;  
    state.V = 0.0;
    state.W = 1.0;
    state.X = 0.0;
    state.Y = 0.0;

    double dphiz;
    if(idz == 0.0)
      dphiz = dphi;
    else
      dphiz = dphi*dz;
    for(long int k = 0; k < nz; ++k)
    {
      state.Z = zmin + dz*(static_cast<double>(k) + 0.5);
      for(long int j = 0; j < nphi; ++j){
	double phi = (static_cast<double>(j)+0.5)*dphi;
	double sphi = sin(phi);
	double cphi = cos(phi);
        for(long int i = 0; i < nr; ++i)
	  {
	    double rho1 = rmin + dr*static_cast<double>(i);
	    double rho  = rho1 + dr*0.5;
            //This is to locate a point and find its material
            state.X = rho*cphi;
            state.Y = rho*sphi;
            long int bin = nr*(k*nphi + j) + i;
            geometry.locate(state);
            if(state.MAT > 0)
	      {
                //Local mass density
                double localdens = materials[state.MAT-1]->readDens();
		double rho2 = rho1+dr;
		double volume = 0.5*(rho2*rho2-rho1*rho1)*dphiz;
		double mass = volume*localdens;
                if (mass > 0.0){imass[bin] = 1.0/mass;}
	      }
	  }
      }
    }

   if(verbose > 1){
    printf("Number of radial bins: %ld\n",nr);
    printf("Minimum r value (cm): %12.5E\n",rmin);
    printf("Bin width radial (cm): %12.5E\n",dr);
    printf("Number of depth bins: %ld\n",nz);
    printf("Minimum z value (cm): %12.5E\n",zmin);
    printf("Bin width depth (cm): %12.5E\n",dz);
    printf("Number of angular bins: %ld\n",nphi);
    printf("Bin width angular (rad): %12.5E\n",dphi);
  }  
    
  //Register data to dump
  dump.toDump(imass,nbins);
  dump.toDump(edep,nbins);
  dump.toDump(edep2,nbins);
  dump.toDump(&nbins,1);

  //Register data to create images
  unsigned elements[] = {
    static_cast<unsigned>(nr),
    static_cast<unsigned>(nphi),
    static_cast<unsigned>(nz)};
  
  float delements[] = {
    static_cast<float>(dr),
    static_cast<float>(dphi),
    static_cast<float>(dz)};

  // ** Calculate the origin
  double origin[3] = {0.0,0.0,0.0};
  // Add the mesh origin
  origin[0] = rmin;
  origin[1] = 0.0;
  origin[2] = zmin;

  addImage<double>("cylindricalDose",3,elements,delements,origin,
		   [=](unsigned long long nhist,
		       size_t i, double& sigma) -> double{

		     const double dhists = static_cast<double>(nhist);
		     const double fact = imass[i];
		     const double q = edep[i]/dhists;
		     sigma = edep2[i]/dhists-q*q;
		     sigma = sigma > 0.0 ? fact*sqrt(sigma/dhists) : 0.0;
		     
		     return q*fact;
		   });  
  
  return 0;   
    
    
}

void pen_CylindricalDoseDistrib::tally_endSim(const unsigned long long /*nhist*/){
        
  flush();
}

void pen_CylindricalDoseDistrib::saveData(const unsigned long long nhist) const{
   
  char buffer[81];
  FILE*out = 0;
  int nzef;
  const double twothird=2.0/3.0;

  //Prepare z factors and filename
  if(idz == 0.0)   // No z bins
    {
      nzef = 0;
      strcpy(buffer,"radialDoseDistrib.dat");
    }
  else
    {
      nzef = nz;
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
  fprintf(out,"# %4ld %12.5E %12.5E\n",nr,rmin,dr);
  fprintf(out,"# Depth (z) info -> Number of bins, minimum value, bin width (cm):\n");
  fprintf(out,"# %4ld %12.5E %12.5E\n",nz,zmin,dz);
  fprintf(out,"# Angular (phi) info -> Number of bins, bin width (rad):\n");
  fprintf(out,"# %4ld %12.5E\n",nphi,dphi);
    
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
    fprintf(out,"phiBinIndex : phiLow(DEG) : phiMiddle(DEG) : ");
    if (nzef != 0)
      {
	fprintf(out,"zBinIndex : zLow(cm) : zMiddle(cm) : ");
      }
  }
  fprintf(out,"   dose   : +-2sigma : energy(eV) : +-2sigma :\n");
                
  // Write data:    
  double invn = 1.0/static_cast<double>(nhist);  
  const double rad2deg = 180.0/M_PI;
  for(long int k = 0; k < nz; ++k){
    double z = zmin + dz*static_cast<double>(k);
    double zmiddle = z + dz*0.5;
    // Since z is not written, give at least a summary
    if(nzef != 0 && prtxyz != 1){
      fprintf(out,"# zBinIndex=%-3ld zMiddle(cm)=%12.5E\n",k,zmiddle);
    }
    for(long int j = 0; j < nphi; ++j){
      double phi       = dphi*(static_cast<double>(j));
      double phimiddle = phi + dphi*0.5;
      if(nphi > 1 && prtxyz != 1){
	fprintf(out,"# phiBinIndex=%-3ld phiMiddle(DEG)=%12.5E\n",j,phimiddle);
      }
      for(long int i = 0; i < nr; ++i){
	double r = rmin + dr*static_cast<double>(i);
	double rave;
	// Note: dr2 = (r+dr)^2-r^2
	double invdr2 = 1.0/(dr*(dr+2.0*r));   
	if(prtxyz == 1){
	  // Average with weight(r)~r
	  rave = twothird*((r+dr)*(r+dr)*(r+dr)-r*r*r)*invdr2; 
	  fprintf(out,"   %5ld    %12.5E %12.5E"
		  "    %5ld       %12.5E     %12.5E ",
		  i,r,rave,j,phi*rad2deg,phimiddle*rad2deg);
	}
	if(nzef != 0 && prtxyz == 1){
	  fprintf(out,"  %5ld    %12.5E %12.5E",k,z,zmiddle);
	}
	long int bin = nr*(k*nphi + j) + i;
	double fact = imass[bin];
	double q = edep[bin]*invn;
	double sigma = sqrt((edep2[bin]*invn-q*q)*invn > 0.0 ? (edep2[bin]*invn-q*q)*invn : 0.0);
	fprintf(out," %12.5E  %8.1E %12.5E  %8.1E\n",
		q*fact,2.0*sigma*fact,q,2.0*sigma);
      }
      // Separate data blocks
      if (nphi > 1){ fprintf(out,"  \n");}	
    }
    // Separate data blocks
    if (nr > 1){ fprintf(out,"  \n");}
  }
        
  fclose(out); 
}

int pen_CylindricalDoseDistrib::sumTally(const pen_CylindricalDoseDistrib& tally){
    
    
    if(nbins != tally.nbins)
        return -1;
    
    for(int i = 0; i < nbins; ++i){
        edep[i] += tally.edep[i];
    }
    for(int i = 0; i < nbins; ++i){
        edep2[i] += tally.edep2[i];
    }

  return 0;
    
    
}


REGISTER_COMMON_TALLY(pen_CylindricalDoseDistrib, CYLINDRICAL_DOSE_DISTRIB)
