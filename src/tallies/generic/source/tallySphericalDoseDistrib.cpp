
//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

 
#include "tallySphericalDoseDistrib.hh"

pen_SphericalDoseDistrib::pen_SphericalDoseDistrib() :
  pen_genericTally( USE_LOCALEDEP |
		    USE_BEGINPART |
		    USE_SAMPLEDPART |
		    USE_STEP |
		    USE_ENDSIM |
		    USE_MOVE2GEO),
  nbins(0)

{
  setResultsGenerator<0>
    ([this](const unsigned long long nhists) -> penred::measurements::results<double, 3>{

      const double invn = 1.0/static_cast<double>(nhists);

      //Create results
      penred::measurements::results<double, 3> results;
      results.initFromLists
	({static_cast<unsigned long>(nr),
	   static_cast<unsigned long>(ntheta),
	   static_cast<unsigned long>(nphi)},
	  {penred::measurements::limitsType(rmin, rmin + static_cast<double>(nr)*dr),
	   penred::measurements::limitsType(0.0, 180.0),
	   penred::measurements::limitsType(0.0, 360.0)});
	  
      results.description =
	"PenRed: Spherical dose distribution report\n"
	"Dose units are: eV/g per history\n";
  
      results.setDimHeader(0, "r (cm)");
      results.setDimHeader(1, "polar (deg)");
      results.setDimHeader(2, "azimuth (deg)");
      results.setValueHeader("Dose (eV/g hist)");

      for(unsigned long i = 0; i < nbins; ++i){
	
	const double fact = imass[i];
	const double q = edep[i]*invn;
	double sigma = (edep2[i]*invn-(q*q))*invn;
	if(sigma > 0.0){
	  sigma = sqrt(sigma);
	}
	else{
	  sigma = 0.0;
	}

	results.data[i] = q*fact;
	results.sigma[i] = sigma*fact;
      }

      return results;
    });
}

void pen_SphericalDoseDistrib::updateEdepCounters(const double dE,
						  const unsigned long long nhist,
						  const double X,
						  const double Y,
						  const double Z,
						  const double WGHT){

  double rho2 = X*X + Y*Y;
  double r2 = rho2 + Z*Z;
  //Not in ROI
  if (r2 < r2min || r2 >= r2max){return;}
  //Bin index where score will be added
  long int ir = static_cast<long int>((sqrt(r2) - rmin)*idr);

  long int itheta = 0;
  if(ntheta > 1){
    if(rho2 > 1.0e-15 || fabs(Z) > 1.0e-12){
      double theta = atan2(sqrt(rho2),Z);
      if(std::signbit(theta))
	theta = fabs(theta);
      itheta = theta*idtheta;
      if(itheta >= ntheta) itheta--;
    }
  }
  
  long int iphi = 0;
  if(nphi > 1){
    if(rho2 > 1.0e-15){
      double phi = atan2(Y,X);
      if(std::signbit(phi))
	phi += 2.0*M_PI;
      iphi = phi*idphi;
      if(iphi >= nphi) iphi--;
    }
  }

  long unsigned ibin = nr*(iphi*ntheta + itheta) + ir;
  
  //Transfer partial to totals when a new history visits
  if(nhist > nlast[ibin])
    {
      edep[ibin]  += edptmp[ibin];
      edep2[ibin] += edptmp[ibin]*edptmp[ibin];
      edptmp[ibin] = dE*WGHT;
      // Add 1/2 to avoid roundoff errors
      nlast[ibin]  = nhist;             
    }
  else
    {
      //Same history as before, keep adding to temp counter
      edptmp[ibin] += dE*WGHT;       
    }
}

void pen_SphericalDoseDistrib::flush(){
  // Dump temp counters and obtain max score:
  for(long unsigned i = 0; i < nbins; ++i)
    {
      // Skip empty bins
      if (nlast[i] == 0){ continue;}      
      // Transfer temp counter
      edep[i]  += edptmp[i];   
      edep2[i] += edptmp[i]*edptmp[i];
      // Reset counter
      edptmp[i]= 0.0;                  
      // Reset last visited to avoid recounting in next report
      nlast[i] = 0;                  
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

void pen_SphericalDoseDistrib::tally_sampledPart(const unsigned long long nhist,
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
					const unsigned verbose){

  double rmax;
  // A bit more than one
  const double oneplus = 1.0+1.0E-12;  
  // A bit less than one
  const double oneminus= 1.0-1.0E-12;  
  int err; 
    
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
  int auxInt;
  err = config.read("nr", auxInt);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SphericalDoseDistrib:configure:unable to read 'nr' in configuration. Integrer expected\n");
    }
    return -3;
  }
  nr = auxInt;

  err = config.read("ntheta", auxInt);
  if(err != INTDATA_SUCCESS){
    auxInt = 1;
  }
  if(auxInt <= 0){
    if(verbose > 0){
      printf("pen_SphericalDoseDistrib:configure: Error: "
	     "number of theta bins must be, as least 1\n"
	     "    Specified ntheta: %d\n",auxInt);
    }
    return -3;      
  }
  ntheta = auxInt;

  err = config.read("nphi", auxInt);
  if(err != INTDATA_SUCCESS){
    auxInt = 1;
  }
  if(auxInt <= 0){
    if(verbose > 0){
      printf("pen_SphericalDoseDistrib:configure: Error: "
	     "number of phi bins must be, as least 1\n"
	     "    Specified nphi: %d\n",auxInt);
    }
    return -3;      
  }
  nphi = auxInt;
    
    
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

  if(nr < 1){
    if(verbose > 0){
      printf("pen_SphericalDoseDistrib:configure:Invalid entry, nr must be at least 1.\n");
    }
    return -6;
  }

  //Obtain the total number of bins
  nbins = nr*ntheta*nphi;

  if(verbose > 1)
    printf("Number of bins:\n"
	   "   radial: %ld\n"
	   "    theta: %ld\n"
	   "      phi: %ld\n"
	   ,nr,ntheta,nphi);

    
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
    
  dr = (rmax - rmin)/static_cast<double>(nr);
  idr = 1.0/dr;
  // (1+-eps) tolerances avoid out of range bin indexes in tally routine
  r2min = rmin*rmin*oneplus;  
  r2max = rmax*rmax*oneminus;

  dtheta  = M_PI/static_cast<double>(ntheta);
  dphi    = 2.0*M_PI/static_cast<double>(nphi);
  idtheta = 1.0/dtheta;
  idphi   = 1.0/dphi;
    
  // Init arrays:
  nlast.resize(nbins);
  edptmp.resize(nbins);
  edep.resize(nbins);
  edep2.resize(nbins);
  imass.resize(nbins);

  std::fill(nlast.begin() , nlast.end() , 0.0);
  std::fill(edptmp.begin(), edptmp.end(), 0.0);
  std::fill(edep.begin()  , edep.end()  , 0.0);
  std::fill(edep2.begin() , edep2.end() , 0.0);
  std::fill(imass.begin() , imass.end() , 0.0);
    
  pen_particleState state;
    
  // For locate() to operate properly
  state.U = 0.0;  
  state.V = 0.0;
  state.W = 1.0;
  state.X = 0.0;
  state.Y = 0.0;
    

  for(long int k = 0; k < nphi; ++k){
    double phi = (static_cast<double>(k)+0.5)*dphi;
    double sphi = sin(phi);
    double cphi = cos(phi);
    for(long int j = 0; j < ntheta; ++j){
      double theta1 = static_cast<double>(j)*dtheta; 
      double theta  = theta1+0.5*dtheta;
      double theta2 = theta1+dtheta;

      double stheta = sin(theta);
      double ctheta = cos(theta);
	
      double ctheta1 = cos(theta1);
      double ctheta2 = cos(theta2);
      for(long int i = 0; i < nr; ++i){
	//This is to locate a point and find its material
	double r1 = rmin + dr*static_cast<double>(i);
	double r = r1 + dr*0.5;
	double rstheta = r*stheta;
	state.X = rstheta*cphi;
	state.Y = rstheta*sphi;
	state.Z = r*ctheta;
	geometry.locate(state);
	if(state.MAT > 0)
	  {
	    //Local mass density
	    double localdens = materials[state.MAT-1]->readDens(); 
	    if (localdens > 0.0){
	      double volume =
		dphi*(ctheta1-ctheta2)*(pow(r1+dr,3)-pow(r1,3))/3.0;
	      long unsigned ibin = nr*(k*ntheta + j) + i;
	      imass[ibin] = 1.0/(volume*localdens);
	    }
	  }
      }
    }
  }

  if(verbose > 1){
    printf("Number of radial bins: %ld\n",nr);
    printf("Minimum r value (cm): %12.5E\n",rmin);
    printf("Bin width (cm): %12.5E\n",dr);
  }  
    
  //Register data to dump
  dump.toDump(nlast.data(),nbins);
  dump.toDump(imass.data(),nbins);
  dump.toDump(edep.data(),nbins);
  dump.toDump(edep2.data(),nbins);
  dump.toDump(&nbins,1);
    
    
  //Register data to create images
  unsigned elements[] = {
    static_cast<unsigned>(nr),
    static_cast<unsigned>(ntheta),
    static_cast<unsigned>(nphi)};
  
  float delements[] = {
    static_cast<float>(dr),
    static_cast<float>(dtheta),
    static_cast<float>(dphi)};

  // ** Calculate the origin
  double origin[3];
  // Add the mesh origin
  origin[0] = rmin;
  origin[1] = 0.0;
  origin[2] = 0.0;

  addImage<double>("spatialDose",3,elements,delements,origin,
		   [=](unsigned long long nhist,
		       size_t i, double& sigma) -> double{

		     const double dhists = static_cast<double>(nhist);
		     const double fact = imass[i];             

		     const double q = edep[i]/dhists;
		     sigma = edep2[i]/dhists-q*q;
		     if(sigma > 0.0){
		       sigma = fact*sqrt(sigma/dhists);
		     }
		     else{
		       sigma = 0.0;
		     }
		     
		     return q*fact;
		   });
    
  return 0;   
    
    
}

void pen_SphericalDoseDistrib::tally_endSim(const unsigned long long /*nhist*/){
        
  flush();
}

void pen_SphericalDoseDistrib::saveData(const unsigned long long nhist) const{
   
  //Prepare output files:
  FILE*out = 0;
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
  fprintf(out,"# No. of bins (r,theta,phi):\n");
  fprintf(out,"#  %ld %ld %ld\n",nr,ntheta,nphi);
  fprintf(out,"# Min r and bin width (cm):\n");
  fprintf(out,"#  %12.5E %12.5E\n", rmin, dr);
  if (prtxyz == 1){
    fprintf(out,"#\n");
    fprintf(out,"# Coordinate print enabled. Two radial values for each bin will be printed, the low and an average radius. the average radial value is a weighted proportionally to r^2.\n");
  }
    
  fprintf(out,"#\n");
  fprintf(out,"# ");
  if (prtxyz == 1){
    fprintf(out,"rBinIndex :  rLow(cm)  :  rAve(cm)  :"
	    " thetaBin : theta(DEG) :"
	    "  phiBin  :  phi(DEG)  : ");
  }
  fprintf(out,"dose(eV/g) :+-2sigma: edep(eV/cm) :+-2sigma\n");
                
  // Write data
  double invn = 1.0/static_cast<double>(nhist);  
    

  const double threefour = 3.0/4.0;
  for(long int i = 0; i < nr; ++i){
    double r = rmin+dr*static_cast<double>(i);
    double rnext = r+dr;
    double invdr3 = 1.0/(rnext*rnext*rnext-r*r*r);
    // Average with weight(r)~r^2
    double rave = threefour*(pow(rnext,4)-pow(r,4))*invdr3; 

    for(long int j = 0; j < ntheta; ++j){
      double theta  = static_cast<double>(j)*dtheta;
      for(long int k = 0; k < nphi; ++k){
	double phi = static_cast<double>(k)*dphi;
	long unsigned ibin = nr*(k*ntheta + j) + i;
	// This is 1/Delta_mass
	double fact = imass[ibin];             
	if (prtxyz == 1){
	  fprintf(out," %6ld      %12.5E %12.5E %6ld    %12.5E  %6ld    %12.5E",
		  i,r,rave,j,theta,k,phi);
	}
	double q = edep[ibin]*invn;
	double sigma = (edep2[ibin]*invn-(q*q))*invn;
	if(sigma > 0.0){
	  sigma = sqrt(sigma);
	}
	else{
	  sigma = 0.0;
	}
        
	double sigma2 = 2.0*sigma;
	fprintf(out," %12.5E %8.1E %12.5E  %8.1E\n",
		q*fact,sigma2*fact,q,sigma2);
      }
      fprintf(out,"\n");
    }
    fprintf(out,"\n");
  }
  fclose(out); 
}

int pen_SphericalDoseDistrib::sumTally(const pen_SphericalDoseDistrib& tally){
    
  if(nbins != tally.nbins)
    return -1;
    
  for(long unsigned i = 0; i < nbins; ++i){
    edep[i] += tally.edep[i];
  }
  for(long unsigned i = 0; i < nbins; ++i){
    edep2[i] += tally.edep2[i];
  }

  return 0;
    
    
}


REGISTER_COMMON_TALLY(pen_SphericalDoseDistrib)
