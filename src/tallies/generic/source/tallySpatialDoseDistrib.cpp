
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

 
#include "tallySpatialDoseDistrib.hh"

void pen_SpatialDoseDistrib::updateEdepCounters(const double dE,
						const unsigned long long nhist,
						const double X,
						const double Y,
						const double Z,
						const double WGHT,
						const unsigned MAT){
    
  int bin;

  //Avoid count energy deposition in range (zmin-dz,zmin)
  if(Z < zmin) return;
  
  // Check if particle is inside tally region:
  int k = (Z - zmin)*idz;
  if(k < 0 || k >= nz){return;}

  //***************************
  //* DEPTH ENERGY DEPOSITION
  //***************************

  if(nhist > nlastdepth[k]){
    edepth[k]  += edepthtmp[k];
    edepth2[k] += edepthtmp[k]*edepthtmp[k];
    edepthtmp[k] = dE*WGHT*imatDens[MAT-1];
    // Add 1/2 to avoid roundoff errors
    nlastdepth[k]  = nhist;
  }
  else{
    edepthtmp[k] += dE*WGHT*imatDens[MAT-1];
  }

  //***************************
  //* 3D ENERGY DEPOSITION
  //***************************

  //Avoid count energy deposition in ranges (xmin-dx,xmin) and (ymin-dy,ymin)
  if(X < xmin || Y < ymin) return;
  
  int i = (X - xmin)*idx;
  if(i < 0 || i >= nx){return;}
  int j = (Y - ymin)*idy;
  if(j < 0 || j >= ny){return;}

  //Map i,j,k into a single index
  bin = i + j*nx + k*nxy;
  
  //Transfer partial to totals when a new history visits
  if(nhist > nlast[bin])
    {
      edep[bin]  += edptmp[bin];
      edep2[bin] += edptmp[bin]*edptmp[bin];
      edptmp[bin] = dE*WGHT;
      // Add 1/2 to avoid roundoff errors
      nlast[bin]  = nhist;             
    }
  else
    {
      //Same history as before, keep adding to temp counter
      edptmp[bin] += dE*WGHT;       
    }
}


void pen_SpatialDoseDistrib::flush(){
  // Dump temp counters and obtain max score:
  for(long int i = 0; i < nbin; ++i)
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
  
  for(int i = 0; i < nz; ++i)
    {
      // Skip empty bins
      if (nlastdepth[i] == 0){ continue;}      
      // Transfer temp counter
      edepth[i]  += edepthtmp[i];   
      edepth2[i] += edepthtmp[i]*edepthtmp[i];
      // Reset counter
      edepthtmp[i]= 0.0;                  
      // Reset last visited to avoid recounting in next report
      nlastdepth[i] = 0;                  
    }
  
}

void pen_SpatialDoseDistrib::tally_localEdep(const unsigned long long nhist,
					     const pen_KPAR /*kpar*/,
					     const pen_particleState& state,
					     const double dE){
    
  if(dE == 0.0){return;}  
  //Nothing to deposit.
  //Energy deposited at material
  updateEdepCounters(dE, nhist, state.X, state.Y, state.Z, state.WGHT, state.MAT);
}

void pen_SpatialDoseDistrib::tally_beginPart(const unsigned long long nhist,
					     const unsigned /*kdet*/,
					     const pen_KPAR /*kpar*/,
					     const pen_particleState& state){
 
  //Extract energy from material to create new particle
  updateEdepCounters(-state.E, nhist, state.X, state.Y, state.Z, state.WGHT, state.MAT);
}

void pen_SpatialDoseDistrib::tally_beginHist(const unsigned long long nhist,
					     const unsigned /*kdet*/,
					     const pen_KPAR /*kpar*/,
					     const pen_particleState& state){

  if(state.MAT > 0){
    //Particle created at non void volume. Add particle energy to compensate
    //substracted one when beginPart will be called.
    updateEdepCounters(state.E, nhist, state.X, state.Y, state.Z, state.WGHT, state.MAT);
  }
}

void pen_SpatialDoseDistrib::tally_step(const unsigned long long nhist,
					const pen_KPAR /*kpar*/,
					const pen_particleState& state,
					const tally_StepData& stepData){

  if(stepData.softDE == 0.0){return;}  //Nothing to deposit.
  //Energy deposited
  updateEdepCounters(stepData.softDE, nhist,
		     stepData.softX, stepData.softY, stepData.softZ,
		     state.WGHT, stepData.originMAT);  
}


void pen_SpatialDoseDistrib::tally_move2geo(const unsigned long long nhist,
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
    updateEdepCounters(state.E, nhist, state.X, state.Y, state.Z, state.WGHT, state.MAT);
  }
}



void pen_SpatialDoseDistrib::clear()
{
  if(nlast != nullptr){
    free(nlast);
    nlast = nullptr;
  }
     
  if(edptmp != nullptr){
    free(edptmp);
    edptmp = nullptr;
  }
     
  if(edep != nullptr){
    free(edep);
    edep = nullptr;
  }
     
  if(edep2 != nullptr){
    free(edep2);
    edep2 = nullptr;
  }
     
  if(ivoxMass != nullptr){
    free(ivoxMass);
    ivoxMass = nullptr;
  }

  if(nlastdepth != nullptr){
    free(nlastdepth);
    nlastdepth = nullptr;
  }
  
  if(edepthtmp != nullptr){
    free(edepthtmp);
    edepthtmp = nullptr;
  }

  if(edepth != nullptr){
    free(edepth);
    edepth = nullptr;
  }

  if(edepth2 != nullptr){
    free(edepth2);
    edepth2 = nullptr;
  }
  
  nx = ny = nz = nxy = nbin = 0;
  dx = dy = dz = 0.0;
  idx = idy = idz = 1.0e35;
  xmin = ymin = zmin = 0.0;
     
}



int pen_SpatialDoseDistrib::configure(const wrapper_geometry& geometry,
				      const abc_material* const materials[constants::MAXMAT],
				      const pen_parserSection& config,
				      const unsigned verbose
				      ){
  
  //Clear previous configuration
  clear();

  //**********************
  //* Material densities *
  //**********************

  for(unsigned i = 0; i < constants::MAXMAT; i++){
    if(materials[i] != nullptr)
      imatDens[i] = 1.0/materials[i]->readDens();
  }
  
  int err;
  double xmax,ymax,zmax;

  //*******
  //*  X  *
  //*******
  
  err = config.read("xmin", xmin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'xmin' in configuration. Double expected\n");
    }
    return -1;
  }
  
  err = config.read("xmax", xmax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'xmax' in configuration. Double expected\n");
    }
    return -2;
  }
  
  err = config.read("nx", nx);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'nx' in configuration. Double expected\n");
    }
    return -3;
  }
      
  if(nx <= 0){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Invalid entry, 'nx' must be at least 1\n");
    }
    return -5;
  }

  if(xmin >= xmax){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Invalid entry, 'xmax' must be greater than 'xmin'\n");
    }
    return -6;
  }  
    
  //Bin width
  dx = (xmax-xmin)/double(nx);
  idx = 1.0/dx;

  //*******
  //*  Y  *
  //*******    
    
  err = config.read("ymin", ymin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'ymin' in configuration. Double expected\n");
    }
    return -7;
  }
  
  err = config.read("ymax", ymax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'ymax' in configuration. Double expected\n");
    }
    return -8;
  }
  
  err = config.read("ny", ny);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'ny' in configuration. Double expected\n");
    }
    return -9;
  }
      
  if(ny <= 0){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Invalid entry, 'ny' must be at least 1\n");
    }
    return -10;
  }

  if(ymin >= ymax){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Invalid entry, 'ymax' must be greater than 'ymin'\n");
    }
    return -11;
  }  
    
  //Bin width
  dy = (ymax-ymin)/double(ny);
  idy = 1.0/dy;


  //*******
  //*  Z  *
  //*******    
    
  err = config.read("zmin", zmin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'zmin' in configuration. Double expected\n");
    }
    return -12;
  }
  
  err = config.read("zmax", zmax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'zmax' in configuration. Double expected\n");
    }
    return -13;
  }
  
  err = config.read("nz", nz);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Unable to read 'nz' in configuration. Double expected\n");
    }
    return -14;
  }
      
  if(nz <= 0){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Invalid entry, 'nz' must be at least 1\n");
    }
    return -15;
  }

  if(zmin >= zmax){
    if(verbose > 0){
      printf("pen_SpatialDoseDistrib:configure: Error: Invalid entry, 'zmax' must be greater than 'zmin'\n");
    }
    return -16;
  }  
    
  //Bin width
  dz = (zmax-zmin)/double(nz);
  idz = 1.0/dz;
  
  //**********************
  //* Density calculation
  //**********************
  
  nxy = nx*ny;
  nbin = nxy*nz;
  
    
  voxVol = dx*dy*dz;
    


  //**************
  //*  Allocate  *
  //**************
    
  //Allocate memory for spatial 3D distribution
  edptmp    = (double*) calloc(nbin,sizeof(double));
  edep      = (double*) calloc(nbin,sizeof(double));
  edep2     = (double*) calloc(nbin,sizeof(double));
  nlast     = (unsigned long long*) calloc(nbin,sizeof(unsigned long long));
  ivoxMass  = (double*) calloc(nbin,sizeof(double));

  //Allocate memory for depth dose distribution
  nlastdepth= (unsigned long long*) calloc(nz,sizeof(unsigned long long));
  edepthtmp = (double*) calloc(nz,sizeof(double));
  edepth    = (double*) calloc(nz,sizeof(double));
  edepth2   = (double*) calloc(nz,sizeof(double));
  
  pen_particleState state;
    
  const int nDivisions = 5;
  double inDiv = 1.0/double(nDivisions);
  double ddx = dx*inDiv;
  double ddy = dy*inDiv;
  double ddz = dz*inDiv;
  double subVoxVol = ddx*ddy*ddz;
    
    
  for(int k = 0; k < nz; k++)
    {
      //This is to locate a point and find its material
      double binZpos = zmin + dz*k;
      int ibinZ = k*nxy;
	
      for(int j = 0; j < ny; j++)
	{
	  double binYpos = ymin + dy*j;
	  int ibinY = ibinZ + j*nx;
	    
	  for(int i = 0; i < nx; i++)
            {

	      double binXpos = xmin + dx*i;
	      int bin = ibinY + i;
	                      
	      double localdens = 0.0;
               
	      for (int kk = 0; kk < nDivisions; kk++)
                {
		  state.Z = binZpos + ddz*((double)kk+0.5);

		  //Ensure direction to bin center
		  if(kk < nDivisions/2)
		    state.W = 1.0;
		  else
		    state.W = -1.0;
		    
		  for (int jj = 0; jj < nDivisions; jj++)
                    {
		      state.Y = binYpos + ddy*((double)jj+0.5);

		      //Ensure direction to bin center
		      if(jj < nDivisions/2)
			state.V = 1.0;
		      else
			state.V = -1.0;
		      
		      for(int ii = 0; ii < nDivisions; ii++)
                        {
			  state.X = binXpos + ddx*((double)ii+0.5);
			  
			  //Ensure direction to bin center
			  if(ii < nDivisions/2)
			    state.U = 1.0;
			  else
			    state.U = -1.0;
			  
			  geometry.locate(state);
                            
			  if(state.MAT > 0)
                            {
			      localdens += materials[state.MAT-1]->readDens();
                            }
                        }
                    }
                }
	      if(localdens > 0.0){
		double voxMass = localdens*subVoxVol;
		ivoxMass[bin] = 1.0/voxMass;
	      }
	      else
		ivoxMass[bin] = 1.0e35;
            }
        }
    }
    
  if(verbose > 1){
    printf("Number of x bins: %d\n",nx);
    printf("Minimum x value (cm): %12.5E\n",xmin);
    printf("x Bin width (cm): %12.5E\n",dx);
    printf("Number of y bins: %d\n",ny);
    printf("Minimum y value (cm): %12.5E\n",ymin);
    printf("y Bin width (cm): %12.5E\n",dy);
    printf("Number of z bins: %d\n",nz);
    printf("Minimum z value (cm): %12.5E\n",zmin);
    printf("z Bin width (cm): %12.5E\n",dz);
    printf("Number of bins %ld\n", nbin);
  }
    
  //Register data to dump
  dump.toDump(ivoxMass,nbin);
  dump.toDump(edep,nbin);
  dump.toDump(edep2,nbin);
  dump.toDump(edepth,nz);
  dump.toDump(edepth2,nz);
    
    
    
  return 0; 
    
            
}

void pen_SpatialDoseDistrib::tally_endSim(const unsigned long long /*nhist*/){        
  flush();
}
    
    
 
void pen_SpatialDoseDistrib::saveData(const unsigned long long nhist) const{
  
  char buffer[81];
  FILE*out;
  double q, sigma, fact, x, y, z;
  double xmiddle, ymiddle, zmiddle;

  double invn = 1.0/static_cast<double>(nhist);

  //*************************
  //*  3D DOSE DISTRIBUTION *
  //*************************
  
  strcpy(buffer,"spatialDoseDistrib.dat");
     
  out = fopen(buffer,"w");
  if (out == NULL)
    {
      printf("*********************************************\n");
      printf("pen_SpatialDoseDistrib: Error: Cannot open output data file\n");
      printf("*********************************************\n");
      fclose(out);                            
      return;                      
    }
    
    
  // Write header 
  fprintf(out,"#------------------------------------------------------------\n");
  fprintf(out,"# PenRed: Spatial dose distribution report\n");
  fprintf(out,"# Dose units are: eV/g per history\n");

  fprintf(out,"#\n");
  fprintf(out,"# No. of bins in x,y,z directions and total:\n");
  fprintf(out,"#  %d %d %d %ld\n",nx,ny,nz,nbin);
  fprintf(out,"# Min values for x,y,z(cm):\n");
  fprintf(out,"#  %12.5E %12.5E %12.5E\n",xmin,ymin,zmin);
  fprintf(out,"# Bin widths for x,y,z(cm):\n");
  fprintf(out,"#  %12.5E %12.5E %12.5E\n",dx,dy,dz);
  fprintf(out,"#\n");
  fprintf(out,"# For plotting purposes, two values per bin coordinate are given, namely, the low end and the middle point of each bin.\n");
  fprintf(out,"#\n");
  fprintf(out,"# ");
  fprintf(out,"xBinIndex : xLow(cm) : xMiddle(cm) : ");
  fprintf(out,"yBinIndex : yLow(cm) : yMiddle(cm) : ");
  fprintf(out,"zBinIndex : zLow(cm) : zMiddle(cm) : ");
  fprintf(out,"dose : +-2sigma\n");
    
    
  //Write data
    
  for(int k = 0; k < nz; k++)
    {
      z = zmin + dz*k;
      zmiddle = z + dz*0.5;
      int kbin = k*nxy;
      
      fprintf(out,"# zBinIndex=%d zMiddle(cm)=%12.5E\n",k,zmiddle);

      for(int j = 0; j < ny; j++)
        {
	  y = ymin+dy*j;
	  ymiddle = y+dy*0.5;
	  int jkbin = j*nx+kbin;

	  fprintf(out,"# yBinIndex=%d yMiddle(cm)=%12.5E\n",j,ymiddle);
        
	  for(int i = 0; i < nx; i++)
            {
	      x = xmin+dx*i;
	      xmiddle = x+dx*0.5;
	      fprintf(out," %5d %12.5E %12.5E",i,x,xmiddle);
	      fprintf(out," %5d %12.5E %12.5E",j,y,ymiddle);
	      fprintf(out," %5d %12.5E %12.5E",k,z,zmiddle);
	      // Map i,j,k into a single index
	      int bin = i + jkbin; 
                
	      fact = ivoxMass[bin];
	      q = edep[bin]*invn;
	      sigma = edep2[bin]*invn - q*q;
	      if(sigma > 0.0)
                {
		  sigma = sqrt(sigma*invn)*fact;
                }
	      else
                {
		  sigma = 0.0;
                }
	      q = q*fact;
	      fprintf(out," %12.5E %7.1E\n",q,2.0*sigma);
            }
        }
      fprintf(out,"  \n  \n"); 
    }

  fclose(out);
  out = nullptr;
  
  //***************************
  //* DEPTH DOSE DISTRIBUTION *
  //***************************

  out = fopen("depth-dose.dat","w");

  fprintf(out,"#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  fprintf(out,"# [SECTION REPORT DEPTH DOSE DISTRIB]\n");
  fprintf(out,"# Dose units are: eV/(g/cm^2) per history\n");

  fprintf(out,"#\n");
  fprintf(out,"# No. of Z bins:\n");
  fprintf(out,"#  %d\n",nz);
  fprintf(out,"# Min value for z(cm):\n");
  fprintf(out,"#  %12.5E\n",zmin);
  fprintf(out,"# Bin widths for z(cm):\n");
  fprintf(out,"#  %12.5E\n",dz);
  fprintf(out,"#\n");
  fprintf(out,"#\n");
  fprintf(out,"# ");
  fprintf(out,"zBinIndex : zLow(cm) : zMiddle(cm) : ");
  fprintf(out,"dose : +-2sigma\n");
  
  for(int k = 0; k < nz; k++){
    
    z = zmin+dz*k;
    zmiddle = z+dz*0.5;
    fprintf(out," %5d %12.5E %12.5E",k,z,zmiddle);
    
    q = edepth[k]*invn;
    sigma = edepth2[k]*invn - q*q;
    if(sigma > 0.0)
      {
	sigma = sqrt(sigma*invn)*idz;
      }
    else
      {
	sigma = 0.0;
      }
    q *= idz;
    fprintf(out," %12.5E %7.1E\n",q,2.0*sigma);    
  }

  fclose(out);
}
 
 
int pen_SpatialDoseDistrib::sumTally(const pen_SpatialDoseDistrib& tally){
    
  if(nbin != tally.nbin)
    return -1;
  if(nz != tally.nz)
    return -2;
  
  for(long int i = 0; i < nbin; ++i){
    edep[i] += tally.edep[i];
  }
  for(long int i = 0; i < nbin; ++i){
    edep2[i] += tally.edep2[i];
  }

  for(long int i = 0; i < nz; ++i){
    edepth[i] += tally.edepth[i];
  }
  for(long int i = 0; i < nz; ++i){
    edepth2[i] += tally.edepth2[i];
  }
  
  return 0;
    
    
}

REGISTER_COMMON_TALLY(pen_SpatialDoseDistrib, SPATIAL_DOSE_DISTRIB)

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
