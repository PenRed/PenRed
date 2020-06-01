
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

#ifdef _PEN_USE_DICOM_
 
#include "tallyDICOMDoseDistrib.hh"

void pen_DICOMDoseDistrib::updateEdepCounters(const double dE,
					      const unsigned long long nhist,
					      const double X,
					      const double Y,
					      const double Z,
					      const double WGHT){

  //Avoid count energy deposition in range (zmin-dz,zmin)
  if(Z < zmin) return;
  
  // Check if particle is inside tally region:
  long int k = (Z - zmin)*idz;
  if(k < 0 || k >= nz){return;}

  //Avoid count energy deposition in ranges (xmin-dx,xmin) and (ymin-dy,ymin)
  if(X < xmin || Y < ymin) return;
  
  long int i = (X - xmin)*idx;
  if(i < 0 || i >= nx){return;}
  long int j = (Y - ymin)*idy;
  if(j < 0 || j >= ny){return;}

  //Map i,j,k into a single index
  long int bin = i + j*static_cast<long int>(nx) + k*nxy;
  
  //Transfer partial to totals when a new history visits
  double wde = dE*WGHT;
  if(nhist > nlast[bin])
    {
      edep[bin]  += edptmp[bin];
      edep2[bin] += edptmp[bin]*edptmp[bin];
      edptmp[bin] = wde;
      // Add 1/2 to avoid roundoff errors
      nlast[bin]  = nhist;             
    }
  else
    {
      //Same history as before, keep adding to temp counter
      edptmp[bin] += wde;       
    }

  //Record energy deposition in the corresponding contour
  if(ncontours > 0){
    int icont = contourVox[bin];
    if(icont >= 0){
      contEdptmp[icont] += wde*ivoxMass[bin];
    }
  }
}


void pen_DICOMDoseDistrib::flush(){
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
}

void pen_DICOMDoseDistrib::tally_localEdep(const unsigned long long nhist,
					   const pen_KPAR /*kpar*/,
					   const pen_particleState& state,
					   const double dE){
    
  if(dE == 0.0){return;}  
  //Nothing to deposit.
  //Energy deposited at material
  updateEdepCounters(dE, nhist, state.X, state.Y, state.Z, state.WGHT);
}

void pen_DICOMDoseDistrib::tally_beginPart(const unsigned long long nhist,
					   const unsigned /*kdet*/,
					   const pen_KPAR /*kpar*/,
					   const pen_particleState& state){
 
  //Extract energy from material to create new particle
  updateEdepCounters(-state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
}

void pen_DICOMDoseDistrib::tally_beginHist(const unsigned long long nhist,
					   const unsigned /*kdet*/,
					   const pen_KPAR /*kpar*/,
					   const pen_particleState& state){

  if(state.MAT > 0){
    //Particle created at non void volume. Add particle energy to compensate
    //substracted one when beginPart will be called.
    updateEdepCounters(state.E, nhist, state.X, state.Y, state.Z, state.WGHT);
  }
}

void pen_DICOMDoseDistrib::tally_step(const unsigned long long nhist,
				      const pen_KPAR /*kpar*/,
				      const pen_particleState& state,
				      const tally_StepData& stepData){

  if(stepData.softDE == 0.0){return;}  //Nothing to deposit.
  //Energy deposited
  updateEdepCounters(stepData.softDE, nhist,
		     stepData.softX, stepData.softY, stepData.softZ,
		     state.WGHT);  
}


void pen_DICOMDoseDistrib::tally_move2geo(const unsigned long long nhist,
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

void pen_DICOMDoseDistrib::tally_endHist(const unsigned long long /*nhist*/){

  for(long int i = 0; i < ncontours; ++i)
    {
      // Transfer temp counter
      contEdep[i]  += contEdptmp[i];   
      contEdep2[i] += contEdptmp[i]*contEdptmp[i];
      // Reset counter
      contEdptmp[i]= 0.0;                  
    }  
}

void pen_DICOMDoseDistrib::clear()
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

  if(contEdptmp != nullptr){
    free(contEdptmp);
    contEdptmp = nullptr;
  }
     
  if(contEdep != nullptr){
    free(contEdep);
    contEdep = nullptr;
  }
     
  if(contEdep2 != nullptr){
    free(contEdep2);
    contEdep2 = nullptr;
  }  

  nx = ny = nz = nxy = nbin = 0;
  dx = dy = dz = 0.0;
  idx = idy = idz = 1.0e35;
  xmin = ymin = zmin = 0.0;

  contourVox = nullptr;
  ncontours = 0;
  contNames.clear();
}



int pen_DICOMDoseDistrib::configure(const wrapper_geometry& geometry,
				    const abc_material* const materials[constants::MAXMAT],
				    const pen_parserSection& /*config*/,
				    const unsigned verbose
				      ){
  
  //Clear previous configuration
  clear();

  //**********************
  //* Material densities *
  //**********************

  for(unsigned i = 0; i < constants::MAXMAT; i++){
    if(materials[i] != nullptr)
      matDens[i] = materials[i]->readDens();
  }

  //Check if this geometry is a DICOM based geometry
  const pen_dicomGeo* pDICOMgeo = dynamic_cast<const pen_dicomGeo*>(&geometry);
  if(pDICOMgeo == nullptr){
    if(verbose > 0)
      printf("pen_DICOMDoseDistrib::configure: Error: DICOM dose distribution "
	     "tally requires a DICOM based geometry.\n");
    return 1;
  }

  //Get DICOM
  const pen_dicom& dicom = pDICOMgeo->readDicom();
  nx = dicom.getNX();
  ny = dicom.getNY();
  nz = dicom.getNZ();

  nxy = nx*ny;
  nbin = dicom.getNVox();

  dx = dicom.getDX();
  dy = dicom.getDY();
  dz = dicom.getDZ();

  idx = 1.0/dx;
  idy = 1.0/dy;
  idz = 1.0/dz;

  xmin = ymin = zmin = 0.0;
  
  voxVol = dicom.getVoxVol();

  //Get number of contours
  ncontours = dicom.nContours();
  
  //**************
  //*  Allocate  *
  //**************
    
  //Allocate memory for spatial 3D distribution
  edptmp    = static_cast<double*>(malloc(sizeof(double)*nbin));
  edep      = static_cast<double*>(malloc(sizeof(double)*nbin));
  edep2     = static_cast<double*>(malloc(sizeof(double)*nbin));
  nlast     =
    static_cast<unsigned long long*>(malloc(sizeof(unsigned long long)*nbin));
  ivoxMass  = static_cast<double*>(malloc(sizeof(double)*nbin));

  if(ncontours > 0){
    contEdptmp    = static_cast<double*>(malloc(sizeof(double)*ncontours));
    contEdep      = static_cast<double*>(malloc(sizeof(double)*ncontours));
    contEdep2     = static_cast<double*>(malloc(sizeof(double)*ncontours));

    //Get voxel contour information
    contourVox = dicom.readContour();

    for(long int i = 0; i < ncontours; ++i)
      contEdptmp[i] = 0.0;
    for(long int i = 0; i < ncontours; ++i)
      contEdep[i] = 0.0;
    for(long int i = 0; i < ncontours; ++i)
      contEdep2[i] = 0.0;

    //Save contour names
    for(long int i = 0; i < ncontours; ++i){
      contNames.push_back(dicom.contour(i).name);
    }
  }
  
  for(long int i = 0; i < nbin; ++i)
    edptmp[i] = 0.0;
  for(long int i = 0; i < nbin; ++i)
    edep[i] = 0.0;
  for(long int i = 0; i < nbin; ++i)
    edep2[i] = 0.0;
  for(long int i = 0; i < nbin; ++i)
    nlast[i] = 0;

  //Get voxels mass
  
  for(long int i = 0; i < nbin; ++i){
    const pen_voxel& element = pDICOMgeo->readElement(i);
    if(element.MATER == 0)
      ivoxMass[i] = 0.0;
    else
      ivoxMass[i] = 1.0/(voxVol*matDens[element.MATER-1]*element.densityFact);
  }
    
  if(verbose > 1){
    printf("Number of x bins: %ld\n",nx);
    printf("x Bin width (cm): %12.5E\n",dx);
    printf("Number of y bins: %ld\n",ny);
    printf("y Bin width (cm): %12.5E\n",dy);
    printf("Number of z bins: %ld\n",nz);
    printf("z Bin width (cm): %12.5E\n",dz);
    printf("Number of bins  : %ld\n\n", nbin);

    printf("Number of contours: %ld\n",ncontours);
    if(ncontours > 0){
      printf("Contour names:\n");
      for(const std::string& contName : contNames)
	printf(" %s\n",contName.c_str());
    }
  }
    
  //Register data to dump
  dump.toDump(edep,nbin);
  dump.toDump(edep2,nbin);

  if(ncontours > 0){
    dump.toDump(contEdep,ncontours);
    dump.toDump(contEdep2,ncontours);
  }
    
    
  return 0; 
    
            
}

void pen_DICOMDoseDistrib::tally_endSim(const unsigned long long /*nhist*/){
  flush();
}
    
    
 
void pen_DICOMDoseDistrib::saveData(const unsigned long long nhist) const{
  
  char buffer[81];
  FILE* out;
  double invn = 1.0/static_cast<double>(nhist);

  const double ev2Gy = 1.0/1.60217662E-16;
  
  strcpy(buffer,"dicomDoseDistrib.dat");
     
  out = fopen(buffer,"w");
  if (out == NULL)
    {
      printf("*********************************************\n");
      printf("pen_DICOMDoseDistrib: Error: Cannot open output data file\n");
      printf("*********************************************\n");
      fclose(out);                            
      return;                      
    }
    
    
  // Write header 
  fprintf(out,"#------------------------------------------------------------\n");
  fprintf(out,"# PenRed: DICOM dose distribution report\n");
  fprintf(out,"# Dose units are: Gy per history\n");

  fprintf(out,"#\n");
  fprintf(out,"# No. of bins in x,y,z directions and total:\n");
  fprintf(out,"#  %ld %ld %ld %ld\n",nx,ny,nz,nbin);
  fprintf(out,"# Min values for x,y,z(cm):\n");
  fprintf(out,"#  %12.5E %12.5E %12.5E\n",xmin,ymin,zmin);
  fprintf(out,"# Bin widths for x,y,z(cm):\n");
  fprintf(out,"#  %12.5E %12.5E %12.5E\n",dx,dy,dz);
  fprintf(out,"#\n");
  fprintf(out,"# ");
  fprintf(out,"zBinIndex : ");
  fprintf(out,"yBinIndex : ");
  fprintf(out,"xBinIndex : ");
  fprintf(out,"dose : +-2sigma\n");
    
    
  //Write data    
  for(long int k = 0; k < nz; k++)
    {
      long int kbin = k*nxy;

      for(long int j = 0; j < ny; j++)
        {
	  long int jkbin = j*static_cast<long int>(nx)+kbin;
        
	  for(long int i = 0; i < nx; i++)
            {
	      fprintf(out," %5ld  %5ld  %5ld ",k,j,i);
	      long int bin = i + jkbin; 
                
	      double fact = ivoxMass[bin]*ev2Gy;
	      double q = edep[bin]*invn;
	      double sigma = edep2[bin]*invn - q*q;
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

  if(ncontours == 0)
    return; //No more data to save
  
  //*****************************
  //* CONTOUR DOSE DISTRIBUTION *
  //*****************************

  out = fopen("contourDose.dat","w");

  fprintf(out,"#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  fprintf(out,"# [SECTION REPORT DICOM CONTOUR DOSE]\n");
  fprintf(out,"# Dose units are: Gy per history\n");

  fprintf(out,"#\n");
  fprintf(out,"# No. of contours:\n");
  fprintf(out,"#  %ld\n",ncontours);
  fprintf(out,"#\n");
  fprintf(out,"#\n");
  fprintf(out,"# ");
  fprintf(out,"     contour     : ");
  fprintf(out,"dose : +-2sigma\n");
  
  for(long int icont = 0; icont < ncontours; icont++){

    fprintf(out," %15s ",contNames[icont].c_str());
    
    double q = contEdep[icont]*invn;
    double sigma = contEdep2[icont]*invn - q*q;
    if(sigma > 0.0)
      {
	sigma = sqrt(sigma*invn);
      }
    else
      {
	sigma = 0.0;
      }
    fprintf(out," %12.5E %7.1E\n",q*ev2Gy,2.0*sigma*ev2Gy);
  }

  fclose(out);
}
 
 
int pen_DICOMDoseDistrib::sumTally(const pen_DICOMDoseDistrib& tally){
    
  if(nbin != tally.nbin)
    return -1;
  if(ncontours != tally.ncontours)
    return -2;
  
  for(long int i = 0; i < nbin; ++i){
    edep[i] += tally.edep[i];
  }
  for(long int i = 0; i < nbin; ++i){
    edep2[i] += tally.edep2[i];
  }

  for(long int i = 0; i < ncontours; ++i){
    contEdep[i] += tally.contEdep[i];
  }
  for(long int i = 0; i < ncontours; ++i){
    contEdep2[i] += tally.contEdep2[i];
  }
  
  return 0;
    
    
}

REGISTER_COMMON_TALLY(pen_DICOMDoseDistrib, DICOM_DOSE_DISTRIB)

#endif
 
 
 
 
 
 
 
 
 
