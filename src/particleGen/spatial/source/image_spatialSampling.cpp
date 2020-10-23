
//
//
//    Copyright (C) 2020 Universitat de València - UV
//    Copyright (C) 2020 Universitat Politècnica de València - UPV
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

#include "image_spatialSampling.hh"

void image_spatialSampling::geoSampling(double pos[3], pen_rand& random) const{

  //Obtain voxel index randomly using Walker's aliasing algorithm
  long int ivox = IRND(F,K,nvox,random);

  long int ix = ivox % nx;
  long int iy = (ivox % nxy)/nx;
  long int iz = ivox / nxy;
  
  pos[0] = dx*(static_cast<double>(ix) + (1.0-2.0*random.rand())*0.5)-imageCx;
  pos[1] = dy*(static_cast<double>(iy) + (1.0-2.0*random.rand())*0.5)-imageCy;
  pos[2] = dz*(static_cast<double>(iz) + (1.0-2.0*random.rand())*0.5)-imageCz;

  /* printf("tallyimage_spatialSampling: voxel index   :   %ld\n",ivox);
  printf("tallyimage_spatialSampling: source chosen :  (%3ld,%3ld,%3ld)\n",ix,iy,iz);
  printf("tallyimage_spatialSampling: center (image source coordinate system):  (%9.3lf,%9.3lf,%9.3lf)\n",
    dx*static_cast<double>(ix)+Ox-dx/2-isocenter[0],
    dy*static_cast<double>(iy)+Oy-dy/2-isocenter[1],
    dz*static_cast<double>(iz)+Oz-dz/2-isocenter[2]);
  printf("tallyimage_spatialSampling: Position (phantom coordinate system): (%9.3lf,%9.3lf,%9.3lf)\n",
    pos[0]+translation[0],pos[1]+translation[1],pos[2]+translation[2]);fflush(stdout); */
  
}

int image_spatialSampling::configure(const pen_parserSection& config,
				     const unsigned verbose){


  int err = 0;

  // Read DICOM directory path
  //****************************
  std::string directoryPath;
  if(config.read("directory",directoryPath) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("image_spatialSampling:configure:Error: Unable to read field "
	     "'directory'. String spected.\n");
    }
    return 1;
  }

  err = config.read("isocenter/x",isocenter[0]);
  if(err != INTDATA_SUCCESS){
    isocenter[0] = 0.0;
  }

  err = config.read("isocenter/y",isocenter[1]);
  if(err != INTDATA_SUCCESS){
    isocenter[1] = 0.0;
  }

  err = config.read("isocenter/z",isocenter[2]);
  if(err != INTDATA_SUCCESS){
    isocenter[2] = 0.0;    
  }

  err = config.read("position/x",translation[0]);
  if(err != INTDATA_SUCCESS){
    translation[0] = 0.0;
  }

  err = config.read("position/y",translation[1]);
  if(err != INTDATA_SUCCESS){
    translation[1] = 0.0;
  }

  err = config.read("position/z",translation[2]);
  if(err != INTDATA_SUCCESS){
    translation[2] = 0.0;    
  }

  bool toRotate = false;
  double omega, theta, phi;
  err = config.read("euler/omega",omega);
  if(err != INTDATA_SUCCESS){
    omega = 0.0;
  }else toRotate = true;

  err = config.read("euler/theta",theta);
  if(err != INTDATA_SUCCESS){
    theta = 0.0;
  }else toRotate = true;

  err = config.read("euler/phi",phi);
  if(err != INTDATA_SUCCESS){
    phi = 0.0;
  }else toRotate = true;

  //Create the rotation
  if(toRotate)
    setRotationZYZ(omega,theta,phi);

  // Try to load the DICOM
  //************************
  pen_dicom dicom;
  err = dicom.loadDicom(directoryPath.c_str(),verbose);
  if(err != PEN_DICOM_SUCCESS){
    if(verbose > 0){
      printf("image_spatialSampling:configure: "
	     "Error loading DICOM '%s'\n",directoryPath.c_str());
      printf("                 Error code: %d\n",err);
    }
    return 2;
  }  

  //Store image dimensions
  nx = dicom.getNX();
  ny = dicom.getNY();
  nz = dicom.getNZ();
  nxy = nx*ny;
  nvox = nxy*nz;

  dx = dicom.getDX();
  dy = dicom.getDY();
  dz = dicom.getDZ();

  imageCx = static_cast<double>(nx)*dx*0.5;
  imageCy = static_cast<double>(ny)*dy*0.5;
  imageCz = static_cast<double>(nz)*dz*0.5;

  Ox = dicom.getOriginX();
  Oy = dicom.getOriginY();
  Oz = dicom.getOriginZ();

  //Get raw image information
  const double* image = dicom.readImage();

  F = static_cast<double*>(malloc(sizeof(double)*nvox));
  K = static_cast<long int*>(malloc(sizeof(long int)*nvox));

  if(F == nullptr){
    if(verbose > 0)
      printf("image_spatialSampling:configure: "
	     "Unable to allocate memory for 'F' array\n");
    return 3;
  }

  if(K == nullptr){
    if(verbose > 0)
      printf("image_spatialSampling:configure: "
	     "Unable to allocate memory for 'K' array\n");
    return 4;
  }

  //Init walker algorithm
  printf("Initializing Walker algorithm ... ");fflush(stdout);
  IRND0(image,F,K,nvox);
  printf("done\n");fflush(stdout);

  char buffer[81];
  char BLINE[100];  
  /* FILE* in;
  
  strcpy(buffer,"WalkersCutoffAlias.dat");
     
  in = fopen(buffer,"r");
  if (in == NULL)
    {
      printf("*********************************************\n");
      printf("image_spatialSampling: Error: Cannot open input data file\n");
      printf("*********************************************\n");
      fclose(in);                            
      return 5;                      
    }

  for(int i=1; i<=8; i++)
    {
      BLINE[0]='\0';
      fscanf(in,"%[^\r\n]%*[^\n]",BLINE);
      if(fgetc(in)=='\r'){fgetc(in);}
    }
  for(long int i = 0; i < nvox; i++)
    {
      //int Num_Elements_Llegits = sscanf(BLINE, "%8ld  %16.8E", &K[i], &F[i]);
      //fprintf(out," %8ld  %16.8E\n",K[i],F[i]);
      //fscanf(in,"%ld %lf%*[^\n]", &K[i], &F[i]);getc(in);
      fscanf(in,"%ld %lf", &K[i], &F[i]);
    }

  fclose(in);
  in = nullptr; */

  FILE* out;
  
  strcpy(buffer,"WalkersCutoffAlias_test.dat");
     
  out = fopen(buffer,"w");
  if (out == NULL)
    {
      printf("*********************************************\n");
      printf("image_spatialSampling: Error: Cannot open output data file\n");
      printf("*********************************************\n");
      fclose(out);                            
      return 6;                      
    }
    
    
  // Write header 
  fprintf(out,"#------------------------------------------------------------\n");
  fprintf(out,"# PenRed: Walker aliasing algorithm\n");
  fprintf(out,"# Cutoff and alias values\n");

  fprintf(out,"#\n");
  fprintf(out,"# No. of voxels in x,y,z directions and total:\n");
  fprintf(out,"#  %ld %ld %ld %ld\n",nx,ny,nz,nvox);
  fprintf(out,"#\n");
  fprintf(out,"# ");
  fprintf(out,"Alias    : ");
  fprintf(out,"Cutoff   \n");
        
  //Write data    
  for(long int i = 0; i < nvox; i++)
    {
      fprintf(out," %8ld  %16.8E\n",K[i],F[i]);
    }

  fclose(out);
  out = nullptr;
  
  if(verbose > 1){
    printf("Image Origin (Ox,Oy,Oz):\n %12.4E %12.4E %12.4E\n",
     Ox,Oy,Oz);fflush(stdout);
    printf("Image Isocenter (Isox,Isoy,Isoz):\n %12.4E %12.4E %12.4E\n",
     isocenter[0],isocenter[1],isocenter[2]);fflush(stdout);
    printf("Image center (x,y,z):\n %12.4E %12.4E %12.4E\n",
	   translation[0],translation[1],translation[2]);fflush(stdout);
    printf("Image voxels (nx,ny,nz):\n %ld %ld %ld\n",
	   nx,ny,nz);fflush(stdout);
    printf("Voxel sizes (dx,dy,dz):\n %12.4E %12.4E %12.4E\n",
	   dx,dy,dz);fflush(stdout);
    if(toRotate)
      printf("Image rotation (omega,theta,phi):\n %12.4E %12.4E %12.4E\n",
	     omega,theta,phi);
    else
      printf("No rotation applied.\n");
  }
  
  return 0;
}

image_spatialSampling::~image_spatialSampling(){
  if(F != nullptr)
    free(F);
  if(K != nullptr)
    free(K);
  F = nullptr;
  K = nullptr;
}

REGISTER_SAMPLER(image_spatialSampling,IMAGE)

#endif
