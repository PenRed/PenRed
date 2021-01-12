 
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
//
//    contact emails:
//
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//
 
#include "pen_geometries.hh"

int main(int argc, char** argv){

  if(argc < 3){
    printf("usage: %s configuration-file kernelSize material\n",argv[0]);
    return 0;
  }

  unsigned verbose = 3;
  int err;
  
  //Get kernel size and material ID
  const int ikernelSize = atoi(argv[2]);
  const int imaterial = atoi(argv[3]);
  
  if(ikernelSize < 1){
      printf("Error: Kernel size must be 1 or greater. Specified kernel size: %d\n",ikernelSize);
      return -1;
  }
  if(imaterial < 1){
      printf("Error: Invalid material index: %d\n",imaterial);
      return -2;
  }
  
  const unsigned material = static_cast<unsigned>(imaterial);
  const unsigned kernelSize = static_cast<unsigned>(ikernelSize);
  const unsigned kernel05 = kernelSize/2;
  
  //**************************
  // Parse configuration
  //**************************

  //Parse configuration file
  pen_parserSection config;
  std::string errorLine;
  unsigned long errorLineNum;  
  err = parseFile(argv[1],config,errorLine,errorLineNum);

  printf("Configuration:\n");
  printf("%s\n", config.stringify().c_str());
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }

  //**************************
  // Configure geometry
  //**************************  
  
  //Create dicom voxelized geometry
  pen_dicomGeo geometry;
  
  err = geometry.configure(config,verbose);
  if(err != 0){
    printf("Error at geometry configuration.\n");
    printf("                 Error code: %d\n",err);
    return -2;
  }

  const pen_voxel* mesh = geometry.readMesh();
  const unsigned long nvox = geometry.getElements();
  const unsigned long nx = geometry.xVox();
  const unsigned long ny = geometry.yVox();
  const unsigned long nz = geometry.zVox();
  
  const unsigned long nxy = nx*ny;
  
  //Create and init material and density arrays
  unsigned nkvox[3] = {kernelSize,kernelSize,kernelSize};
  double dsvox[3] = {geometry.xSize(),geometry.ySize(),geometry.zSize()};
  unsigned long kernelVox = kernelSize*kernelSize*kernelSize;
  unsigned* mats = static_cast<unsigned*>(malloc(sizeof(unsigned)*kernelVox));
  double* densFact = static_cast<double*>(malloc(sizeof(double)*kernelVox));
  
  for(unsigned imat = 0; imat < kernelVox; ++imat)
      mats[imat] = 0;
  for(unsigned ifact = 0; ifact < kernelVox; ++ifact)
      densFact[ifact] = 1.0;

  //Create global kernel information file
  FILE* fout = fopen("kernel-info.dat","w");
  if(fout == nullptr){
      printf("Error: Unable to open kernel information file.\n");
      return -3;
  }
  
  //Extract kernels
  size_t nkernels = 0;
  for(unsigned long ivox = 0; ivox < nvox; ++ivox){
      //Check voxel material
      if(mesh[ivox].MATER == material){
          //Get index i,j,k
          unsigned long k = ivox/nxy;
          unsigned long planeVox = ivox%nxy;
          
          unsigned long j = planeVox/nx;
          unsigned long i = planeVox%nx;
          
          //Check if the kernel fit inside the image
          if(kernel05 > i || kernel05 > j || kernel05 > k)
              continue;
          if(i+kernel05 >= nx || j+kernel05 >= ny || k+kernel05 >= nz)
              continue;
          
          unsigned long iinit = i - kernel05;
          unsigned long jinit = j - kernel05;
          unsigned long kinit = k - kernel05;
          
          size_t index = 0;
          for(unsigned long kk = 0; kk < kernelSize; ++kk){
              unsigned long kkvox = kinit + k;
              for(unsigned long jj = 0; jj < kernelSize; ++jj){
                  unsigned long jjvox = kkvox*nxy + (jinit+jj)*nx;
                  for(unsigned long ii = 0; ii < kernelSize; ++ii){
                      unsigned long iivox = jjvox + iinit + ii;
                      mats[index] = mesh[iivox].MATER;
                      densFact[index] = mesh[iivox].densityFact;
                      ++index;
                  }
              }
          }
          pen_voxelGeo kernelGeo;
          kernelGeo.setVoxels(nkvox,dsvox,mats,densFact,verbose);
          char filename[100];
          sprintf(filename,"kernel-%lu-%ux%ux%u.dump",
                  nkernels,kernelSize,kernelSize,kernelSize);
          kernelGeo.dump2File(filename);
          
          //Print kernel information
          fprintf(fout,"%lu ",nkernels);
          for(unsigned long ik = 0; ik < kernelVox; ++ik){
            fprintf(fout,"%u %le ",mats[ik],densFact[ik]);
          }
          fprintf(fout,"\n");
          ++nkernels;
      }
  }
  
  fclose(fout);
  
  return 0;
}
