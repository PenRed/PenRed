
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



#include <cstdio>
#include "pen_geometries.hh"

bool diff(const double val1,
	  const double val2,
	  const double tol = 1.0e-14){
  
  if(val1 == 0.0 || val2 == 0.0){
    if(val1 != val2){
      return true;
    }
  }
  else{
    if(fabs(1.0-val1/val2) > tol){
      return true;
    }
  }
  
  return false;
}

int voxGeocmp(const pen_voxelGeo& geo1,
	      const pen_voxelGeo& geo2,
	      const double tol = 1.0e-14){

  printf("  nx = %u\n",geo2.xVox());
  printf("  ny = %u\n",geo2.yVox());
  printf("  nz = %u\n",geo2.zVox());
  printf("nvox = %u\n",geo2.nVox());

  if(geo1.xVox() != geo2.xVox() ||
     geo1.yVox() != geo2.yVox() ||
     geo1.zVox() != geo2.zVox() ||
     geo1.nVox() != geo2.nVox()){
    printf("Error: Number of voxels mismatch.\n");
    printf("  nx = %u\n",geo1.xVox());
    printf("  ny = %u\n",geo1.yVox());
    printf("  nz = %u\n",geo1.zVox());
    printf("nvox = %u\n",geo1.nVox());    
    return -3;
  }

  printf("dx = %12.4E\n",geo2.xSize());
  printf("dy = %12.4E\n",geo2.ySize());
  printf("dz = %12.4E\n",geo2.zSize());

  if(geo1.xSize() != geo2.xSize() ||
     geo1.ySize() != geo2.ySize() ||
     geo1.zSize() != geo2.zSize()){
    printf("Error: Voxel sizes mismatch.\n");
    printf("dx = %12.4E\n",geo1.xSize());
    printf("dy = %12.4E\n",geo1.ySize());
    printf("dz = %12.4E\n",geo1.zSize());
    return -4;
  }

  printf("\n\n");
  unsigned err2 = 0;
  for(unsigned i = 0; i < geo2.nVox(); i++){

    const pen_voxel voxel1 = geo1.readElement(i);
    const pen_voxel voxel2 = geo2.readElement(i);
    if(voxel1.MATER != voxel2.MATER ||
       diff(voxel1.densityFact,voxel2.densityFact, tol)){
      printf("Error: voxel %u mismatch material or density factor:\n",i);
      printf("     voxel mat = %u\n",voxel2.MATER);
      printf("  expected mat = %u\n",voxel1.MATER);
      printf("    voxel dens = %15.6E\n",voxel2.densityFact);
      printf(" expected dens = %15.6E\n",voxel1.densityFact);
      err2++;
    }
  }

  if(err2 > 0){
    return -5;
  }

  return 0;
}

int main(){

  int err;
  
  pen_voxelGeo voxelgeo1;
  pen_voxelGeo voxelgeo2;
  pen_voxelGeo voxelgeo3;
  pen_voxelGeo voxelgeo4;

  pen_rand random;
  random.rand0(2); //Set initial random seeds

  const unsigned nvox[3] = {312,524,208};
  const size_t tvox = nvox[0]*nvox[1]*nvox[2];
  const double sizes[3] = {3.2,0.2,1.6};

  unsigned* voxMats = nullptr;
  double*   voxDens = nullptr;
  
  //Create random voxelized geometry
  voxMats = (unsigned*) malloc(tvox*sizeof(unsigned));
  voxDens = (double*)   malloc(tvox*sizeof(double));

  if(voxMats == nullptr || voxDens == nullptr){
    printf("Error: Unable to allocate memory for voxel data\n");
    return -1;
  }

  //Sample random values for materials and density factors
  for(size_t i = 0; i < tvox; i++){
    //Generate random integer [1-16]    
    voxMats[i] = unsigned(random.rand()*15.5)+1;

    //Generate random density factor
    voxDens[i] = random.rand()*1.0e6+1.0e-3;
  }

  //******************************
  // Initialize voxel geometry 1
  //******************************
  err = voxelgeo1.setVoxels(nvox,sizes,voxMats,voxDens,3);
  if(err != 0){
    printf("Error using 'setVoxels': %d\n",err);
    free(voxMats);
    free(voxDens);
    return -2;
  }

  printf("\nGeometry 1 initialized from data.\n\n");
  
  //**************************
  // Check initialized data
  //**************************

  printf("  nx = %u\n",voxelgeo1.xVox());
  printf("  ny = %u\n",voxelgeo1.yVox());
  printf("  nz = %u\n",voxelgeo1.zVox());
  printf("nvox = %u\n",voxelgeo1.nVox());

  if(nvox[0] != voxelgeo1.xVox() ||
     nvox[1] != voxelgeo1.yVox() ||
     nvox[2] != voxelgeo1.zVox() ||
     tvox    != voxelgeo1.nVox()){
    printf("Error: Initialized number of voxels mismatch.\n");
    printf("  nx = %u\n",nvox[0]);
    printf("  ny = %u\n",nvox[1]);
    printf("  nz = %u\n",nvox[2]);
    printf("nvox = %lu\n",tvox);
    free(voxMats);
    free(voxDens);
    return -3;
  }

  printf("dx = %12.4E\n",voxelgeo1.xSize());
  printf("dy = %12.4E\n",voxelgeo1.ySize());
  printf("dz = %12.4E\n",voxelgeo1.zSize());

  if(sizes[0] != voxelgeo1.xSize() ||
     sizes[1] != voxelgeo1.ySize() ||
     sizes[2] != voxelgeo1.zSize()){
    printf("Error: Initialized voxel sizes mismatch.\n");
    printf("dx = %12.4E\n",sizes[0]);
    printf("dy = %12.4E\n",sizes[1]);
    printf("dz = %12.4E\n",sizes[2]);
    free(voxMats);
    free(voxDens);
    return -4;
  }

  printf("\n\n");
  unsigned err2 = 0;
  for(unsigned i = 0; i < tvox; i++){

    const pen_voxel voxel = voxelgeo1.readElement(i);
    if(voxel.MATER != voxMats[i] || voxel.densityFact != voxDens[i]){
      printf("Error: voxel %u mismatch material or density factor:\n",i);
      printf("     voxel mat = %u\n",voxel.MATER);
      printf("  expected mat = %u\n",voxMats[i]);
      printf("    voxel dens = %15.6E\n",voxel.densityFact);
      printf(" expected dens = %15.6E\n",voxDens[i]);
      err2++;
    }
  }

  if(err2 > 0){
    return -5;
  }

  //Write voxel data in ASCII format
  printf("Write voxel data in ASCII format\n");
  FILE* fascii = fopen("asciiVox.geo", "w");  
  fprintf(fascii, "%u %u %u\n", nvox[0], nvox[1], nvox[2]);
  fprintf(fascii, "%E %E %E\n", sizes[0], sizes[1], sizes[2]);
  for(size_t i = 0; i < tvox; ++i){
    fprintf(fascii, "%u %15.6E\n", voxMats[i], voxDens[i]);
  }
  fclose(fascii);
  printf("Voxel ASCII file completed, load it.\n");  

  free(voxMats);
  free(voxDens);
  
  printf("\nLoaded data in geometry 1 cheked.\n");
  
  //*********************************
  // Dump voxel data from geometry 1
  //*********************************

  unsigned char* pdata;
  size_t dataSize;
  err = voxelgeo1.dump(pdata,dataSize,3);
  if(err != 0){
    printf("Error dumping voxel data: %d\n",err);
    return -6;
  }

  printf("\nVoxel data dumped (%lu B).\n",dataSize);

  //*********************************
  // Read voxel data to geometry 2
  //*********************************

  size_t pos = 0;
  err = voxelgeo2.loadData(pdata,pos,3);
  if(err != 0){
    printf("Error reading dumped data to geometry 2: %d\n",err);
    free(pdata);
    return -7;
  }

  if(pos != dataSize){
    printf("Error: read and dumped data size mismatch!");
    printf("        dumped: %lu\n",dataSize);
    printf("          read: %lu\n",pos);
    free(pdata);
    return -8;
  }

  printf("Compare loaded data in geometry 2 with data contained in geometry 1:\n");
  if(voxGeocmp(voxelgeo1,voxelgeo2) != 0){
    printf("Loaded voxel data at geometry 1 and 2 misatch.\n");
    return -9;
  }

  printf("Voxel data in Geometry 1 and 2 matches.\n");
  
  //*****************************************
  // Read voxel data to geometry 3 from file
  //*****************************************

  FILE* fdump = nullptr;
  fdump = fopen("dumpVox.dump","wb");
  if(fdump == nullptr){
    printf("Unable to open file 'dumpVox.dump'\n");
    free(pdata);
    return -10;
  }

  size_t written = fwrite(pdata,1,dataSize,fdump);
  free(pdata);
  fclose(fdump);
  if(written != dataSize){
    printf("Unable to write all data:\n");
    printf("      written: %lu/%lu\n",written,dataSize);
    return -11;
  }

  err = voxelgeo3.loadFile("dumpVox.dump",3);
  if(err != 0){
    printf("Error loading voxel geometry from file: %d\n",err);
    return -12;
  }

  printf("Compare loaded data in geometry 3 with data contained in geometry 1:\n");
  if(voxGeocmp(voxelgeo1,voxelgeo3) != 0){
    printf("Loaded voxel data at geometry 1 and 3 misatch.\n");
    return -13;
  }

  printf("Voxel data in Geometry 1 and 3 matches.\n");


  printf("Load voxel data in ASCII format in geometry 4.\n");

  err = voxelgeo4.loadASCII("asciiVox.geo", 3);
  if(err != 0){
    printf("Error loading voxel geometry from file: %d\n",err);
    return -14;
  }

  printf("Compare loaded data in geometry 4 with data contained in geometry 1:\n");
  if(voxGeocmp(voxelgeo1,voxelgeo4, 1.0e-4) != 0){
    printf("Loaded voxel data at geometry 1 and 4 misatch.\n");
    return -15;
  }
  
  printf("Voxel data in Geometry 1 and 4 matches.\n");
  
  printf("\n\nTest pased!\n\n");
  
  return 0;
}
