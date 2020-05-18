
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


#include "voxel_geo.hh"


const double pen_voxelGeo::pZero =  1.0e-17;
const double pen_voxelGeo::nZero = -1.0e-17;
const double pen_voxelGeo::inf = 1.0e35;


pen_voxelGeo::pen_voxelGeo() : nx(0), ny(0), nz(0), nxy(0),
			       dx(0.0), dy(0.0), dz(0.0),
			       idx(inf), idy(inf), idz(inf),
			       Mdx(0.0), Mdy(0.0), Mdz(0.0)
{}


int pen_voxelGeo::configure(const pen_parserSection& config,
			     const unsigned verbose){

  int err = 0;

  // Read voxels filename
  //////////////////////////
  std::string filename;
  if(config.read("filename",filename) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_voxelGeo:configure:Error: Unable to read field 'filename'. String spected.\n");
    }
    err++;
  }

  int err2 = loadFile(filename.c_str(),verbose);  
  if(err2 != 0){
    if(verbose > 0){
      printf("pen_voxelGeo:configure: Error reading and loading file '%s': %d\n", filename.c_str(),err2);
    }
    err++;
    configStatus = err;
    return err;
  }
  
  // Read number of voxels
  //////////////////////////
  int aux;
  if(config.read("nvoxels/nx",aux) == INTDATA_SUCCESS){
    if(aux != (int)nx){
      if(verbose > 0){
	printf("pen_voxelGeo:configure:Error: Read value 'nvoxels/nx' mismatch with data stored in %s.\n",filename.c_str());
	printf("                       %d != %u\n",aux,nx);
      }
      err++;
    }
    else if(verbose > 1)
      printf("Number of voxels in x axis match!\n");
  }

  if(config.read("nvoxels/ny",aux) == INTDATA_SUCCESS){
    if(aux != (int)ny){
      if(verbose > 0){
	printf("pen_voxelGeo:configure:Error: Read value 'nvoxels/ny' mismatch with data stored in %s.\n",filename.c_str());
	printf("                       %d != %u\n",aux,ny);
      }
      err++;
    }
    else if(verbose > 1)
      printf("Number of voxels in y axis match!\n");
  }
  
  if(config.read("nvoxels/nz",aux) == INTDATA_SUCCESS){
    if(aux != (int)nz){
      if(verbose > 0){
	printf("pen_voxelGeo:configure:Error: Read value 'nvoxels/nz' mismatch with data stored in %s.\n",filename.c_str());
	printf("                       %d != %u\n",aux,nz);
      }
      err++;
    }
    else if(verbose > 1)
      printf("Number of voxels in z axis match!\n");
  }
  
  // Read voxels size
  //////////////////////////
  double auxd;
  if(config.read("voxel-size/dx",auxd) == INTDATA_SUCCESS){
    if(auxd != dx){
      if(verbose > 0){
	printf("pen_voxelGeo:configure:Error: Read value 'voxe-size/dx' mismatch with data stored in %s.\n",filename.c_str());
	printf("                       %12.4E != %12.4E\n",auxd,dx);
      }
      err++;
    }
    else if(verbose > 1)
      printf("Voxels size in x axis match!\n");
  }

  if(config.read("voxel-size/dy",auxd) == INTDATA_SUCCESS){
    if(auxd != dy){
      if(verbose > 0){
	printf("pen_voxelGeo:configure:Error: Read value 'voxe-size/dy' mismatch with data stored in %s.\n",filename.c_str());
	printf("                       %12.4E != %12.4E\n",auxd,dy);
      }
      err++;
    }
    else if(verbose > 1)
      printf("Voxels size in y axis match!\n");
  }
  
  if(config.read("voxel-size/dz",auxd) == INTDATA_SUCCESS){
    if(auxd != dz){
      if(verbose > 0){
	printf("pen_voxelGeo:configure:Error: Read value 'voxe-size/dz' mismatch with data stored in %s.\n",filename.c_str());
	printf("                       %12.4E != %12.4E\n",auxd,dz);
      }
      err++;
    }
    else if(verbose > 1)
      printf("Voxels size in z axis match!\n");
  }

  // Read DSMAX
  //////////////////////////

  if(config.isSection("dsmax")){
    //Get bodies alias
    std::vector<std::string> bodyAliasVect;
    config.ls("dsmax",bodyAliasVect);

    if(bodyAliasVect.size() > 0){
      for(unsigned iname = 0; iname < bodyAliasVect.size(); iname++){
	unsigned bodyIndex = getIBody(bodyAliasVect[iname].c_str());
	if(bodyIndex > nBodies){
	  if(verbose > 0){
	    printf("pen_voxelGeo:configure:Error: Body alias '%s' not found. Can't set DSMAX.\n",bodyAliasVect[iname].c_str());
	  }
	  err++;
	}
	else{
	  //Read dsmax
	  std::string key = "dsmax/" + bodyAliasVect[iname];
	  double auxDSmax; 
	  if(config.read(key,auxDSmax) != INTDATA_SUCCESS){
	    if(verbose > 0){
	      printf("pen_voxelGeo:configure:Error: Unable to read DSMAX for body alias '%s'. Double expected.\n",bodyAliasVect[iname].c_str());
	    }
	    err++;
	  }
	  else{
	    if(auxDSmax > 0.0){
	      DSMAX[bodyIndex] = auxDSmax;
	    }
	    else{
	      if(verbose > 0){
		printf("pen_voxelGeo:configure:Error: Invalid value of DSMAX for body alias '%s'.\n",bodyAliasVect[iname].c_str());
	      }
	      err++;
	    }
	  }
	}
      }
    }
  }

  // Check print image option
  //***************************
  bool toASCII = true;
  if(config.read("print-ASCII",toASCII) != INTDATA_SUCCESS){
    toASCII = false;
  }

  if(toASCII){
    printImage("voxelsASCII.rep");    
  }
  
  //Print report
  if(verbose > 1){
    printf("Number of voxels in the mesh (x,y,z):\n");
    printf(" %u %u %u\n",nx,ny,nz);
    
    printf("Voxels sizes (dx,dy,dz):\n");
    printf(" %12.4E %12.4E %12.4E\n",dx,dy,dz);

    printf("Mesh dimensions (Dx,Dy,Dz):\n");
    printf(" %12.4E %12.4E %12.4E\n\n",Mdx,Mdy,Mdz);

    printf("Printing to ASCII %s\n",toASCII ? "enabled" : "disabled");
    
    for(unsigned i = 0; i < nBodies; i++){
      printf("Body %u:\n",i);
      printf("   DSMAX: %12.5E\n",DSMAX[i]);
      printf("    KDET: %u\n",KDET[i]);
      printf("-------------------------------\n\n");
    }
  }
  
  configStatus = err;
  return err;
}


void pen_voxelGeo::locate(pen_particleState& state) const{

  long int ix, iy, iz;

  ix = state.X/dx;
  iy = state.Y/dy;
  iz = state.Z/dz;

  if(ix < 0 || (long unsigned)ix >= nx ||
     iy < 0 || (long unsigned)iy >= ny ||
     iz < 0 || (long unsigned)iz >= nz){
    //Particle scapes from geometry mesh
    state.IBODY = constants::MAXMAT;
    state.MAT = 0;
  }
  else{
    //Particle is into geometry mesh, calculate vox
    //index and its material
    long int index = iz*nxy + iy*nx + ix;
    state.MAT = mesh[index].MATER;
    state.IBODY = state.MAT-1;
  }
  
}

void pen_voxelGeo::step(pen_particleState& state,
			double DS,
			double &DSEF,
			double &DSTOT,
			int &NCROSS) const{
  
  //Check if particle is out of geometry mesh
  if(state.MAT == 0){
    
    //Check if the particle reaches the mesh in
    //each axis
    double dsxIn, dsyIn, dszIn;
    double dsxOut, dsyOut, dszOut;
    if(!moveIn(state.X, state.U, dsxIn, dsxOut, Mdx) ||
       !moveIn(state.Y, state.V, dsyIn, dsyOut, Mdy) ||
       !moveIn(state.Z, state.W, dszIn, dszOut, Mdz)){
    
      state.IBODY = constants::MAXMAT;
      state.MAT = 0;
      DSEF = 1.0e35;
      DSTOT = 1.0e35;
      NCROSS = 0; //No interface crossed

      move(inf,state);
      return;
    }

    //Move the particle the required distance to reach
    //geometry mesh on all axis.
    double ds2allIn = std::max(std::max(dsxIn,dsyIn),dszIn);
    //Check if the particle go out of the mesh when this distance is traveled
    if(ds2allIn > 0.0){
      double ds2someOut = std::min(std::min(dsxOut,dsyOut),dszOut);
      if(ds2allIn >= ds2someOut){
	//The particle doesn't reaches the mesh
	state.IBODY = constants::MAXMAT;
	state.MAT = 0;
	DSEF = 1.0e35;
	DSTOT = 1.0e35;
	NCROSS = 0; //No interface crossed

	move(inf,state);	
	return;
      }
    }    
    move(ds2allIn,state);

    //Get voxel index and material 
    locate(state);
    if(state.MAT == 0){
      DSEF = 1.0e35;
      DSTOT = 1.0e35;
      NCROSS = 0; //No interface crossed

      move(inf,state);
      return;      
    }
    
    NCROSS = 1; //Interface crossed
    DSEF = ds2allIn; //Distance traveled in void
    DSTOT = ds2allIn; //Total distance traveled
    
    return;
  }
  
  //Get voxel index for x, y and z axis
  long int ix, iy, iz;
  unsigned long ivox;

  ix = state.X/dx;
  iy = state.Y/dy;
  iz = state.Z/dz;

  ivox = iz*nxy + iy*nx + ix;
  
  //Particle direction cosinus inverse
  double invU, invV, invW;
  //Required distance to cross a voxel in each direction
  double vox_dx, vox_dy, vox_dz;
  //Distance until next voxel wall
  double ds_x, ds_y, ds_z;
  //Voxel index increment
  int voxInc_x, voxInc_y, voxInc_z;
  int voxIncGlob_x, voxIncGlob_y, voxIncGlob_z;
  //Get initial material from actual voxel
  unsigned imat = state.MAT;
  
  //Calculate the distance to travel to cross a voxel in each direction
  if(state.U != 0.0E0){
    invU = 1.0/state.U;
    vox_dx = fabs(dx*invU);  //Distance to travel to cross a voxel on X axis
  }
  else{
    invU = inf;
    vox_dx = inf;
  }
  
  if(state.V != 0.0E0){
    invV = 1.0/state.V;
    vox_dy = fabs(dy*invV);  //Distance to travel to cross a voxel on Y axis
  }
  else{
    invV = inf;
    vox_dy = inf;
  }

  if(state.W != 0.0E0){
    invW = 1.0/state.W;
    vox_dz = fabs(dz*invW);  //Distance to travel to cross a voxel on Z axis
  }
  else{
    invW = inf;
    vox_dz = inf;
  }

  // Calculate the distance to the nex three voxel walls
  if(std::signbit(invU)){
    // Moves backward on X
    ds_x = (ix*dx-state.X)*invU; //Notice the sign of invU
    voxInc_x = -1;
    voxIncGlob_x = -1;    
  }
  else{
    //Moves forward on X
    ds_x = ((ix+1)*dx-state.X)*invU;
    voxInc_x = +1;
    voxIncGlob_x = +1;    
  }

  if(std::signbit(invV)){
    // Moves backward on Y
    ds_y = (iy*dy-state.Y)*invV;
    voxInc_y = -1;    
    voxIncGlob_y = -nx;    
  }
  else{
    //Moves forward on Y
    ds_y = ((iy+1)*dy-state.Y)*invV;
    voxInc_y = +1;
    voxIncGlob_y = +nx;    
  }

  if(std::signbit(invW)){
    // Moves backward on Z
    ds_z = (iz*dz-state.Z)*invW;
    voxInc_z = -1;    
    voxIncGlob_z = -nxy;    
  }
  else{
    //Moves forward on Z
    ds_z = ((iz+1)*dz-state.Z)*invW;
    voxInc_z = +1;
    voxIncGlob_z = +nxy;    
  }
  
  NCROSS = 0;
  DSEF = 0.0;
  DSTOT = 0.0;
  for(;;){
    if(ds_x < ds_y && ds_x < ds_z){

      ds_x = std::max(ds_x,0.0);
      if(crossVox(ds_x,imat,nx,voxInc_x,voxIncGlob_x,
		  DS,ivox,ix,DSEF,DSTOT,NCROSS,state)){
	//Particle cross an interface or consumed the step
	return;
      }

      //Update distances until next walls
      ds_y -= ds_x;
      ds_z -= ds_x;
      ds_x = vox_dx;      
    }
    else if(ds_y < ds_z){

      ds_y = std::max(ds_y,0.0);      
      if(crossVox(ds_y,imat,ny,voxInc_y,voxIncGlob_y,
		  DS,ivox,iy,DSEF,DSTOT,NCROSS,state)){
	//Particle cross an interface or consumed the step
	return;
      }
      //Update distances until next walls
      ds_x -= ds_y;            
      ds_z -= ds_y;
      ds_y = vox_dy;      
    }
    else{

      ds_z = std::max(ds_z,0.0);      
      if(crossVox(ds_z,imat,nz,voxInc_z,voxIncGlob_z,
		  DS,ivox,iz,DSEF,DSTOT,NCROSS,state)){
	//Particle cross an interface or consumed the step
	return;
      }
      
      //Update distances until next walls
      ds_x -= ds_z;            
      ds_y -= ds_z;            
      ds_z = vox_dz;
    }
  }
}


int pen_voxelGeo::setVoxels(const unsigned nvox[3],
			    const double size[3],
			    const unsigned* mats,
			    const double* dens,
			    const unsigned verbose){

  unsigned err = 0;

  //Check number of voxels
  if(nvox[0] <= 0 || nvox[1] <= 0 || nvox[2] <= 0){
    if(verbose > 0){
      printf("pen_voxelGeo:setVoxels:Error: Number of voxels in each axis must be greater than 0.\n");
    }
    err++;
  }

  //Check voxel size
  if(size[0] <= 0.0 || size[1] <= 0.0 || size[2] <= 0.0){
    if(verbose > 0){
      printf("pen_voxelGeo:setVoxels:Error: Voxels size in each axis must be greater than 0.\n");
    }
    err++;
  }

  if(err > 0) return err;

  //Resize mesh
  long unsigned ntvox = 
  static_cast<long unsigned>(nvox[0])*static_cast<long unsigned>(nvox[1])*
  static_cast<long unsigned>(nvox[2]);
  resizeMesh(ntvox);
  if(getStatus() != PEN_MESH_INITIALIZED){
    if(verbose > 0){
      printf("pen_voxelGeo:setVoxels:Error: Unable to initialize the mesh with %u voxels.\n",nvox[0]*nvox[1]*nvox[2]);
    }
    err++;
    return err; 
  }

  //Store voxel number and size
  nx = nvox[0]; ny = nvox[1]; nz = nvox[2];
  nxy = static_cast<unsigned long>(nx)*static_cast<unsigned long>(ny);
  nElements = nxy*static_cast<unsigned long>(nz);
    
  dx = size[0]; dy = size[1]; dz = size[2];

  //Calculate inverses of voxel sizes
  idx = 1.0/dx;
  idy = 1.0/dy;
  idz = 1.0/dz;

  //Calculate mesh size
  Mdx = static_cast<double>(nx)*dx;
  Mdy = static_cast<double>(ny)*dy;
  Mdz = static_cast<double>(nz)*dz;
  
  
  //Fill the mesh
  for(unsigned long i = 0; i < getElements(); i++){
    if(mats[i] == 0 || mats[i] >= constants::MAXMAT|| dens[i] <= 0.0){
      if(verbose > 0){
	unsigned ix, iy, iz;
	iz = i / nxy;
	iy = (i-iz*nxy)/nx;
	ix = i % nx;	
	printf("pen_voxelGeo:setVoxels:Error: Voxels can't be filled with Void nor null density material: voxel %lu (%u,%u,%u).\n",i, ix, iy, iz);
	printf("                             mat: %u\n",mats[i]);
	printf("        density factor (vox/mat): %12.4E\n",dens[i]);
      }
      err++;
    }
    else{
      //The number of bodies is equal to
      //the greater material index
      if(nBodies < mats[i]){
	nBodies = mats[i];
      }
      mesh[i].MATER = mats[i];
      mesh[i].densityFact = dens[i];
    }
  }
  
  return err;
  
}

unsigned pen_voxelGeo::getIBody(const char* bname) const {

  int index = atoi(bname);
  if(index < 1 || index > (int)nBodies)
    return constants::MAXMAT;
  else
    return unsigned(index-1);
}

int pen_voxelGeo::loadData(const unsigned char* data,
			   size_t& pos,
			   const unsigned verbose){

  int err;
  
  //Read number of voxels
  uint32_t nvox[3];
  unsigned totalVox = 0;
  
  memcpy(nvox,&data[pos],3*sizeof(uint32_t));
  pos += 3*sizeof(uint32_t);

  totalVox = nvox[0]*nvox[1]*nvox[2];
  if(nvox[0] < 1 || nvox[1] < 1 || nvox[2] < 1){
    if(verbose > 0){
      printf("pen_voxelGeo:loadFile: Error: Number of voxels must be greater than 0 on each axis.\n");
      printf("                             nx: %u\n",nvox[0]);
      printf("                             ny: %u\n",nvox[1]);
      printf("                             nz: %u\n",nvox[2]);
    }
    return -1;
  }
  
  //Allocate memory to extract data
  double sizes[3];
  unsigned* materials = nullptr;
  double* densityFacts = nullptr;
  materials    = (unsigned*) malloc(sizeof(unsigned)*totalVox);
  densityFacts = (double*)   malloc(sizeof(double)*totalVox);

  if(materials == nullptr || densityFacts == nullptr){
    if(verbose > 0){
      printf("pen_voxelGeo:loadFile: Error: Unable to allocate memory to read voxel data.\n");
    }
    if(materials != nullptr)
      free(materials);
    if(densityFacts != nullptr)
      free(densityFacts);
    return -2;
  }
  
  //Create a dump reader
  pen_dump reader;
  reader.toDump(sizes,3);
  reader.toDump(materials,totalVox);
  reader.toDump(densityFacts,totalVox);

  err = reader.read(data,pos,verbose);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("pen_voxelGeo:loadFile: Error reading voxel data: %d.\n",err);
    }
    free(materials);
    free(densityFacts);    
    return -3;
  }

  err = setVoxels(nvox,sizes,materials,densityFacts,verbose);
  if(err != 0){
    if(verbose > 0){
      printf("pen_voxelGeo:loadFile: Error storing voxel data: %d.\n",err);
    }
    free(materials);
    free(densityFacts);
    return -4;
  }
  
  free(materials);
  free(densityFacts);
  return 0;
}

int pen_voxelGeo::dump(unsigned char*& data,
		       size_t& pos,
		       const unsigned verbose) const {

  if(getElements() == 0lu){
    if(verbose > 0){
      printf("pen_voxelGeo:dump: Error: No elements to dump\n");
    }
    return -1;
  }
  
  int err;
  pos = 0;
  size_t dataSize;

  //Allocate memory to write data
  double sizes[3] = {dx,dy,dz};
  uint32_t nvox[3] = {nx,ny,nz};
  unsigned long totalVox = getElements();
  unsigned* materials = nullptr;
  double* densityFacts = nullptr;
  materials    = (unsigned*) malloc(sizeof(unsigned)*totalVox);
  densityFacts = (double*)   malloc(sizeof(double)*totalVox);

  if(materials == nullptr || densityFacts == nullptr){
    if(verbose > 0){
      printf("pen_voxelGeo:dump: Error: Unable to allocate memory to write voxel data.\n");
    }
    if(materials != nullptr)
      free(materials);
    if(densityFacts != nullptr)
      free(densityFacts);
    return -2;
  }

  //Fill materials and densities
  for(unsigned long i = 0; i < totalVox; i++){
    materials[i]    = mesh[i].MATER;
    densityFacts[i] = mesh[i].densityFact;
  }
  
  //Create a dump writer
  pen_dump writer;
  writer.toDump(sizes,3);
  writer.toDump(materials,totalVox);
  writer.toDump(densityFacts,totalVox);

  //Get required data size
  dataSize = writer.memory() + 3*sizeof(uint32_t);

  //Allocate buffer
  data = nullptr;
  data = (unsigned char*) malloc(dataSize);
  if(data == nullptr){
    if(verbose > 0){
      printf("pen_voxelGeo:dump: Error: Unable to allocate memory to write voxel data.\n");
    }
    free(materials);
    free(densityFacts);
    return -3;
  }
  
  //Write number of voxels
  memcpy(&data[pos],nvox,3*sizeof(uint32_t));
  pos += 3*sizeof(uint32_t);

  //Dump voxel data
  unsigned char* data2;
  size_t pos2;
  err = writer.dump(data2,pos2,0,verbose);
  //Free material and density arrays
  free(materials);
  free(densityFacts);
  
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("pen_voxelGeo:dump: Error dumping voxel data: %d.\n",err);
    }
    free(data);
    return -4;
  }

  //Copy dumped voxel data to 'data' 
  memcpy(&data[pos],data2,pos2);  
  pos += pos2;
  //Free auxiliar memory
  free(data2);
  return 0;
}

int pen_voxelGeo::dump2File(const char* filename){

  if(filename == nullptr)
    return -1;

  //Open output file
  FILE* fout = nullptr;
  fout = fopen(filename,"wb");
  if(fout == nullptr)
    return -2;

  //Dump data
  unsigned char* data = nullptr;
  size_t size;
  dump(data,size,3);

  //Write data
  fwrite(data,sizeof(unsigned char),size,fout);

  fclose(fout);

  return 0;
}

int pen_voxelGeo::loadFile(const char* filename,
			   const unsigned verbose){

  if(filename == nullptr){
    if(verbose > 0)
      printf("pen_voxelGeo:loadFile: Error: filename is null.\n");
    return -1;
  }

  int err;
  size_t fsize;
  unsigned char* pdata = nullptr;
  FILE* fin = nullptr;
  
  //Open file
  fin = fopen(filename,"rb");
  if(fin == nullptr){
    printf("pen_voxelGeo:loadFile: Error: Unable to open file '%s'\n",filename);
    return -2;
  }

  //Get file size
  fseek(fin, 0L, SEEK_END);
  fsize = ftell(fin);
  rewind(fin);

  //Allocate memory to read the entire file
  pdata = (unsigned char*) malloc(fsize);
  if(pdata == nullptr){
    if(verbose > 0)
      printf("pen_voxelGeo:loadFile: Error: Unable to allocate memory to read voxel file.\n");
    
    fclose(fin);
    return -3;
  }

  size_t read = fread(pdata,
			sizeof(unsigned char),
			fsize,
			fin);
  //Close file
  fclose(fin);
  
  if(read != fsize){
    if(verbose > 0)
      printf("pen_voxelGeo:loadFile: Error: file size and data read mismatch.\n");
    free(pdata);
    return -4;    
  }

  size_t pos = 0;
  err = loadData(pdata,pos,verbose);
  if(err != 0){
    if(verbose > 0){
      printf("pen_voxelGeo:loadFile: Error loading data: %d\n",err);
    }
    free(pdata);
    return -5;
  }
  
  free(pdata);  
  return 0;
  
}


int pen_voxelGeo::printImage(const char* filename) const{

  if(filename == nullptr)
    return -1;
  
  //Create a file to store contours data
  FILE* OutVox = nullptr;
  OutVox = fopen(filename,"w");
  if(OutVox == nullptr){
    return -2;
  }

  fprintf(OutVox,"# \n");
  fprintf(OutVox,"# Voxel geometry file\n");
  fprintf(OutVox,"# Nº of voxels (nx,ny,nz):\n");  
  fprintf(OutVox,"# %5u %5u %5u\n",nx,ny,nz);
  fprintf(OutVox,"# Voxel sizes (dx,dy,dz):\n");  
  fprintf(OutVox,"# %8.5E %8.5E %8.5E\n",dx,dy,dz);
  fprintf(OutVox,"# Voxel data:\n");
  fprintf(OutVox,"#    X(cm)   |    Y(cm)   | MAT | density(vox)/density(mat)\n");

  //Iterate over Z planes
  for(unsigned k = 0; k < nz; k++){
    size_t indexZ = nxy*k;
    fprintf(OutVox,"# Index Z = %4d\n",k);

    //Iterate over rows
    for(unsigned j = 0; j < ny; j++){
      size_t indexYZ = indexZ + j*nx;
      fprintf(OutVox,"# Index Y = %4d\n",j);

      //Iterate over columns
      for(unsigned i = 0; i < nx; i++){
	size_t ivoxel = indexYZ + i;

	//Save voxel X Y and intensity
	fprintf(OutVox," %12.5E %12.5E %4u   %12.5E\n", i*dx, j*dy, mesh[ivoxel].MATER, mesh[ivoxel].densityFact);
      }
      
    }
    //Set a space between planes
    fprintf(OutVox,"\n\n\n");    
  }

  return 0;
}


bool moveIn(double pos,
	    const double dir,
	    double& ds2in,
	    double& ds2out,
	    const double max){

  //Check if a particle at position 'pos' and direction 'dir'
  //is in the geometry limits or will enter to it.  

  // This function asumes that geometry is in the interval [0,max)
  // Also, stores in 'ds' the distance until geometry is reached

  const double inf = 1.0e35;
  const double eps = 1.0e-8;

  if(std::signbit(pos) || pos >= max){
    // position component is out of geometry,
    // check if is moving to it
    if(dir == 0.0){
      //Particle is not moving on this axy, so never will reaches the mesh
      ds2in = ds2out = inf;	 
      return false;
    }

    if(std::signbit(dir) == std::signbit(pos)){
      //Particle is moving away from geometry
      ds2in = ds2out = inf;	 
      return false;	
    }

    //Particle will reach the geometry
    if(std::signbit(dir)){
      //Is moving on negative direction
      ds2in  = (max-pos)/dir+eps;
      ds2out = pos/fabs(dir)-eps;
    }
    else{
      //Particle is moving on positive direction
      ds2in  = fabs(pos)/dir+eps;
      ds2out = (max-pos)/dir-eps;
    }
    return true;
  }
  ds2in = 0.0;
  if(std::signbit(dir)){
    //Negative direction
    ds2out = pos/fabs(dir)-eps;
  }
  else{
    //Positive direction
    ds2out = (max-pos)/dir-eps;
  }
  return true; //Particle is in geometry limits
}


REGISTER_GEOMETRY(pen_voxelGeo,VOXEL)
