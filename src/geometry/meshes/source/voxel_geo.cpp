
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

  //Check file type, ASCII or binary
  bool ASCII = false;
  if(config.read("ascii",ASCII) != INTDATA_SUCCESS){
    ASCII = false;
  }
  if(verbose > 1){
    printf("Expected voxelized file in %s format\n\n", ASCII ? "ASCII" : "binary" );
  }

  if(ASCII){    
    int err2 = loadASCII(filename.c_str(),verbose);  
    if(err2 != 0){
      if(verbose > 0){
	printf("pen_voxelGeo:configure: Error reading and "
	       "loading ASCII voxel geometry file '%s': Corrupted file\n",
	       filename.c_str());
	
	switch(err2){
	case -1:
	  printf("Unexpected error: No filename provided. "
		 "Please, report this error.\n");
	  break;
	case -2:
	  printf("Unable to open geometry file.\n");
	  break;
	case -3:
	  printf("Number of voxels must be positive.\n");	  
	  break;
	case -4:
	  printf("Voxels sizes must be positive.\n");
	  break;
	case -5:
	  printf("Unexpected end of file reached.\n");
	  break;
	case -6:
	  printf("Invalid material index. Material index "
		 "must be greater than 0.\n");
	  break;
	}
      }
      err++;
      configStatus = err;
      return err;
    }
  }
  else{
    int err2 = loadFile(filename.c_str(),verbose);  
    if(err2 != 0){
      if(verbose > 0){
	printf("pen_voxelGeo:configure: Error reading and "
	       "loading binary voxel geometry file '%s': %d\n",
	       filename.c_str(),err2);
      }
      switch(err2){
      case -1:
	printf("Unexpected error: No filename provided. "
	       "Please, report this error.\n");
	break;
      case -2:
	printf("Unable to open geometry file.\n");
	break;
      }
      
      err++;
      configStatus = err;
      return err;
    }
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
    if(fabs((auxd-dx)/dx) > 1.0E-6){
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
    if(fabs((auxd-dy)/dy) > 1.0E-6){
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
    if(fabs((auxd-dz)/dz) > 1.0E-6){
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
  
  // Read enclosure information
  ///////////////////////////////
  double dr;
  if(config.read("enclosure-margin",dr) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_voxelGeo:configure:Error: Enclosure margin value"
	     " 'enclosure-margin' not found. "
	     "Double expected.\n");
    }
    err++;
  }else{
    if(dr < 1.0e-6){
      if(verbose > 0){
	printf("pen_voxelGeo:configure:Error: Enclosure "
	       "margin value must be greater than zero\n");
      }
      err++;
    }
      
    enclosureMargin = dr;
    
    //Precalculate enclosure high margins (+x,+y,+z)
    enclosureXlimit = Mdx + dr;
    enclosureYlimit = Mdy + dr;
    enclosureZlimit = Mdz + dr;

    //Precalculate the enclosure limits x,y,z moving it to the origin (0,0,0)
    enclosureXlimit0 = enclosureXlimit + dr;
    enclosureYlimit0 = enclosureYlimit + dr;
    enclosureZlimit0 = enclosureZlimit + dr;
  }
  
  // Material enclosure
  //**********
	
  //Read material ID
  int auxMat;
  if(config.read("enclosure-material",auxMat) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: Error: Unable to read enclosure material ID for material. Integer expecteed");
    }
    err++;
  }
  //Check material ID
  if(auxMat < 1 || auxMat > (int)constants::MAXMAT){
    if(verbose > 0){
      printf("pen_dicomGeo:configure: Error: Invalid ID specified for enclosure material.");
      printf("                         ID: %d\n",aux);
      printf("Maximum number of materials: %d\n",constants::MAXMAT);
    }
    err++;
  }
  else
    enclosureMat=auxMat;
  //Create a interface between the enclosure and the mesh
  //assigning a detector to the enclosure
  KDET[0] = 1;

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
    
  //Check if it is in the voxel mesh
  locateInMesh(state);
  if(state.MAT == 0){
    //It is not in the mesh, check if it is in the enclosure
    if(state.X > -enclosureMargin && state.X < enclosureXlimit &&
       state.Y > -enclosureMargin && state.Y < enclosureYlimit &&
       state.Z > -enclosureMargin && state.Z < enclosureZlimit){

      //It is in the enclosure
      state.MAT = enclosureMat;
      state.IBODY = 0;
    }
    else{
      state.MAT = 0;
      state.IBODY = constants::MAXMAT+1;
    }
  }
    
}


void pen_voxelGeo::locateInMesh(pen_particleState& state) const{

  long int ix, iy, iz;

  ix = state.X/dx;
  iy = state.Y/dy;
  iz = state.Z/dz;

  if(ix < 0 || (long unsigned)ix >= nx ||
     iy < 0 || (long unsigned)iy >= ny ||
     iz < 0 || (long unsigned)iz >= nz){
    //Particle scapes from geometry mesh
    state.IBODY = constants::MAXMAT+1;
    state.MAT = 0;
  }
  else{
    //Particle is into geometry mesh, calculate vox
    //index and its material
    long int index = iz*nxy + iy*nx + ix;
    state.MAT = mesh[index].MATER;
    state.IBODY = state.MAT;
  }
  
}

void pen_voxelGeo::step(pen_particleState& state,
			double DS,
			double &DSEF,
			double &DSTOT,
			int &NCROSS) const{
    
    //Check if the particle is outside the enclosure
    if(state.MAT == 0){

      //The particle is out, check if it enters the enclosure
      double ds2enclosure;
      if(!enterEnclosure(state.X,state.Y,state.Z,
			 state.U,state.V,state.W,
			 ds2enclosure)){
	
	//Do not enters the enclosure
	state.MAT = 0;
	state.IBODY = constants::MAXMAT+1;
	
	DSEF = 1.0e35;
	DSTOT = 1.0e35;
	NCROSS = 0; //No interface crossed
	
	move(inf,state);
	return;
      }
              
      //The particle reaches the enclosure
      state.MAT = enclosureMat;
      state.IBODY = 0;
      DSEF = DSTOT = ds2enclosure;
      NCROSS = 1;

      move(ds2enclosure,state);
        
      return;        
    }
    
    //Check if the particle is inside the enclosure,
    //but not in the mesh
    if(state.IBODY == 0){
        //Is inside the enclosure, get distance to cross the mesh
        DSEF = DSTOT = 0.0;
        double ds2mesh;
        if(crossMesh(state.X,state.Y,state.Z,
                       state.U,state.V,state.W,
                       ds2mesh)){
            //Aim to the mesh
            if(DS < ds2mesh){
                
                //No reaches the mesh
                move(DS,state);
                DSEF = DS;
                DSTOT = DS;
                NCROSS = 0; //No interface crossed
                return;
            }
            
            //Reaches the mesh
            move(ds2mesh,state);
            
            //Locate the particle in the mesh
            locate(state);
            if(state.IBODY != 0){
                //Has been located in the mesh
                DSEF = ds2mesh;
                DSTOT = ds2mesh;
                NCROSS = 1;
                return;      
            }
            
            //Still remains in the enclosure, move ahead
            DS -= ds2mesh; //Remove traveled distance
            DSEF += ds2mesh;
            DSTOT += ds2mesh;
        }
        
        //The mesh has not been crossed, move inside the enclosure

        //Calculate the distance to exit the enclosure
        double ds2enclosure;
        exitEnclosure(state.X,state.Y,state.Z,
		      state.U,state.V,state.W,
		      ds2enclosure);            
        
        //Check the enclosure limit
        if(ds2enclosure > DS){
            //The enclosure limit has not been reached, the
            //particle remains inside
            move(DS,state);
            DSEF += DS;
            DSTOT += DS;
            NCROSS = 0;            
        } else{
            //The enclosure limit has been reached, the particle
            //scapes to void
            move(inf,state);
            DSEF += ds2enclosure;
            DSTOT = 1.0e35;
            NCROSS = 1;
            
            state.MAT = 0;
            state.IBODY = constants::MAXMAT+1;
        }
        return;      
    }
    
    //The particle is moving inside the mesh
    stepInMesh(state,DS,DSEF,DSTOT,NCROSS);
}

void pen_voxelGeo::stepInMesh(pen_particleState& state,
			double DS,
			double &DSEF,
			double &DSTOT,
			int &NCROSS) const{
  
                
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
  long int voxIncGlob_x, voxIncGlob_y, voxIncGlob_z;
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
    voxIncGlob_y = -static_cast<long int>(nx);    
  }
  else{
    //Moves forward on Y
    ds_y = ((iy+1)*dy-state.Y)*invV;
    voxInc_y = +1;
    voxIncGlob_y = +static_cast<long int>(nx);    
  }

  if(std::signbit(invW)){
    // Moves backward on Z
    ds_z = (iz*dz-state.Z)*invW;
    voxInc_z = -1;    
    voxIncGlob_z = -static_cast<long int>(nxy);    
  }
  else{
    //Moves forward on Z
    ds_z = ((iz+1)*dz-state.Z)*invW;
    voxInc_z = +1;
    voxIncGlob_z = +static_cast<long int>(nxy);    
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
          if(state.IBODY == 0) //The particle escapes to the enclosure
              state.X = ix < 0 ? -1.0e-6 : Mdx + 1.0e-6;
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
          if(state.IBODY == 0) //The particle escapes to the enclosure
              state.Y = iy < 0 ? -1.0e-6 : Mdy + 1.0e-6;          
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
          if(state.IBODY == 0) //The particle escapes to the enclosure
              state.Z = iz < 0 ? -1.0e-6 : Mdz + 1.0e-6;                    
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
      //the greater material index plus one for the enclosure
      if(nBodies < mats[i]){
	nBodies = mats[i];
      }
      mesh[i].MATER = mats[i];
      mesh[i].densityFact = dens[i];
    }
  }
  
  //Add an extra body to allocate the enclosure
  nBodies += 1;
  
  return err;
  
}

unsigned pen_voxelGeo::getIBody(const char* bname) const {

  int index = atoi(bname);
  if(index >= (int)nBodies)
    return constants::MAXMAT+1;
  else
    return unsigned(index);
}

std::string pen_voxelGeo::getBodyName(const unsigned ibody) const{
  if(ibody < nBodies){
    return std::to_string(ibody);
  }else{
    return std::string("NONE");
  }
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

int pen_voxelGeo::loadASCII(const char* filename,
			    const unsigned verbose) {

  if(filename == nullptr)
    return -1;
  
  FILE* input = nullptr;
  input = fopen(filename, "r");
  if(input == nullptr){
    return -2;
  }
  
  //Read number of voxels
  int nvoxAux[3];
  fscanf(input, " %d %d %d", &nvoxAux[0], &nvoxAux[1], &nvoxAux[2]);

  if(nvoxAux[0] < 1 ||
     nvoxAux[1] < 1 ||
     nvoxAux[2] < 1){
    fclose(input);
    return -3;
  }
  
  unsigned nvox[3] = {
    static_cast<unsigned>(nvoxAux[0]),
    static_cast<unsigned>(nvoxAux[1]),
    static_cast<unsigned>(nvoxAux[2])
  };

  //Read voxel sizes
  double size[3];
  fscanf(input, " %lE %lE %lE", &size[0], &size[1], &size[2]);  

  if(size[0] <= 0.0 ||
     size[1] <= 0.0 ||
     size[2] <= 0.0){
    fclose(input);
    return -4;
  }
  
  //Create vectors to store materials and density factors
  size_t nbins = nvox[0]*nvox[1]*nvox[2];
  std::vector<unsigned> materials(nbins);
  std::vector<double> densityFacts(nbins);
  
  for(size_t i = 0; i < nbins; ++i){
    int auxMat;
    int nread = fscanf(input, " %d %lE ", &auxMat, &densityFacts[i]);
    if(nread != 2){
      fclose(input);
      return -5;
    }

    if(auxMat <= 0){
      fclose(input);
      return -6;
    }
    materials[i] = static_cast<unsigned>(auxMat);
  }

  fclose(input);
  
  int err = setVoxels(nvox,size,materials.data(),densityFacts.data(),verbose);

  return err;
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
  
  fprintf(OutVox," %5u %5u %5u\n",nx,ny,nz);
  fprintf(OutVox,"%8.5E %8.5E %8.5E\n",dx,dy,dz);
  
  //Iterate over Z planes
  for(unsigned k = 0; k < nz; k++){
    size_t indexZ = nxy*k;

    //Iterate over rows
    for(unsigned j = 0; j < ny; j++){
      size_t indexYZ = indexZ + j*nx;

      //Iterate over columns
      for(unsigned i = 0; i < nx; i++){
	size_t ivoxel = indexYZ + i;

	//Save voxel X Y and intensity
	fprintf(OutVox,"%4u   %12.5E\n", mesh[ivoxel].MATER, mesh[ivoxel].densityFact);
      }      
    }
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
      //Particle is not moving on this axis, so never will reaches the mesh
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
