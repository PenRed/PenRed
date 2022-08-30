 
//
//
//    Copyright (C) 2021-2022 Universitat de València - UV
//    Copyright (C) 2021-2022 Universitat Politècnica de València - UPV
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
//        vicente.gimenez@uv.es (Vicent Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Olver Gil)
//

#include "pen_geoView.hh"

bool pen_geoView::checkStep(pen_particleState& state,
			    double toTravel,
			    geoError& err) const{

  double dsef,dstot;
  int ncross;

  //Get particle material and body
  geometry->locate(state);

  //Get particle material and body
  //in the objective position
  pen_particleState objState = state;
  objState.X += objState.U*toTravel;
  objState.Y += objState.V*toTravel;
  objState.Z += objState.W*toTravel;

  geometry->locate(objState);

  unsigned int BODY = objState.IBODY;
  unsigned int MAT  = objState.MAT;

  //Save initial position and index
  err.from[0] = state.X;
  err.from[1] = state.Y;
  err.from[2] = state.Z;

  err.iIBODY = state.IBODY;
  err.iMAT   = state.MAT;
    
  err.eIBODY = BODY;
  err.eMAT = MAT;
  /////////////////////////////////
    
  //Move the particle
  do{
    geometry->step(state,toTravel,dsef,dstot,ncross);
    toTravel -= dstot;
  }while(toTravel > 1.0e-6);

  //Save final position
  err.to[0] = state.X;
  err.to[1] = state.Y;
  err.to[2] = state.Z;    

  err.fIBODY = state.IBODY;
  err.fMAT = state.MAT;
  /////////////////////////////////

  if(toTravel < -1.0e-6){
    //A void region has been crossed in the final step

    //Flags the end in a void region
    err.endInVoid = true;
    /////////////////////////////////
      
    if(MAT != 0){
      //locate and step do not match
      return false;
    }
  }else{
    //Non void region has been crossed in the final step

    //Flags the end in a non void region
    err.endInVoid = false;
    /////////////////////////////////
      
    if(BODY != state.IBODY){
      return false;
    }    
    if(MAT != state.MAT){
      return false;
    }
  }
  //The particle information match with the expected indexes
  return true;
}

int pen_geoView::init(const pen_parserSection& configIn,
		      const unsigned verbose){
  pen_parserSection config = configIn;

  //Append dummy material information to geometry section
  for(unsigned imat = 0; imat < constants::MAXMAT; ++imat){
    char key[400];
    sprintf(key,"materials/mat%03d/ID",imat+1);
    config.set(key,(int)imat+1);
    sprintf(key,"materials/mat%03d/density",imat+1);
    config.set(key,1.0);
  }

  //Get geometry type
  std::string geoType;
  if(config.read("type",geoType) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_geoView: init: Error: field 'type' not specified. "
	     "String expected.\n");
    }
    return -1;
  }

  //Remove a possible previous geometry
  if(geometry != nullptr){
    delete geometry;
    geometry = nullptr;
  }

  //Create geometry
  geometry = penGeoRegister_create(geoType.c_str());
  if(geometry == nullptr){
    if(verbose > 0){
      printf("pen_geoView: init: Error creating geometry instance "
	     "of type '%s'\n", geoType.c_str());
    }
    return -2;
  }
    
  //Configure geometry  
  geometry->name.assign("geometry");    
  geometry->configure(config,verbose);

  //Check errors
  if(geometry->configureStatus() != 0){
    if(verbose > 0)
      printf("pen_geoView: init: Error: Fail on geometry configuration.\n");
    return -3;
  }

  if(verbose > 2){
    printf("pen_geoView: init: Geometry load successfully\n");
  }
  
  return 0;
}

void pen_geoView::testX(std::vector<geoError>& errors,
			const float x, const float y, const float z,
			const float dy, const float dz,
			const unsigned ny, const unsigned nz) const{
  
  //Calculate the image origin (top left corner)
  float oz, oy;
  oy = y-dy*static_cast<float>(ny)/2.0;
  oz = z+dz*static_cast<float>(nz)/2.0;

  pen_particleState state;

  state.X = x;
  
  for(unsigned j = 0; j < nz; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;
    for(unsigned i = 0; i < ny; ++i){
      state.Y = oy + (static_cast<float>(i)+0.5)*dy;

      pen_particleState stateYp = state;
      pen_particleState stateYm = state;
      pen_particleState stateZp = state;
      pen_particleState stateZm = state;
      
      stateYp.U = 0.0; stateYp.V = +1.0; stateYp.W = +0.0;
      stateYm.U = 0.0; stateYm.V = -1.0; stateYm.W = +0.0;
      stateZp.U = 0.0; stateZp.V = +0.0; stateZp.W = +1.0;
      stateZm.U = 0.0; stateZm.V = +0.0; stateZm.W = -1.0;
      
      geoError err;
      
      if(!checkStep(stateYp,dy,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateYm,dy,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateZp,dz,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateZm,dz,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
    }
  }
}

void pen_geoView::renderX(unsigned char* renderMat,unsigned int* renderBody,
			  const float x, const float y, const float z,
			  const float dy, const float dz,
			  const unsigned ny, const unsigned nz,
			  const unsigned nthreads) const{

  //Calculate the image origin (top left corner)
  const float oy = y-dy*static_cast<float>(ny)/2.0;
  const float oz = z+dz*static_cast<float>(nz)/2.0;

  if(nthreads < 2){
    pen_particleState state;

    state.X = x;
  
    for(unsigned j = 0; j < nz; ++j){
      state.Z = oz - (static_cast<float>(j)+0.5)*dz;
      for(unsigned i = 0; i < ny; ++i){
	state.Y = oy + (static_cast<float>(i)+0.5)*dy;

	geometry->locate(state);

	renderMat[j*ny+i] = static_cast<unsigned char>(state.MAT);
	renderBody[j*ny+i] = static_cast<unsigned int>(state.IBODY);
      }
    }
  }else{

    std::vector<std::thread> threads;

    for(size_t ith = 0; ith < nthreads; ++ith){
    
      threads.push_back(std::thread([&,ith,nthreads](){

	pen_particleState state;

	state.X = x;
  
	for(unsigned j = ith; j < nz; j += nthreads){
	  state.Z = oz - (static_cast<float>(j)+0.5)*dz;
	  for(unsigned i = 0; i < ny; ++i){
	    state.Y = oy + (static_cast<float>(i)+0.5)*dy;

	    geometry->locate(state);

	    renderMat[j*ny+i] = static_cast<unsigned char>(state.MAT);
	    renderBody[j*ny+i] = static_cast<unsigned int>(state.IBODY);
	  }
	}
	
      }));
    }

    for(auto& thread : threads){
      thread.join();
    }
  }
  
}

void pen_geoView::renderXtoLeft(unsigned char* renderMat,
				unsigned int* renderBody,
				const unsigned nPixels,				
				const float x, const float y, const float z,
				const float dy, const float dz,
				const unsigned ny, const unsigned nz) const{

  //Calculate the image origin (top left corner)
  float oz, oy;
  oy = y-dy*(static_cast<float>(ny)/2.0 + nPixels);
  oz = z+dz*static_cast<float>(nz)/2.0;

  pen_particleState state;

  state.X = x;
  
  for(unsigned j = 0; j < nz; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;

    //Copy the already rendered image zone
    for(unsigned i = ny-1; i >= nPixels; --i){
      renderMat[j*ny+i] = renderMat[j*ny+i-nPixels];
      renderBody[j*ny+i] = renderBody[j*ny+i-nPixels];
    }
    
    //Render the new zone
    for(unsigned i = 0; i < nPixels; ++i){
      state.Y = oy + (static_cast<float>(i)+0.5)*dy;

      geometry->locate(state);

      renderMat[j*ny+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*ny+i] = static_cast<unsigned int>(state.IBODY);
    }

  }
  
}

void pen_geoView::renderXtoRight(unsigned char* renderMat,
				 unsigned int* renderBody,
				 const unsigned nPixels,
				 const float x, const float y, const float z,
				 const float dy, const float dz,
				 const unsigned ny, const unsigned nz) const{

  //Calculate the image origin (top left corner)
  float oz, oy;
  oy = y-dy*(static_cast<float>(ny)/2.0 - nPixels);
  oz = z+dz*static_cast<float>(nz)/2.0;

  pen_particleState state;

  state.X = x;
  
  for(unsigned j = 0; j < nz; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;

    //Copy the already rendered image zone
    for(unsigned i = 0; i < ny-nPixels; ++i){
      renderMat[j*ny+i] = renderMat[j*ny+i+nPixels];
      renderBody[j*ny+i] = renderBody[j*ny+i+nPixels];
    }
    
    //Render the new zone
    for(unsigned i = ny-nPixels; i < ny; ++i){
      state.Y = oy + (static_cast<float>(i)+0.5)*dy;

      geometry->locate(state);

      renderMat[j*ny+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*ny+i] = static_cast<unsigned int>(state.IBODY);
    }

  }
  
}

void pen_geoView::renderXtoUp(unsigned char* renderMat,
				unsigned int* renderBody,
				const unsigned nPixels,				
				const float x, const float y, const float z,
				const float dy, const float dz,
				const unsigned ny, const unsigned nz) const{

  //Calculate the image origin (top left corner)
  float oz, oy;
  oy = y-dy*static_cast<float>(ny)/2.0;
  oz = z+dz*(static_cast<float>(nz)/2.0 + nPixels);

  pen_particleState state;

  state.X = x;
  
  //Copy the already rendered image zone
  for(unsigned j = nz-1; j >= nPixels; --j){

    for(unsigned i = 0; i < ny; ++i){
      renderMat[j*ny+i] = renderMat[(j-nPixels)*ny+i];
      renderBody[j*ny+i] = renderBody[(j-nPixels)*ny+i];
    }
  }

  //Render the new zone
  for(unsigned j = 0; j < nPixels; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;

    for(unsigned i = 0; i < ny; ++i){
      state.Y = oy + (static_cast<float>(i)+0.5)*dy;

      geometry->locate(state);

      renderMat[j*ny+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*ny+i] = static_cast<unsigned int>(state.IBODY);
    }
  }
  
}

void pen_geoView::renderXtoDown(unsigned char* renderMat,
				 unsigned int* renderBody,
				 const unsigned nPixels,
				 const float x, const float y, const float z,
				 const float dy, const float dz,
				 const unsigned ny, const unsigned nz) const{

  //Calculate the image origin (top left corner)
  float oz, oy;
  oy = y-dy*static_cast<float>(ny)/2.0;
  oz = z+dz*(static_cast<float>(nz)/2.0 - nPixels);

  pen_particleState state;

  state.X = x;

  //Copy already rendered zone
  for(unsigned j = 0; j < nz-nPixels; ++j){

    for(unsigned i = 0; i < ny; ++i){
      renderMat[j*ny+i] = renderMat[(j+nPixels)*ny+i];
      renderBody[j*ny+i] = renderBody[(j+nPixels)*ny+i];      
    }
  }

  for(unsigned j = nz-nPixels; j < nz; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;

    for(unsigned i = 0; i < ny; ++i){
      state.Y = oy + (static_cast<float>(i)+0.5)*dy;

      geometry->locate(state);

      renderMat[j*ny+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*ny+i] = static_cast<unsigned int>(state.IBODY);
    }    
  }  
}

void pen_geoView::testY(std::vector<geoError>& errors,
			const float x, const float y, const float z,
			const float dx, const float dz,
			const unsigned nx, const unsigned nz) const{
  
  //Calculate the image origin (top left corner)
  float ox, oz;
  ox = x-dx*static_cast<float>(nx)/2.0;
  oz = z+dz*static_cast<float>(nz)/2.0;
  
  pen_particleState state;

  state.Y = y;  
  
  for(unsigned j = 0; j < nz; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;
    for(unsigned i = 0; i < nx; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      pen_particleState stateXp = state;
      pen_particleState stateXm = state;
      pen_particleState stateZp = state;
      pen_particleState stateZm = state;
      
      stateXp.U = +1.0; stateXp.V = +0.0; stateXp.W = +0.0;
      stateXm.U = -1.0; stateXm.V = +0.0; stateXm.W = +0.0;
      stateZp.U = +0.0; stateZp.V = +0.0; stateZp.W = +1.0;
      stateZm.U = +0.0; stateZm.V = +0.0; stateZm.W = -1.0;

      geoError err;
      
      if(!checkStep(stateXp,dx,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateXm,dx,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateZp,dz,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateZm,dz,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
    }
  }
}

void pen_geoView::renderY(unsigned char* renderMat,unsigned int* renderBody,
			  const float x, const float y, const float z,
			  const float dx, const float dz,
			  const unsigned nx, const unsigned nz,
			  const unsigned nthreads) const{

  //Calculate the image origin (top left corner)
  float ox, oz;
  ox = x-dx*static_cast<float>(nx)/2.0;
  oz = z+dz*static_cast<float>(nz)/2.0;

  if(nthreads < 2){
    pen_particleState state;

    state.Y = y;
  
    for(unsigned j = 0; j < nz; ++j){
      state.Z = oz - (static_cast<float>(j)+0.5)*dz;
      for(unsigned i = 0; i < nx; ++i){
	state.X = ox + (static_cast<float>(i)+0.5)*dx;

	geometry->locate(state);

	renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
	renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
      }
    }
  }else{

    std::vector<std::thread> threads;

    for(size_t ith = 0; ith < nthreads; ++ith){


      threads.push_back(std::thread([&,ith,nthreads](){

	pen_particleState state;

	state.Y = y;
      
	for(unsigned j = ith; j < nz; j += nthreads){
	  state.Z = oz - (static_cast<float>(j)+0.5)*dz;
	  for(unsigned i = 0; i < nx; ++i){
	    state.X = ox + (static_cast<float>(i)+0.5)*dx;

	    geometry->locate(state);

	    renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
	    renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
	  }
	}
      }));
    }

    for(auto& thread : threads){
      thread.join();
    }
  }
  
}

void pen_geoView::renderYtoLeft(unsigned char* renderMat,
				unsigned int* renderBody,
				const unsigned nPixels,				
				const float x, const float y, const float z,
				const float dx, const float dz,
				const unsigned nx, const unsigned nz) const{

  //Calculate the image origin (top left corner)
  float oz, ox;
  ox = x-dx*(static_cast<float>(nx)/2.0 + nPixels);
  oz = z+dz*static_cast<float>(nz)/2.0;

  pen_particleState state;

  state.Y = y;
  
  for(unsigned j = 0; j < nz; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;

    //Copy the already rendered image zone
    for(unsigned i = nx-1; i >= nPixels; --i){
      renderMat[j*nx+i] = renderMat[j*nx+i-nPixels];
      renderBody[j*nx+i] = renderBody[j*nx+i-nPixels];
    }
    
    //Render the new zone
    for(unsigned i = 0; i < nPixels; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      geometry->locate(state);

      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
    }

  }
  
}

void pen_geoView::renderYtoRight(unsigned char* renderMat,
				 unsigned int* renderBody,
				 const unsigned nPixels,
				 const float x, const float y, const float z,
				 const float dx, const float dz,
				 const unsigned nx, const unsigned nz) const{

  //Calculate the image origin (top left corner)
  float oz, ox;
  ox = x-dx*(static_cast<float>(nx)/2.0 - nPixels);
  oz = z+dz*static_cast<float>(nz)/2.0;

  pen_particleState state;

  state.Y = y;
  
  for(unsigned j = 0; j < nz; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;

    //Copy the already rendered image zone
    for(unsigned i = 0; i < nx-nPixels; ++i){
      renderMat[j*nx+i] = renderMat[j*nx+i+nPixels];
      renderBody[j*nx+i] = renderBody[j*nx+i+nPixels];
    }
    
    //Render the new zone
    for(unsigned i = nx-nPixels; i < nx; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      geometry->locate(state);

      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
    }

  }
  
}

void pen_geoView::renderYtoUp(unsigned char* renderMat,
				unsigned int* renderBody,
				const unsigned nPixels,				
				const float x, const float y, const float z,
				const float dx, const float dz,
				const unsigned nx, const unsigned nz) const{

  //Calculate the image origin (top left corner)
  float oz, ox;
  ox = x-dx*static_cast<float>(nx)/2.0;
  oz = z+dz*(static_cast<float>(nz)/2.0 + nPixels);

  pen_particleState state;

  state.Y = y;
  
  //Copy the already rendered image zone
  for(unsigned j = nz-1; j >= nPixels; --j){

    for(unsigned i = 0; i < nx; ++i){
      renderMat[j*nx+i] = renderMat[(j-nPixels)*nx+i];
      renderBody[j*nx+i] = renderBody[(j-nPixels)*nx+i];
    }
  }

  //Render the new zone
  for(unsigned j = 0; j < nPixels; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;

    for(unsigned i = 0; i < nx; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      geometry->locate(state);

      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
    }
  }
  
}

void pen_geoView::renderYtoDown(unsigned char* renderMat,
				 unsigned int* renderBody,
				 const unsigned nPixels,
				 const float x, const float y, const float z,
				 const float dx, const float dz,
				 const unsigned nx, const unsigned nz) const{

  //Calculate the image origin (top left corner)
  float oz, ox;
  ox = x-dx*static_cast<float>(nx)/2.0;
  oz = z+dz*(static_cast<float>(nz)/2.0 - nPixels);

  pen_particleState state;

  state.Y = y;

  //Copy already rendered zone
  for(unsigned j = 0; j < nz-nPixels; ++j){

    for(unsigned i = 0; i < nx; ++i){
      renderMat[j*nx+i] = renderMat[(j+nPixels)*nx+i];
      renderBody[j*nx+i] = renderBody[(j+nPixels)*nx+i];      
    }
  }

  for(unsigned j = nz-nPixels; j < nz; ++j){
    state.Z = oz - (static_cast<float>(j)+0.5)*dz;

    for(unsigned i = 0; i < nx; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      geometry->locate(state);

      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
    }
  }
}

void pen_geoView::testZ(std::vector<geoError>& errors,
			const float x, const float y, const float z,
			const float dx, const float dy,
			const unsigned nx, const unsigned ny) const{
  
  //Calculate the image origin (top left corner)
  float ox, oy;
  ox = x-dx*static_cast<float>(nx)/2.0;
  oy = y+dy*static_cast<float>(ny)/2.0;

  pen_particleState state;

  state.Z = z;
  
  for(unsigned j = 0; j < ny; ++j){
    state.Y = oy - (static_cast<float>(j)+0.5)*dy;
    for(unsigned i = 0; i < nx; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      pen_particleState stateXp = state;
      pen_particleState stateXm = state;
      pen_particleState stateYp = state;
      pen_particleState stateYm = state;
      
      stateXp.U = +1.0; stateXp.V = +0.0; stateXp.W = +0.0;
      stateXm.U = -1.0; stateXm.V = +0.0; stateXm.W = +0.0;
      stateYp.U = +0.0; stateYp.V = +1.0; stateYp.W = +0.0;
      stateYm.U = +0.0; stateYm.V = -1.0; stateYm.W = +0.0;

      geoError err;
      
      if(!checkStep(stateXp,dx,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateXm,dx,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateYp,dy,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
      if(!checkStep(stateYm,dy,err)){
	//Inconsistency detected
	errors.push_back(err);
      }
    }
  }
}

void pen_geoView::renderZ(unsigned char* renderMat,unsigned int* renderBody,
			  const float x, const float y, const float z,
			  const float dx, const float dy,
			  const unsigned nx, const unsigned ny,
			  const unsigned nthreads) const{

  //Calculate the image origin (top left corner)
  float ox, oy;
  ox = x-dx*static_cast<float>(nx)/2.0;
  oy = y+dy*static_cast<float>(ny)/2.0;

  
  if(nthreads < 2){
  
    pen_particleState state;

    state.Z = z;
  
    for(unsigned j = 0; j < ny; ++j){
      state.Y = oy - (static_cast<float>(j)+0.5)*dy;
      for(unsigned i = 0; i < nx; ++i){
	state.X = ox + (static_cast<float>(i)+0.5)*dx;

	geometry->locate(state);

	renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
	renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
      }
    }

  }else{

    std::vector<std::thread> threads;

    for(size_t ith = 0; ith < nthreads; ++ith){
    
      threads.push_back(std::thread([&,ith,nthreads](){

	pen_particleState state;

	state.Z = z;
  
	for(unsigned j = ith; j < ny; j += nthreads){
	  state.Y = oy - (static_cast<float>(j)+0.5)*dy;
	  for(unsigned i = 0; i < nx; ++i){
	    state.X = ox + (static_cast<float>(i)+0.5)*dx;

	    geometry->locate(state);

	    renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
	    renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
	  }
	}
	
      }));
    }

    for(auto& thread : threads){
      thread.join();
    }
  }
  
}

void pen_geoView::renderZtoLeft(unsigned char* renderMat,
				unsigned int* renderBody,
				const unsigned nPixels,				
				const float x, const float y, const float z,
				const float dx, const float dy,
				const unsigned nx, const unsigned ny) const{

  //Calculate the image origin (top left corner)
  float oy, ox;
  ox = x-dx*(static_cast<float>(nx)/2.0 + nPixels);
  oy = y+dy*static_cast<float>(ny)/2.0;

  pen_particleState state;

  state.Z = z;
  
  for(unsigned j = 0; j < ny; ++j){
    state.Y = oy - (static_cast<float>(j)+0.5)*dy;

    //Copy the already rendered image zone
    for(unsigned i = nx-1; i >= nPixels; --i){
      renderMat[j*nx+i] = renderMat[j*nx+i-nPixels];
      renderBody[j*nx+i] = renderBody[j*nx+i-nPixels];
    }
    
    //Render the new zone
    for(unsigned i = 0; i < nPixels; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      geometry->locate(state);

      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
    }

  }
  
}

void pen_geoView::renderZtoRight(unsigned char* renderMat,
				 unsigned int* renderBody,
				 const unsigned nPixels,
				 const float x, const float y, const float z,
				 const float dx, const float dy,
				 const unsigned nx, const unsigned ny) const{

  //Calculate the image origin (top left corner)
  float oy, ox;
  ox = x-dx*(static_cast<float>(nx)/2.0 - nPixels);
  oy = y+dy*static_cast<float>(ny)/2.0;

  pen_particleState state;

  state.Z = z;
  
  for(unsigned j = 0; j < ny; ++j){
    state.Y = oy - (static_cast<float>(j)+0.5)*dy;

    //Copy the already rendered image zone
    for(unsigned i = 0; i < nx-nPixels; ++i){
      renderMat[j*nx+i] = renderMat[j*nx+i+nPixels];
      renderBody[j*nx+i] = renderBody[j*nx+i+nPixels];
    }
    
    //Render the new zone
    for(unsigned i = nx-nPixels; i < nx; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      geometry->locate(state);

      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
    }

  }
  
}

void pen_geoView::renderZtoUp(unsigned char* renderMat,
				unsigned int* renderBody,
				const unsigned nPixels,				
				const float x, const float y, const float z,
				const float dx, const float dy,
				const unsigned nx, const unsigned ny) const{

  //Calculate the image origin (top left corner)
  float oy, ox;
  ox = x-dx*static_cast<float>(nx)/2.0;
  oy = y+dy*(static_cast<float>(ny)/2.0 + nPixels);

  pen_particleState state;

  state.Z = z;
  
  //Copy the already rendered image zone
  for(unsigned j = ny-1; j >= nPixels; --j){

    for(unsigned i = 0; i < nx; ++i){
      renderMat[j*nx+i] = renderMat[(j-nPixels)*nx+i];
      renderBody[j*nx+i] = renderBody[(j-nPixels)*nx+i];
    }
  }

  //Render the new zone
  for(unsigned j = 0; j < nPixels; ++j){
    state.Y = oy - (static_cast<float>(j)+0.5)*dy;

    for(unsigned i = 0; i < nx; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      geometry->locate(state);

      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
    }
  }
  
}

void pen_geoView::renderZtoDown(unsigned char* renderMat,
				 unsigned int* renderBody,
				 const unsigned nPixels,
				 const float x, const float y, const float z,
				 const float dx, const float dy,
				 const unsigned nx, const unsigned ny) const{

  //Calculate the image origin (top left corner)
  float oy, ox;
  ox = x-dx*static_cast<float>(nx)/2.0;
  oy = y+dy*(static_cast<float>(ny)/2.0 - nPixels);

  pen_particleState state;

  state.Z = z;

  //Copy already rendered zone
  for(unsigned j = 0; j < ny-nPixels; ++j){

    for(unsigned i = 0; i < nx; ++i){
      renderMat[j*nx+i] = renderMat[(j+nPixels)*nx+i];
      renderBody[j*nx+i] = renderBody[(j+nPixels)*nx+i];      
    }
  }

  for(unsigned j = ny-nPixels; j < ny; ++j){
    state.Y = oy - (static_cast<float>(j)+0.5)*dy;

    for(unsigned i = 0; i < nx; ++i){
      state.X = ox + (static_cast<float>(i)+0.5)*dx;

      geometry->locate(state);

      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
    }
  }
}

int pen_geoView::render3Dortho(unsigned char* renderMat,unsigned int* renderBody,
			       const float x, const float y, const float z,
			       const float u, const float v, const float w,
			       const float roll, float& phi,
			       const float dx, const float dy,
			       const unsigned nx, const unsigned ny,
			       float* distances,
			       float& minDistance, float& maxDistance,
			       const float threshold) const{
  
  float norm = sqrt(u*u + v*v + w*w);
  if(norm < 1.0e-6){
    return -1;
  }

  float unorm = u/norm;
  float vnorm = v/norm;
  float wnorm = w/norm;

  float maxD = -1.0;
  float minD = 1.0e35;
  
  //Get the rotation to set the camera
  float rotation[9];
  phi = rollAlignf(unorm,vnorm,wnorm,roll,rotation,phi,threshold);
  
  //Calculate the image origin
  float ox, oy;
  ox = -dx*static_cast<float>(nx)/2.0;
  oy = -dy*static_cast<float>(ny)/2.0;

  pen_particleState state;
  
  state.U = unorm;
  state.V = vnorm;
  state.W = wnorm;

  float px, py;
  for(unsigned j = 0; j < ny; ++j){
    py = oy + (static_cast<float>(j)+0.5)*dy;
    for(unsigned i = 0; i < nx; ++i){
      px = ox + (static_cast<float>(i)+0.5)*dx;

      float pos[3] = {px,py,0.0};
      matmul3D(rotation,pos);

      state.X = pos[0] + x;
      state.Y = pos[1] + y;
      state.Z = pos[2] + z;
      
      geometry->locate(state);

      double dsef,dstot;
      int ncross;
      geometry->step(state,1.0e35,dsef,dstot,ncross);
      
      renderMat[j*nx+i] = static_cast<unsigned char>(state.MAT);
      renderBody[j*nx+i] = static_cast<unsigned int>(state.IBODY);
      const float fdstot = static_cast<float>(dstot);
      distances[j*nx+i] = fdstot;
      if(state.MAT != 0 && maxD < fdstot) maxD = fdstot;
      if(minD > fdstot) minD = fdstot;      
    }
  }

  maxDistance = maxD;
  minDistance = minD;
  
  return 0;
  
}

void pen_geoView::set3DResolution(const unsigned nx, const unsigned ny,
				  const float dx, const float dy,
				  const float perspective){

  perspective3D = perspective;
  
  nx3D = nx;
  ny3D = ny;
  nxy3D = static_cast<unsigned long>(nx3D)*static_cast<unsigned long>(ny3D);

  dx3D = dx;
  dy3D = dy;

  float dxtot = dx3D*static_cast<float>(nx3D);
  float dytot = dy3D*static_cast<float>(ny3D);

  float anglePerDs = perspective3D/std::max(dxtot,dytot);
  
  perspectiveVect.resize(nxy3D*3);

  //Calculate the image origin
  ox3D = -dx3D*static_cast<float>(nx3D)/2.0;
  oy3D = -dy3D*static_cast<float>(ny3D)/2.0;

  
  
  for(unsigned j = 0; j < ny3D; ++j){
    float y = oy3D + (static_cast<float>(j)+0.5)*dy3D;
    //Rotation on X axis
    float gamma = -anglePerDs*y;
    float cgamma = cos(gamma);
    float sgamma = sin(gamma);
          
    for(unsigned i = 0; i < nx3D; ++i){

      float x = ox3D + (static_cast<float>(i)+0.5)*dx3D;
      //Rotation on Y axis
      float beta = anglePerDs*x;
      float cbeta = cos(beta);
      float sbeta = sin(beta);

      unsigned long index = (j*nx3D+i)*3;

      //Obtain rotation matrix for the normal vector
      float rotation[9];
      rotation[0] = cbeta;
      rotation[1] = sbeta*sgamma;
      rotation[2] = sbeta*cgamma;
      rotation[3] = 0.0;
      rotation[4] = cgamma;
      rotation[5] = -sgamma;
      rotation[6] = -sbeta;
      rotation[7] = cbeta*sgamma;
      rotation[8] = cbeta*cgamma;

      //Rotate normal
      perspectiveVect[index]   = 0.0;
      perspectiveVect[index+1] = 0.0;
      perspectiveVect[index+2] = 1.0;
      matmul3D(rotation,&perspectiveVect[index]);
    }
  }
  
}

int pen_geoView::render3D(unsigned char* renderMat,unsigned int* renderBody,
			  const float x, const float y, const float z,
			  const float u, const float v, const float w,
			  const float roll, float& phi,
			  float* distances,
			  float& minDistance, float& maxDistance,
			  const float threshold) const{
  
  float norm = sqrt(u*u + v*v + w*w);
  if(norm < 1.0e-6){
    return -1;
  }

  float unorm = u/norm;
  float vnorm = v/norm;
  float wnorm = w/norm;

  float maxD = -1.0;
  float minD = 1.0e35;
  
  //Get the rotation to set the camera
  float rotation[9];
  phi = rollAlignf(unorm,vnorm,wnorm,roll,rotation,phi,threshold);
  
  pen_particleState state;
  
  float px, py;
  for(unsigned j = 0; j < ny3D; ++j){
    py = oy3D + (static_cast<float>(j)+0.5)*dy3D;
    for(unsigned i = 0; i < nx3D; ++i){
      px = ox3D + (static_cast<float>(i)+0.5)*dx3D;

      float pos[3] = {px,py,0.0};
      matmul3D(rotation,pos);

      state.X = pos[0] + x;
      state.Y = pos[1] + y;
      state.Z = pos[2] + z;

      //Set perspective on direction
      unsigned long index = (j*nx3D+i)*3;
      float dir[3] = {perspectiveVect[index],
	perspectiveVect[index+1],
	perspectiveVect[index+2]};
      
      matmul3D(rotation,dir);
      
      state.U = dir[0];
      state.V = dir[1];
      state.W = dir[2];
      
      geometry->locate(state);
      unsigned int initMat = state.MAT;
      unsigned int initBody = state.IBODY;
      
      double dsef,dstot;
      int ncross;
      geometry->step(state,1.0e35,dsef,dstot,ncross);
      
      renderMat[j*nx3D+i] =
	static_cast<unsigned char>(state.MAT == 0 ? initMat : state.MAT);
      renderBody[j*nx3D+i] =
	static_cast<unsigned int>(state.MAT == 0 ? initBody : state.IBODY);
      const float fdsef = static_cast<float>(dsef);
      distances[j*nx3D+i] = fdsef;
      if(state.MAT != 0 && maxD < fdsef) maxD = fdsef;
      if(minD > fdsef) minD = fdsef;
    }
  }

  maxDistance = maxD;
  minDistance = minD;
  
  return 0;
  
}
