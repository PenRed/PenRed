
//
//
//    Copyright (C) 2022-2023 Universitat de València - UV
//    Copyright (C) 2022-2023 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es
//        vicente.gimenez@uv.es
//


#include "meshBody_geo.hh"

bool pen_meshBody::inside(const v3D pos) const{

  if(!boundary.in(pos, crossThreshold))
    return false;
  v3D origin(pos);
  origin.x = boundary.minx()-1.0;
  v3D dir(1.0,0.0,0.0);
  double dsMax = pos.x - origin.x;

  //Flag if is outside or inside
  bool back = false; //Outside requires frontface check

  while(dsMax > 0.0){
    double ds;
    if(cross(origin,dir,ds,back,dsMax)){
      origin.x += ds;
      dsMax -= ds;
      back = !back;
    }else{
      return back;
    }
  }

  return back;
}

    
bool pen_meshBody::cross(const v3D pos,
			 const v3D dir,
			 double& ds,
			 const bool back,
			 const double maxDs) const{
  
  const double eps = crossThreshold;
  const double maxDsThres = maxDs + crossThreshold;
  double lowestDS = maxDsThres;

  //If the position is outside, check if the ray reaches the object boundary
  if(!back){
    double ds2In;
    bool goesIn = boundary.toIn(pos.x,pos.y,pos.z,dir.x,dir.y,dir.z,ds2In);
    if(!goesIn || ds2In > lowestDS)
      return false;
  }

  //Iterate over mesh body super regions
  for(const superRegion& supRegion : regions){

    //Check if this super region is crossed
    double ds2In;
    bool goesIn = supRegion.toIn(pos.x,pos.y,pos.z,dir.x,dir.y,dir.z,ds2In);
    if(!goesIn || ds2In > lowestDS)
      continue;
    
    //Iterate over super regions sub regions
    for(const triangleRegion& region : supRegion.elements){
    
      //Check if this region is crossed
      goesIn = region.toIn(pos.x,pos.y,pos.z,dir.x,dir.y,dir.z,ds2In);
      if(goesIn && ds2In <= lowestDS){

	//Check region triangles
	for(const meshBodyTriangle& triangle : region.elements){      

	  if(triangle.canCross(pos, dir, lowestDS)){
	    
	    double t0;
	    if(triangle.intersect(pos, dir, t0, back)){
	      //If distance t from origin to intersection point
	      //is lower than ds, update ds as the minimum distance
	      if(!std::signbit(t0)){
		if(t0 < lowestDS){
		  lowestDS = t0;
		}
	      }
	    }
            
	  }
	}
      }
    }

  }
  
    
  if(lowestDS < maxDsThres){
    ds = lowestDS + eps;
    return true;
  } else { //No intersection with any triangle, ds has not been updated
    ds = 1.0e35;
    return false;
  }
  
}

bool meshBodyTriangle::intersect(const v3D pos,
				 const v3D dir,
                                 double& t,
				 const bool back) const{
    
  const double eps = crossThreshold;

  //Vectorial product and determinant
  v3D pvec         = dir^edge2;
  const double det = edge1*pvec;

  if(back){ //Testing only backface triangles
    if(det > -eps)
      return false;
  }else{ //Testing only frontface triangles
    if(det < eps)
      return false;
  }

  //Calculate distance between v1 and ray orig
  //v3D aux = pos-v1;
  v3D tvec = pos - v1;
    
  //Calculate barycentric coordinates and check it
  //const double u = (aux*pvec)/det;
  const double u = (tvec*pvec)/det;
        
  if(u < 0.0 || u > 1.0)
    return false;
  
  //aux.crossProd(edge1);
  v3D qvec = tvec^edge1;
  //const double v = (dir*aux)/det;
  const double v = (dir*qvec)/det;

  if(v < 0.0 || u + v > 1.0)
    return false;
        
  //Ray intersects the triangle
  //Calculate t
  //t = (edge2*aux)/det;
  t = (edge2*qvec)/det;

  if(std::signbit(t))
    return false;  
  else
    return true;
}

void meshBodyTriangle::fill(const v3D v1in, const v3D v2in, const v3D v3in){
        
  double fact = 1.0/3.0;
        
  v1 = v1in;
  v2 = v2in;
  v3 = v3in;

  edge1 = v2 - v1;
  edge2 = v3 - v1;  
        
  c.x = fact*(v1in.x + v2in.x + v3in.x);
  c.y = fact*(v1in.y + v2in.y + v3in.y);
  c.z = fact*(v1in.z + v2in.z + v3in.z);
            
  refresh();
}
    
void meshBodyTriangle::refresh(){
  //This function obtain the sphere radious which
  //enclose the triangle.This radious will be set to the
  //larger distance from the center to each vertex
  
  double d1 = (v1 - c).mod2();
  double d2 = (v2 - c).mod2();
  double d3 = (v3 - c).mod2();
            
  r2 = std::max(std::max(d1,d2),d3) + crossThreshold;
            
  r = sqrt(r2);        
}


//
//    pen_meshBodyGeo 
// 


void pen_meshBodyGeo::locate(pen_particleState& state) const{
    
  v3D pos(state.X, state.Y, state.Z);
    
  //Check if is inside the world
  if(!bodies[iworld].inside(pos)){
    //The particle is outside the world
    state.IBODY = getBodies();
    state.MAT = 0;
    return;
  }
    
  //The particle is inside the world. Check daughters
  unsigned currentBody;
  unsigned nextBody = iworld;    
  do{
    currentBody = nextBody;
    //Check if it is inside of some daughter
    for(unsigned i = 0; i < bodies[nextBody].nDaughters; ++i){
      unsigned daugIndex = bodies[nextBody].daughters[i];
      if(bodies[daugIndex].inside(pos)){
	nextBody = daugIndex;
	break;
      }
    }
        
  }while(currentBody != nextBody);
    
  state.MAT = bodies[currentBody].MATER;
  state.IBODY = currentBody;

}

void pen_meshBodyGeo::step(pen_particleState& state, 
                           double DS, 
                           double &DSEF, double &DSTOT, 
                           int &NCROSS) const{

  const double inf = 1.0e36;
    
  v3D pos(state.X, state.Y, state.Z);
  v3D dir(state.U, state.V, state.W);

  //Create local variables for DSEF, DSTOT and NCROSS
  double dsef, dstot;
  dsef = dstot = 0.0;
  int ncross = 0;
                               
  //Check if it is outside the geometry system
  if(state.IBODY >= getBodies()){
    //Is outside. Check if aims to the world
    double dsIn;
    if(bodies[iworld].cross(pos,dir,dsIn,false)){
      //The particle enters the world
      state.IBODY = iworld;
      state.MAT = bodies[iworld].MATER;
      dstot = dsIn;
      ncross = 1;

      move(dsIn,state);

      //If the world is not void, stop the particle
      if(state.MAT != 0){
	DSEF = dsIn;
	DSTOT = dsIn;
	NCROSS = 1;
	return;
      }
    }else{
      //The particle escapes
      state.IBODY = getBodies();
      state.MAT = 0;
      DSEF = inf;
      DSTOT = inf;
      NCROSS = 0; //No interface crossed

      move(inf,state);
      return;            
    }
  }
  

  //The particle is inside the geometry system. It could be in a void region
    
  const unsigned MAT0 = state.MAT;
  bool inVoid = state.MAT == 0 ? true : false;
    
  //The particle is inside the geometry system.
  unsigned currentBody;
  unsigned MATNext = state.MAT;
  unsigned nextBody = state.IBODY;
  double toTravel = DS;
  for(;;){
        
    currentBody = nextBody;
    //Check if the current body or some daughters can be crossed
    const pen_meshBody& body = bodies[currentBody];

    double travel = MATNext == 0 ? inf : toTravel;
    int travelType = 0; // 0 -> Self body, 1 -> To Parent, 2 -> To Daughter

    //Check if the parent can be crossed within the maximum distance
    double ds2Up;
    if(body.cross(pos, dir, ds2Up, true, travel)){
      travelType = 1; //Flag travel as go to parent
      travel = ds2Up;
      nextBody = body.parent;
    }

    //Check if any daughter is closer
    for(unsigned i = 0; i< body.nDaughters; ++i){
      const unsigned iDaugh = body.daughters[i];

      double dsDaugh;
      if(bodies[iDaugh].cross(pos,dir,dsDaugh,false,travel)){
	travelType = 2; //Flag travel as go to daugther
	travel = dsDaugh;
	nextBody = iDaugh;
      }
    }

    //Check the closest type value
    if(travelType == 0){ //Remains in the same body.
      //Move the particle and stop tracking      
      dsef += travel;
      move(travel,dir,pos);
      break;
    }
    else if(travelType == 1){ //Cross the actual body boundary

      if(currentBody == iworld){ //Particle escapes the world
	state.IBODY = getBodies();
	state.MAT = 0;
	if(MAT0 == 0) //Crossed only void regions
	  DSEF = inf;
	else if(MATNext == 0) //World material is void. The travel must not be added to dsef
	  DSEF = dsef;
	else //World material is not void, add the travel to dsef
	  DSEF = dsef + travel;
	DSTOT = inf;
	NCROSS = ncross+1;

	move(inf,state);
	return;
      }

      //Remains in the geometry system. Check overlaps
      solveOverlapsUp(travel,pos,dir,currentBody,nextBody);      
	
    }else{ //Cross some daughter body. Check overlaps
      solveOverlapsDown(travel,pos,dir,nextBody,nextBody);
    }

    //Update maximum remaining distance to travel
    if(MATNext == 0){
      if(MAT0 == 0){
	dsef += travel;
      }else{
	dstot += travel;
      }
    }else{
      dsef += travel;
      toTravel -= travel;
    }

    //Move the particle
    move(travel,dir,pos);

    //Get material of the next body
    MATNext = bodies[nextBody].MATER;
            
    if(MATNext == 0){ //Entering in a void region
      if(!inVoid){
	++ncross;
	inVoid = true;
      }
    } else if(inVoid){
      //Particle comes from a void region to a non void region. Stop it
      ++ncross;
      break;
    }else if(MATNext == MAT0){
      //Entering in a new body with the same material without crossing void regions
      if(bodies[nextBody].KDET != body.KDET){
	//The body belongs to another detector.
	//Stop the tracking
	++ncross;
	break;
      }
    }else{ //Entering in a different material
      //Stop the tracking
      ++ncross;
      break;
    }  
    
  }
    
  //Update particle final state
  DSTOT = dstot + dsef;
  if(MAT0 == 0)
    DSEF = DSTOT;
  else
    DSEF = dsef;
  NCROSS = ncross;
                
  state.X = pos.x;
  state.Y = pos.y;
  state.Z = pos.z;
    
  state.MAT = MATNext;
  state.IBODY = nextBody;

}


int pen_meshBodyGeo::configure(const pen_parserSection& config,
			       const unsigned verbose){
  
  int err;
  //Read input file from configuration
  std::string infilename;
  if(config.read("input-file",infilename) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("meshBodyGeo:configure:Error: 'input-file' field missing at configuration section.\n");
    }
    configStatus = PEN_MESHBODY_GEO_INPUT_SECTION;
    return configStatus;
  }  
    
  //Open input file
  //*****************
  FILE* in = nullptr;
  in = fopen(infilename.c_str(),"r");
  if(in == nullptr){
    if(verbose > 0){
      printf("pen_meshBodyGeo:configure:Error: unable to open input file '%s'.\n", infilename.c_str());
    }
    configStatus = PEN_MESHBODY_GEO_INPUT_SECTION;
    return configStatus;
  }

  
  if(verbose > 0){
    printf("pen_meshBodyGeo:configure:\n");
    printf("                 input filename : %s\n",infilename.c_str());
  }
  
  //Load geometry
  //*****************
  err = GEOMESH(in,verbose);
  if(err != PEN_MESHBODY_GEO_SUCCESS){
    if(verbose > 0){
      printf("pen_meshBodyGeo:configure: Error loading geometry.\n");
      printf("                          Error code: %d\n",err);
    }
    fclose(in);
    configStatus = err;
    return configStatus;
  }
  fclose(in);
  
  //Load dsmax
  //*****************
  std::vector<std::string> bodiesAlias;
  err = config.ls("dsmax",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No dsmax specified for any body\n");
    }
  }
  else{
    printf("dsmax specified for %lu bodies:\n\n",bodiesAlias.size());    
    for(unsigned i = 0; i < bodiesAlias.size(); i++){
      bool found = false;
      std::string key("dsmax/");      
      key += bodiesAlias[i];
      
      for(unsigned j = 0; j < getElements(); j++){
	//Check if body alias is the expected one
	if(bodiesAlias[i].compare(bodies[j].BALIAS) == 0){
	  //Get DSMAX
	  double auxDSmax;
	  err = config.read(key,auxDSmax);
	  if(err != INTDATA_SUCCESS){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error reading 'dsmax' of body %s\n",bodies[j].BALIAS);
	      printf("                     key: %s\n",key.c_str());
	    }
	    configStatus = PEN_MESHBODY_GEO_BAD_READ_DSMAX;
	    return configStatus;
	  }
	  
	  //Check dsmax value
	  if(auxDSmax <= 0.0){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error: 'DSMAX' must be greater than zero.\n");
	      printf("            Specified for body %s: %12.4E\n",bodies[j].BALIAS,auxDSmax);
	    }
	    configStatus = PEN_QUAD_GEO_INVALID_DSMAX;
	    return PEN_QUAD_GEO_INVALID_DSMAX;
	  }
	  //Assign dsmax
	  bodies[j].DSMAX = auxDSmax;
	  if(verbose > 1){
	    printf("Set DSMAX %12.4E for body %s\n",bodies[j].DSMAX,bodies[j].BALIAS);
	  }
	  found = true;
	  break;
	}
      }
      if(!found && verbose > 1){
	printf("Warning: body '%s' not found (key %s).\n",
	       bodiesAlias[i].c_str(),key.c_str());
      }
    }
  }

  if(verbose > 1)
    printf("\n\n");
  
  // Load Detectors
  //*****************
  bodiesAlias.clear();
  err = config.ls("kdet",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No detector specified for any body\n");
    }
  }
  else{
    printf("Detector specified for %lu bodies:\n\n",bodiesAlias.size());    
    for(unsigned i = 0; i < bodiesAlias.size(); i++){
      std::string key("kdet/");      
      key += bodiesAlias[i];

      //Get bodi index
      unsigned bIndex = getIBody(bodiesAlias[i].c_str());
      if(bIndex < getElements()){
	
	//Get KDET
	int auxKDET;
	err = config.read(key,auxKDET);
	if(err != INTDATA_SUCCESS){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error reading 'kdet' for body %s\n",bodies[bIndex].BALIAS);
	    printf("                     key: %s\n",key.c_str());
	  }
	  configStatus = PEN_MESHBODY_GEO_BAD_READ_KDET;
	  return configStatus;
	}

	//Check kdet value
	if(auxKDET <= 0){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: 'KDET' must be greater than zero.\n");
	    printf("            Specified for body %s: %d\n",bodies[bIndex].BALIAS,auxKDET);
	  }
	  configStatus = PEN_QUAD_GEO_INVALID_KDET;
	  return PEN_QUAD_GEO_INVALID_KDET;
	}
	//Assign dsmax
	bodies[bIndex].KDET = (unsigned)auxKDET;
	if(verbose > 1){
	  printf("Set KDET %u for body '%s'\n",bodies[bIndex].KDET,bodies[bIndex].BALIAS);
	}
      }
      else if(verbose > 1){
	printf("Warning: body '%s' not found (key %s).\n",
	       bodiesAlias[i].c_str(),key.c_str());
      }
    }
  }

  // Load Region elements
  //***********************

  // Default value
  int defaultRegionElements;
  err = config.read("DefaultRegionElements",defaultRegionElements);
  if(err != INTDATA_SUCCESS){
    if(verbose > 2){
      printf("No default region size specified\n");
    }
    defaultRegionElements = 40;
  }else{
    if(defaultRegionElements < 1){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error: 'DefaultRegionElements' must "
	       "be greater than zero.\n");
      }
      configStatus = PEN_MESHBODY_GEO_INVALID_REGIONSIZE;
      return PEN_MESHBODY_GEO_INVALID_REGIONSIZE;
    }
  }

  // Specific values
  bodiesAlias.clear();
  err = config.ls("RegionElements",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No region size specified for any body\n");
    }
  }
  else{
    printf("Region size specified for %lu bodies:\n\n",bodiesAlias.size());    
    for(unsigned i = 0; i < bodiesAlias.size(); i++){
      std::string key("RegionElements/");      
      key += bodiesAlias[i];

      //Get bodi index
      unsigned bIndex = getIBody(bodiesAlias[i].c_str());
      if(bIndex < getElements()){
	
	//Get region size
	int regionElements = 0;
	err = config.read(key,regionElements);
	if(err != INTDATA_SUCCESS){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error reading 'RegionElements' for body %s\n",
		   bodies[bIndex].BALIAS);
	    printf("                     key: %s\n",key.c_str());
	  }
	  configStatus = PEN_MESHBODY_GEO_BAD_READ_REGIONSIZE;
	  return configStatus;
	}

	//Check region size value
	if(regionElements < 0){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: 'regionElements' must "
		   "be greater than zero.\n");
	    printf("            Specified for body %s: %d\n",
		   bodies[bIndex].BALIAS,regionElements);
	  }
	  configStatus = PEN_MESHBODY_GEO_INVALID_REGIONSIZE;
	  return PEN_MESHBODY_GEO_INVALID_REGIONSIZE;
	}
	//Assign region size
	bodies[bIndex].meanTrianglesRegion =
	  static_cast<unsigned>(regionElements);
	if(verbose > 1){
	  printf("Set region size %lu for body '%s'\n",
		 bodies[bIndex].meanTrianglesRegion,bodies[bIndex].BALIAS);
	}
      }
      else if(verbose > 1){
	printf("Warning: body '%s' not found (key %s).\n",
	       bodiesAlias[i].c_str(),key.c_str());
      }
    }
  }

  // Load Super Region Elements
  //*****************************

  // Default value
  int defaultSuperRegionElements;
  err = config.read("DefaultSuperRegionElements",defaultSuperRegionElements);
  if(err != INTDATA_SUCCESS){
    if(verbose > 2){
      printf("No default super region size specified\n");
    }
    defaultSuperRegionElements = 20;
  }else{
    if(defaultSuperRegionElements < 1){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error: 'DefaultSuperRegionElements' "
	       "must be greater than zero.\n");
      }
      configStatus = PEN_MESHBODY_GEO_INVALID_REGIONSIZE;
      return PEN_MESHBODY_GEO_INVALID_REGIONSIZE;
    }
  }

  // Specific values  
  bodiesAlias.clear();
  err = config.ls("SuperRegionElements",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No super region size specified for any body\n");
    }
  }
  else{
    printf("Super region size specified for %lu bodies:\n\n",
	   bodiesAlias.size());    
    for(unsigned i = 0; i < bodiesAlias.size(); i++){
      std::string key("SuperRegions/");
      key += bodiesAlias[i];

      //Get bodi index
      unsigned bIndex = getIBody(bodiesAlias[i].c_str());
      if(bIndex < getElements()){
	
	//Get region size
	int superRegionElements = 0;
	err = config.read(key,superRegionElements);
	if(err != INTDATA_SUCCESS){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error reading 'SuperRegions' "
		   "for body %s\n",
		   bodies[bIndex].BALIAS);
	    printf("                     key: %s\n",key.c_str());
	  }
	  configStatus = PEN_MESHBODY_GEO_BAD_READ_REGIONSIZE;
	  return configStatus;
	}

	//Check region size value
	if(superRegionElements < 0){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: 'SuperRegions' must "
		   "be greater than zero.\n");
	    printf("            Specified for body %s: %d\n",
		   bodies[bIndex].BALIAS,superRegionElements);
	  }
	  configStatus = PEN_MESHBODY_GEO_INVALID_REGIONSIZE;
	  return PEN_MESHBODY_GEO_INVALID_REGIONSIZE;
	}
	//Assign region size
	bodies[bIndex].meanRegionsSuperRegion =
	  static_cast<unsigned>(superRegionElements);
	if(verbose > 1){
	  printf("Set super region size %lu for body '%s'\n",
		 bodies[bIndex].meanRegionsSuperRegion,bodies[bIndex].BALIAS);
	}
      }
      else if(verbose > 1){
	printf("Warning: body '%s' not found (key %s).\n",
	       bodiesAlias[i].c_str(),key.c_str());
      }
    }
  }

  //Split regions

#ifdef _PEN_USE_THREADS_
  
  unsigned int nSplitThreads =
    std::max(static_cast<unsigned int>(2),
	     std::thread::hardware_concurrency());
  std::vector<std::thread> splitThreads;

  std::atomic<unsigned int> atomicCount{0};

  if(verbose > 2){
    printf(" * Using %u threads to split bodies in regions\n", nSplitThreads);
  }
  
  for(size_t ith = 0; ith < nSplitThreads; ++ith){
    splitThreads.push_back(std::thread([&,ith](){
      
      unsigned int ibody = atomicCount++;
      while(ibody < this->getElements()){

	if(verbose > 2){
	  printf("   - Body %u split begins in split thread %lu\n",
		 ibody, static_cast<unsigned long>(ith));
	  fflush(stdout);
	}

	//Split the body regions
	pen_meshBody& body = this->bodies[ibody];
    
	//Check if a region size has been specified for this body
	if(body.meanTrianglesRegion == 0){
	  body.meanTrianglesRegion =
	    std::max(static_cast<unsigned long>(defaultRegionElements),
		     body.nTriangles/static_cast<unsigned long>(pen_meshBody::MAX_REGIONS/2));
	}

	if(body.meanRegionsSuperRegion == 0){
	  body.meanRegionsSuperRegion = defaultSuperRegionElements;
	}

	//Check if this body requires more regions
	pen_meshBody::superRegion& superRegion = body.regions[0];

	pen_meshBody::triangleRegion::splitUntil(body.meanTrianglesRegion,
						 pen_meshBody::crossThreshold,
						 superRegion.elements,
						 pen_meshBody::MAX_REGIONS);

	//Enlarge regions to avoid rounding errors
	//for(pen_meshBody::triangleRegion& region : superRegion.elements){
	//  region.enlarge(pen_meshBody::crossThreshold);
	//}

	//Check total number of triangles
	unsigned long nTrianglesAux = 0;
	for(const pen_meshBody::triangleRegion& region : superRegion.elements){
	  nTrianglesAux += region.nElements();
	}
	if(nTrianglesAux != body.nTriangles){
	  printf("pen_meshBodyGeo:configure: Error: Triangles lost on first "
		 "region split in body '%s.'\n"
		 "      Expected triangles : %lu\n"
		 "      Remaining triangles: %lu\n"
		 " Please, report this issue.\n",
		 body.BALIAS, body.nTriangles, nTrianglesAux);
	    
	  fflush(stdout);
	  throw std::range_error("Lost triangles");
	}
	    
	pen_meshBody::superRegion::splitUntil(body.meanRegionsSuperRegion,
					      pen_meshBody::crossThreshold,
					      body.regions,
					      pen_meshBody::MAX_SUP_REGIONS);

	//Enlarge super regions to avoid rounding errors
	//for(pen_meshBody::superRegion& supRegion : body.regions){
	//  supRegion.enlarge(pen_meshBody::crossThreshold);
	//}

	//Check total number of triangles
	nTrianglesAux = 0;
	for(const pen_meshBody::superRegion& supRegion : body.regions){
	  for(const pen_meshBody::triangleRegion& region : supRegion.elements){
	    nTrianglesAux += region.nElements();
	  }
	}
	if(nTrianglesAux != body.nTriangles){
	  printf("pen_meshBodyGeo:configure: Error: Triangles lost on "
		 "super-region split in body '%s.'\n"
		 "      Expected triangles : %lu\n"
		 "      Remaining triangles: %lu\n"
		 " Please, report this issue.\n",
		 body.BALIAS, body.nTriangles, nTrianglesAux);
	  fflush(stdout);
	  throw std::range_error("Lost triangles");
	}

	if(verbose > 2){
	  printf("   + Body %u split completed in thread %lu\n",
		 ibody, static_cast<unsigned long>(ith));
	  fflush(stdout);	  
	}
	
	ibody = atomicCount++;	
      }
      
    }));
  }

  //Wait until all threas have been finished
  for(std::thread& t : splitThreads){
    t.join();
  }

#else

  for(std::size_t ibody = 0; ibody < getElements(); ++ibody){

	if(verbose > 2){
	  printf("   - Body %u split begins\n",ibody);
	  fflush(stdout);
	}

	//Split the body regions
	pen_meshBody& body = bodies[ibody];
    
	//Check if a region size has been specified for this body
	if(body.meanTrianglesRegion == 0){
	  body.meanTrianglesRegion =
	    std::max(static_cast<unsigned long>(40),
		     body.nTriangles/static_cast<unsigned long>(pen_meshBody::MAX_REGIONS/2));
	}

	if(body.meanRegionsSuperRegion == 0){
	  body.meanRegionsSuperRegion = 20;
	}

	//Check if this body requires more regions
	pen_meshBody::superRegion& superRegion = body.regions[0];

	pen_meshBody::triangleRegion::splitUntil(body.meanTrianglesRegion,
						 pen_meshBody::crossThreshold,
						 superRegion.elements,
						 pen_meshBody::MAX_REGIONS);

	pen_meshBody::superRegion::splitUntil(body.meanRegionsSuperRegion,
					      pen_meshBody::crossThreshold,
					      body.regions,
					      pen_meshBody::MAX_SUP_REGIONS);

	if(verbose > 2){
	  printf("   + Body %u split completed\n",ibody);
	  fflush(stdout);
	}
    
  }
  
#endif
  
  // Load Absortion Energies
  //*************************
  bodiesAlias.clear();
  err = config.ls("eabs",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No absortion energies specified for any body\n");
    }
  }
  else{
    printf("Absortion energies specified for %lu bodies:\n\n",bodiesAlias.size());
    for(unsigned i = 0; i < bodiesAlias.size(); i++){
      std::string key("eabs/");      
      key += bodiesAlias[i];

      //Get body index
      unsigned bIndex = getIBody(bodiesAlias[i].c_str());
      if(bIndex < getElements()){
	
	//Get particle names
	std::vector<std::string> particleNames;
	config.ls(key.c_str(),particleNames);

	//Get particle eabs
	double bodyEABS[constants::nParTypes];
	for(unsigned j = 0; j < constants::nParTypes; j++)
	  bodyEABS[j] = -1.0;
	  
	for(unsigned j = 0; j < particleNames.size(); j++){

	  unsigned kpar = particleID(particleNames[j].c_str());
	  if(kpar >= ALWAYS_AT_END){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error on 'eabs' field, unknown particle '%s' on body '%s'.\n",particleNames[j].c_str(),bodiesAlias[i].c_str());
	    }
	    return PEN_MESHBODY_GEO_UNKNOWN_PARTICLE;
	  }
	  
	  std::string key2 = key + std::string("/") + particleNames[j];
	  double eabs;
	  err = config.read(key2,eabs);
	  if(err != INTDATA_SUCCESS){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error reading energy absortion at field '%s'. Double expected.\n",key2.c_str());
	    }
	    return PEN_MESHBODY_GEO_BAD_READ_EABS;
	  }

	  if(eabs <= 0.0){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error: Invalid energy absortion %12.4E for body '%s' particle '%s'. Must be greater than zero.\n",eabs,bodiesAlias[i].c_str(),particleNames[j].c_str());
	    }
	    return PEN_MESHBODY_GEO_INVALID_EABS;
	  }
	  
	  bodyEABS[kpar] = eabs;
	}

	//Set body eabs for each specified particle
	if(setBodyEabs(bIndex,bodyEABS) != 0){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error on 'eabs' field, unknown body '%s'\n",bodiesAlias[i].c_str());
	  }
	  return PEN_MESHBODY_GEO_UNDEF_BODY_LABEL;	  
	}
	
	if(verbose > 1){
	  printf("Absotion energies (eV) for body '%s':\n",bodies[bIndex].BALIAS);
	  for(unsigned j = 0; j < constants::nParTypes; j++){
	    printf(" %20.20s: %14.5E\n",particleName(j),bodies[bIndex].localEABS[j]);
	  }
	}
      }
      else if(verbose > 1){
	printf("Warning: body '%s' not found (key %s).\n",
	       bodiesAlias[i].c_str(),key.c_str());
      }
    }
  }
  
  if(verbose > 1){
    printf("\nBodies information:\n\n");

    for(unsigned j = 0; j < getElements(); j++){
      printf("  Body '%s':\n",bodies[j].BALIAS);
      printf("    Material: %u\n",bodies[j].MATER);
      printf("    DSMAX   : %10.4E\n",bodies[j].DSMAX);
      printf("    KDET    : %u\n",bodies[j].KDET);
      printf("    EABS    :\n");
      for(unsigned k = 0; k < constants::nParTypes; k++)
	printf("  %20.20s: %14.5E\n",particleName(k),bodies[j].localEABS[k]);
    }
    printf("\n");
  }

  //Read report geometry file
  std::string reportFilename;
  if(config.read("report-file",reportFilename) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No report file selected.\n");
    }
  }else{
    FILE* freport = nullptr;
    freport = fopen(reportFilename.c_str(), "w");

    if(freport == nullptr){
      printf("Unable to open report file '%s'\n", reportFilename.c_str());
    }else{

      //Print bodies information
      for(unsigned ib = 0; ib < getElements(); ++ib){

	const pen_meshBody& body = bodies[ib];

	//Print generic information
	fprintf(freport,"  Body '%s':\n",body.BALIAS);
	fprintf(freport,"    Material: %u\n",body.MATER);
	fprintf(freport,"    DSMAX   : %10.4E\n",body.DSMAX);
	fprintf(freport,"    KDET    : %u\n",body.KDET);
	fprintf(freport,"    EABS    :\n");
	for(unsigned k = 0; k < constants::nParTypes; k++)
	  fprintf(freport,"  %20.20s: %14.5E\n",
		  particleName(k),body.localEABS[k]);

	//Print daughters information
	fprintf(freport,"\n\n    Parent: %s %s\n",
		body.PALIAS,
		body.canOverlapParent ? "(Overlaps)" : "(No overlaps)");
	fprintf(freport,"    Body daughters: %u\n",
		body.nDaughters);
	for(unsigned id = 0; id < body.nDaughters; ++id){
	  const pen_meshBody& daughter = bodies[body.daughters[id]];
	  fprintf(freport,"     - %s: %s\n",
		  daughter.BALIAS,
		  daughter.canOverlapParent ? "Overlaps" : "No Overlaps");
	}
	
	//Print sister overlaps
	fprintf(freport,"\n    Body sister overlaps: %u\n",body.nOverlap);
	for(unsigned is = 0; is < body.nOverlap; ++is){
	  fprintf(freport,"     - %s\n",bodies[body.overlapedBodies[is]].BALIAS);
	}

	unsigned long nRegions = 0;
	for(unsigned isr = 0; isr < body.nSupRegions(); ++isr){
	  const pen_meshBody::superRegion& supReg = body.regions[isr];
	  nRegions += supReg.nElements();
	}
	
	//Print boundary and regions information
	fprintf(freport,"\n\n    Body boundary: %s\n",
		body.boundary.stringify().c_str());
	fprintf(freport,"      Number of super regions          : %lu\n",
		static_cast<unsigned long>(body.nSupRegions()));
	fprintf(freport,"      Number of regions                : %lu\n",
		nRegions);
	fprintf(freport,"      Number of triangles              : %lu\n",
		body.nTriangles);	
	fprintf(freport,"      Expected triangles per region    : %lu\n",
		body.meanTrianglesRegion);
	fprintf(freport,"      Expected regions per super region: %lu\n",
		body.meanRegionsSuperRegion);

	for(unsigned isr = 0; isr < body.nSupRegions(); ++isr){

	  const pen_meshBody::superRegion& supReg = body.regions[isr];
	  
	  fprintf(freport,"\n      * Super Region %u: %s\n", isr,
		  supReg.stringify().c_str());
	  fprintf(freport,"           Number regions  : %lu\n",
		  static_cast<unsigned long>(supReg.nElements()));

	  for(unsigned ir = 0; ir < supReg.nElements(); ++ir){
	  
	    const pen_meshBody::triangleRegion& region = supReg.elements[ir];
	  
	    fprintf(freport,"     - Region %u: %s\n", ir,
		    region.stringify().c_str());
	    fprintf(freport,"        Triangles: %lu\n", region.nElements());
	    for(const auto& t : region.elements){
	      fprintf(freport,"         %s\n",t.stringify().c_str());
	    }
	  }
	}
      }
      
      fclose(freport);
    }
    
  }

  //Check if the user has requested testing the mesh
  bool toTest = false;
  if(config.read("test-mesh",toTest) == INTDATA_SUCCESS){

    if(toTest){

      if(verbose > 1){
	printf(" ** Testing mesh.\n\n");

	printf("  + Body to body intersections:\n");
      }

      bool intersectionFound = false;
      for(unsigned ibody = 0; ibody < getElements(); ++ibody){
	  
	const pen_meshBody& body = bodies[ibody];

	if(verbose > 1)
	  printf("    - Body %s (%u):\n",body.BALIAS,ibody);
	  
	bool localIntersectionFound = false;
	  
	//Check only overlapping bodies

	//Parent
	if(body.canOverlapParent && body.parent != ibody){
	  const pen_meshBody& parent = bodies[body.parent];

	  //All points should be inside the parent, check it
	  for(const pen_meshBody::superRegion& supRegion : body.regions){
	    for(const pen_meshBody::triangleRegion& region : supRegion.elements){
	      for(const meshBodyTriangle& triangle : region.elements){

		//Check if the three triangle vertex are inside the parent
		if(!parent.inside(triangle.readV1()) ||
		   !parent.inside(triangle.readV2()) ||
		   !parent.inside(triangle.readV3())){
		  //Intersection found
		  if(verbose > 2){
		    printf("      - Intersection with parent '%s' (%u): %s\n",
			   body.PALIAS,body.parent,triangle.stringify().c_str());
		  }
		  else if(verbose > 1){
		    printf("      - Intersection with parent '%s' (%u)\n",
			   body.PALIAS,body.parent);
		  }
		  localIntersectionFound = true;
		  intersectionFound = true;
		  if(verbose <= 2)
		    break;
		}
		  
	      }
	      if(verbose <= 2 && localIntersectionFound) break;
	    }
	    if(verbose <= 2 && localIntersectionFound) break;
	  }
	}

#ifdef _PEN_USE_THREADS_
  
	unsigned int nOverlapThreads =
	  std::max(static_cast<unsigned int>(2),
		   std::thread::hardware_concurrency());
	std::vector<std::thread> overlapThreads;

	atomicCount = 0;
	std::atomic<bool> sharedIntersect{false};

	for(size_t ith = 0; ith < nOverlapThreads; ++ith){
	  overlapThreads.push_back(std::thread([&,ith](){

	    unsigned int iover = atomicCount++;
	    while(iover < body.nOverlap){
	      const unsigned overlapBodyIndex = body.overlapedBodies[iover];
	      const pen_meshBody& overlapBody = bodies[overlapBodyIndex];

	      //All points should be inside the parent, check it
	      bool overlapIntersectionFound = false;
	      for(const pen_meshBody::superRegion& supRegion : body.regions){
		for(const pen_meshBody::triangleRegion& region : supRegion.elements){
		  for(const meshBodyTriangle& triangle : region.elements){

		    //Check if the three triangle vertex are outside the overlapping body
		    if(overlapBody.inside(triangle.readV1()) ||
		       overlapBody.inside(triangle.readV2()) ||
		       overlapBody.inside(triangle.readV3())){
		      //Intersection found
		      if(verbose > 2){
			printf("      - Intersection with overlapping "
			       "body '%s' (%u): %s\n",
			       overlapBody.BALIAS,
			       overlapBodyIndex,
			       triangle.stringify().c_str());			
		      }
		      else if(verbose > 1){
			printf("      - Intersection with overlapping body '%s' (%u)\n",
			       overlapBody.BALIAS,overlapBodyIndex);
		      }
		      sharedIntersect = true;
		      overlapIntersectionFound = true;
		      if(verbose <= 2)
			break;
		    }
		  }
		  if(verbose <= 2 && overlapIntersectionFound) break;
		}
		if(verbose <= 2 && overlapIntersectionFound) break;
	      }

	      iover = atomicCount++;	      
	    }
      
	  }));
	}

	//Wait until all threas have been finished
	for(std::thread& t : overlapThreads){
	  t.join();
	}	

	if(sharedIntersect){
	  localIntersectionFound = true;
	  intersectionFound = true;    
	}
      
  
#else
	
	//Overlap bodies
	for(unsigned iover = 0; iover < body.nOverlap; ++iover){
	  const unsigned overlapBodyIndex = body.overlapedBodies[iover];
	  const pen_meshBody& overlapBody = bodies[overlapBodyIndex];

	  //All points should be inside the parent, check it
	  bool overlapIntersectionFound = false;
	  for(const pen_meshBody::superRegion& supRegion : body.regions){
	    for(const pen_meshBody::triangleRegion& region : supRegion.elements){
	      for(const meshBodyTriangle& triangle : region.elements){

		//Check if the three triangle vertex are outside the overlapping body
		if(overlapBody.inside(triangle.readV1()) ||
		   overlapBody.inside(triangle.readV2()) ||
		   overlapBody.inside(triangle.readV3())){
		  //Intersection found
		  if(verbose > 2){
		    printf("      - Intersection with overlapping "
			   "body '%s' (%u): %s\n",
			   overlapBody.BALIAS,
			   overlapBodyIndex,
			   triangle.c_str());			
		  }
		  else if(verbose > 1){
		    printf("      - Intersection with overlapping body '%s' (%u)\n",
			   overlapBody.BALIAS,overlapBodyIndex);
		  }
		  localIntersectionFound = true;
		  intersectionFound = true;
		  overlapIntersectionFound = true;
		  if(verbose <= 2)
		    break;
		}
	      }
	      if(verbose <= 2 && overlapIntersectionFound) break;
	    }
	    if(verbose <= 2 && overlapIntersectionFound) break;
	  }
	}
	
#endif

	if(!localIntersectionFound){
	  if(verbose > 1)
	    printf("      + No intersections found\n");
	}
      }

      if(intersectionFound){
	configStatus = PEN_MESHBODY_GEO_BODY_INTERSECTIONS_FOUND;
	return configStatus;
      }

    }
  }
  
  configStatus = PEN_MESHBODY_GEO_SUCCESS;
  return configStatus;
}

int pen_meshBodyGeo::GEOMESH(FILE* in,const unsigned verbose){
  //Read comment lines (start with #)
  char line[50000];
  unsigned long nlines;
  unsigned long nRead = 0;
  int nBodies;
    
  if(pen_getLine(in,50000,line,nlines) == 0)
    {
      nRead += nlines;
      //Read number of bodies in the geometry file
      if(sscanf(line," %d",&nBodies) != 1){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error reading number of objects."
		 "Unexpected format in line %lu: \n %s\n", nRead, line);
	}
	return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
      }

      if(nBodies <= 0){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: number of bodies must be greater than zero.\n");
	  printf("        Number of bodies read is %d. \n",nBodies);
	}
	return PEN_MESHBODY_GEO_INVALID_NBODIES;
      }
        
      NBODYS = static_cast<unsigned int>(nBodies);
        
      if(NBODYS > NB){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: number of bodies read from the geometry file: %u, is greater than the maximum number of bodies allowed: %u\n.",NBODYS, NB);
	}
	return PEN_MESHBODY_GEO_INVALID_NBODIES;
      }

            
      if(verbose > 1)
        {
	  printf("Number of bodies read: %u\n",NBODYS);
	  printf("\n");
        }
    }
  
    
    
  //Create vector to save bodies vertex info
  std::vector<v3D>vertex;
  unsigned long nvertex;
    
  if(verbose > 1){
    printf("Bodies mesh information:\n");
  }
    
  for(size_t i=0; i<NBODYS; ++i){
    if(pen_getLine(in,50000,line,nlines) == 0){
      nRead += nlines;
      //Read characteristics of each body in the geometry file
      //Create auxiliar variable to ensure material is greater than zero
      int mat;
      if(sscanf(line," %d  %lu  %lu  %s   %s", &mat, &bodies[i].nTriangles, &nvertex, bodies[i].BALIAS, bodies[i].PALIAS) != 5){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error reading object header information.\n"
		 "Unexpected format in line %lu: \n %s\n", nRead, line);
	  printf("Expected format is #MAT      #NFACES     #NVERTEX     #NAME        #PARENT NAME\n");
	}
	return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
      }
      if(mat < 0){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: 'material' must be greater than zero.\n");
	  printf("            Specified for body %ld: %d\n",i,mat);
	}
	return PEN_MESHBODY_GEO_INVALID_MAT;
      }
      else
	bodies[i].MATER = static_cast<unsigned>(mat);
                
      if(bodies[i].MATER > constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: material number read from the geometry file: %u, is greater than the maximum number of materials allowed: %u\n", bodies[i].MATER, constants::MAXMAT);
	}
	return PEN_MESHBODY_GEO_INVALID_MAT;
      }
        
            
      //Resize vector of vertex with the number of vertex
      //read for each body from the geometry file
      vertex.resize(nvertex);
            
      if(vertex.size() != nvertex){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error allocating memory for vertex reading. Can't allocate memory for %lu vertex.\n", nvertex);
	}
	return PEN_MESHBODY_BAD_MEMORY_ALLOCATION;
      }

      for(size_t j=0; j<nvertex; ++j){
	if(pen_getLine(in,50000,line,nlines) == 0){
	  int index;
	  double x,y,z;
	  nRead += nlines;
	  //Read vertex of each body
	  if(sscanf(line," %d  %le  %le  %le", &index, &x, &y, &z) != 4){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error reading vertex information."
		     "Unexpected format in line %lu: \n %s\n", nRead, line);
	    }
	    return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
	  }
	  
	  if(index < 0){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error: vertex index must be positive.\n");
	    }
	    return PEN_MESHBODY_GEO_INVALID_VERTEX_INDEX;
	  }

	  //Save vertex
	  vertex[index].x = x;
	  vertex[index].y = y;
	  vertex[index].z = z;	  
                    
	  //Save minimum and maximum value of each coordenate
	  if(j==0){
	    //Set mins and maxs as the first read vertex coordenates
	    bodies[i].boundary.set(x,y,z,x,y,z);
	  } else {
	    //Enlarge boundary region
	    bodies[i].boundary.enlarge(vertex[index]);
	  }
	}
      }

      //Enlarge boundary to avoid rounding errors
      bodies[i].boundary.enlarge(threshold);

      //Create a initial region for this body
      pen_meshBody::triangleRegion region(bodies[i].boundary);
      //Resize triangle vector
      std::vector<meshBodyTriangle>& triangles = region.elements;
      triangles.resize(bodies[i].nTriangles);
            
      for(size_t k=0; k < bodies[i].nTriangles; ++k){
	unsigned int index[3];
	if(pen_getLine(in,50000,line,nlines) == 0){
	  nRead += nlines;
	  //Read each triangle or face of each body 
	      if(sscanf(line," %u  %u  %u", &index[0], &index[1], &index[2]) != 3){
		if(verbose > 0){
		  printf("pen_meshBodyGeo:configure: Error reading triangle faces."
			 "Unexpected format in line %lu: \n %s\n", nRead, line);
		}
		return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
	      }
                    
	  //Fill triangle
	  triangles[k].fill(vertex[index[0]],vertex[index[1]],vertex[index[2]]);
	}
      }

      //Add the region to the regions vector of the first super region
      bodies[i].regions.emplace_back(region);
      bodies[i].regions[0].elements.push_back(region);
    }
        
    if(verbose > 1){
      printf("Body number %lu\n"
	     "Body alias: %s\n"
	     "Parent alias: %s\n"
	     "Material read: %u\n"
	     "Number of triangles read: %lu\n"
	     "Number of vertex read: %lu\n",i, bodies[i].BALIAS, bodies[i].PALIAS, bodies[i].MATER, bodies[i].nTriangles, nvertex);

      printf("\n");
      //Print minimum and maximum coordenates of the read body
           
      printf("Body limits:\n"
	     "xmin = %15.8E    xmax = %15.8E\n"
	     "ymin = %15.8E    ymax = %15.8E\n"
	     "zmin = %15.8E    zmax = %15.8E\n", 
	     bodies[i].boundary.minx(),
	     bodies[i].boundary.maxx(),
	     bodies[i].boundary.miny(),
	     bodies[i].boundary.maxy(),
	     bodies[i].boundary.minz(),
	     bodies[i].boundary.maxz());
    
      printf("\n");
            
           
    }
  }
    
  //Assign index to each parent alias
  for(size_t i=0; i<NBODYS; ++i){
    if(strcmp(bodies[i].PALIAS,"VOID") == 0 ||
       strcmp(bodies[i].PALIAS,"void") == 0){
      if(worldFound){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: only a single world is allowed.\n"
		 "Defined worlds are '%s' and '%s'\n",bodies[iworld].BALIAS,bodies[i].BALIAS);
	}
	return PEN_MESHBODY_MULTIPLE_WORLDS;
      }
      else{
	bodies[i].parent = i;
	iworld = i;
	worldFound = true;
      }
    }
    else{
      bodies[i].parent = getIBody(bodies[i].PALIAS);
    }
  }
    
  //Check if world is found
  if(!worldFound){
    if(verbose > 0){
      printf("pen_meshBodyGeo:configure: Error: world not found.\n Expected one object with 'VOID' or 'void' parent name.\n");
    }
    return PEN_MESHBODY_WORLD_NOT_FOUND;
  }
    
  //Assign daughter information to each parent body
  for(size_t i=0; i<NBODYS; ++i){
    if(i != iworld){ //Avoid world to assign itself as child
      bodies[bodies[i].parent].addDaughter(i);
    }
  }
    
  //Construct possible collisions
  checkCross(iworld);

  if(verbose > 1){
    for(size_t i=0; i<NBODYS; ++i){
      printf("Body %lu (%s)\n", i, bodies[i].BALIAS);
      if(bodies[i].canOverlapParent){
	printf("overlaps with her parent %u (%s).\n", bodies[i].parent, bodies[i].PALIAS);
      }
      else{
	printf("no overlaps with her parent %u (%s).\n", bodies[i].parent, bodies[i].PALIAS);
      }
            
      if(bodies[i].nOverlap > 0){
	for(size_t j=0; j<bodies[i].nOverlap; ++j){
	  long unsigned int iOverlap = bodies[i].overlapedBodies[j];
	  printf("overlaps with her sister: %lu (%s).\n",
		 iOverlap, bodies[iOverlap].BALIAS);
	}
  
      }
      else{
	printf("no overlaps with any sister.\n");
      }

      printf("\n");
    }
  }
  if(verbose > 2){
    printf(" * Geometry tree:\n");
        
    bodies[iworld].printDaughters(1,bodies);
  }
    
  return PEN_MESHBODY_GEO_SUCCESS;
}

bool pen_meshBodyGeo::canOverlapParent(const unsigned ibody) const {
    
  //Check if the limits of the two bodies overlap
  const pen_meshBody& body = bodies[ibody];
  return body.boundary.commonBoundary(bodies[body.parent].boundary, threshold);
}

bool pen_meshBodyGeo::canOverlap(const unsigned ibody1,
				 const unsigned ibody2)const {

  //Get bodies references
  const pen_meshBody& body1 = bodies[ibody1];
  const pen_meshBody& body2 = bodies[ibody2];
    
  //Calculate lateral semi-sizes of boxes in each axis
  double dx1 = (body1.boundary.dx())/2.0;
  double dy1 = (body1.boundary.dy())/2.0;
  double dz1 = (body1.boundary.dz())/2.0;

  double dx2 = (body2.boundary.dx())/2.0;
  double dy2 = (body2.boundary.dy())/2.0;
  double dz2 = (body2.boundary.dz())/2.0;
    
  //Calculate boxes centers
  double cx1 = body1.xmin() + dx1;
  double cy1 = body1.ymin() + dy1;
  double cz1 = body1.zmin() + dz1;

  double cx2 = body2.xmin() + dx2;
  double cy2 = body2.ymin() + dy2;
  double cz2 = body2.zmin() + dz2;

  //Check if can overlap
  if( fabs(cx1-cx2) - (dx1+dx2) > threshold || 
      fabs(cy1-cy2) - (dy1+dy2) > threshold ||
      fabs(cz1-cz2) - (dz1+dz2) > threshold){
    return false;
  }

  return true;
}

void pen_meshBodyGeo::checkCross(const unsigned iparent){
    
  //Get parent reference
  const pen_meshBody& parent = bodies[iparent];
    
  for(unsigned i = 0; i < parent.nDaughters; ++i){
        
    //Save daughter index
    const unsigned idaugh1 = parent.daughters[i];
    //Get body reference
    pen_meshBody& body1 = bodies[idaugh1];
        
    //Check if this daughter can cross with their parent
    if(canOverlapParent(idaugh1)){
      body1.canOverlapParent = true;
    }else{
      body1.canOverlapParent = false;
    }
        
    //Check crosses with other daughters (sisters)
    for(unsigned j = i+1; j < parent.nDaughters; ++j){
            
      //Save second daughter index
      const unsigned idaugh2 = parent.daughters[j];
      //Get body reference
      pen_meshBody& body2 = bodies[idaugh2];
            
      if(canOverlap(idaugh1,idaugh2)){
	body1.overlapedBodies[body1.nOverlap++] = idaugh2;
	body2.overlapedBodies[body2.nOverlap++] = idaugh1;
      }
    }
        
    //Propagate cross checks to daughters
    checkCross(idaugh1);
  }    
}

REGISTER_GEOMETRY(pen_meshBodyGeo,MESH_BODY)
