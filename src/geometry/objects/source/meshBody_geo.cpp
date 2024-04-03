
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

  //Load vertex transformations 
  //******************************
  std::map<std::string, std::vector<pen_meshTransform::group>> transMap;

  //Get bodies with defined transformations
  std::vector<std::string> transBodyNames;
  config.ls("transforms/",transBodyNames);

  //Iterate over each body with used vertex group
  for(size_t itransbody = 0; itransbody < transBodyNames.size(); ++itransbody){

    //Read transform section for this body
    std::string transBodyKey = std::string("transforms/") + transBodyNames[itransbody];
    pen_parserSection transBodySec;
    if(config.readSubsection(transBodyKey,transBodySec) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error reading '%s',"
	       " section expected.\n", transBodyKey.c_str()); 
      }
      return PEN_MESHBODY_GEO_INPUT_SECTION;
    }

    //Create an entry in the map for this body and save the reference 
    std::vector<pen_meshTransform::group>& bodyTransformGroups =
      transMap[transBodyNames[itransbody]];
    
    // * Read vertex groups names
    std::vector<std::string> transGroupNames;
    transBodySec.ls(transGroupNames);

    //Resize the number of vertex groups according to the number of transformation group names
    bodyTransformGroups.resize(transGroupNames.size());
    
    //Create a vector to flag the used transformations positions
    std::vector<bool> usedTransGroup(transGroupNames.size(), false);

    //Iterate over transformation groups to configure them
    for(size_t iTransGroup = 0; iTransGroup < transGroupNames.size(); ++iTransGroup){

      //Read group section
      pen_parserSection transGroupSec;
      if(transBodySec.readSubsection(transGroupNames[iTransGroup],
				     transGroupSec) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error reading '%s/%s',"
		 " section expected.\n",
		 transBodyKey.c_str(),transGroupNames[iTransGroup].c_str()); 
	}
	return PEN_MESHBODY_GEO_INPUT_SECTION;
      }

      //Read the group index
      int transGroupPos;
      if(transGroupSec.read("index", transGroupPos) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error reading '%s/%s/index'."
		 " Integer expected.\n", transBodyKey.c_str(),
		 transGroupNames[iTransGroup].c_str());
	}
	return PEN_MESHBODY_GEO_INVALID_VERTEX_GROUP_INDEX;
      }

      //Check it
      if(transGroupPos < 0 ||
	 transGroupPos > static_cast<int>(bodyTransformGroups.size())){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: Invalid transformation index"
		 " for group '%s' in body '%s'.\n"
		 "   Expected in range [0,%u), provided %d\n",
		 transGroupNames[iTransGroup].c_str(), transBodyNames[itransbody].c_str(),
		 static_cast<unsigned>(bodyTransformGroups.size()), transGroupPos);
	}
	return PEN_MESHBODY_GEO_INVALID_VERTEX_GROUP_INDEX;
      }

      //Ensure this position is unused
      if(usedTransGroup[transGroupPos]){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: Assigned vertex group index"
		 " %d for group '%s' in body '%s' already used.\n",
		 transGroupPos, transGroupNames[iTransGroup].c_str(),
		 transBodyNames[itransbody].c_str());
	}
	return PEN_MESHBODY_GEO_INVALID_VERTEX_GROUP_INDEX;	  
      }

      //Read the vertex group where the transforms are applied
      std::string vertexGroup;
      if(transGroupSec.read("vertex-group", vertexGroup) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error reading '%s/%s/vertex-group'."
		 " String expected.\n", transBodyKey.c_str(),
		 transGroupNames[iTransGroup].c_str());
	}
	return PEN_MESHBODY_GEO_INVALID_VERTEX_GROUP_INDEX;
      }

      //Set this position as used and save the group name
      usedTransGroup[transGroupPos] = true;
      bodyTransformGroups[transGroupPos].name.assign(vertexGroup);

      //Get a reference to this group
      pen_meshTransform::group& transG = bodyTransformGroups[transGroupPos];

      //Read transform names in the group
      std::vector<std::string> transformNames;
      transGroupSec.ls("transforms",transformNames);

      //Resize transformation vectors
      transG.resize(transformNames.size());

      //Create a vector to flag used transformation positions
      std::vector<bool> usedTrans(transformNames.size(),false);

      for(size_t it = 0; it < transformNames.size(); ++it){

	//Get transformation section
	std::string transKey = std::string("transforms/") + transformNames[it];
	pen_parserSection transSec;
	if(transGroupSec.readSubsection(transKey, transSec) != INTDATA_SUCCESS){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error reading transform section '%s/%s/%s',"
		   " section expected.\n", transBodyKey.c_str(),
		   transGroupNames[iTransGroup].c_str(),
		   transKey.c_str()); 
	  }
	  return PEN_MESHBODY_GEO_INPUT_SECTION;	  
	}

	//Read transformation index
	int transIndex;
	if(transSec.read("index", transIndex) != INTDATA_SUCCESS){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error reading transformation %s 'index'."
		   " Integer expected.\n", transformNames[it].c_str());
	  }
	  return PEN_MESHBODY_GEO_INVALID_TRANSFORMATION_INDEX;
	}

	//Check it
	if(transIndex < 0 ||
	   transIndex > static_cast<int>(transformNames.size())){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: Invalid transformation '%s' index"
		   " for vertex group '%s' in body '%s'.\n"
		   "   Expected in range [0,%u), provided %d\n",
		   transformNames[it].c_str(), transGroupNames[iTransGroup].c_str(),
		   transBodyNames[itransbody].c_str(),
		   static_cast<unsigned>(transformNames.size()),transIndex);
	  }
	  return PEN_MESHBODY_GEO_INVALID_TRANSFORMATION_INDEX;
	}

	//Ensure this position is unused
	if(usedTrans[transIndex]){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: Assigned index %d"
		   " for transformation %s, in vertex group '%s', in body '%s' already used.\n",
		   transIndex, transformNames[it].c_str(),
		   transGroupNames[iTransGroup].c_str(),
		   transBodyNames[itransbody].c_str());
	  }
	  return PEN_MESHBODY_GEO_INVALID_TRANSFORMATION_INDEX;	  
	}

	//Set this position as used
	usedTrans[transIndex] = true;

	//Read the transformation type
	std::string transType;
	if(transSec.read("type", transType) != INTDATA_SUCCESS){
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: Unable to read type for "
		   "transformation %s, in vertex group '%s', in body '%s'. String expected.\n",
		   transformNames[it].c_str(), transGroupNames[iTransGroup].c_str(),
		   transBodyNames[itransbody].c_str());
	  }
	  return PEN_MESHBODY_GEO_INVALID_TRANSFORMATION_TYPE;	  
	}

	if(transType.compare("TRANSLATION")   == 0 ||
	   transType.compare("TRANSLATION_X") == 0 ||
	   transType.compare("TRANSLATION_Y") == 0 ||
	   transType.compare("TRANSLATION_Z") == 0){

	  //Is a translation, read the translated distance

	  double ds;
	  if(transSec.read("ds", ds) != INTDATA_SUCCESS){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error: Unable to read translation "
		     "distance 'ds' for translation %s, in vertex group '%s', in body '%s'."
		     " Double expected.\n", transformNames[it].c_str(),
		     transGroupNames[iTransGroup].c_str(),
		     transBodyNames[itransbody].c_str());	      
	    }
	    return PEN_MESHBODY_GEO_INVALID_DS;
	  }

	  //Check which translation is
	  if(transType.compare("TRANSLATION_X") == 0){
	    if(transG.setTranslationX(transIndex, ds) != 0){
	      printf("UNEXPECTED ERROR: Unable to create translation on X axis."
		     " Please, report this.");
	      return PEN_MESHBODY_GEO_UNEXPECTED_ERROR;
	    }
	  }
	  else if(transType.compare("TRANSLATION_Y") == 0){
	    if(transG.setTranslationY(transIndex, ds) != 0){
	      printf("UNEXPECTED ERROR: Unable to create translation on Y axis."
		     " Please, report this.");
	      return PEN_MESHBODY_GEO_UNEXPECTED_ERROR;
	    }
	  }
	  else if(transType.compare("TRANSLATION_Z") == 0){
	    if(transG.setTranslationZ(transIndex, ds) != 0){
	      printf("UNEXPECTED ERROR: Unable to create translation on Z axis."
		     " Please, report this.");
	      return PEN_MESHBODY_GEO_UNEXPECTED_ERROR;
	    }
	  }
	  else{
	    //Is a generic translation, read direction
	    v3D dir;

	    // dir.x
	    if(transSec.read("u", dir.x) != INTDATA_SUCCESS){
	      if(verbose > 0){
		printf("pen_meshBodyGeo:configure: Error: Unable to read 'u' direction for "
		       "translation %s, in vertex group '%s', in body '%s'. Double expected.\n",
		       transformNames[it].c_str(), transGroupNames[iTransGroup].c_str(),
		       transBodyNames[itransbody].c_str());	      
	      }
	      return PEN_MESHBODY_GEO_INVALID_DIR;
	    }

	    // dir.y
	    if(transSec.read("v", dir.y) != INTDATA_SUCCESS){
	      if(verbose > 0){
		printf("pen_meshBodyGeo:configure: Error: Unable to read 'v' direction for "
		       "translation %s, in vertex group '%s', in body '%s'. Double expected.\n",
		       transformNames[it].c_str(), transGroupNames[iTransGroup].c_str(),
		       transBodyNames[itransbody].c_str());	      
	      }
	      return PEN_MESHBODY_GEO_INVALID_DIR;
	    }

	    // dir.z
	    if(transSec.read("w", dir.z) != INTDATA_SUCCESS){
	      if(verbose > 0){
		printf("pen_meshBodyGeo:configure: Error: Unable to read 'w' direction for "
		       "translation %s, in vertex group '%s', in body '%s'. Double expected.\n",
		       transformNames[it].c_str(), transGroupNames[iTransGroup].c_str(),
		       transBodyNames[itransbody].c_str());	      
	      }
	      return PEN_MESHBODY_GEO_INVALID_DIR;
	    }

	    //Create the translation
	    if(transG.setTranslation(transIndex, dir, ds) != 0){
	      printf("UNEXPECTED ERROR: Unable to create generic translation."
		     " Please, report this.");
	      return PEN_MESHBODY_GEO_UNEXPECTED_ERROR;
	    }	    
	  }
	}
	else if(transType.compare("SCALE")   == 0 ||
		transType.compare("SCALE_XY") == 0 ||
		transType.compare("SCALE_XZ") == 0 ||
		transType.compare("SCALE_YZ") == 0){

	  //Is a scale transform, read the scale factor

	  double f;
	  if(transSec.read("factor", f) != INTDATA_SUCCESS){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error: Unable to read scale "
		     "factor 'factor' for transform %s, in vertex group '%s', in body '%s'."
		     " Double expected.\n", transformNames[it].c_str(),
		     transGroupNames[iTransGroup].c_str(),
		     transBodyNames[itransbody].c_str());
	    }
	    return PEN_MESHBODY_GEO_INVALID_SCALE;
	  }

	  //Check which scale transform is
	  if(transType.compare("SCALE_XY") == 0){
	    if(transG.setScaleXY(transIndex, f) != 0){
	      printf("UNEXPECTED ERROR: Unable to create scale transform on XY plane."
		     " Please, report this.");
	      return PEN_MESHBODY_GEO_UNEXPECTED_ERROR;
	    }
	  }
	  else if(transType.compare("SCALE_XZ") == 0){
	    if(transG.setScaleXZ(transIndex, f) != 0){
	      printf("UNEXPECTED ERROR: Unable to create scale transform on XZ plane."
		     " Please, report this.");
	      return PEN_MESHBODY_GEO_UNEXPECTED_ERROR;
	    }
	  }
	  else if(transType.compare("SCALE_YZ") == 0){
	    if(transG.setScaleYZ(transIndex, f) != 0){
	      printf("UNEXPECTED ERROR: Unable to create scale transform on YZ plane."
		     " Please, report this.");
	      return PEN_MESHBODY_GEO_UNEXPECTED_ERROR;
	    }
	  }
	  else{
	    //Generic scale transform
	    if(transG.setScale(transIndex, f) != 0){
	      printf("UNEXPECTED ERROR: Unable to create scale transform."
		     " Please, report this.");
	      return PEN_MESHBODY_GEO_UNEXPECTED_ERROR;
	    }
	  }
	}
	else{
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: Unknown type for "
		   "transformation %s, in vertex group '%s', in body '%s'."
		   " Provided type: %s\n",
		   transformNames[it].c_str(), transGroupNames[iTransGroup].c_str(),
		   transBodyNames[itransbody].c_str(), transType.c_str());
	  }
	  return PEN_MESHBODY_GEO_INVALID_TRANSFORMATION_TYPE;
	}
	
      }
    }

  }

  //Set geometry stream
  //**********************

  //Check if the geometry must be read from memory or a file
  bool inMemoryGeo;
  if(config.read("memory-file",inMemoryGeo) == INTDATA_SUCCESS){
    if(verbose > 1){
      if(inMemoryGeo)
	printf("Geometry will be read from memory buffer.\n");
      else
	printf("Geometry will be read from external file.\n");
    }
  }else{
    inMemoryGeo = false;
  }

  if(inMemoryGeo){
    //Geometry read from memory (variable "preloadGeo")
    std::istringstream in(preloadGeo);
    if(in.good()){
      //Load geometry
      //*****************
      err = GEOMESH(in,transMap,verbose);
      if(err != PEN_MESHBODY_GEO_SUCCESS){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error loading geometry.\n");
	  printf("                          Error code: %d\n",err);
	}
	configStatus = err;
	return configStatus;
      }
    }else{
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure:Error: unable to read "
	       "memory stream as input file.\n");
      }
      configStatus = PEN_MESHBODY_GEO_INPUT_SECTION;
      return configStatus;
    }
  }else{
    //Geometry read from external file

    //Read input file from configuration
    std::string infilename;
    if(config.read("input-file",infilename) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("meshBodyGeo:configure:Error: 'input-file' field missing at "
	       "configuration section.\n");
      }
      configStatus = PEN_MESHBODY_GEO_INPUT_SECTION;
      return configStatus;
    }
    
    //Open input file
    //*****************
    std::ifstream in(infilename, std::ifstream::in);
    if(!in.good()){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure:Error: unable to open "
	       "input file '%s'.\n", infilename.c_str());
      }
      configStatus = PEN_MESHBODY_GEO_INPUT_SECTION;
      return configStatus;
    }

  
    if(verbose > 0){
      printf(" input filename : %s\n",infilename.c_str());
    }  
  
    //Load geometry
    //*****************
    err = GEOMESH(in,transMap,verbose);
    if(err != PEN_MESHBODY_GEO_SUCCESS){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error loading geometry.\n");
	printf("                          Error code: %d\n",err);
      }
      in.close();
      configStatus = err;
      return configStatus;
    }
    in.close();    
  }
  
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
    if(verbose > 1)
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
    if(verbose > 1)
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
    if(verbose > 1)
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
    if(verbose > 1)
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
	  printf("   - Body %lu split begins\n",static_cast<unsigned long>(ibody));
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
	  printf("   + Body %lu split completed\n",static_cast<unsigned long>(ibody));
	  fflush(stdout);
	}
    
  }
  
#endif
  
  // Load absorption Energies
  //*************************
  bodiesAlias.clear();
  err = config.ls("eabs",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No absorption energies specified for any body\n");
    }
  }
  else{
    if(verbose > 1)
      printf("Absorption energies specified"
	     " for %lu bodies:\n\n",bodiesAlias.size());
    
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
	      printf("pen_meshBodyGeo:configure: Error reading"
		     " energy absorption at field '%s'. "
		     "Double expected.\n",key2.c_str());
	    }
	    return PEN_MESHBODY_GEO_BAD_READ_EABS;
	  }

	  if(eabs <= 0.0){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error: Invalid energy "
		     "absorption %12.4E for body '%s' particle '%s'. "
		     "Must be greater than zero.\n",
		     eabs,bodiesAlias[i].c_str(),particleNames[j].c_str());
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
			   triangle.stringify().c_str());			
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

int pen_meshBodyGeo::GEOMESH(std::istream& in,
			     std::map<std::string, std::vector<pen_meshTransform::group>>& transMap,
			     const unsigned verbose){

  //Create a vector for inlcuded files
  std::vector<std::ifstream> includes;
  
  //Read comment lines (start with #)
  std::string line;
  unsigned long nlines;
  unsigned long nRead = 0;
  int nBodies;
    
  if(meshGetLine(includes,in,line,nlines) == 0){
    nRead += nlines;
    //Read number of bodies in the geometry file
    if(sscanf(line.c_str()," %d",&nBodies) != 1){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error reading number of objects."
	       "Unexpected format in line %lu: \n %s\n", nRead, line.c_str());
      }
      return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
    }

    if(nBodies <= 0){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error: number of bodies must "
	       "be greater than zero.\n");
	printf("        Number of bodies read is %d. \n",nBodies);
      }
      return PEN_MESHBODY_GEO_INVALID_NBODIES;
    }
        
    NBODYS = static_cast<unsigned int>(nBodies);
        
    if(NBODYS > NB){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error: number of bodies read "
	       "from the geometry file: %u, is greater than the maximum "
	       "number of bodies allowed: %u\n.",NBODYS, NB);
      }
      return PEN_MESHBODY_GEO_INVALID_NBODIES;
    }

            
    if(verbose > 1)
      {
	printf("Number of bodies read: %u\n",NBODYS);
	printf("\n");
      }
  }
  else{
    if(verbose > 0){
      printf("pen_meshBodyGeo:configure: Error: Invalid geometry file.\n"
	     "Unable to read a single data line. Possible end of file reached.");
    }
    return PEN_MESHBODY_GEO_INVALID_FILE;
  }
    
    
  //Create vector to save bodies vertex info
  std::vector<v3D> vertex;
  unsigned long nvertex;
    
  if(verbose > 1){
    printf("Bodies mesh information:\n");
  }
    
  for(size_t i = 0; i < NBODYS; ++i){

    //Create a map to store vertex groups
    std::map<std::string,std::vector<unsigned>> vgMap;
    
    if(meshGetLine(includes,in,line,nlines) == 0){
      nRead += nlines;
      //Read characteristics of each body in the geometry file
      //Create auxiliar variable to ensure material is greater than zero
      int mat;
      int nVertexGroups;
      long int nTrianglesAux;
      long int nVertexAux;
      int nParamRead = sscanf(line.c_str()," %d  %ld  %ld  %s   %s %d",
			      &mat, &nTrianglesAux, &nVertexAux,
			      bodies[i].BALIAS, bodies[i].PALIAS, &nVertexGroups);
      if(nParamRead != 5 && nParamRead != 6){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error reading object header information.\n"
		 "Unexpected format in line %lu: \n %s\n", nRead, line.c_str());
	  printf("Expected format is #MAT      #NFACES     #NVERTEX     "
		 "#NAME        #PARENT NAME    (#NVGROUPS)\n");
	}
	return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
      }

      if(nParamRead == 5){
	//Assume no vertex group
	nVertexGroups = 0;
      }

      //Check material
      if(mat < 0){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: 'material' must be greater "
		 "or equal to zero.\n");
	  printf("            Specified for body %ld: %d\n",i,mat);
	}
	return PEN_MESHBODY_GEO_INVALID_MAT;
      }
      else{
	bodies[i].MATER = static_cast<unsigned>(mat);
      }

      //Check vertex groups
      if(nVertexGroups < 0){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: number of vertex groups must "
		 "be greater or equal to zero.\n");
	  printf("            Specified for body %ld: %d\n",i,nVertexGroups);
	}
	return PEN_MESHBODY_GEO_INVALID_N_VERTEX_GROUP;
      }
                
      if(bodies[i].MATER > constants::MAXMAT){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: material number read from the "
		 "geometry file: %u, is greater than the maximum number of "
		 "materials allowed: %u\n", bodies[i].MATER, constants::MAXMAT);
	}
	return PEN_MESHBODY_GEO_INVALID_MAT;
      }

      if(nTrianglesAux < 0){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: Invalid number of triangles. "
		 "Must be greater than zero.\n");
	  printf("            Specified for body %ld: %ld\n",i,nTrianglesAux);
	}
	return PEN_MESHBODY_GEO_INVALID_TRIANGLES_NUMBER;	
      }

      bodies[i].nTriangles = static_cast<unsigned long>(nTrianglesAux);
      
      //Check number of vertex
      if(nVertexAux < 0){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error: Invalid number of vertex. "
		 "Must be greater than zero.\n");
	  printf("            Specified for body %ld: %ld\n",i,nVertexAux);
	}
	return PEN_MESHBODY_GEO_INVALID_VERTEX_NUMBER;	
      }else{
	nvertex = static_cast<unsigned long>(nVertexAux);
      }
            
      //Resize vector of vertex with the number of vertex
      //read for each body from the geometry file
      vertex.resize(nvertex);
            
      if(vertex.size() != nvertex){
	if(verbose > 0){
	  printf("pen_meshBodyGeo:configure: Error allocating memory for vertex"
		 " reading. Can't allocate memory for %lu vertex.\n", nvertex);
	}
	return PEN_MESHBODY_BAD_MEMORY_ALLOCATION;
      }
      
      //Read vertex groups data from geometry file
      for(int j = 0; j < nVertexGroups; ++j){
	//Read vertex group header
	if(meshGetLine(includes,in,line,nlines) == 0){
	  nRead += nlines;
	  char groupName[100];
	  long int nGroupVertex;
	  if(sscanf(line.c_str(), " %s %ld ", groupName, &nGroupVertex) != 2){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error reading vertex group information."
		     "Unexpected format in line %lu: \n %s\n", nRead, line.c_str());
	    }
	    return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
	  }

	  //Check vertex number
	  if(nGroupVertex <= 0){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error: The number of vertex in a group"
		     " must be greater than 0.\n"
		     " Vertex in group '%s': %ld\n", groupName, nGroupVertex);
	    }
	    return PEN_MESHBODY_GEO_INVALID_N_VERTEX_GROUP;	    
	  }

	  //Create a entry for this vertex group and resize it
	  std::vector<unsigned>& vgIndex = vgMap[std::string(groupName)];
	  vgIndex.resize(nGroupVertex);

	  //Read vertex indexes belonging this group
	  for(long int iv = 0; iv < nGroupVertex; ++iv){
	    if(meshGetLine(includes,in,line,nlines) == 0){
	      nRead += nlines;
	      long int vIndex;
	      if(sscanf(line.c_str(), " %ld ", &vIndex) != 1){
		if(verbose > 0){
		  printf("pen_meshBodyGeo:configure: Error reading vertex group "
			 "information. Unexpected format in line %lu: \n %s\n",
			 nRead, line.c_str());
		}
		return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
	      }

	      //Check vertex index
	      if(vIndex < 0 ||
		 vIndex >= static_cast<long int>(nvertex)){
		if(verbose > 0){
		  printf("pen_meshBodyGeo:configure: Invalid vertex index in group '%s'."
			 "Index %ld is not in the interval [0,%lu).\n",
			 groupName, vIndex, nvertex);
		}
		return PEN_MESHBODY_GEO_INVALID_VERTEX_INDEX;		
	      }

	      //Save index
	      vgIndex[iv] = static_cast<unsigned>(vIndex);
	    }
	    else{
	      if(verbose > 0){
		printf("pen_meshBodyGeo:configure: Error: Invalid geometry file.\n"
		       "Unable to read vertex number %ld for vertex group %d "
		       "in body %lu. Possible end of file reached.",
		       iv,j,static_cast<unsigned long>(i));
	      }
	      return PEN_MESHBODY_GEO_INVALID_FILE;
	    }
	  }
	  
	}
	else{
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: Invalid geometry file.\n"
		   "Unable to read vertex group %d data for body %lu. "
		   "Possible end of file reached.",
		   j, static_cast<unsigned long>(i));
	  }
	  return PEN_MESHBODY_GEO_INVALID_FILE;
	}	
      }

      //** Ensure that defined vertex groups during the
      //configuration are defined also in the geometry file

      //First, check if some transformation has been defined for this body

      const auto search = transMap.find(bodies[i].BALIAS);
      if(search != transMap.end()){
	//Transformations have been defined for this body.
	//Ensure vertex groups are defined in the geometry file
	for(const pen_meshTransform::group& g : search->second){
	  if(vgMap.count(g.name) == 0){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error: Vertex group '%s'"
		     " not found in body '%s'.\n",
		     g.name.c_str(), bodies[i].BALIAS);
	    }
	    return PEN_MESHBODY_GEO_VG_NOT_FOUND;
	  }
	}
      }

      //Read all vertex for this body
      for(size_t j = 0; j < nvertex; ++j){
	if(meshGetLine(includes,in,line,nlines) == 0){
	  long int index;
	  double x,y,z;
	  nRead += nlines;
	  if(sscanf(line.c_str()," %ld  %le  %le  %le", &index, &x, &y, &z) != 4){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Error reading vertex information."
		     "Unexpected format in line %lu: \n %s\n", nRead, line.c_str());
	    }
	    return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
	  }
	  
	  if(index < 0 ||
	     index > static_cast<long int>(nvertex)){
	    if(verbose > 0){
	      printf("pen_meshBodyGeo:configure: Invalid vertex index in body '%s'."
		     "Index %ld is not in the interval [0,%lu).\n",
		     bodies[i].BALIAS, index, nvertex);
	    }
	    return PEN_MESHBODY_GEO_INVALID_VERTEX_INDEX;
	  }

	  //Save vertex data
	  vertex[index].x = x;
	  vertex[index].y = y;
	  vertex[index].z = z;

	}
	else{
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: Invalid geometry file.\n"
		   "Unable to read vertex %lu data for body %lu. "
		   "Possible end of file reached.",
		   static_cast<unsigned long>(j), static_cast<unsigned long>(i));	
	  }
	  return PEN_MESHBODY_GEO_INVALID_FILE;
	}
      }

      //Apply transformations to vertex groups, if required
      if(transMap.count(bodies[i].BALIAS) != 0){
	const std::vector<pen_meshTransform::group>& transGroups =
	  transMap.at(bodies[i].BALIAS);
	for(size_t ivg = 0; ivg < transGroups.size(); ++ivg){
	  transGroups[ivg].apply(vgMap.at(transGroups[ivg].name), vertex);
	}
      }

      //Create the boundary box
      //Set mins and maxs as the first read vertex coordenates
      bodies[i].boundary.set(vertex[0].x,vertex[0].y,vertex[0].z,
			     vertex[0].x,vertex[0].y,vertex[0].z);      
      for(size_t j = 1; j < nvertex; ++j){
	//Enlarge boundary region
	bodies[i].boundary.enlarge(vertex[j]);	
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
	if(meshGetLine(includes,in,line,nlines) == 0){
	  nRead += nlines;
	  //Read each triangle or face of each body 
	      if(sscanf(line.c_str()," %u  %u  %u",
			&index[0], &index[1], &index[2]) != 3){
		if(verbose > 0){
		  printf("pen_meshBodyGeo:configure: Error reading triangle faces."
			 "Unexpected format in line %lu: \n %s\n", nRead, line.c_str());
		}
		return PEN_MESHBODY_GEO_UNEXPECTED_LINE_FORMAT;
	      }
                    
	  //Fill triangle
	  triangles[k].fill(vertex[index[0]],vertex[index[1]],vertex[index[2]]);
	}
	else{
	  if(verbose > 0){
	    printf("pen_meshBodyGeo:configure: Error: Invalid geometry file.\n"
		   "Unable to read triangle %lu data for body %lu. "
		   "Possible end of file reached.",
		   static_cast<unsigned long>(k), static_cast<unsigned long>(i));
	  }
	  return PEN_MESHBODY_GEO_INVALID_FILE;
	}
      }

      //Add the region to the regions vector of the first super region
      bodies[i].regions.emplace_back(region);
      bodies[i].regions[0].elements.push_back(region);
    }
    else{
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error: Invalid geometry file.\n"
	       "Unable to read data for body %lu. Possible end of file reached.",
	       static_cast<unsigned long>(i));	
      }
      return PEN_MESHBODY_GEO_INVALID_FILE;
    }
    
    if(verbose > 1){
      printf("Body number %lu\n"
	     "Body alias: %s\n"
	     "Parent alias: %s\n"
	     "Material read: %u\n"
	     "Number of triangles read: %lu\n"
	     "Number of vertex read: %lu\n"
	     "Number of vertex groups: %lu\n",
	     i, bodies[i].BALIAS, bodies[i].PALIAS,
	     bodies[i].MATER, bodies[i].nTriangles, nvertex,
	     static_cast<unsigned long>(vgMap.size()));

      if(vgMap.size() > 0){

	printf("Vertex groups information: \n");
	
	//Search transformation groups for this body
	const auto& search = transMap.find(bodies[i].BALIAS);
	if(search == transMap.end()){
	  //No transformations applied to VG
	  for(const auto& vg : vgMap){
	    printf("  %20.20s : %8lu vertex, no transformation applied\n",
		   vg.first.c_str(),
		   static_cast<unsigned long>(vg.second.size()));
	  }
	}else{
	  //Trasformations have been applied to this body
	  const std::vector<pen_meshTransform::group>& transGroups = search->second;
	
	  for(const auto& vg : vgMap){
	    printf("  %20.20s : %8lu vertex",
		   vg.first.c_str(),
		   static_cast<unsigned long>(vg.second.size()));

	    bool found = false;
	    for(const pen_meshTransform::group& g : transGroups){
	      if(g.name.compare(vg.first) == 0){
		printf(", %4lu transformations applied:\n\n %s\n",
		       static_cast<unsigned long>(g.size()),
		       g.stringify().c_str());
		
		found = true;
		break;
	      }
	    }
	    if(!found){
	      printf(",   no transformations applied\n");
	    }
	  }
	}
      }

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
    
  //Ensure all transformation bodies have been found
  for(const auto& pair : transMap){
    //Find the body for this transformation groups
    bool found = false;
    for(size_t i = 0; i < NBODYS; ++i){
      if(pair.first.compare(bodies[i].BALIAS) == 0){
	found = true;
	break;
      }
    }
    if(!found){
      if(verbose > 0){
	printf("pen_meshBodyGeo:configure: Error: Transformations for body '%s' "
	       "have been defined in the configuration file, but it is not found"
	       " in the geometry file.\n",pair.first.c_str());
      }
      return PEN_MESHBODY_GEO_BODY_NOT_FOUND;
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

int pen_meshBodyGeo::meshGetLine(std::vector<std::ifstream>& included,
				 std::istream& root,
				 std::string&line,
				 unsigned long& nRead){

  //Clear line
  line.clear();

  //Read until non empty line

  while(line.empty()){
    //Check if we must read from main istream or from a included file
    if(included.size() > 0){
      //Read from last included file
      std::ifstream& in = included.back();
      int err = pen_getLine(in, line, nRead); 
      if(err != 0){
	//Error while reading, check if end of file has been reached
	if(in.eof()){
	  //End of file reached, close the file and remove the stream
	  in.close();
	  included.pop_back();
	  return meshGetLine(included, root, line, nRead);
	}else{
	  //Error reading data, return the error
	  return err;
	}
      }
    }else{
      //Read from root stream
      int err = pen_getLine(root, line, nRead); 
      if(err != 0){
	return err;
      }
    }

    //After successful read, check if line stores a "include" instruction
      
    //Find non white characters
    const char* whiteChars = " \n\t\r";
    const std::string::size_type firstCharPos = line.find_first_not_of(whiteChars);
    if(firstCharPos == std::string::npos){
      //Empty line
      line.clear();
    }else{
      //Line with data, try to get two words
      char word1[200], word2[200];
      word1[0] = '\0';
      word2[0] = '\0';
      sscanf(line.c_str(), " %s %s ", word1, word2);

      //Check if it is an include instruction
      if(std::strcmp(word1,"include") == 0){
	//Open the file and read from it
	included.emplace_back(word2);
	if(!included.back().is_open()){
	  printf("pen_meshBodyGeo:configure: Error: Unable to "
		 "open included file '%s'\n", word2);
	  return -1;
	}
	return meshGetLine(included, root, line, nRead);
      }
    }
  }

  return 0;
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
