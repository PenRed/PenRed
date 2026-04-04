
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
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
//    
//


#include "filter_geo.hh"

penred::errors::Error pen_filterGeo::specificConfigure(const pen_parserSection& config, unsigned verbose){

  penred::errors::SpecificError<pen_filterGeo> error;
  
  //Read origin
  int err = config.read("origin", origin);
  if(err != INTDATA_SUCCESS){
    origin = 0.0;
  }

  if(verbose > 1){
    printf("Origin (cm): %E\n\n", origin);
  }
  
  //Read filters names
  std::vector<std::string> filters;
  err = config.ls("filters", filters);
  if(err != INTDATA_SUCCESS){
    error.code = PEN_FILTER_GEO_NO_FILTERS;
    error.description = "pen_filterGeo:configure:Error: No section 'filter' specified";
    return error;
  }

  NBODYS = filters.size(); //The geometry has as many bodies as filters
  if(NBODYS < 1){
    error.code = PEN_FILTER_GEO_INVALID_NUMBER_OF_FILTERS;
    error.description = "pen_filterGeo:configure:Error: No filters configured";
    return error;
  }

  if(NBODYS > NB){
    error.code = PEN_FILTER_GEO_INVALID_NUMBER_OF_FILTERS;
    error.description = "pen_filterGeo:configure: Error: Maximum number of filters exceeded: " +
      std::to_string(NB);
    return error;    
  }

  std::array<bool,NB> usedPosition;
  usedPosition.fill(false);
  
  for(const std::string& filter : filters){

    //Process each configured filter
    
    std::string filterPath = std::string("filters/") + filter;
    pen_parserSection filterSection;
    err = config.readSubsection(filterPath, filterSection);
    if(err != INTDATA_SUCCESS){
      error.code = PEN_FILTER_GEO_MISSING_PARAMETER;
      error.description = "pen_filterGeo:configure:Error: Unable to read '"
	+ filter + "' filter configuration";
      return error;
    }

    // ** Filter position
    
    int position;
    err = filterSection.read("position", position);
    if(err != INTDATA_SUCCESS){
      error.code = PEN_FILTER_GEO_MISSING_PARAMETER;
      error.description = "pen_filterGeo:configure:Error: Unable to read '"
	+ filter + "' filter 'position'";      
      return error;
    }

    if(position < 0 || position >= static_cast<int>(NBODYS)){
      error.code = PEN_FILTER_GEO_BAD_VALUE;
      error.description = "pen_filterGeo:configure:Error: Bad position for filter '"
	+ filter + "', must be positive and lesser than the number of filters";
      return error;
    }

    if(usedPosition[position]){
      error.code = PEN_FILTER_GEO_BAD_VALUE;
      error.description = "pen_filterGeo:configure:Error: Bad position for filter '"
	+ filter + "', already used by filter '" + bodies[position].name + "'";
      return error;      
    }

    //Get body corresponding to this position
    pen_filterBody& body = bodies[position];

    //Flag this position as used
    usedPosition[position] = true;

    //Set filter name
    body.name = filter;

    // Filter width
    double width;
    err = filterSection.read("width", width);
    if(err != INTDATA_SUCCESS){
      error.code = PEN_FILTER_GEO_MISSING_PARAMETER;
      error.description = "pen_filterGeo:configure:Error: Unable to read '"
	+ filter + "' filter 'width'";      
      return error;          
    }

    if(width <= 0.0){
      error.code = PEN_FILTER_GEO_BAD_VALUE;
      error.description = "pen_filterGeo:configure:Error: Invalid width for filter '"
	+ filter + "'. Must be greater than 0";
      return error;          
    }
    
    //Set width
    body.width = width;
    
    // ** Filter material index
    int matIndex;
    err = filterSection.read("material", matIndex);
    if(err != INTDATA_SUCCESS){
      error.code = PEN_FILTER_GEO_MISSING_PARAMETER;
      error.description = "pen_filterGeo:configure:Error: Unable to read '"
	+ filter + "' filter 'material'. Integer expected";            
      return error;          
    }

    if(matIndex < 0){
      error.code = PEN_FILTER_GEO_BAD_VALUE;
      error.description = "pen_filterGeo:configure:Error: Invalid material index for filter '"
	+ filter + "'. Must be greater or equal to 0";
      return error;
    }

    body.MATER = static_cast<unsigned>(matIndex);

    // ** Filter detector index
    int detector;
    err = filterSection.read("detector", detector);
    if(err == INTDATA_SUCCESS){

      if(detector < 0){
	error.code = PEN_FILTER_GEO_BAD_VALUE;
	error.description = "pen_filterGeo:configure:Error: Invalid detector index for filter '"
	+ filter + "'. Must be greater or equal to 0";
	return error;
      }
      
      body.KDET = detector;
    }

    // ** DSMAX
    body.DSMAX = width/10.0;
    double dsmax;
    err = filterSection.read("dsmax", dsmax);
    if(err == INTDATA_SUCCESS){
      if(dsmax > 0.0){
	body.DSMAX = dsmax;
      }
    }

    // ** Local EABS


    //Get particle names
    std::vector<std::string> particleNames;
    filterSection.ls("eabs",particleNames);

    //Get particle eabs
    double filterEABS[constants::nParTypes];
    for(unsigned j = 0; j < constants::nParTypes; j++)
      filterEABS[j] = -1.0;
	  
    for(unsigned j = 0; j < particleNames.size(); j++){

      unsigned kpar = particleID(particleNames[j].c_str());
      if(kpar >= ALWAYS_AT_END){
	error.code = PEN_FILTER_GEO_BAD_VALUE;
	error.description = "pen_filterGeo:configure:Error: Error on 'eabs' field for filter '"
	  + filter + "'. Unknown particle type: " + particleNames[j];
	return error;
      }
      
      std::string eabsPath = std::string("eabs/") + particleNames[j];
      double eabs;
      err = filterSection.read(eabsPath,eabs);
      if(err != INTDATA_SUCCESS){
	error.code = PEN_FILTER_GEO_MISSING_PARAMETER;
	error.description = "pen_filterGeo:configure:Error: Unable to read '"
	  + filter + "' filter absorption energy from path: " + eabsPath;            
	return error;
      }

      if(eabs <= 0.0){
	error.code = PEN_FILTER_GEO_BAD_VALUE;
	error.description = "pen_filterGeo:configure: Error on 'eabs' field for filter '"
	  + filter + "' and particle '" + particleNames[j] + "'. Absorption energy must be "
	  "greater than 0";
	return error;
      }
	  
      filterEABS[kpar] = eabs;
    }

    //Set filter eabs for each specified particle
    if(setBodyEabs(position,filterEABS) != 0){
      error.code = PEN_FILTER_GEO_BAD_VALUE;
      error.description = "pen_filterGeo:configure:Error: Error on 'eabs' field, "
	"unknown body '" + filter + "'";            
      return error;
    }
    
    if(verbose > 1){
      printf("%s:\n", body.name.c_str());
      printf("  Position: %d\n",position);
      printf("  Width   : %f\n",body.width);
      printf("  Material: %d\n",body.MATER);
      printf("  DSMAX   : %E\n",body.DSMAX);
      printf("  KDET    : %d\n",body.KDET);
      printf("  EABS    :\n");
      for(unsigned k = 0; k < constants::nParTypes; k++)
	printf("  %20.20s: %14.5E\n",particleName(k),body.localEABS[k]);
    }
    
  }

  //Complete filter bodies configuration
  bodies[0].origin = origin;
  //Filters grow through negative Z direction
  bodies[0].limit = origin - bodies[0].width;

  if(verbose > 1){
    printf("\n\n Geometry schema: \n\n");
    printf("------------------------- %15.5E\n", origin);
    printf("          %3d            \n", bodies[0].MATER);
    printf("------------------------- %15.5E\n", bodies[0].limit);
  }
  
  for(unsigned i = 1; i < NBODYS; ++i){

    const double localOrigin = bodies[i-1].limit;
    bodies[i].origin = localOrigin;
    bodies[i].limit = localOrigin - bodies[i].width;

    if(verbose > 1){
      printf("          %3d            \n", bodies[i].MATER);
      printf("------------------------- %15.5E\n", bodies[i].limit);
    }
    
  }
  
  return error;
}

void pen_filterGeo::locate(pen_particleState& state) const{

  for(unsigned i = 0; i < getBodies(); ++i)
    if(bodies[i].isIn(state)){
      state.MAT = bodies[i].MATER;
      state.IBODY = i;
      return;
    }

  state.IBODY = getBodies();
  state.MAT = 0;
  
}
  
 

void pen_filterGeo::step(pen_particleState& state, double DS, double &DSEF, double &DSTOT, int &NCROSS) const{

  double dsef = 0.0;
  double dstot = 0.0;
  
  //Check if the particle is outside the geometry system
  if(state.IBODY >= getBodies()){

    //Check if the particle aims to the geoemtry system

    if(state.Z >= origin){
      //Particle is over the geometry system
      if(state.W < 0.0){
	//Move the particle to reach the first filter
	double ds2filter = (bodies[0].origin - state.Z)/state.W; //Notice W sign
	state.X += state.U * ds2filter;
	state.Y += state.V * ds2filter;
	state.Z += state.W * ds2filter;
	state.IBODY = 0;
	state.MAT = bodies[0].MATER;

	dsef = ds2filter;
	dstot = ds2filter;
      }else{
	//Can not reach the geometry system
	state.X += state.U * 1.0e35;
	state.Y += state.V * 1.0e35;
	state.Z += state.W * 1.0e35;
	    
	state.IBODY = getBodies();
	state.MAT = 0;

	DSEF = DSTOT = 1.0e35;
	NCROSS = 0;
	return;	
      }
      
    }else{

      //Particle is under the geometry system
      if(state.W > 0.0){
	//Move the particle to reach the last filter
	const unsigned last = getBodies()-1;
	double ds2filter = (bodies[last].limit - state.Z)/state.W;
	state.X += state.U * ds2filter;
	state.Y += state.V * ds2filter;
	state.Z += state.W * ds2filter;
	state.IBODY = last;
	state.MAT = bodies[last].MATER;

	dsef = ds2filter;
	dstot = ds2filter;
      }else{
	//Can not reach the geometry system
	state.X += state.U * 1.0e35;
	state.Y += state.V * 1.0e35;
	state.Z += state.W * 1.0e35;
	    
	state.IBODY = getBodies();
	state.MAT = 0;

	DSEF = DSTOT = 1.0e35;
	NCROSS = 0;
	return;	
      }
      
    }

    //Check if the reached filter is a void region
    if(state.MAT != 0){
      //Is a non void region, finish
      NCROSS = 1;
      DSEF = dsef;
      DSTOT = dstot;
      return;
    }
  }

  // ** The particle is in the geometry system
  
  //Check if the particle is moving on Z axis
  if(state.W == 0.0){

    if(state.MAT == 0){
      //Can not reach a non void region
      state.X += state.U * 1.0e35;
      state.Y += state.V * 1.0e35;
      state.Z += state.W * 1.0e35;
	    
      state.IBODY = getBodies();
      state.MAT = 0;

      DSEF = DSTOT = 1.0e35;
      NCROSS = 0;
      return;	
      
    }
    
    //Unable to reach any filter limit
    state.X += DS*state.U;
    state.Y += DS*state.V;
    state.Z += DS*state.W;

    DSEF = DS;
    DSTOT = DS;
    NCROSS = 0;
    
    return;
  }

  //Set tracking variables
  vector3D<double> pos(state.X, state.Y, state.Z);
  const vector3D<double> dir(state.U, state.V, state.W);
  int ncross = 0;
  double toTravel = DS;
  const unsigned initMat = state.MAT;
  bool inVoid = initMat == 0;
  unsigned filterIndex = state.IBODY;
  const int filterIncrement = state.W > 0.0 ? -1 : 1;
  const unsigned lastFilter = filterIncrement < 0 ? 0 : getBodies()-1;

  for(;;){
    
    //Get the actual filter
    const pen_filterBody& filter = bodies[filterIndex];

    //Calculate distance to next filter
    //                                 ^
    // ----------------------- origin  |
    //                                 |
    //            x particle           |
    // ----------------------- limit   Z
    
    double ds2limit;
    if(filterIncrement < 0){
      //Moving up
      ds2limit = (filter.origin - pos.z)/dir.z;
	
    }else{
      //Moving down
      ds2limit = (filter.limit - pos.z)/dir.z;  //Notice W sign
      
    }

    //Travel distance does not affect void regions    
    if(filter.MATER != 0){

      //Check if the travel distance reaches the limit
      if(toTravel < ds2limit){
	//The particle is unable to reach the next interface
	pos = pos + dir*toTravel;

	if(filter.MATER == initMat)
	  dsef += toTravel;
	dstot += toTravel;
	break;
      }

      //Update distance to travel
      toTravel -= ds2limit;
    }

    //The particle will reach the next filter, move it
    pos = pos + dir*ds2limit;

    if(filter.MATER == initMat)
      dsef += ds2limit;
    dstot += ds2limit;
    
    
    //Check if the particle escapes the geometry
    if(filterIndex == lastFilter){
      //Particle escapes
      state.X += state.U * 1.0e35;
      state.Y += state.V * 1.0e35;
      state.Z += state.W * 1.0e35;
	    
      state.IBODY = getBodies();
      state.MAT = 0;

      if(initMat == 0)
	DSEF = 1.0e35;
      else
	DSEF = dsef;
      DSTOT = 1.0e35;

      if(!inVoid)
	++ncross;
      NCROSS = ncross;
      return;      
    }

    //Obtain next filter index
    filterIndex += filterIncrement;

    //Check if a interface is crossed
    const unsigned MATNext = bodies[filterIndex].MATER;
    if(MATNext == 0){ //Entering in a void region
      if(!inVoid){
	++ncross;
	inVoid = true;
      }
    } else if(inVoid){
      //Particle comes from a void region to a non void region. Stop it
      ++ncross;
      break;
    }else if(MATNext == initMat){
      //Entering in a new filter with the same material without crossing void regions
      if(bodies[filterIndex].KDET != filter.KDET){
	//The filter belongs to another detector.
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
  DSTOT = dstot;
  DSEF = dsef;
  NCROSS = ncross;
                
  state.X = pos.x;
  state.Y = pos.y;
  state.Z = pos.z;
    
  state.MAT = bodies[filterIndex].MATER;
  state.IBODY = filterIndex;  
}

REGISTER_GEOMETRY(pen_filterGeo,FILTER)
