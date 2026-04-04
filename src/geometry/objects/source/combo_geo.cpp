
//
//
//    Copyright (C) 2023-2024 Universitat de València - UV
//    Copyright (C) 2023-2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2025-2026 Vicent Giménez Alventosa
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
//        sanolgi@upvnet.upv.es              (Sandra Oliver Gil)
//    
//


#include "combo_geo.hh"

penred::errors::Error pen_comboGeo::specificConfigure(const pen_parserSection& config,
						      unsigned verbose){

  penred::errors::SpecificError<pen_meshBodyGeo> error;
  
  //Read the materials section
  pen_parserSection materialsSection;
  if(config.readSubsection("materials",materialsSection) != INTDATA_SUCCESS){
    error.code = SECTION_READ_FAIL;
    error.description = "pen_comboGeo:configure:Error: Unable to read 'materials' section";
    return error;
  }
  
  //Read the geometries section
  pen_parserSection geometriesSection;
  std::vector<std::string> geometriesNames;
  if(config.readSubsection("geometries",geometriesSection) != INTDATA_SUCCESS){
    error.code = SECTION_READ_FAIL;
    error.description = "pen_comboGeo:configure:Error: No 'geometries' section provided";
    return error;
  }
  
  //Extract material names
  geometriesSection.ls(geometriesNames);

  //Resize the geometry vector
  geometries.resize(geometriesNames.size(), nullptr);
  firstIBody.resize(geometriesNames.size());
  
  //Iterate over each geometry
  for(const std::string& geoName : geometriesNames){

    if(verbose > 1){
      printf(" - Creating geometry '%s'\n", geoName.c_str());
    }
    
    //Get geometry section
    pen_parserSection geometrySection;
    std::string geometrySectionPath("geometries/");
    geometrySectionPath += geoName;
    if(config.readSubsection(geometrySectionPath,
			     geometrySection) != INTDATA_SUCCESS){
      error.code = SECTION_READ_FAIL;
      error.description = "pen_comboGeo:configure:Error: Unable to read 'geometries/" +
	geoName + "' section";
      
      return error;
    }

    // ** Priority
    //Read geometry priority
    int priorityAux;
    if(geometrySection.read("priority", priorityAux) != INTDATA_SUCCESS){
      error.code = MISSING_PARAMETER;
      error.description = "pen_comboGeo:configure:Error: Unable to read 'priority' "
	"for geometry '" + geoName + "'.";
      
      return error;
    }

    //Check the priority value
    if(priorityAux < 0){
      error.code = BAD_VALUE;
      error.description = "pen_comboGeo:configure:Error: Invalid 'priority' value "
	"for geometry '" + geoName + "'. Priority must be positive.";
      
      return error;      
    }
    unsigned priority = static_cast<unsigned>(priorityAux);
    if(priority >= geometriesNames.size()){
      error.code = BAD_VALUE;
      error.description = "pen_comboGeo:configure:Error: Invalid 'priority' value "
	"for geometry '" + geoName + "'. Priority must be lesser than the number of geometries.";
      
      return error;
    }

    if(verbose > 1){
      printf("    + Priority: %d\n", priority);
    }
    
    //Check if this priority has already been used
    if(geometries[priority] != nullptr){
      error.code = BAD_VALUE;
      error.description = "pen_comboGeo:configure:Error: Priority value " + std::to_string(priority) +
	" assigned twice, in geometries '" + geoName + "' and '" + geometries[priority]->name + "'.";
      
      return error;
    }

    //Read the configuration section for this geometry
    pen_parserSection geometryConfigSection;
    if(geometrySection.readSubsection("config",
				      geometryConfigSection) != INTDATA_SUCCESS){
      error.code = SECTION_READ_FAIL;
      error.description = "pen_comboGeo:configure:Error: Unable to read 'config' "
	"section for geometry '" + geoName + "'";
      
      return error;
    }

    //Add the material section to the configuration
    geometryConfigSection.addSubsection("materials", materialsSection);
    
    //Try to instantiate the geometry and configure it
    //Get geometry type
    std::string geoType;
    if(geometryConfigSection.read("type",geoType) != INTDATA_SUCCESS){
      error.code = MISSING_PARAMETER;
      error.description = "pen_comboGeo:configure:Error: Field 'type' not specified "
	"for geometry '" + geoName + "'.";
      
      return error;
    }

    if(verbose > 1){
      printf("    + Type    : %s\n", geoType.c_str());
    }
    
    //Instantiate the geometry
    geometries[priority] = penGeoRegister_create(geoType.c_str());
    if(geometries[priority] == nullptr){
      error.code = BAD_VALUE;
      error.description = "pen_comboGeo:configure: Error creating geometry "
	"'" + geoName + "' of type '" + geoType + "'";
      
      return error;
    }

    //Set name
    geometries[priority]->name.assign(geoName);

    //Configure the geometry
    penred::errors::Error subGeoError = geometries[priority]->configure(geometryConfigSection, verbose);
    if(subGeoError){
      error.code = ERROR_CONFIGURING_SUBGEOEMTRY;
      error.description = "pen_comboGeo:configure: Error configuring geometry "
	"'" + geoName + "' of type '" + geoType + "'";
      error.setTrace(subGeoError);
      
      return error;
    }
  }

  //Get body info for each geometry and fill the global body array
  if(verbose > 1){
    printf("\n\n");
  }
  unsigned nextFirstIBody = 0;
  for(size_t igeo = 0; igeo < geometries.size(); ++igeo){
    
    //Set the first IBody for this geometry
    firstIBody[igeo] = nextFirstIBody;

    //Check body limit
    if(nextFirstIBody + geometries[igeo]->getBodies() > NB){
      error.code = BODY_LIMIT_REACHED;
      error.description = "pen_comboGeo:configure:Error: Global maximum number of bodies "
	"(" + std::to_string(NB) + ") reached at geometry '" + geometries[igeo]->name +
	"' (Priority " + std::to_string(igeo) + ").";
      
      return error;
    }
    
    //Iterate over geometry bodies
    for(size_t ibody = 0; ibody < geometries[igeo]->getBodies(); ++ibody){

      //Get global index
      unsigned iBodyIndex = nextFirstIBody+ibody;

      //Get body
      pen_comboBody& body = bodies[iBodyIndex];
      
      //Set body name
      body.name = geometries[igeo]->name + "_" + geometries[igeo]->getBodyName(ibody);

      
      //Set body material
      body.MATER = geometries[igeo]->getMat(ibody);

      //Set body detector
      body.KDET = geometries[igeo]->getDET(ibody);

      //Set body DSMAX
      body.DSMAX = geometries[igeo]->getDSMAX(ibody);

      //Set absorption energies
      for(unsigned ipar = 0; ipar < constants::nParTypes; ++ipar){
	body.localEABS[ipar] = geometries[igeo]->getEabs(ibody,ipar);
      }
      
      //Set geometry index
      body.geoIndex = igeo;
      
    }

    //Update next first IBody
    nextFirstIBody += geometries[igeo]->getBodies();
    
  }

  //Save the total number of bodies
  NBODYS = nextFirstIBody;

  return error;
}

unsigned pen_comboGeo::getIBody(const char* bodyName) const{

  //Extract the geometry index from name
  for(size_t i = 0; i < getBodies(); ++i){
    if(bodies[i].name.compare(bodyName) == 0){
      return i;
    }
  }
  return getBodies();
}

std::string pen_comboGeo::getBodyName(const unsigned ibody) const{
  if(ibody < getBodies()){
    return std::string(bodies[ibody].name);
  }else{
    return std::string("NONE");
  }
}

void pen_comboGeo::locate(pen_particleState& state) const{

  //Iterate through all geometries to find the body where
  //the particle is located
  for(size_t igeo = 0; igeo < geometries.size(); ++igeo){

    geometries[igeo]->locate(state);
    if(state.MAT != 0){
      //The particle is located in this geometry. Get the global body index
      unsigned auxBodyIndex = firstIBody[igeo] + state.IBODY;
      state.IBODY = auxBodyIndex;
      return;
    }
  }

  //The particle is located in a void/overwrite region in all geometries
  state.MAT = 0;
  state.IBODY = getBodies();
}

void pen_comboGeo::step(pen_particleState& state,
			double DS, double &DSEF, double &DSTOT,
			int &NCROSS) const{

  // ** Check if the particle is in a void region
  if(state.MAT == 0){

    pen_particleState finalState(state);
    double closestDS = 1.0e35;
    int finalNCross = 0;
    
    //Find the closest distance to a non void region
    for(size_t i = 0; i < geometries.size(); ++i){

      unsigned firstBodyIndex = firstIBody[i];

      //Copy the particle state
      pen_particleState auxState(state);

      //Locate the particle
      geometries[i]->locate(auxState);      
      
      double dsef, dstot;
      int totalNCross;
      geometries[i]->step(auxState, 1.0e35, dsef, dstot, totalNCross);

      //Check if the particle is not in the void after the step.
      if(auxState.MAT != 0){

	//The particle reaches a non void region in this geometry.
	//Check if it is the closest one
	if(closestDS > dstot){ //If it is the closest geometry, update it
	  auxState.IBODY += firstBodyIndex;
	  finalState = auxState;
	  finalNCross = totalNCross;
	  closestDS = dstot;
	}
	
      }
    }

    //Check if some non void region has been
    //reached in any geometry
    if(finalState.MAT != 0){
      state = finalState;
      DSEF = closestDS;
      DSTOT = closestDS;
      NCROSS = finalNCross;
    }else{
      //The particle does not reach any non void region
      //in any geometry
      state.X += state.U*1.0e35;
      state.Y += state.V*1.0e35;
      state.Z += state.W*1.0e35;

      state.MAT = 0;
      state.IBODY = getBodies();
      DSEF = 1.0e35;
      DSTOT = 1.0e35;
      NCROSS = 0;
    }
    
    return;
  }

  // ** The particle is NOT in a void region

  // ** Transport in the current geometry
  //****************************************
  
  //Get actual body
  const pen_comboBody& actualBody = bodies[state.IBODY];
  //Get the first global body index for this geometry
  const unsigned currentGeo = actualBody.geoIndex;
  unsigned firstCurrentBodyIndex = firstIBody[currentGeo];


  //First, calculate the step in the current geometry
  pen_particleState currentGeoState(state);
  currentGeoState.IBODY -= firstCurrentBodyIndex;
  double dsefCurrent, dstotCurrent;
  int ncrossCurrent;
  geometries[currentGeo]->step(currentGeoState,
			       DS, dsefCurrent, dstotCurrent,
			       ncrossCurrent);

  //Set the global body index
  if(currentGeoState.MAT == 0){
    currentGeoState.IBODY = getBodies();
  }else{
    currentGeoState.IBODY += firstCurrentBodyIndex;
  }
  

  // ** Transport in geometries with higher priority
  //************************************************
  
  //Save the maximum distance to reach a higher priority geometry
  double maxToTravel = dstotCurrent;

  //Save the closest higher priority geometry information
  unsigned finalBodyHighPriority;
  unsigned finalMatHighPriority = 0;

  //Notice that the particle must be in the void for all
  //geometries with a higher priority
  const pen_particleState origStateHighPriority = state;
    
  //Iterate over higher priority geometries
  for(size_t i = 0; i < actualBody.geoIndex; ++i){

    //Create an auxiliary state and set the local body index
    pen_particleState auxState(origStateHighPriority);

    //Locate the state in the geometry
    geometries[i]->locate(auxState);

    //Check if it is inside a non void region
    if(auxState.MAT != 0){
      //The particle reach a region, update particle state
      auxState.IBODY += firstIBody[i];
      state = auxState;

      DSEF = 0.0;
      DSTOT = 0.0;
      NCROSS = 1;
      return;
    }
    
    //Try to reach a non void region
    double dsef, dstot;
    int ncross;
    geometries[i]->step(auxState, 1.0e35, dsef, dstot, ncross);

    //Check if a non void region has been reached
    if(auxState.MAT != 0){
      //Check if it is the closest one
      if(dstot < maxToTravel){
	//Update maximum distance to travel
	maxToTravel = dstot;
	finalBodyHighPriority = auxState.IBODY + firstIBody[i];
	finalMatHighPriority = auxState.MAT;
      }
    }
  }

  // ** Transport in geometry with lower priority
  //************************************************
  
  //Check if the particle crosses an interface in the actual geometry.
  //If no interface is crossed, a lower priority geometry can't overwritte
  //the actual geometry. In addition, the traveled distance inside
  //the actual material must be lesser than the required to reach a
  //higher priority geometry 
  if(ncrossCurrent != 0 && dsefCurrent < maxToTravel){
    //An interface is crossed and the particle reach void region.
    //In addition, high priority geometries are not reached before
    //the void region.

    //Check if a non void region of a lower priority geometry is reached.
    //Calculate the position when the particle just enter the void region
    pen_particleState stateOnVoid(state);
    //Move the state to the frontier with the void region
    double dsef2Void = dsefCurrent+1.0e-8;
    stateOnVoid.X += stateOnVoid.U*dsef2Void;
    stateOnVoid.Y += stateOnVoid.V*dsef2Void;
    stateOnVoid.Z += stateOnVoid.W*dsef2Void;

    double maxToTravelFromVoid = maxToTravel-dsefCurrent;

    //Save the closest lower priority geometry information
    unsigned finalBodyLowPriority;
    unsigned finalMatLowPriority = 0;
      
    //Iterate over lower priority geometries
    for(size_t i = actualBody.geoIndex+1; i < geometries.size(); ++i){
	
      //Create an auxiliary state
      pen_particleState auxState(stateOnVoid);

      //Locate the state in the geometry
      geometries[i]->locate(auxState);

      //Check if it is insde a non void region
      if(auxState.MAT != 0){
	//The particle reach a region, update particle state
	auxState.IBODY += firstIBody[i];
	state = auxState;

	DSEF = dsefCurrent;
	DSTOT = dsefCurrent;
	NCROSS = ncrossCurrent;
	  
	return;
      }

      //Try to reach a non void region
      double dsef, dstot;
      int ncross;
      geometries[i]->step(auxState, 1.0e35, dsef, dstot, ncross);

      //Check if a non void region has been reached
      if(auxState.MAT != 0){
	//Check if it is the closest one
	if(dstot < maxToTravelFromVoid){
	  //Update maximum distance to travel
	  maxToTravelFromVoid = dstot;
	  finalBodyLowPriority = auxState.IBODY + firstIBody[i];
	  finalMatLowPriority = auxState.MAT;
	}
      }
    }
    
    //Check if a lower priority geometry has been reached before
    //the high priority geometries
    if(finalMatLowPriority != 0){
      //The particle reach a region, update the particle state
      stateOnVoid.IBODY = finalBodyLowPriority;
      stateOnVoid.MAT = finalMatLowPriority;

      stateOnVoid.X += stateOnVoid.U*maxToTravelFromVoid;
      stateOnVoid.Y += stateOnVoid.V*maxToTravelFromVoid;
      stateOnVoid.Z += stateOnVoid.W*maxToTravelFromVoid;

      state = stateOnVoid;

      DSEF = dsefCurrent;
      DSTOT = dsefCurrent + maxToTravelFromVoid;
      NCROSS = ncrossCurrent + 1;

      return;
    }
  }

  // ** Final state 
  //*****************
  
  //Check if a high priority geometry has been reached
  if(finalMatHighPriority != 0){
      //The particle reach a region, update particle state
      state.IBODY = finalBodyHighPriority;
      state.MAT = finalMatHighPriority;

      state.X += state.U*maxToTravel;
      state.Y += state.V*maxToTravel;
      state.Z += state.W*maxToTravel;
      
      DSEF = maxToTravel < dsefCurrent ? maxToTravel : dsefCurrent;
      DSTOT = maxToTravel;
      NCROSS = 1;
	
      return;    
  }

  //No higher nor lower priority geometries have been reached
  DSEF = dsefCurrent;
  DSTOT = dstotCurrent;
  NCROSS = ncrossCurrent;

  state = currentGeoState;
}

REGISTER_GEOMETRY(pen_comboGeo,COMBO)
