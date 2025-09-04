
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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


#include "pen_samplers.hh"


abc_spatialSampler::abc_spatialSampler(){
  //Set rotation to identity matrix
  rotation[0] = 1.0; rotation[1] = 0.0; rotation[2] = 0.0;
  rotation[3] = 0.0; rotation[4] = 1.0; rotation[5] = 0.0;
  rotation[6] = 0.0; rotation[7] = 0.0; rotation[8] = 1.0;
  
  //Set translation to zero
  translation[0] = 0.0; translation[1] = 0.0; translation[2] = 0.0;

  rotate = true;
}

void abc_spatialSampler::sample(pen_particleState& state, pen_rand& random) const{

  double pos[3] = {0.0,0.0,0.0};
  // Sample position
  geoSampling(pos,random);
  // Apply rotation
  if(rotate)
    matmul3D(rotation,pos);
  // Apply translation
  pos[0] += translation[0];
  pos[1] += translation[1];
  pos[2] += translation[2];

  state.X = pos[0];
  state.Y = pos[1];
  state.Z = pos[2];
}

abc_directionSampler::abc_directionSampler()
{}

void abc_directionSampler::sample(pen_particleState& state, pen_rand& random) const{
  double dir[3] = {0.0,0.0,0.0};
  directionSampling(dir,random);

  state.U = dir[0];
  state.V = dir[1];
  state.W = dir[2];
}

abc_energySampler::abc_energySampler() : maxEnergy(1.0e9), minimumEnergy(50.0)
{}

void abc_energySampler::sample(pen_particleState& state, pen_rand& random) const{

  double E;
  // Sample energy
  energySampling(E,random);

  //Check energy interval
  if(E > maxEnergy  || E < minimumEnergy){
    char error[400];
    sprintf(error,"Sampled energy %12.5E out of range (%12.5E-%12.5E)",E,minimumEnergy, maxEnergy);
    throw std::out_of_range(error);
  }

  state.E = E;
}

abc_timeSampler::abc_timeSampler(){}

void abc_timeSampler::sample(pen_particleState& state, pen_rand& random) const{

  double time;
  // Sample energy
  timeSampling(time,random);

  state.PAGE = time;
}


pen_genericStateGen::pen_genericStateGen() : spatialSampler(nullptr),
					     directionSampler(nullptr),
					     energySampler(nullptr),
					     timeSampler(nullptr),
					     geometry(nullptr),
					     configStatus(0),
					     Emax(0.0),
					     rotate(false),
					     translation{0.0,0.0,0.0},
					     sourceBody(-1),
					     sourceMat(0),
					     name("unamed"),
					     LAGE(false),
					     kpar(PEN_PHOTON)
{}

#ifdef _PEN_USE_LB_
//Define handlefinish for load balance compiled source
bool pen_genericStateGen::handleFinish(const unsigned iw,
				       const unsigned long long nDone,
				       unsigned long long& assigned,
				       const unsigned verbose){
  //End of simulation, try to end the task

  //Sleep few seconds to ensure minimum time between reports
  std::this_thread::sleep_for (std::chrono::seconds(2));
  //Perform a report
  int errReport;
  report(iw,nDone,&errReport,verbose);
  int why;
  if(workerFinish(iw,why,verbose)){
    //This worker can finish
    return true;
  }
  else{
    //This worker can't finish yet, check the reason
    switch(why){
    case 0:
      {
	//Update assigned histories
	const unsigned long long newAssign = toDo(iw);
	const unsigned long long lastAssign = assigned;
	assigned = newAssign;
	if(verbose > 1){
	  printf("Source '%s': Worker %u: Assigned histories "
		 "updated from %llu to %llu\n",
		 name.c_str(),iw,lastAssign,newAssign);
	}
      }break;
    case 1:
      {
	//Perform a checkpoint
	checkPoint(verbose);
	//Update histories to do
	const unsigned long long lastAssign = assigned;
	assigned = toDo(iw);
	if(lastAssign != assigned && verbose > 1)
	  printf("Source '%s': Worker %u: Assigned histories "
		 "updated from %llu to %llu\n",
		 name.c_str(),iw,lastAssign,assigned);
	
      }break;
    case 2:
      {
	//Rank 0 has not sent permission to finish the task
	std::this_thread::sleep_for (std::chrono::seconds(10));
	//Add only 10000 iterations and ask again
	const unsigned long long lastAssign = assigned;
	assigned += std::min(10000ull,nDone/1000ull);
	if(verbose > 1){	
	  printf("Source '%s': Worker %u: Assigned histories "
		 "updated from %llu to %llu\n",
		 name.c_str(),iw,lastAssign,assigned);
	}
      }break;
    default:
      {
	if(verbose > 1)
	  printf("Source '%s': Worker %u: Permission denied, "
		 "unexpected reason: %d\n",
		 name.c_str(),iw,why);
	assigned += 10000ull; //Add only 10000 iterations and ask again
      }
    }
    return false;
  }
}
#endif
      

std::string pen_genericStateGen::samplersList(){
  std::string aux;
  aux += " --- Spatial:\n";
  aux += spatialSamplers().typesList();
  aux += " --- Direction:\n";
  aux += directionSamplers().typesList();
  aux += " --- Energy:\n";
  aux += energySamplers().typesList();
  aux += " --- Time:\n";
  aux += timeSamplers().typesList();
  return aux;
}
std::string pen_genericStateGen::samplersList(std::vector<std::string>& spatial,
					      std::vector<std::string>& direction,
					      std::vector<std::string>& energy,
					      std::vector<std::string>& time){
  std::string aux;
  aux += " --- Spatial:\n";
  aux += spatialSamplers().typesList(spatial);
  aux += " --- Direction:\n";
  aux += directionSamplers().typesList(direction);
  aux += " --- Energy:\n";
  aux += energySamplers().typesList(energy);
  aux += " --- Time:\n";
  aux += timeSamplers().typesList(time);
  return aux;
}

int pen_genericStateGen::setGeometry(const wrapper_geometry* geometryIn){
  if(geometryIn == nullptr)
    return -1;
  geometry = geometryIn;

  //Update samplers geometry
  if(spatialSampler != nullptr)
    spatialSampler->updateGeometry(geometry);
  
  if(directionSampler != nullptr)
    directionSampler->updateGeometry(geometry);
  
  if(energySampler != nullptr)
    energySampler->updateGeometry(geometry);
  
  if(timeSampler != nullptr)
    timeSampler->updateGeometry(geometry);
  
  return 0;
}

int pen_genericStateGen::selectSpatialSampler(const char* ID,
					      const pen_parserSection& config,
					      const unsigned verbose){


  if(spatialSampler != nullptr){
    delete spatialSampler;
    spatialSampler = nullptr;    
  }
  
  spatialSampler = spatialSamplers().createInstance(ID);
  if(spatialSampler == nullptr){
    if(verbose > 0){
      printf("SelectSpatialSampler: Error: Unable to create "
	     "spatial sampler '%s'\n",ID);
    }
    return -1;
  }
  
  int errConfig = spatialSampler->configure(config, verbose);
  if(errConfig != 0){
    if(verbose > 0){
      printf("SelectSpatialSampler: Error: Unable to "
	     "configurate spatial sampler '%s'\n",ID);
    }
    delete spatialSampler;
    spatialSampler = nullptr;
    return -2;
  }

  if(geometry != nullptr)
    spatialSampler->updateGeometry(geometry);
  
  return 0;
}

int pen_genericStateGen::selectDirectionSampler(const char* ID,
						const pen_parserSection& config,
						const unsigned verbose){
  
  if(directionSampler != nullptr){
    delete directionSampler;
    directionSampler = nullptr;
  }

  directionSampler = directionSamplers().createInstance(ID);
  if(directionSampler == nullptr){
    if(verbose > 0){
      printf("SelectDirectionSampler: Error: Unable to "
	     "create direction sampler '%s'\n",ID);
    }
    return -1;
  }

  int errConfig = directionSampler->configure(config,verbose);
  if(errConfig != 0){
    if(verbose > 0){
      printf("SelectDirectionSampler: Error: Unable to "
	     "configure direction sampler '%s'\n",ID);
    }
    delete directionSampler;
    directionSampler = nullptr;
    return -2;
  }

  if(geometry != nullptr)
    directionSampler->updateGeometry(geometry);  

  return 0;
}

int pen_genericStateGen::selectEnergySampler(const char* ID,
					     const pen_parserSection& config,
					     const unsigned verbose){
  if(energySampler != nullptr){
    delete energySampler;
    energySampler = nullptr;
  }

  energySampler = energySamplers().createInstance(ID);
  if(energySampler == nullptr){
    if(verbose > 0){
      printf("SelectEnergySampler: Error: Unable to create "
	     "energy sampler '%s'\n",ID);
    }
    return -1;
  }

  int errConfig = energySampler->configure(Emax,config,verbose);
  if(errConfig != 0){
    if(verbose > 0){
      printf("SelectEnergySampler: Error: Unable to "
	     "configure energy sampler '%s'\n",ID);
    }
    Emax = 0.0;
    delete energySampler;
    energySampler = nullptr;
    return -2;
  }

  if(geometry != nullptr)
    energySampler->updateGeometry(geometry);  
  
  return 0;
}

int pen_genericStateGen::selectTimeSampler(const char* ID,
					   const pen_parserSection& config,
					   const unsigned verbose){
  if(timeSampler != nullptr){
    delete timeSampler;
    timeSampler = nullptr;
  }

  timeSampler = timeSamplers().createInstance(ID);
  if(timeSampler == nullptr){
    if(verbose > 0){
      printf("SelectTimeSampler: Error: Unable to create "
	     "time sampler '%s'\n",ID);
    }
    return -1;
  }

  int errConfig = timeSampler->configure(config,verbose);
  if(errConfig != 0){
    if(verbose > 0){
      printf("SelectTimeSampler: Error: Unable to "
	     "configure time sampler '%s'\n",ID);
    }
    delete timeSampler;
    timeSampler = nullptr;
    return -2;
  }

  if(geometry != nullptr)
    timeSampler->updateGeometry(geometry);  
  
  return 0;
}

void pen_genericStateGen::sample(pen_particleState& state, pen_rand& random) const{

  //Reset state
  state.reset();

  //Perform time sampling if exists
  if(timeSampler != nullptr){
    timeSampler->sample(state,random);
  } else {
    state.PAGE = 0.0;
  }

  state.LAGE=LAGE;

  //Perform direction sampling
  directionSampler->sample(state,random);
  if(rotate){
    //Rotate sampled direction
    double dir[3] = {state.U,state.V,state.W};
    matmul3D(rotation, dir);
    state.U = dir[0];
    state.V = dir[1];
    state.W = dir[2];
  }
  
  //Perform spatial sampling and locate the particle
  spatialSampler->sample(state,random);
  if(rotate){
    //Rotate and move sampled position
    double pos[3] = {state.X,state.Y,state.Z};
    matmul3D(rotation, pos);
    state.X = pos[0] + translation[0];
    state.Y = pos[1] + translation[1];
    state.Z = pos[2] + translation[2];
  }else{
    //Apply post-translation
    state.X += translation[0];
    state.Y += translation[1];
    state.Z += translation[2];
  }
  
  geometry->locate(state);

  //Check if source is restricted to specified body
  //or material
  if(sourceBody >= 0){
    while(state.IBODY != (unsigned)sourceBody){	
      spatialSampler->sample(state,random);
      geometry->locate(state);
    }
  }
  else if(sourceMat > 0){
    while(state.MAT != sourceMat){	
      spatialSampler->sample(state,random);
      geometry->locate(state);
    }      
  }
  
  //Perform energy sampling
  energySampler->sample(state,random);

  //Its a primary (source) particle, set ILB[0] = 1
  state.ILB[0] = 1;
}

void pen_genericStateGen::configure(const pen_parserSection& config,
				    const unsigned verbose){

  int err;

  //Clear generator
  clear();

  //Check registered types to ensure static library linking of the register variable
  if(!penred::sampler::checkRegistered<pen_particleState>(verbose)){
    if(verbose > 0){
      printf("Warning: Some generic sampler types are not properly registered\n");
    }
  }
  
  //*******************************
  // Check if age must be recorded
  //*******************************

  err = config.read("record-time",LAGE);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("genericStateGen: configure: Field ('record-time') "
	     "not found. Age recording disabled.\n");
    }
    LAGE = false;
  }

  if(verbose > 1){
    printf("Age recording: %s\n", LAGE ? "enabled" : "disabled");
  }
  
  //*******************
  // Samplers sections
  //*******************
    
  //Extract 'spatial' sampler field
  pen_parserSection spatialSection;
  err = config.readSubsection("spatial",spatialSection);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("genericStateGen: configure: Error: '%s': Unable "
	     "to read field 'spatial'.\n",name.c_str());
    }
    configStatus = -1;
    return;
  }

  //Extract 'direction' sampler field
  pen_parserSection directionSection;
  err = config.readSubsection("direction",directionSection);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("genericStateGen: configure: Error: '%s': Unable "
	     "to read field 'direction'.\n",name.c_str());
    }
    configStatus = -2;
    return;
  }

  //Extract 'energy' sampler field
  pen_parserSection energySection;
  err = config.readSubsection("energy",energySection);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){    
      printf("genericStateGen: configure: Error: '%s': Unable "
	     "to read field 'energy'.\n",name.c_str());
    }
    configStatus = -3;
    return;
  }
  
  //Extract 'time' sampler field
  bool UseTimeSampling = true;
  pen_parserSection timeSection;
  err = config.readSubsection("time",timeSection);
  if(err != INTDATA_SUCCESS){
    //No time sampler
    UseTimeSampling = false;
  }


  //**************
  // Samplers IDs
  //**************    
    
  //Get spatial sampler ID
  std::string spatialType;
  err = spatialSection.read("type",spatialType);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("genericStateGen:configure: Error: '%s': Unable "
	     "to read field 'spatial/type'. String expected.\n",name.c_str());
    }
    configStatus = -1;
    return;
  }

  //Get direction sampler ID
  std::string directionType;
  err = directionSection.read("type",directionType);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("genericStateGen:configure: Error: '%s': Unable to "
	     "read field 'direction/type'. String expected.\n",name.c_str());
    }
    configStatus = -2;
    return;
  }

  //Get energy sampler ID
  std::string energyType;
  err = energySection.read("type",energyType);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("genericStateGen:configure: Error: '%s': Unable to "
	     "read field 'energy/type'. String expected.\n",name.c_str());
    }
    configStatus = -3;
    return;
  }
    
  std::string timeType;
  if(UseTimeSampling){
    //Get time sampler ID
    err = timeSection.read("type",timeType);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("genericStateGen:configure: Error: '%s': Unable to "
	       "read field 'time/type'. String expected.\n",name.c_str());
      }
      configStatus = -4;
      return;
    }
  }

  
  //********************************
  // Create and configure samplers
  //********************************

  if(verbose > 1){
    printf("\n------------------------------------\n");
    printf("\n **** Source '%s'\n",name.c_str());
    printf("\n------------------------------------\n\n");
    printf("\nSpatial sampler '%s':\n\n",spatialType.c_str());
  }
  err = selectSpatialSampler(spatialType.c_str(),spatialSection,verbose);
  if(verbose > 1){printf("\n------------------------------------\n\n");}     
  if(err != INTDATA_SUCCESS){
    clear();
    configStatus = -1;
    return;
  }

  if(verbose > 1){printf("\nDirection sampler '%s':\n\n",directionType.c_str());}
  err = selectDirectionSampler(directionType.c_str(),directionSection,verbose);
  if(verbose > 1){printf("\n------------------------------------\n\n");} 
  if(err != INTDATA_SUCCESS){
    clear();
    configStatus = -2;
    return;
  }

  if(verbose > 1){printf("\nEnergy sampler '%s':\n\n",energyType.c_str());}    
  err = selectEnergySampler(energyType.c_str(),energySection,verbose);
  if(verbose > 1){printf("\n------------------------------------\n\n");}
  if(err != INTDATA_SUCCESS){
    clear();
    configStatus = -3;
    return;
  }

  if(UseTimeSampling){
    if(verbose > 1){printf("\nTime sampler '%s':\n\n",timeType.c_str());}
    err = selectTimeSampler(timeType.c_str(),timeSection,verbose);
    if(verbose > 1){printf("\n------------------------------------\n\n");}
    if(err != INTDATA_SUCCESS){
      clear();
      configStatus = -4;
      return;
    }
  }

  //********************************
  // Optional parameters
  //********************************

  //Get source body
  int sbody;
  sourceBody = -1;
  err = config.read("source-body",sbody);
  if(err == INTDATA_SUCCESS && sbody >= 0){
    if(sbody >= int(pen_geoconst::NB)){
      if(verbose > 0){
	printf("genericStateGen: configure: Warning: '%s': Selected "
	       "source body (%d) out of range [0-%u)\n",
	       name.c_str(),sbody,pen_geoconst::NB);	
      }
    }
    else{
      sourceBody = sbody;
      if(verbose > 1){
	printf("Selected source body: %d\n",sourceBody);
      }
    }
  } else if(verbose > 1){
    printf("No source body selected\n");
  }

  //Get source material
  int smat;
  sourceMat = 0;
  err = config.read("source-material",smat);
  if(err == INTDATA_SUCCESS){
    if(smat < 1 || smat > int(constants::MAXMAT)){
      if(verbose > 0){
	printf("genericStateGen: configure: Warning: '%s': Selected "
	       "source material (%d) out of range (0-%u]\n",
	       name.c_str(),smat,constants::MAXMAT);	
      }
    }
    else{
      sourceMat = unsigned(smat);
      if(verbose > 1){
	printf("Selected source material: %u\n",sourceMat);
      }
    }
  } else if(verbose > 1){
    printf("No source material selected\n");
  }
  
  if(verbose > 1){

    printf("\n");
    if(!UseTimeSampling){
      printf("No time sampler selected\n\n");
    }
    printf("Generator '%s' configured with generic samplers:\n",name.c_str());
    printf("Spatial   -> %s\n", spatialID());
    printf("Direction -> %s\n", directionID());
    printf("Energy    -> %s\n", energyID());
    printf("Time      -> %s\n", timeID());
  }

  configStatus = 0;  
}

void pen_genericStateGen::clear(){
  if(spatialSampler != nullptr){
      delete spatialSampler;
      spatialSampler = nullptr;
  }
  if(directionSampler != nullptr){
      delete directionSampler;
      directionSampler = nullptr;
  }
  if(energySampler != nullptr){
      delete energySampler;
      energySampler = nullptr;
  }
  if(timeSampler != nullptr){
      delete timeSampler;
      timeSampler = nullptr;
  }
  geometry = nullptr;
  configStatus = 0;
  sourceBody = -1;
  sourceMat = 0;
  Emax = 0.0;
}

//Instantiators construct on first use
instantiator<abc_spatialSampler>& pen_genericStateGen::spatialSamplers(){
  static instantiator<abc_spatialSampler>* ans =
    new instantiator<abc_spatialSampler>;
  return *ans;
}
instantiator<abc_directionSampler>& pen_genericStateGen::directionSamplers(){
  static instantiator<abc_directionSampler>* ans =
    new instantiator<abc_directionSampler>;
  return *ans;
}
instantiator<abc_energySampler>& pen_genericStateGen::energySamplers(){
  static instantiator<abc_energySampler>* ans =
    new instantiator<abc_energySampler>;
  return *ans;
}
instantiator<abc_timeSampler>& pen_genericStateGen::timeSamplers(){
  static instantiator<abc_timeSampler>* ans =
    new instantiator<abc_timeSampler>;
  return *ans;
}

//Include defined samplers
#include "specificSamplers.cpp"
#include "spatialSamplers.cpp"
#include "directionSamplers.cpp"
#include "energySamplers.cpp"
#include "timeSamplers.cpp"


namespace penred{
  namespace sampler{

    template<>
    bool checkRegistered<pen_particleState>(const unsigned verbose){
      return checkRegistersSpatial<0>(verbose) &&
	checkRegistersDirection<0>(verbose) &&
	checkRegistersEnergy<0>(verbose) &&
	checkRegistersTime<0>(verbose) &&
	checkRegistersSpecificCommon<0>(verbose);
    }

    template<>
    bool checkRegistered<pen_state_gPol>(const unsigned verbose){
      return checkRegistersSpecificGPol<0>(verbose);
    }
  }
}
