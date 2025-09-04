
//
//
//    Copyright (C) 2023-2024 Universitat de València - UV
//    Copyright (C) 2023-2024 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

 
#include "tallyPSS.hh"

const std::array<std::string,11> pen_PSS::compatibleTallies = {"ANGULAR_DET",
							       "CYLINDRICAL_DOSE_DISTRIB",
							       "DICOM_DOSE_DISTRIB",
							       "DICOM_KERMA_TRACK_LENGTH",
							       "EDEP_BODY",
							       "EDEP_MAT",
							       "EMERGING_PART_DISTRIB",
							       "KERMA_TRACK_LENGTH",
							       "SPATIAL_DOSE_DISTRIB",
							       "SPHERICAL_DOSE_DISTRIB"};

const std::array<std::string,6> pen_PSS::noCompatibleTallies = {"CT_SINOGRAM",
								"IMPACT_DET",
								"PSF",
								"PSS",
								"SECONDARY_GEN",
								"TRACK"};

void pen_PSS::clear(){

  //Clear counters
  for(particleBKP& bkp : bkps)
    bkp.clear();

  if(primary != nullptr){
    delete primary;
    primary = nullptr;
  }
  if(scatter != nullptr){
    delete scatter;
    scatter = nullptr;
  }
  if(multiScatter != nullptr){
    delete multiScatter;
    multiScatter = nullptr;
  }

  std::fill(sourceMats.begin(),
	    sourceMats.end(),
	    false);

  lastKnock = false;
}

void pen_PSS::tally_beginSim(){

  //Clear counters
  for(particleBKP& bkp : bkps)
    bkp.clear();

  if(primary->has_beginSim()){
    primary->tally_beginSim();
    scatter->tally_beginSim();
    multiScatter->tally_beginSim();
  }
}

void pen_PSS::tally_endSim(const unsigned long long nhist){

  //Clear counters
  for(particleBKP& bkp : bkps)
    bkp.clear();
  
  if(primary->has_endSim()){
    primary->tally_endSim(nhist);
    scatter->tally_endSim(nhist);
    multiScatter->tally_endSim(nhist);
  }
}

void pen_PSS::tally_sampledPart(const unsigned long long nhist,
				const unsigned long long dhist,
				const unsigned kdet,
				const pen_KPAR kpar,
				const pen_particleState& state){

  //Clear counters
  for(particleBKP& bkp : bkps)
    bkp.clear();
  
  if(primary->has_sampledPart()){
    primary->tally_sampledPart(nhist,dhist,kdet,kpar,state);
  }

  //Init bkp
  actualBKP = 0;

  //Init last knock flag
  lastKnock = false;
}

void pen_PSS::tally_endHist(const unsigned long long nhist){

  //Flag end history in all tallies
  if(primary->has_endHist()){
    primary->tally_endHist(nhist);
    scatter->tally_endHist(nhist);
    multiScatter->tally_endHist(nhist);    
  }
  
}

void pen_PSS::tally_move2geo(const unsigned long long nhist,
			     const unsigned kdet,
			     const pen_KPAR kpar,
			     const pen_particleState& state,
			     const double dsef,
			     const double dstot){

  if(primary->has_move2geo()){
    primary->tally_move2geo(nhist,kdet,kpar,state,dsef,dstot);
  }
  
}

void pen_PSS::tally_beginPart(const unsigned long long nhist,
			      const unsigned kdet,
			      const pen_KPAR kpar,
			      const pen_particleState& state){

  //Check if this particle has been retrieved from the stack
  //or has been sampled
  unsigned bkp2Extract = actualBKP;
  if(bkps[kpar].stored() > 0){
    //Particle read from the stack
    const abc_particleStack* stack = readStack(kpar);

    //Retrieve the bkp of the next particle
    bkpInfo info = bkps[kpar].get();
    actualBKP = info.val;
    bkp2Extract = info.bkp2Extract;

    //Get actual number of stored particles in the bkp store
    const unsigned int lastStored = bkps[kpar].stored();

    //Compare with the actual number of particles in the stack
    const unsigned int actualStored = stack->getNSec();

    // The expected number of particles stored in the particle stack is
    //
    //  lastStored
    //
    // However, due a call to VR techniques, some particles could
    // be added since the end of the previous particle simulation
    //
  
    //Check if particles has been added to the stack
    if(lastStored < actualStored){
      //Particles created with VR techniques inherits the last stored
      //particle in the stack. Therefore, them must get the new particle bpk

      for(size_t i = 0; i < actualStored-lastStored; ++i)
	bkps[kpar].store(actualBKP,bkp2Extract);
    }
  }

  //Call the corresponding tally
  if(primary->has_beginPart()){
    if(bkp2Extract == 0)
      primary->tally_beginPart(nhist,kdet,kpar,state);
    else if(bkp2Extract == 1)
      scatter->tally_beginPart(nhist,kdet,kpar,state);
    else
      multiScatter->tally_beginPart(nhist,kdet,kpar,state);    
  }

  //Init last knock flag
  lastKnock = false;
  
}

void pen_PSS::tally_endPart(const unsigned long long nhist,
			    const pen_KPAR kpar,
			    const pen_particleState& state){

  // ** Check if some particles have been created by annihilation
  for(unsigned ikpar = 0; ikpar < constants::nParTypes; ++ikpar){

    //Avoid particles of the same type
    if(static_cast<pen_KPAR>(ikpar) == kpar)
      continue;
    
    //Read particle stack
    const abc_particleStack* stack = readStack(static_cast<pen_KPAR>(ikpar));

    //Get actual number of stored particles in the bkp store
    const unsigned int lastStored = bkps[ikpar].stored();

    //Compare with the actual number of particles in the stack
    const unsigned int actualStored = stack->getNSec();

    if(actualStored > lastStored){
      //Particles produced by the knock call must inherite the parent bkp
      //if are non enabled particles or increased by 1 otherwise

      if(ikpar == enabledKPAR){
	for(size_t i = 0; i < actualStored-lastStored; ++i)
	  bkps[ikpar].store(actualBKP+1, actualBKP);
      }else{
	for(size_t i = 0; i < actualStored-lastStored; ++i)
	  bkps[ikpar].store(actualBKP);
      }
    }
  }

  // ** Check also if VR techniques have been added particles to the stack
  
  //Read particle stack
  const abc_particleStack* stack = readStack(kpar);

  //Get actual number of stored particles in the bkp store
  const unsigned int lastStored = bkps[kpar].stored();

  //Compare with the actual number of particles in the stack
  const unsigned int actualStored = stack->getNSec();

  // The expected number of particles stored in the particle stack is
  //
  //  lastStored
  //
  // However, due a call to some VR technique, some particles could
  // be added since the knock of the previous particle simulation
  //
  
  //Check if particles has been added to the stack
  if(lastStored < actualStored){
    //Particles created with VR techniques inherits the last stored
    //particle in the stack. Therefore, them must get the actual bpk

    for(size_t i = 0; i < actualStored-lastStored; ++i)
      bkps[kpar].store(actualBKP);
  }
  
  //Call the corresponding tally
  if(primary->has_endPart()){
    if(actualBKP == 0)
      primary->tally_endPart(nhist,kpar,state);
    else if(actualBKP == 1)
      scatter->tally_endPart(nhist,kpar,state);
    else
      multiScatter->tally_endPart(nhist,kpar,state);    
  }
  
}

void pen_PSS::tally_localEdep(const unsigned long long nhist,
			      const pen_KPAR kpar,
			      const pen_particleState& state,
			      const double dE){

  //Call the corresponding tally
  if(primary->has_localEdep()){
    if(actualBKP == 0)
      primary->tally_localEdep(nhist,kpar,state,dE);
    else if(actualBKP == 1)
      scatter->tally_localEdep(nhist,kpar,state,dE);
    else
      multiScatter->tally_localEdep(nhist,kpar,state,dE);    
  }

  //Check if this local Edep follows a knock call
  if(lastKnock){

    //Increase the bkp if this particle is enabled
    if(kpar == enabledKPAR){

      //Avoid counting interactions in source materials for primary radiation
      if(actualBKP > 0 || !sourceMats[state.MAT]){
	++actualBKP;
      }
    }

    lastKnock = false;
  }
  
}

void pen_PSS::tally_step(const unsigned long long nhist,
			 const pen_KPAR kpar,
			 const pen_particleState& state,
			 const tally_StepData& stepData){

  //Call the corresponding tally
  if(primary->has_step()){
    if(actualBKP == 0)
      primary->tally_step(nhist,kpar,state,stepData);
    else if(actualBKP == 1)
      scatter->tally_step(nhist,kpar,state,stepData);
    else
      multiScatter->tally_step(nhist,kpar,state,stepData);    
  }  
  
}

void pen_PSS::tally_interfCross(const unsigned long long nhist,
				const unsigned kdet,
				const pen_KPAR kpar,
				const pen_particleState& state){

  //Call the corresponding tally
  if(primary->has_interfCross()){
    if(actualBKP == 0)
      primary->tally_interfCross(nhist,kdet,kpar,state);
    else if(actualBKP == 1)
      scatter->tally_interfCross(nhist,kdet,kpar,state);
    else
      multiScatter->tally_interfCross(nhist,kdet,kpar,state);    
  }
}

void pen_PSS::tally_matChange(const unsigned long long nhist,
			      const pen_KPAR kpar,
			      const pen_particleState& state,
			      const unsigned prevMat){

  //Call the corresponding tally
  if(primary->has_matChange()){
    if(actualBKP == 0)
      primary->tally_matChange(nhist,kpar,state,prevMat);
    else if(actualBKP == 1)
      scatter->tally_matChange(nhist,kpar,state,prevMat);
    else
      multiScatter->tally_matChange(nhist,kpar,state,prevMat);    
  }
}

void pen_PSS::tally_jump(const unsigned long long nhist,
			 const pen_KPAR kpar,
			 const pen_particleState& state,
			 const double ds){

  //Read particle stack
  const abc_particleStack* stack = readStack(kpar);

  //Get actual number of stored particles in the bkp store
  const unsigned int lastStored = bkps[kpar].stored();

  //Compare with the actual number of particles in the stack
  const unsigned int actualStored = stack->getNSec();

  // The expected number of particles stored in the particle stack is
  //
  //  lastStored
  //
  // However, due a call to some VR technique, some particles could
  // be added since the knock of the previous particle simulation
  //
  
  //Check if particles has been added to the stack
  if(lastStored < actualStored){
    //Particles created with VR techniques inherits the last stored
    //particle in the stack. Therefore, them must get the actual bpk

    for(size_t i = 0; i < actualStored-lastStored; ++i)
      bkps[kpar].store(actualBKP);
  }
  
  //Call the corresponding tally
  if(primary->has_jump()){
    if(actualBKP == 0)
      primary->tally_jump(nhist,kpar,state,ds);
    else if(actualBKP == 1)
      scatter->tally_jump(nhist,kpar,state,ds);
    else
      multiScatter->tally_jump(nhist,kpar,state,ds);    
  }
}

void pen_PSS::tally_knock(const unsigned long long nhist,
			  const pen_KPAR kpar,
			  const pen_particleState& state,
			  const int icol){

  //After knock, all stacks could be added particles
  //due the interaction. We must check them.
  for(unsigned ikpar = 0; ikpar < constants::nParTypes; ++ikpar){
    //Read particle stack
    const abc_particleStack* stack = readStack(static_cast<pen_KPAR>(ikpar));

    //Get actual number of stored particles in the bkp store
    const unsigned int lastStored = bkps[ikpar].stored();

    //Compare with the actual number of particles in the stack
    const unsigned int actualStored = stack->getNSec();

    if(actualStored > lastStored){
      //Particles produced by the knock call must inherite the parent bkp
      //if are non enabled particles or increased by 1 otherwise

      if(ikpar == enabledKPAR){
	for(size_t i = 0; i < actualStored-lastStored; ++i)
	  bkps[ikpar].store(actualBKP+1, actualBKP);
      }else{
	for(size_t i = 0; i < actualStored-lastStored; ++i)
	  bkps[ikpar].store(actualBKP);
      }
    }
  }

  //Call the corresponding tally
  if(primary->has_knock()){
    if(actualBKP == 0)
      primary->tally_knock(nhist,kpar,state,icol);
    else if(actualBKP == 1)
      scatter->tally_knock(nhist,kpar,state,icol);
    else
      multiScatter->tally_knock(nhist,kpar,state,icol);    
  }

  //Ensure this is not a delta interaction
  if(kpar == PEN_PHOTON && icol == 4)
    return;
  if(kpar == PEN_ELECTRON && icol == 4)
    return;
  if(kpar == PEN_POSITRON && icol == 5)
    return;
  
  //Flag knock as last call
  lastKnock = true;
  
}

void pen_PSS::tally_lastHist(const unsigned long long lasthist){

  //Call last history in all tallies
  if(primary->has_lastHist()){
    primary->tally_lastHist(lasthist);
    scatter->tally_lastHist(lasthist);
    multiScatter->tally_lastHist(lasthist);    
  }
}

int pen_PSS::configure(const wrapper_geometry& geometry,
		       const abc_material* const materials[constants::MAXMAT],
		       const pen_parserSection& config,
		       const unsigned verbose){

  clear();
  
  //Get the active particle
  std::string particleType;
  enabledKPAR = PEN_PHOTON; //Photon by default
  if(config.read("PSS-particle", particleType) == INTDATA_SUCCESS){
    //Get the particle kpar
    enabledKPAR = static_cast<pen_KPAR>(particleID(particleType.c_str()));
    if(enabledKPAR == ALWAYS_AT_END){
      if(verbose > 0){
	printf("pen_PSS:configure: Error: Unknown particle type '%s'\n",
	       particleType.c_str());
      }
      return -1;
    }    
  }

  //Read source material list
  pen_parserArray sourceMatsList;
  if(config.read("PSS-source-mats", sourceMatsList) == INTDATA_SUCCESS){

    for(size_t i = 0; i < sourceMatsList.size(); ++i){
      int auxMat;
      int err = sourceMatsList[i].read(auxMat);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("pen_PSS:configure: Error: Unable to read material "
		 "source at position %lu. Integer expected\n",
		 static_cast<unsigned long>(i));
	}
	return -1;
      }

      //Check the material number
      if(auxMat <= 0 || auxMat >= static_cast<int>(constants::MAXMAT)){
	if(verbose > 0){
	  printf("pen_PSS:configure: Error: Source material must be in the interval (0,%u).\n"
		 "                   Invalid value: %d\n",
		 constants::MAXMAT,auxMat);
	}
	return -1;	
      }

      sourceMats[auxMat] = true;
    }
    
  }

  //Read subtally type
  std::string subTallyType;
  if(config.read("PSS-tally", subTallyType) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pen_PSS:configure: Error reading 'PSS-tally' field. A string "
	     "with the tally to be used must be provided.\n");
    }
    return -2;
    
  }

  if(verbose > 1)
    printf(" PSS particle type: %s\n", particleName(enabledKPAR));

  
  //Create the tallies with the specified type and configure them
  primary      = pen_commonTallyCluster::createInstance(subTallyType);
  scatter      = pen_commonTallyCluster::createInstance(subTallyType);
  multiScatter = pen_commonTallyCluster::createInstance(subTallyType);
  if(primary == nullptr || scatter == nullptr || multiScatter == nullptr){
    if(verbose > 0){
      printf("pen_PSS:configure: Error: unable to create a tally "
	     "of type '%s'.\n",subTallyType.c_str());
    }
    return -3;
  }

  if(verbose > 1)
    printf(" PSS tally type: %s\n", subTallyType.c_str());

  //Check if the selected tally is compatible
  if(std::find(noCompatibleTallies.cbegin(),
	       noCompatibleTallies.cend(),
	       subTallyType) != noCompatibleTallies.cend()){
    if(verbose > 0){
      printf("pen_PSS:configure: Error: Tally type '%s' is not compatible with PSS.\n",
	     subTallyType.c_str());
    }
    return -3;
  }  
  if(std::find(compatibleTallies.cbegin(),
	       compatibleTallies.cend(),
	       subTallyType) == compatibleTallies.cend()){
    if(verbose > 1){
      printf("pen_PSS:configure: Warning: It is not known if tally type '%s' "
	     "is compatible  with PSS. The behavioud could be unexpected.\n",
	     subTallyType.c_str());
    }
  }
  

  //Set tally name
  std::string subTallyName = readName() + std::string("-primary");
  primary->setName(subTallyName);
  subTallyName = readName() + std::string("-scatter");
  scatter->setName(subTallyName);
  subTallyName = readName() + std::string("-multiScatter");
  multiScatter->setName(subTallyName);

  //Set thread
  primary->setThread(getThread());
  scatter->setThread(getThread());
  multiScatter->setThread(getThread());

  //Set OutDir
  primary->setOutputDirPath(readOutputDirPath());  
  scatter->setOutputDirPath(readOutputDirPath());  
  multiScatter->setOutputDirPath(readOutputDirPath());  

  //Set stacks
  for(unsigned is = 0; is < constants::nParTypes; ++is){
    const pen_KPAR iskpar = static_cast<pen_KPAR>(is);
    primary->setStack(iskpar, readStack(iskpar));
    scatter->setStack(iskpar, readStack(iskpar));
    multiScatter->setStack(iskpar, readStack(iskpar));
  }
  
  //Configure tallies
  int errConfig = primary->configure(geometry,materials,config,verbose);
  if(errConfig != 0){
    clear();
    if(verbose > 0){
      printf("pen_PSS:configure: Error configuring a tally of type '%s'\n",
	     subTallyType.c_str());
    }
    return -4;
  }
  errConfig = scatter->configure(geometry,materials,config,verbose);
  if(errConfig != 0){
    clear();
    if(verbose > 0){
      printf("pen_PSS:configure: Error configuring a tally of type '%s'\n",
	     subTallyType.c_str());
    }
    return -4;
  }
  errConfig = multiScatter->configure(geometry,materials,config,verbose);
  if(errConfig != 0){
    clear();
    if(verbose > 0){
      printf("pen_PSS:configure: Error configuring a tally of type '%s'\n",
	     subTallyType.c_str());
    }
    return -4;
  }

  //Configure dumps
  addSubDump(*primary);
  addSubDump(*scatter);
  addSubDump(*multiScatter);
  
  return 0;
}

int pen_PSS::sharedConfig(const pen_PSS& tally){
  
  int err1 = primary->shareConfig(*tally.primary);
  int err2 = scatter->shareConfig(*tally.scatter);
  int err3 = multiScatter->shareConfig(*tally.multiScatter);

  if(err1 != 0){
    printf("pen_PSS::sharedConfig: Error adding primary radiation tallies:\n"
	   "                       Error code: %d\n", err1);
  }
  if(err2 != 0){
    printf("pen_PSS::sharedConfig: Error adding scatter radiation tallies:\n"
	   "                       Error code: %d\n", err2);
  }
  if(err3 != 0){
    printf("pen_PSS::sharedConfig: Error adding multi-scatter radiation tallies:\n"
	   "                       Error code: %d\n", err3);
  }
  
  return err1 + err2 + err3;
}


void pen_PSS::saveData(const unsigned long long nhist) const{

  primary->saveData(nhist);
  scatter->saveData(nhist);
  multiScatter->saveData(nhist);
  
}

void pen_PSS::flush(void){
  primary->flush();
  scatter->flush();
  multiScatter->flush();
}

int pen_PSS::sumTally(const pen_PSS& tally){

  //Sum tallies
  int err1 = primary->sum(*(tally.primary));
  int err2 = scatter->sum(*(tally.scatter));
  int err3 = multiScatter->sum(*(tally.multiScatter));

  if(err1 != 0){
    printf("pen_PSS::sumTally: Error adding primary radiation tallies:\n"
	   "                   Error code: %d\n", err1);
  }
  if(err2 != 0){
    printf("pen_PSS::sumTally: Error adding scatter radiation tallies:\n"
	   "                   Error code: %d\n", err2);
  }
  if(err3 != 0){
    printf("pen_PSS::sumTally: Error adding multi-scatter radiation tallies:\n"
	   "                   Error code: %d\n", err3);
  }
  
  return err1 + err2 + err3;
}

REGISTER_COMMON_TALLY(pen_PSS)
