
//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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


#include "tallyDetectionEDep.hh"

//Reader functions
int tallyReader_DetectionEDep::storeElement(const std::string& pathInSection,
					    const pen_parserData& element,
					    const unsigned){

  if(pathInSection.compare("spatial/xmin") == 0){
    xmin = element;
  }
  else if(pathInSection.compare("spatial/xmax") == 0){
    xmax = element;
  }
  else if(pathInSection.compare("spatial/nx") == 0){
    nx = element;
  }
  else if(pathInSection.compare("spatial/ymin") == 0){
    ymin = element;
  }
  else if(pathInSection.compare("spatial/ymax") == 0){
    ymax = element;
  }
  else if(pathInSection.compare("spatial/ny") == 0){
    ny = element;
  }
  else if(pathInSection.compare("spatial/zmin") == 0){
    zmin = element;
  }
  else if(pathInSection.compare("spatial/zmax") == 0){
    zmax = element;
  }
  else if(pathInSection.compare("spatial/nz") == 0){
    nz = element;
  }
  else if(pathInSection.compare("detector") == 0){
    kdet = element;
  }
  else if(pathInSection.compare("time/nbins") == 0){
    nt = element;
  }
  else if(pathInSection.compare("time/min") == 0){
    tmin = element;
  }
  else if(pathInSection.compare("time/max") == 0){
    tmax = element;
  }
  else if(pathInSection.compare("printBins") == 0){
    printBins = element;
  }
  else if(pathInSection.compare("printCoord") == 0){
    printCoord = element;
  }  
  else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
  
}

int tallyReader_DetectionEDep::storeString(const std::string& pathInSection,
					   const std::string& element,
					   const unsigned verbose){

  if(pathInSection.compare("particle") == 0){

    if(element.compare("all") == 0){
      //All particles are registered
      ipar = constants::nParTypes;
    }
    else{
      //Only the specified particle type is registered
      ipar = particleID(element.c_str());
      if(ipar == constants::nParTypes){
	if(verbose > 0){
	  printf("Unknown particle type '%s'. Available particles are:\n", element.c_str());
	  for(size_t i = 0; i < constants::nParTypes; ++i){
	    printf(" - %s\n", particleName(i));
	  }
	}
	return errors::UNKNOWN_PARTICLE;
      }
    }
  }
  else if(pathInSection.compare("binding") == 0){
    bodyBind = element;    
  }
  else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
}

//Tally functions
void pen_DetectionEDep::tally_step(const unsigned long long nhist,
				   const pen_KPAR /*kpar*/,
				   const pen_particleState& state,
				   const tally_StepData& stepData){
                  

  //Record energy deposition during step if needed
  if(detect){

    if(stepData.softDE == 0.0){return;}  //Nothing to score
    vector3D<double> pos(stepData.softX, stepData.softY, stepData.softZ);
    transformAndAdd(pos, state.PAGE, stepData.softDE*state.WGHT, nhist);
  }
}


void pen_DetectionEDep::tally_localEdep(const unsigned long long nhist,
					const pen_KPAR /*kpar*/,
					const pen_particleState& state,
					const double dE){

  if(detect){
    vector3D<double> pos(state.X, state.Y, state.Z);
    transformAndAdd(pos, state.PAGE, dE*state.WGHT, nhist);
  }
}

void pen_DetectionEDep::tally_interfCross(const unsigned long long /*nhist*/,
					  const unsigned kdet,
					  const pen_KPAR kpar,
					  const pen_particleState& /*state*/){
        
  //Check if the particle is within the detector
  if(idet == kdet && (iPar >= constants::nParTypes || kpar == iPar)){
    detect = true; //Valid particle Inside detector
  }
  else{
    detect = false;
  }
}


void pen_DetectionEDep::tally_move2geo(const unsigned long long nhist,
				       const unsigned kdet,
				       const pen_KPAR kpar,
				       const pen_particleState& state,
				       const double /*dsef*/,
				       const double /*dstot*/){
  if(state.MAT > 0){
    //The particle reach the geometry system after sampled and moved from void region
    if(idet == kdet && (iPar >= constants::nParTypes || kpar == iPar)){
      detect = true; //Valid particle Inside detector
      //Score particle energy to compensate the extracted during beginPart call
      vector3D<double> pos(state.X, state.Y, state.Z);
      transformAndAdd(pos, state.PAGE, state.E*state.WGHT, nhist);
    }
    else{
      detect = false;
    }
  }
}


void pen_DetectionEDep::tally_beginPart(const unsigned long long nhist,
					const unsigned kdet,
					const pen_KPAR kpar,
					const pen_particleState& state){
        
  if(idet == kdet && (iPar >= constants::nParTypes || kpar == iPar)){
    detect = true;   //Valid particle Inside detector
    
    // "Extract" particle energy from material
    vector3D<double> pos(state.X, state.Y, state.Z);
    transformAndAdd(pos, state.PAGE, -state.E*state.WGHT, nhist);
  }
  else {
    detect = false;
  }
}

void pen_DetectionEDep::tally_sampledPart(const unsigned long long nhist,
					  const unsigned long long /*dhist*/,
					  const unsigned kdet,
					  const pen_KPAR kpar,
					  const pen_particleState& state){

  //Check if the particle has been sampled at void volume.
  //If this happens, energy will be added within move2geo call
  if(state.MAT > 0){    

    if(idet == kdet && (iPar >= constants::nParTypes || kpar == iPar)){
      detect = true; //Valid particle Inside detector
      //Score particle energy to compensate the extracted during beginPart call
      vector3D<double> pos(state.X, state.Y, state.Z);
      transformAndAdd(pos, state.PAGE, state.E*state.WGHT, nhist);
    }
    else{
      detect = false;
    }
  }
}

int pen_DetectionEDep::configure(const wrapper_geometry& geometry,
				 const abc_material* const /*materials*/[constants::MAXMAT],
				 const pen_parserSection& config,
				 const unsigned verbose){

  //Read material information from config section
  tallyReader_DetectionEDep reader;
  int err = reader.read(config,verbose);
  if(err != tallyReader_DetectionSpatialDistrib::SUCCESS){
    return err;
  }

  //Check body binding
  iBodyBind = geometry.getBodies();
  if(!reader.bodyBind.empty()){
    //Check if the specified body exists
    iBodyBind = geometry.getIBody(reader.bodyBind.c_str());
    if(iBodyBind >= geometry.getBodies()){
      if(verbose > 0){
	printf("pen_DetectionEDep: configure: Error: Body '%s' doesn't "
	       "exists in the loaded geometry. Unable to bind tally to body transforms\n",
	       reader.bodyBind.c_str());
      }
      return -1;
    }
    //Save geometry pointer
    geo = &geometry;
  }
  else{
    //No binding defined
    geo = nullptr;
  }

  //Save detector
  idet = reader.kdet;

  //Save particle type to be tallied
  iPar = reader.ipar;

  //Save printing options
  printBins = reader.printBins;
  printCoord = reader.printCoord;
  
  //Configure the corresponding measurement
  
  measure.description = "Tally energy deposition in detector " + std::to_string(idet);
  
  measure.setDimHeader(0, "x (cm)");
  measure.setDimHeader(1, "y (cm)");
  measure.setDimHeader(2, "z (cm)");
  measure.setDimHeader(3, "t (s)");
  measure.setValueHeader("E(eV/hist)");
    
  //Check if energy must be tallierd
  measure.initFromLists({reader.nx, reader.ny, reader.nz, reader.nt},
			{penred::measurements::limitsType(reader.xmin, reader.xmax),
			 penred::measurements::limitsType(reader.ymin, reader.ymax),
			 penred::measurements::limitsType(reader.zmin, reader.zmax),
			 penred::measurements::limitsType(reader.tmin, reader.tmax)});

  //Set measure to dump
  toDump(measure);

  if(verbose > 1){
    //Print summary
    printf("\n%s\n", measure.stringifyInfo().c_str());
    //Print binding
    if(geo != nullptr){
      printf("Binded to object '%s'\n",
	     reader.bodyBind.c_str());
    }
  }

  return 0;
}

void pen_DetectionEDep::saveData(const unsigned long long nhist) const{

  std::string filename("spatial-detection-");
  filename += std::to_string(idet);
  filename += ".dat";

  FILE* out = nullptr;
  out = fopen(filename.c_str(), "w");

  measure.print(out, nhist, 2, printCoord, printBins);
  fclose(out);

  //Save cumulative results in time
  filename.assign("spatial-detection-cummulative-");
  filename += std::to_string(idet);
  filename += ".dat";

  out = nullptr;
  out = fopen(filename.c_str(), "w");
  
  measure.cummulative().print(out, nhist, 2, printCoord, printBins);
  fclose(out);
}

int pen_DetectionEDep::sumTally(const pen_DetectionEDep& tally){

  return measure.add(tally.measure);
  
}

REGISTER_COMMON_TALLY(pen_DetectionEDep)
