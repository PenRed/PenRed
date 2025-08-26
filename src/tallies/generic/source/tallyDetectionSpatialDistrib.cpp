
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


#include "tallyDetectionSpatialDistrib.hh"

//Reader functions
int tallyReader_DetectionSpatialDistrib::storeElement(const std::string& pathInSection,
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
  else if(pathInSection.compare("energy/nbins") == 0){
    nEBins = element;
  }
  else if(pathInSection.compare("energy/emin") == 0){
    emin = element;
  }
  else if(pathInSection.compare("energy/emax") == 0){
    emax = element;
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

int tallyReader_DetectionSpatialDistrib::storeString(const std::string& pathInSection,
						     const std::string& element,
						     const unsigned verbose){

  if(pathInSection.compare("particle") == 0){
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
  }else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
}

void pen_DetectionSpatialDistrib::tally_interfCross(const unsigned long long nhist,
						    const unsigned kdet,
						    const pen_KPAR kpar,
						    const pen_particleState& state){

  if(kdet == idet && kpar == iPar){
    results.add({state.E, state.X, state.Y, state.Z}, state.WGHT, nhist);
  }
}

void pen_DetectionSpatialDistrib::tally_move2geo(const unsigned long long nhist,
						 const unsigned kdet,
						 const pen_KPAR kpar,
						 const pen_particleState& state,
						 const double /*dsef*/,
						 const double /*dstot*/){
  if(kdet == idet && kpar == iPar){
    results.add({state.E, state.X, state.Y, state.Z}, state.WGHT, nhist);
  }  
}

int pen_DetectionSpatialDistrib::configure(const wrapper_geometry& /*geometry*/,
					   const abc_material* const /*materials*/[constants::MAXMAT],
					   const pen_parserSection& config,
					   const unsigned verbose){

  //Read material information from config section
  tallyReader_DetectionSpatialDistrib reader;
  int err = reader.read(config,verbose);
  if(err != tallyReader_DetectionSpatialDistrib::SUCCESS){
    return err;
  }

  //Save detector
  idet = reader.kdet;

  //Save particle type to be tallied
  iPar = reader.ipar;

  //Save printing options
  printBins = reader.printBins;
  printCoord = reader.printCoord;
  
  //Configure the corresponding measurement
  
  results.description = "Tally spatial distribution in detector";
  
  results.setDimHeader(0, "E (eV)");
  results.setDimHeader(1, "x (cm)");
  results.setDimHeader(2, "y (cm)");
  results.setDimHeader(3, "z (cm)");
  results.setValueHeader("Prob(1/hist)");
    
  //Check if energy must be tallierd
  results.initFromLists({reader.nEBins, reader.nx, reader.ny, reader.nz},
			{penred::measurements::limitsType(reader.emin, reader.emax),
			 penred::measurements::limitsType(reader.xmin, reader.xmax),
			 penred::measurements::limitsType(reader.ymin, reader.ymax),
			 penred::measurements::limitsType(reader.zmin, reader.zmax)});

  //Set results to dump
  toDump(results);

  if(verbose > 1){
    //Print summary
    printf("\n%s\n", results.stringifyInfo().c_str());
  }
  
  return 0;
}

void pen_DetectionSpatialDistrib::saveData(const unsigned long long nhist) const{

  std::string filename("spatial-detection-");
  filename += std::to_string(idet);
  filename += ".dat";

  FILE* out = nullptr;
  out = fopen(filename.c_str(), "w");

  results.print(out, nhist, 2, printCoord, printBins);
}





int pen_DetectionSpatialDistrib::sumTally(const pen_DetectionSpatialDistrib& tally){

  return results.add(tally.results);
  
}

REGISTER_COMMON_TALLY(pen_DetectionSpatialDistrib)
