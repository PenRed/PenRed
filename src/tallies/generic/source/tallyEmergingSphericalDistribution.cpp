
//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
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

 
#include "tallyEmergingSphericalDistribution.hh"

//Reader functions

int tallyReader_EmergingSphericalDistrib::beginSectionFamily(const std::string& pathInSection,
							     const size_t,
							     const unsigned){
  if(family == -1){ //Root subsection
    if(pathInSection.compare("particle") == 0){
      family = 0; //Move to particle family
      return errors::SUCCESS;
    }
  }

  return errors::UNHANDLED;
}

int tallyReader_EmergingSphericalDistrib::endSectionFamily(const unsigned){

  if(family == 0){ //Particles section
    family = -1; //Return to root
    return errors::SUCCESS;
  }
  
  return errors::UNHANDLED;  
}

int tallyReader_EmergingSphericalDistrib::beginSection(const std::string& name,
						       const unsigned verbose){

  if(family == 0){ //Particles section

    //Get particle ID from particle "name"
    actualPartID = particleID(name.c_str());
    if(actualPartID == ALWAYS_AT_END){
      if(verbose > 0){
	printf("Error: Unknown particle type %s\n"
	       " Known particles:\n", name.c_str());
	for(size_t i = 0; i < constants::nParTypes; ++i){
	  printf("    %s\n", particleName(i));
	}
	
      }
      return errors::UNKNOWN_PARTICLE;
    }
    
    return errors::SUCCESS;
  }
  
  return errors::UNHANDLED;  
}

int tallyReader_EmergingSphericalDistrib::endSection(const unsigned){
  return errors::SUCCESS;
}  


int tallyReader_EmergingSphericalDistrib::storeElement(const std::string& pathInSection,
						       const pen_parserData& element,
						       const unsigned){

  if(family == -1){ //Root elements
  
    if(pathInSection.compare("energy/nbins") == 0){
      ne = element;
    }
    else if(pathInSection.compare("energy/min") == 0){
      emin = element;
    }
    else if(pathInSection.compare("energy/max") == 0){
      emax = element;
    }
    else if(pathInSection.compare("theta/nbins") == 0){
      nt = element;
    }
    else if(pathInSection.compare("theta/min") == 0){
      tmin = element;
      tmin *= constants::PI/180.0;
    }
    else if(pathInSection.compare("theta/max") == 0){
      tmax = element;
      tmax *= constants::PI/180.0;
    }
    else if(pathInSection.compare("phi/nbins") == 0){
      np = element;
    }
    else if(pathInSection.compare("phi/min") == 0){
      pmin = element;
      pmin *= constants::PI/180.0;
    }
    else if(pathInSection.compare("phi/max") == 0){
      pmax = element;
      pmax *= constants::PI/180.0;
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
  }
  else if(family == 0){ //Particle type subsections
    if(!pathInSection.empty())
      return errors::UNHANDLED;
    enabledPart[actualPartID] = element;
  }
  else{
    return errors::UNHANDLED;
  }
  
  return errors::SUCCESS;
  
}


void pen_EmergingSphericalDistrib::tally_move2geo(const unsigned long long nhist,
						  const unsigned /*kdet*/,
						  const pen_KPAR kpar,
						  const pen_particleState& state,
						  const double /*dsef*/,
						  const double /*dstot*/){
     
  if(state.MAT == 0 && enabled[kpar]){
    //Particle scape the geometry
    //The particle doesn't reached the material system

    const double twopi = 2.0*constants::PI;
    
    //Calculate theta and phi
    double theta, phi;
    
    theta = acos(state.W);
    if(fabs(state.V) > 1.0e-16 || fabs(state.U) > 1.0e-16){
      phi = atan2(state.V,state.U);
    }
    else{
      phi = 0.0;
    }
    
    if(phi < 0.0){
      phi = twopi + phi;
    }    
    
    results[kpar].add({state.E, theta, phi}, state.WGHT, nhist);
  }
    
}
 
 
void pen_EmergingSphericalDistrib::tally_matChange(const unsigned long long nhist,
						   const pen_KPAR kpar,
						   const pen_particleState& state,
						   const unsigned /*prevMat*/){

  if(state.MAT == 0 && enabled[kpar]){  
    //Particle scape the geometry
    //The particle scaped from the material system    

    const double twopi = 2.0*constants::PI;
    
    //Calculate theta and phi
    double theta, phi;
    
    theta = acos(state.W);
    if(fabs(state.V) > 1.0e-16 || fabs(state.U) > 1.0e-16){
      phi = atan2(state.V,state.U);
    }
    else{
      phi = 0.0;
    }
    
    if(phi < 0.0){
      phi = twopi + phi;
    }
    
    results[kpar].add({state.E, theta, phi}, state.WGHT, nhist);    
  }
}


void pen_EmergingSphericalDistrib::flush(){
  for(size_t ip = 0; ip < constants::nParTypes; ++ip){
    if(enabled[ip]){
      results[ip].flush();
    }
  }
}


void pen_EmergingSphericalDistrib::tally_endSim(const unsigned long long /*nhist*/){        
  flush();
}

int pen_EmergingSphericalDistrib::configure(const wrapper_geometry& /*geometry*/,
					    const abc_material* const /*materials*/[constants::MAXMAT],
					    const pen_parserSection& config,
					    const unsigned verbose){
                   
  //Read material information from config section
  tallyReader_EmergingSphericalDistrib reader;
  int err = reader.read(config,verbose);
  if(err != tallyReader_DetectionSpatialDistrib::SUCCESS){
    return err;
  }

  //Save enabled particles
  enabled = reader.enabledPart;
  
  //Save printing options
  printBins = reader.printBins;
  printCoord = reader.printCoord;

  //Configure the corresponding measurements
    
  //Init results measurements
  bool printedInfo = false;
  for(size_t ip = 0; ip < constants::nParTypes; ++ip){
    if(enabled[ip]){

      results[ip].description = "Tally emerging spherical distribution";
  
      results[ip].setDimHeader(0, "E (eV)");
      results[ip].setDimHeader(1, "theta (rad)");
      results[ip].setDimHeader(2, "phi (rad)");
      results[ip].setValueHeader("Prob(1/hist)");
  
      results[ip].
	initFromLists({reader.ne, reader.nt, reader.np},
		      {penred::measurements::limitsType(reader.emin, reader.emax),
		       penred::measurements::limitsType(reader.tmin, reader.tmax),
		       penred::measurements::limitsType(reader.pmin, reader.pmax)});

      if(verbose > 1 && !printedInfo){
	//Print summary
	printf("\n%s\n", results[ip].stringifyInfo().c_str());
	printedInfo = true;
      }
    }
  }
  
  //Register dumps
  for(size_t ip = 0; ip < constants::nParTypes; ++ip){
    if(enabled[ip]){  
      toDump(results[ip]);
    }
  }
  
  return 0;  
}
               


void pen_EmergingSphericalDistrib::saveData(const unsigned long long nhist) const{
  
  for(size_t ip = 0; ip < constants::nParTypes; ++ip){
    if(enabled[ip]){
      std::string filename("emerging-spherical-distrib-");
      filename += particleName(ip);
      filename += ".dat";
      FILE* out = nullptr;
      out = fopen(filename.c_str(), "w");
      
      results[ip].print(out, nhist, 2, printCoord, printBins);

      fclose(out);
    }
  }
}


int pen_EmergingSphericalDistrib::sumTally(const pen_EmergingSphericalDistrib& tally){

  //Check enabled particles
  for(size_t ip = 0; ip < constants::nParTypes; ++ip){
    if(enabled[ip] != tally.enabled[ip])
      return -1;
  }
  
  int err = 0;
  for(size_t ip = 0; ip < constants::nParTypes; ++ip){
    if(enabled[ip]){
      err = results[ip].add(tally.results[ip]);
      if(err != 0)
	return err;
    }
  }
  return err;
}


REGISTER_COMMON_TALLY(pen_EmergingSphericalDistrib, EMERGING_SPHERICAL_DISTRIB)
 
        
    
    
    
    
    
    
    
    
    
      


      
      
      
      
      
      
      
      
      
      
      
