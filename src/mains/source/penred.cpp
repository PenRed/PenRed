//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2024-2025 Vicent Giménez Alventosa
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//


#include "pen_simulation.hh"

int main(int argc, char** argv){

  std::string configFilename("config.in");
  if(argc < 2){
    printf("No argument provided, the default filename "
	   "(config.in) will be used to configure the simulation.\n"
	   "To specify a different file use: %s config-filename\n",argv[0]);
  }else{
    configFilename.assign(argv[1]);
  }

  // Create simulator object
  //**************************
  penred::simulation::simulator<pen_context> simula;

  //Print version
  std::cout << simula.versionMessage() << std::endl;

  //**************************
  // Parse configuration
  //**************************

  //Parse configuration file
  pen_parserSection config;
  std::string errorLine;
  unsigned long errorLineNum;
  int err = parseFile(configFilename.c_str(),config,errorLine,errorLineNum);
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }

  //Set log files
  std::string configLogFilename;
  if(config.read("log/configuration",configLogFilename) != INTDATA_SUCCESS){
    configLogFilename.assign("config.log");
  }
  std::string simLogFilename;
  if(config.read("log/simulation",simLogFilename) != INTDATA_SUCCESS){
    simLogFilename.assign("simulation.log");
  }

  penred::logs::logger log;  
  if(!configLogFilename.empty()){
    printf("Configuration log redirected to '%s'\n", configLogFilename.c_str());
    log.setConfigurationLogFile(configLogFilename.c_str());
  }
  if(!simLogFilename.empty()){
    printf("Simulation log redirected to '%s'\n", simLogFilename.c_str());
    log.setSimulationLogFile(simLogFilename.c_str());  
  }

  // Configure simulator
  //***********************
  simula.configure(config);

  // Run simulation
  //***********************
  simula.simulate();
}
