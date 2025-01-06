 

//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
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
//

#include "base_simulation_functions.hh"

namespace penred{

  namespace simulation{

    //End line defition
    constexpr simConfig::Endl simConfig::endl;

    simConfig::simConfig() :
      logger(penred::logs::SIMULATION),
      iSeed1(1), iSeed2(1),
      lSeed1(1), lSeed2(1),
      firstSourceIndex(0),
      simulatedHistsInFirstSource(0),
      actualSource(""),
      maxSimTime(1000000000000000),
      iThread(0),
      dumpTime(1000000000000000),
      dumpFilename("dump"),
      writePartial(false),
      fSimFinish(simConfig::noFinishSim),
      verbose(2){}

    simConfig::simConfig(const unsigned& iThreadIn,
			 const long long int& dumpTimeIn,
			 const std::string& dumpFilenameIn,
			 const bool& writePartialIn,
			 const long long int& maxSimTimeIn,
			 const int& iSeed1In, const int& iSeed2In,
			 const unsigned& verboseIn) :
      logger(penred::logs::SIMULATION),
      iSeed1(iSeed1In), iSeed2(iSeed2In),
      lSeed1(iSeed1In), lSeed2(iSeed2In),
      firstSourceIndex(0),
      simulatedHistsInFirstSource(0),
      actualSource(""),
      maxSimTime(maxSimTimeIn),
      iThread(iThreadIn),
      dumpTime(dumpTimeIn),
      dumpFilename(dumpFilenameIn),
      writePartial(writePartialIn),
      fSimFinish(simConfig::noFinishSim),
      verbose(verboseIn){}


    std::string simConfig::stringifyConfig() const{
      std::string result;

      result += "Maximum time (s): " + std::to_string(maxSimTime/1000) + "\n";      
      result += "Dump time    (s): " + std::to_string(dumpTime/1000) + "\n";
      result += "Dump filename   : " + dumpFilename + "\n";
      result += "Verbose level   : " + std::to_string(verbose) + "\n";
      if(writePartial)
	result += "Partial results: enabled\n";
      else
	result += "Partial results: disabled\n";

      return result;
    }
    
    int simConfig::configure(const pen_parserSection& config){

      // Get verbose level
      //********************
      int auxVerbose;
      if(config.read("verbose",auxVerbose) == INTDATA_SUCCESS){
	auxVerbose = std::max(0,auxVerbose);
	verbose = static_cast<unsigned>(auxVerbose);
      }

      // Get initial seeds
      //*******************************
      if(config.read("seed1",iSeed1) != INTDATA_SUCCESS){
	iSeed1 = 1;
      }
      if(config.read("seed2",iSeed2) != INTDATA_SUCCESS){
	iSeed2 = 1;
      }

      if(iSeed1 <= 0 || iSeed2 <= 0){
	return errors::ERROR_INVALID_SEEDS;
      }

      lSeed1 = iSeed1;
      lSeed2 = iSeed2;

      // Write partial results option
      //*******************************
      if(config.read("partial-results",writePartial) != INTDATA_SUCCESS){
	writePartial = false;
      }

      // Get maximum simulation time
      //*******************************
      double auxMaxSimTime;
      if(config.read("max-time",auxMaxSimTime) == INTDATA_SUCCESS){
	if(auxMaxSimTime <= 0.0 || auxMaxSimTime > 1.0e12){
	  //No maximum simulation time
	  maxSimTime = 1000000000000000;
	}else{
	  maxSimTime = static_cast<long long int>(1000*auxMaxSimTime);
	}
      }else{
	//No maximum simulation time
	maxSimTime = 1000000000000000;
      }

      // Get time between dumps
      //*******************************
      double auxDumpTime;
      if(config.read("dump-interval",auxDumpTime) == INTDATA_SUCCESS){
	if(auxDumpTime <= 0.0 || auxDumpTime > 1.0e12){
	  //Dump "disabled"
	  dumpTime = 1000000000000000;
	}else{
	  dumpTime = static_cast<long long int>(1000*auxDumpTime);
	}
      }else{
	//No dump time specified
	dumpTime = 1000000000000000;
      }

      // Get dump filename
      //*******************************
      if(config.read("dump2write",dumpFilename) != INTDATA_SUCCESS){
	dumpFilename = "dump.dat";
      }

      return errors::SUCCESS;
    }
    
  } // namespace simulation
  
} // namespace penred
