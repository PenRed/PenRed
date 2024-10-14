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
//    
//

#include "phantom.hh"

namespace penred{

  namespace xray{

    int simCylPhantom(const char* matFilename,
		      double eMin,
		      const vector3D<double>& isocenter,		      
		      const double phantomDiam,
		      const double iso2Detector,
		      const double detDx, const double detDy,
		      std::vector<detectedPart>& particlesIn,
		      const bool onlyPhotons,
		      const unsigned verbose,
		      const unsigned threads2Use,
		      const std::string geoFilename){

      // Simulates the transport of the provided particles through a cylindrical
      // phantom, and saves the detecteded ones. The geometry is as follows:
      //
      //          _________                             n  Z
      //         /         \                            |  
      //        /           \                           |
      //       |  isocenter  |    <--- Phantom          |----------> X
      //       |      *      |
      //        \     n     /
      //         \____|____/
      //              |
      //              | <-------------- isocenter to detector
      //        ______v_______
      //       |______________|   <--- Detector
      //
      //
      //
      // Input:
      //
      // matFilename  : Phantom material file
      // eMin         : Minimum energy to be tallied
      // isocenter    : Isocenter position (x,y,z)
      // phantomDiam  : Phantom diameter in cm
      // iso2Detector : Distance isocenter to detector in cm
      // detDx, detDy : Detector size in X and Y axis
      // particlesIn  : Input particles to be simulated
      // onlyPhotons  : Enable/disable only photons transport
      // verbose      : verbose level
      // threads2Use  : Number of threads to use
      //
      //
      // output:
      //
      // particlesIn  : Tallied particles after the simulation
      //
      //
      // Returns error::SUCCESS on success, an error code otherwise
      

      if(phantomDiam <= 0.0){
	return errors::NEGATIVE_DISTANCE;
      }

      if(iso2Detector <= phantomDiam/2.0){
	return errors::INVALID_DISTANCE;
      }

      //Check minimum energy
      eMin = std::max(eMin, 50.0);

      //Ensure some particle must be simulated
      if(particlesIn.size() == 0){
	return errors::SUCCESS;
      }

      //Get maximum energy and total number of histories
      double maxE = 0.0;
      unsigned long long nHists = 1;
      for(const detectedPart& part : particlesIn){
	if(maxE < part.state.E)
	  maxE = part.state.E;
	nHists += part.dh;
      }

      if(verbose > 1){
	printf("\n + Cylindrical phantom simulation:\n"
	       "     Isocenter               : %s cm\n"
	       "     Phantom diameter        : %.4f cm\n"
	       "     Isocenter to detector   : %.4f cm\n"
	       "     Detector X size         : %.4f cm\n"
	       "     Detector Y size         : %.4f cm\n"
	       "     Maximum incident energy : %.4f keV\n"
	       "     Energy threshold        : %.4f keV\n"
	       "     %s\n",
	       isocenter.stringify().c_str(),
	       phantomDiam,
	       iso2Detector,
	       detDx,
	       detDy,
	       maxE/1000.0,
	       eMin/1000.0,
	       onlyPhotons ? "Only photons will be simulated" : "");
      }

      //Get the number of threads to use
      unsigned nThreads = threads2Use;
#ifdef _PEN_USE_THREADS_
      if(nThreads == 0){
	nThreads = std::max(static_cast<unsigned int>(2),
			    std::thread::hardware_concurrency());
      }
#else
      if(nThreads != 0){
	if(verbose > 1)
	  printf("simCylPhantom: Warning: Number of threads has been specified,"
		 " but the code has been compiled with no multithreading support.\n"
		 " Only one thread will be used.\n");
      }
      nThreads = 1;
#endif      
      
      // ** Create simulation context  
      //*******************************
  
      //Create simulation context
      std::shared_ptr<pen_context> pcontext = createContext<pen_context>();
      //Get context reference. Notice that pcontext will
      //not be released until the function ends
      pen_context& contextSim = *pcontext.get();

      //Create context configuration
      pen_parserSection contextConf;
      contextConf.set("context-log", "context-simCylPhantom.rep");

      //Materials
      contextConf.set("materials/phantom/number", 1);
      contextConf.set("materials/phantom/eabs/electron", onlyPhotons ? 1.0e35 : eMin);
      contextConf.set("materials/phantom/eabs/positron", onlyPhotons ? 1.0e35 : eMin);
      contextConf.set("materials/phantom/eabs/gamma", eMin);
      contextConf.set("materials/phantom/C1", 0.05);
      contextConf.set("materials/phantom/C2", 0.05);
      contextConf.set("materials/phantom/WCC", std::min(5e3,maxE/100.0));
      contextConf.set("materials/phantom/WCR", std::min(5e3,maxE/100.0));
      contextConf.set("materials/phantom/filename", matFilename);


      // ** Init context
      //********************
      pen_parserSection matInfo;
      if(contextSim.configure(maxE,
			      contextConf,
			      matInfo,
			      verbose > 3 ? verbose : 1) != pen_context::SUCCESS){
	printf("simAnodesimCylPhantom: Error at simulation context initialization. "
	       "See context report.\n");
	return errors::ERROR_ON_CONTEXT_INITIALIZATION;
      }

      // ** Geometry
      //**************

      //Create a mesh geometry instance
      pen_meshBodyGeo* geometry = new pen_meshBodyGeo;

      //Set geometry file in geometry instance
      geometry->preloadGeo.assign(preloadGeos::phantomCylGeoFile);

      if(!geoFilename.empty()){
	FILE* fout = fopen(geoFilename.c_str(), "w");
	if(fout != nullptr){
	  fprintf(fout,"%s", preloadGeos::phantomCylGeoFile);
	  fclose(fout);
	}
      }
      
      //Set the geometry configuration
      pen_parserSection config;

      //Set inmemory geometry
      config.set("memory-file", true);

      //Set detector object as perfect absorber
      config.set("eabs/detector/electron", 1.0e35);      
      config.set("eabs/detector/positron", 1.0e35);      
      config.set("eabs/detector/gamma"   , 1.0e35);

      //Set detector kdet
      config.set("kdet/detector", 1);
      
      // * Phantom
      
      //Scale phantom to fit the specified diameter and move it
      config.set("transforms/phantom/conform/index", 0);
      config.set("transforms/phantom/conform/vertex-group", "all");

      //Fit the diameter
      config.set("transforms/phantom/conform/transforms/setDia/type" , "SCALE_XZ");
      config.set("transforms/phantom/conform/transforms/setDia/index", 0);
      config.set("transforms/phantom/conform/transforms/setDia/factor", phantomDiam);

      //Fit the detector length in Y axis
      config.set("transforms/phantom/conform/transforms/setLength/type" , "SCALE_Y");
      config.set("transforms/phantom/conform/transforms/setLength/index", 1);
      config.set("transforms/phantom/conform/transforms/setLength/factor", 2.0*detDy);

      //Move to isocenter
      config.set("transforms/phantom/conform/transforms/move/type" , "TRANSLATION");
      config.set("transforms/phantom/conform/transforms/move/index", 2);
      config.set("transforms/phantom/conform/transforms/move/u",  isocenter.x);
      config.set("transforms/phantom/conform/transforms/move/v",  isocenter.y);
      config.set("transforms/phantom/conform/transforms/move/w",  isocenter.z);
      config.set("transforms/phantom/conform/transforms/move/ds", isocenter.mod());

      // * Detector
      
      //Scale detector to fit the specified size and move it
      config.set("transforms/detector/conform/index", 0);
      config.set("transforms/detector/conform/vertex-group", "all");

      // Size in X
      config.set("transforms/detector/conform/transforms/setSizeX/type" , "SCALE_X");
      config.set("transforms/detector/conform/transforms/setSizeX/index", 0);
      config.set("transforms/detector/conform/transforms/setSizeX/factor", detDx);

      // Size in Y
      config.set("transforms/detector/conform/transforms/setSizeY/type" , "SCALE_Y");
      config.set("transforms/detector/conform/transforms/setSizeY/index", 1);
      config.set("transforms/detector/conform/transforms/setSizeY/factor", detDy);
      
      // Translate the detector to isocented
      config.set("transforms/detector/conform/transforms/move/type" , "TRANSLATION");
      config.set("transforms/detector/conform/transforms/move/index", 2);
      config.set("transforms/detector/conform/transforms/move/u",  isocenter.x);
      config.set("transforms/detector/conform/transforms/move/v",  isocenter.y);
      config.set("transforms/detector/conform/transforms/move/w",  isocenter.z);
      config.set("transforms/detector/conform/transforms/move/ds", isocenter.mod());

      // Move the detector down
      config.set("transforms/detector/conform/transforms/shiftZ/type" , "TRANSLATION_Z");
      config.set("transforms/detector/conform/transforms/shiftZ/index", 3);
      config.set("transforms/detector/conform/transforms/shiftZ/ds", -iso2Detector - 0.05);

      // * World

      //Scale the world to fit all components inside and move it
      config.set("transforms/world/conform/index", 0);
      config.set("transforms/world/conform/vertex-group", "all");

      // Size in X
      const double worldDx = std::max(detDx, phantomDiam)*1.10; 
      config.set("transforms/world/conform/transforms/setSizeX/type" , "SCALE_X");
      config.set("transforms/world/conform/transforms/setSizeX/index", 0);
      config.set("transforms/world/conform/transforms/setSizeX/factor", worldDx);

      // Size in Y
      config.set("transforms/world/conform/transforms/setSizeY/type" , "SCALE_Y");
      config.set("transforms/world/conform/transforms/setSizeY/index", 1);
      config.set("transforms/world/conform/transforms/setSizeY/factor", detDy*2.1);

      // Size in Z
      const double worldDz = (phantomDiam + 2.0*iso2Detector + 2.0)*1.10; 
      config.set("transforms/world/conform/transforms/setSizeZ/type" , "SCALE_Z");
      config.set("transforms/world/conform/transforms/setSizeZ/index", 2);
      config.set("transforms/world/conform/transforms/setSizeZ/factor", worldDz);
      
      // Move it to isocented
      config.set("transforms/world/conform/transforms/move/type" , "TRANSLATION");
      config.set("transforms/world/conform/transforms/move/index", 3);
      config.set("transforms/world/conform/transforms/move/u",  isocenter.x);
      config.set("transforms/world/conform/transforms/move/v",  isocenter.y);
      config.set("transforms/world/conform/transforms/move/w",  isocenter.z);
      config.set("transforms/world/conform/transforms/move/ds", isocenter.mod());      


      // ** Configure geometry
      if(geometry->configure(config,verbose > 3 ? verbose : 1) != PEN_MESHBODY_GEO_SUCCESS){
	if(verbose > 0){
	  printf("simCylPhantom: Unexpected Error: Unable to construct "
		 "the geometry. Please, report this error\n");
	}
	delete geometry;
	return errors::ERROR_ON_GEOMETRY_INITIALIZATION;
      }

      if(!geoFilename.empty()){
	FILE* fout = fopen(("config_" + geoFilename).c_str(), "w");
	if(fout != nullptr){
	  config.set("type", "MESH_BODY");
	  config.set("input-file", geoFilename);	  
	  config.set("memory-file", false, true);
	  
	  fprintf(fout,"%s", config.stringify().c_str());
	  fclose(fout);
	}
      }      

      //Set the geometry to the simulation context
      contextSim.setGeometry(geometry);
      
      //Run context configuration step with geometry
      if(contextSim.configureWithGeo(contextConf,
				     verbose > 3 ? verbose : 1) != pen_context::SUCCESS){
	printf("simCylPhantom: Error at simulation context initialization "
	       "with geometry. Report this error.\n");
	delete geometry;
	return errors::ERROR_ON_CONTEXT_INITIALIZATION;
      }

      // ** Simulation configuraiton
      //-------------------------------

      //Create a base configuration
      penred::simulation::simConfig basicConf;
      basicConf.verbose = verbose > 3 ? verbose : 1;

      //Generate configurations for each thread
      std::vector<penred::simulation::simConfig> simConf =
	basicConf.generateThreadConfigs(nThreads);

      //Create vectors to store resulting particles in each thread
      std::vector<std::vector<detectedPart>> result(nThreads);

      //Calculate number of histories per thread
      unsigned long long nHistsPerThread = nHists/nThreads;      

      // ** Simulate
      //--------------

#ifdef _PEN_USE_THREADS_

      std::vector<std::thread> threads;

      for(size_t ith = 0; ith < nThreads; ++ith){

	threads.emplace_back([&, ith, nHists, nThreads, nHistsPerThread]
			     (){	  

	  unsigned long long hists2simulate = nHistsPerThread;
	  if(ith == nThreads-1)
	    hists2simulate += nHists % nThreads;
	  
	  simulateVectorAndDetect(simConf[ith],contextSim,
				  result[ith],
				  particlesIn,
				  ith*nHistsPerThread,
				  hists2simulate,
				  "cyl-phantom",
				  1); //Detector index
	  
	});
      }

      for(size_t ith = 0; ith < nThreads; ++ith){
	threads[ith].join();	
      }
      
#else

      simulateVectorAndDetect(simConf[0], context, result[0],
			      particlesIn,
			      nHists,
			      "cyl-phantom",
			      1); //Detector index      
#endif

      //Clear input particles and store the resulting ones
      particlesIn.clear();
      for(size_t ith = 0; ith < nThreads; ++ith){
	particlesIn.insert(particlesIn.end(),
			   result[ith].begin(), result[ith].end());
      }

      //Delete geometry
      delete geometry;

      return errors::SUCCESS;      
    }

  } // namespace xray
} // namespace penred
