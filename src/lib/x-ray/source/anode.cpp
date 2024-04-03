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


#include "anode.hh"

namespace penred{

  namespace xray{

    //Get anode geometry file at compile time
    const char* preloadGeos::anodeGeoFile = {
#include "baseAnode.geo"
    };

    void runAnodeSimulation(const unsigned long long nHists,
			    const unsigned long long maxParticles,
			    const double Einit,
			    const double beamRad,
			    const pen_context& context,
			    const pen_VRCluster<pen_state_gPol>& photonVR,
			    std::vector<detectedPart>& results,
			    int& seed1, int& seed2,
			    const double colAngle,
			    const bool onlyPhotons){

      //Constants
      constexpr double pi = 3.141592653589793;
      constexpr double pi2 = 2.0 * pi;
      
      //Create random generator
      pen_rand random;
      random.setSeeds(seed1, seed2);

      //Create particle stacks for electrons, positrons and gamma
      pen_particleStack<pen_particleState> stackE;
      pen_particleStack<pen_particleState> stackP;
      pen_particleStack<pen_state_gPol> stackG;
      
      //Create particle instances
      pen_betaE betaE(context,stackE,stackG);
      pen_gamma gamma(context,stackE,stackP,stackG);
      pen_betaP betaP(context,stackE,stackG,stackP);

      //Register x-ray splitting VR
      gamma.registerSpecificVR(photonVR);

      //Create common initial state 
      pen_particleState genState;
      genState.X = 0.0;
      genState.Y = +10.0;
      genState.Z = 0.0;

      genState.U = 0.0;
      genState.V = -1.0;
      genState.W = 0.0;
      
      genState.E = Einit;

      //Locate the particle
      context.readGeometry()->locate(genState);

      //Init history counter
      unsigned long long hist = 0;

      //Init last history registered counter
      unsigned long long lastHist = 0;

      //Calculate maximum deviation. Note that scoring particles
      //come with -Z direction
      const double colAngleRad = colAngle*pi/180.0;
      const double minW = -cos(colAngleRad);
      
      //Define scoring function
      simulation::tallyFuncType f;
      if(onlyPhotons){
	f = [&hist, &lastHist, &results, minW]
	  (const pen_particleState& state,
	   const pen_KPAR kpar,
	   const unsigned long long,
	   const int val){
	  if(val == 1){
	    if(kpar == PEN_PHOTON){

	      if(state.W < minW){ //Note that scoring particles has W < 0	      
		results.emplace_back(state,
				     hist-lastHist,
				     kpar);
		lastHist = hist;
	      }
	    }
	  }
	};
      }else{
	f = [&hist, &lastHist, &results, minW]
	  (const pen_particleState& state,
	   const pen_KPAR kpar,
	   const unsigned long long,
	   const int val){
	  if(val == 1){
	    if(state.W < minW){ //Note that scoring particles has W < 0	      
	    
	      results.emplace_back(state,
				   hist-lastHist,
				   kpar);
	      lastHist = hist;
	    }
	  }
	};	
      }
      
      //Simulation loop
      while(hist < nHists && results.size() < maxParticles){

	//Increase history counter
	++hist;

	//Sample particle position in the beam radius
	double r = beamRad * sqrt(random.rand());
	double theta = random.rand() * pi2; 

	genState.X = r * cos(theta);
	genState.Z = r * sin(theta);

	//Copy generated state to the beta instance
	stateCopy(betaE.getState(), genState);
	
	//Simulate shower
	simulation::simulateShowerKnownCond(hist,
					    random,
					    simulation::finishTypes::DETECTOR_REACHED, 1,
					    f,
					    betaE, gamma, betaP);	
      }

      //Update seeds
      random.getSeeds(seed1,seed2);      
    }
  
    int simAnode(const char* matFilename,
		 const double eEnergy,
		 const double eMin,
		 const double focalSpot,
		 const double angle,
		 const unsigned long long nHists,
		 const unsigned long long maxParticles,
		 double& dReg,
		 std::vector<detectedPart>& results,
		 const double colAngle,
		 const bool onlyPhotons,
		 const unsigned verbose,
		 const unsigned threads2Use){

      //Simulates a monoenergetic electron beam aiming to an anode of the specified material.
      //Returns the particles registered at the specified distance
      //
      // Input:
      // 
      // matFilename : Anode material filename
      // eEnergy     : Energy of the monoenergetic electron beam
      // eMin        : Minimum energy to register
      // focalSpot   : Anode effective focal spot
      // angle       : Anode angle in deg
      // nHists      : Number of histories to simulate
      // colAngle    : Collimation angle in Deg. Particles with a direction angle (compared with -Z) greater than the specified value will be ignored
      // onlyPhotons : Flags if only photons, or all particles, must be registered 
      //
      // Output:
      //
      // dReg  : Distance from the beam center where the particles have been recorded
      // results: Vector with pairs <history increment, particleState> storing produced particles
      //
      // return: Returns 0 on success, non zero values on failure
      //

      constexpr double pi = 3.141592653589793;
      constexpr double pi05 = pi/2.0;
      constexpr double deg2rad = pi/180.0;
      
      // ** Parameters check
      //***********************
  
      //Check specified minimum energy
      if(eMin < 50.0){
	if(verbose > 0)
	  printf("simAnode: Error: Minimum energy must be greater than 50 eV.\n");
	return -1;
      }

      //Check beam energy
      if(eEnergy < 50.0){
	if(verbose > 0)
	  printf("simAnode: Error: Beam energy must be greater than 50 eV.\n");
	return -1;
      }
      if(eEnergy >= 1.0e6){
	if(verbose > 0)
	  printf("simAnode: Error: Beam energy must be lesser than 1 GeV.\n");
	return -1;
      }

      if(eMin >= eEnergy){
	if(verbose > 0)
	  printf("simAnode: Error: Minimum scoring energy greater "
		 "than initial beam energy.\n");
	return -1;
      }

      //Check focal spot
      if(focalSpot <= 0.0){
	if(verbose > 0)
	  printf("simAnode: Error: Invalid focal sport. Must be greater than zero.\n");
	return -1;
      }
      if(focalSpot >= 1.0){
	if(verbose > 0)
	  printf("simAnode: Error: Invalid focal spot. Must be lesser than 1cm");
	return -1;
      }

      //Check angle
      if(angle < 0.0 || angle >= 90.0){
	if(verbose > 0)
	  printf("simAnode: Error: Anode angle must be in the range [0,90) Deg.\n");
	return -1;
      }

      const double angleRad = deg2rad*angle;
      
      //Check number of histories
      if(nHists == 0){
	if(verbose > 0)
	  printf("simAnode: Error: No histories must be simulated.\n");
	return -1;
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
	  printf("simAnode: Warning: Number of threads has been specified,"
		 " but the code has been compiled with no multithreading support.\n"
		 " Only one thread will be used.\n");
      }
      nThreads = 1;
#endif

      //Calculate beam radius
      const double beamRad = focalSpot*tan(pi05-angleRad);

      //Set dreg according to anode geometry
      dReg = 0.6;

      if(verbose > 1){
	printf("          x-ray kvp    : %.2f\n"
	       "Minimum tallied energy : %.2f keV\n"
	       "Efective focal spot    : %.4f\n"
	       "  Anode angle (Deg)    : %.4f\n"
	       "   Histories to sim    : %llu\n",
	       eEnergy/1000.0, eMin/1000.0, focalSpot, angle,nHists);
      }

      // ** Create simulation context  
      //*******************************
  
      //Create a context
      pen_context contextSim;

      //Create context configuration
      pen_parserSection contextConf;
      contextConf.set("context-log", "context-simAnode.rep");

      //Materials
      contextConf.set("materials/anode/number", 1);
      contextConf.set("materials/anode/eabs/electron", eMin);
      contextConf.set("materials/anode/eabs/positron", eMin);
      contextConf.set("materials/anode/eabs/gamma", eMin);
      contextConf.set("materials/anode/C1", 0.05);
      contextConf.set("materials/anode/C2", 0.05);
      contextConf.set("materials/anode/WCC", std::min(5e3,eEnergy/100.0));
      contextConf.set("materials/anode/WCR", std::min(5e3,eEnergy/100.0));
      contextConf.set("materials/anode/filename", matFilename);

      //VR
      contextConf.set("VR/IForcing/bremss/particle", "electron");
      contextConf.set("VR/IForcing/bremss/interaction", BETAe_HARD_BREMSSTRAHLUNG);
      contextConf.set("VR/IForcing/bremss/factor", 400);
      contextConf.set("VR/IForcing/bremss/min-weight", 0.1);
      contextConf.set("VR/IForcing/bremss/max-weight", 2.0);
      contextConf.set("VR/IForcing/bremss/bodies/anode", true);

      contextConf.set("VR/IForcing/innerShell/particle", "electron");
      contextConf.set("VR/IForcing/innerShell/interaction", BETAe_HARD_INNER_SHELL);
      contextConf.set("VR/IForcing/innerShell/factor", 400);
      contextConf.set("VR/IForcing/innerShell/min-weight", 0.1);
      contextConf.set("VR/IForcing/innerShell/max-weight", 2.0);
      contextConf.set("VR/IForcing/innerShell/bodies/anode", true);

      contextConf.set("VR/bremss/split4/splitting", 4);
      contextConf.set("VR/bremss/split4/bodies/anode", true);      

      // ** Init context
      //********************
      pen_parserSection matInfo;
      if(contextSim.configure(eEnergy,
			      contextConf,
			      matInfo,
			      verbose) != pen_context::SUCCESS){
	printf("simAnode: Error at simulation context initialization. "
	       "See context report.\n");
	return -3;
      }

      // ** Geometry
      //**************

      //Create a mesh geometry instance
      pen_meshBodyGeo* geometry = new pen_meshBodyGeo;

      //Set geometry file in geometry instance
      geometry->preloadGeo.assign(preloadGeos::anodeGeoFile);

      //Set the geometry configuration
      pen_parserSection config;

      //Set inmemory geometry
      config.set("memory-file", true);
      
      //Enlarge anode back face 1 cm
      config.set("transforms/anode/enlargeY/index", 0);
      config.set("transforms/anode/enlargeY/vertex-group", "back");
      config.set("transforms/anode/enlargeY/transforms/enlargeY/type" , "TRANSLATION_Y");
      config.set("transforms/anode/enlargeY/transforms/enlargeY/index", 0);
      config.set("transforms/anode/enlargeY/transforms/enlargeY/ds   ", -1.0);

      //Set anode angle
      // 
      //  ------      --------
      //  |    |      |     /
      //  |    |  ->  |    /
      //  ------      -----
      //
      const double ds = std::tan(angleRad)/2.0; //Notice that the anode height is 1cm
      config.set("transforms/anode/setAngleTop/index", 1);
      config.set("transforms/anode/setAngleTop/vertex-group", "front_up");
      config.set("transforms/anode/setAngleTop/transforms/forward/type" , "TRANSLATION_Y");
      config.set("transforms/anode/setAngleTop/transforms/forward/index", 0);
      config.set("transforms/anode/setAngleTop/transforms/forward/ds   ", ds);

      config.set("transforms/anode/setAngleBot/index", 2);
      config.set("transforms/anode/setAngleBot/vertex-group", "front_down");
      config.set("transforms/anode/setAngleBot/transforms/forward/type" , "TRANSLATION_Y");
      config.set("transforms/anode/setAngleBot/transforms/forward/index", 0);
      config.set("transforms/anode/setAngleBot/transforms/forward/ds   ", -ds);

      //Set detector object as perfect absorber
      config.set("eabs/detector/electron", 1.0e35);      
      config.set("eabs/detector/positron", 1.0e35);      
      config.set("eabs/detector/gamma"   , 1.0e35);

      //Assign kdet 1 to detector object
      config.set("kdet/detector", 1);

      //Set dsmax for anode object
      config.set("dsmax/anode", 2.0e-2);

      //Configure geometry
      if(geometry->configure(config,verbose) != PEN_MESHBODY_GEO_SUCCESS){
	if(verbose > 0){
	  printf("Unexpected Error: Unable to construct the geometry. "
		 "Please, report this error\n");
	}
	delete geometry;
	return -4;
      }

      //Set the geometry to the simulation context
      contextSim.setGeometry(geometry);

      //Run context configuration step with geometry
      if(contextSim.configureWithGeo(contextConf,
				     verbose) != pen_context::SUCCESS){
	printf("simAnode: Error at simulation context initialization with geometry. "
	       "Report this error.\n");
	return -5;
      }

      // ** Convigure VR
      //*******************

      //x-ray splitting
      pen_parserSection configVRsplitting;
      configVRsplitting.set("x-ray/type", "XRAY_SPLITTING");
      configVRsplitting.set("x-ray/bodies/anode/splitting", 4);

      pen_VRCluster<pen_state_gPol> photonVR;
      photonVR.name.assign("Photon-VR");
      photonVR.configure(configVRsplitting,*geometry,verbose);
      if(photonVR.configureStatus() != 0){
	printf("Unexpected error: Unable to configure "
	       "photon variance reduction. Please, report this.");
	delete geometry;
	return -5;
      }
      
      // ** Configure seeds
      //**********************

      //Set initial seeds for each thread
      std::vector<int> seeds1(nThreads);
      std::vector<int> seeds2(nThreads);

      for(unsigned ith = 0; ith < nThreads; ++ith){
	rand0(ith, seeds1[ith], seeds2[ith]);
      }

      // ** Tally
      //***********

      //Create a results vector for each thread
      std::vector<std::vector<detectedPart>> localResults(nThreads);
      // ** Run simulations
      //*********************

#ifdef _PEN_USE_THREADS_
      std::vector<std::thread> threads;
      for(size_t ith = 0; ith < nThreads; ++ith){
	threads.emplace_back(runAnodeSimulation,
			     nHists/nThreads,
			     maxParticles/nThreads+1,
			     eEnergy,
			     beamRad,
			     std::ref(contextSim),
			     std::ref(photonVR),
			     std::ref(localResults[ith]),
			     std::ref(seeds1[ith]), std::ref(seeds2[ith]),
			     colAngle,
			     onlyPhotons);
      }

      //Wait until threads finish and join results
      for(size_t ith = 0; ith < nThreads; ++ith){
	threads[ith].join();
	results.insert(results.end(),
		       localResults[ith].begin(), localResults[ith].end());
	localResults[ith].clear();
      }
#else
      //Simulate using a single thread
      runAnodeSimulation(nHists, maxParticles,
			 eEnergy, beamRad, context, photonVR,
			 results, seeds1[0], seeds2[0],
			 colAngle,
			 onlyPhotons);
#endif
      delete geometry;
      return 0;
    }

  };
};
