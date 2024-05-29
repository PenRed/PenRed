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
	f = [&lastHist, &results, minW]
	  (const pen_particleState& state,
	   const pen_KPAR kpar,
	   const unsigned long long histSim,
	   const int val){
	  if(val == 1){
	    if(kpar == PEN_PHOTON){

	      if(state.W < minW){ //Note that scoring particles have W < 0	      
		results.emplace_back(state,
				     histSim-lastHist,
				     kpar);
		lastHist = histSim;
	      }
	    }
	  }
	};
      }else{
	f = [&lastHist, &results, minW]
	  (const pen_particleState& state,
	   const pen_KPAR kpar,
	   const unsigned long long histSim,
	   const int val){
	  if(val == 1){
	    if(state.W < minW){ //Note that scoring particles have W < 0	      
	    
	      results.emplace_back(state,
				   histSim-lastHist,
				   kpar);
	      lastHist = histSim;
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

    void runAnodeDistribSimulation(const unsigned long long nHists,
				   const double Einit,
				   const pen_context& context,
				   const pen_VRCluster<pen_state_gPol>& photonVR,
				   measurements::measurement<double,1>& spectrum,
				   measurements::measurement<double,2>& spatialDistrib,
				   int& seed1, int& seed2,
				   const unsigned verbose){
      
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
      
      //Define scoring function
      simulation::tallyFuncType f = [&spectrum, &spatialDistrib]
	(const pen_particleState& state,
	 const pen_KPAR kpar,
	 const unsigned long long histSim,
	 const int val){
	if(val == 1){
	  if(kpar == PEN_PHOTON){	    
	    spectrum.add({state.E},state.WGHT,histSim);
	    spatialDistrib.add({state.X, state.Y}, state.WGHT, histSim);
	  }
	}
      };

      if(verbose > 1)
	printf("Starting thread simulation (%llu histories) \n", nHists);
      
      //Simulation loop
      while(hist < nHists){

	//Increase history counter
	++hist;

	//Copy initial state to the beta instance
	stateCopy(betaE.getState(), genState);
	
	//Simulate shower
	simulation::simulateShowerKnownCond(hist,
					    random,
					    simulation::finishTypes::DETECTOR_REACHED,
					    1,
					    f,
					    betaE, gamma, betaP);

	if(verbose > 1 && hist % 10000 == 0)
	  printf("Simulated: %llu/%llu \n", hist, nHists);
      }

      if(verbose > 1)
	printf("Simulated: %llu/%llu \n", hist, nHists);

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

    int simAnodeDistrib(const char* matFilename,
			const double eEnergy,
			const double eMin,
			const double pixelSize,
			const double angle,
			const unsigned long long nHists,
			double& dReg,
			measurements::measurement<double,1>& spectrum,
			measurements::measurement<double,2>& spatialDistrib,
			const double colAngle,
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
      // pixelSize   : Spatial pixel size in cm
      // angle       : Anode angle in deg
      // nHists      : Number of histories to simulate
      // colAngle    : Collimation angle in Deg. Particles with a direction angle (compared with -Z) greater than the specified value will be ignored
      //
      // Output:
      //
      // dReg  : Distance from the beam center at the anode surface to the detector
      // spectrum: Tallied photons spectrum from anode simulation
      // spatialDistrib: Tallied XY spatial distribution of photons
      //
      // return: Returns 0 on success, non zero values on failure
      //

      constexpr double pi = 3.141592653589793;
      constexpr double deg2rad = pi/180.0;
      
      // ** Parameters check
      //***********************

      if(pixelSize <= 0.0){
	if(verbose > 0)
	  printf("simAnodeDistrib: Error: Pixel size must be "
		 "greater than 0 cm.\n");
	return -1;
      }
  
      //Check specified minimum energy
      if(eMin < 50.0){
	if(verbose > 0)
	  printf("simAnodeDistrib: Error: Minimum energy must be "
		 "greater than 50 eV.\n");
	return -1;
      }

      //Check beam energy
      if(eEnergy < 50.0){
	if(verbose > 0)
	  printf("simAnodeDistrib: Error: Beam energy must be "
		 "greater than 50 eV.\n");
	return -1;
      }
      if(eEnergy >= 1.0e6){
	if(verbose > 0)
	  printf("simAnodeDistrib: Error: Beam energy must be "
		 "lesser than 1 GeV.\n");
	return -1;
      }

      if(eMin >= eEnergy){
	if(verbose > 0)
	  printf("simAnodeDistrib: Error: Minimum scoring energy greater "
		 "than initial beam energy.\n");
	return -1;
      }

      //Check angle
      if(angle < 0.0 || angle >= 90.0){
	if(verbose > 0)
	  printf("simAnodeDistrib: Error: Anode angle must be "
		 "in the range [0,90) Deg.\n");
	return -1;
      }

      const double angleRad = deg2rad*angle;
      
      //Check number of histories
      if(nHists == 0){
	if(verbose > 0)
	  printf("simAnodeDistrib: Error: No histories must be simulated.\n");
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
	  printf("simAnodeDistrib: Warning: Number of threads has been specified,"
		 " but the code has been compiled with no multithreading support.\n"
		 " Only one thread will be used.\n");
      }
      nThreads = 1;
#endif

      //Set dreg according to anode geometry
      dReg = 0.6;

      //Init spatial distribution

      // X
      double width = tan(std::min(colAngle,89.9)*deg2rad)*dReg;
      unsigned long nPixelsX = static_cast<unsigned long>(2.0*(width+pixelSize)/pixelSize);
      if(nPixelsX % 2 != 0)
	nPixelsX += 1;
      width = static_cast<double>(nPixelsX/2) * pixelSize;

      // Y
      double heightLow = tan(angleRad)*dReg;
      unsigned long nPixelsY = static_cast<unsigned long>(2.0*heightLow/pixelSize);
      heightLow = static_cast<double>(nPixelsY)*pixelSize/2.0;
      nPixelsY += nPixelsX/2;
      
      spatialDistrib.initFromLists({nPixelsX, nPixelsY},
				   {std::pair<double,double>(-width, width),
				    std::pair<double,double>(-2.0*heightLow,
							     width)});

      if(verbose > 1){
	printf("          x-ray kvp    : %.2f\n"
	       "Minimum tallied energy : %.2f keV\n"
	       "  Anode angle (Deg)    : %.4f\n"
	       "   Histories to sim    : %llu\n",
	       eEnergy/1000.0, eMin/1000.0, angle, nHists);
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
			      verbose > 2 ? verbose : 1) != pen_context::SUCCESS){
	printf("simAnodeDistrib: Error at simulation context initialization. "
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
      config.set("transforms/anode/enlargeY/transforms/enlargeY/ds", -1.0);

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
      config.set("transforms/anode/setAngleTop/transforms/forward/ds", ds);

      config.set("transforms/anode/setAngleBot/index", 2);
      config.set("transforms/anode/setAngleBot/vertex-group", "front_down");
      config.set("transforms/anode/setAngleBot/transforms/forward/type" , "TRANSLATION_Y");
      config.set("transforms/anode/setAngleBot/transforms/forward/index", 0);
      config.set("transforms/anode/setAngleBot/transforms/forward/ds", -ds);

      //Set detector object as perfect absorber
      config.set("eabs/detector/electron", 1.0e35);      
      config.set("eabs/detector/positron", 1.0e35);      
      config.set("eabs/detector/gamma"   , 1.0e35);

      //Assign kdet 1 to detector object
      config.set("kdet/detector", 1);

      //Set dsmax for anode object
      config.set("dsmax/anode", 2.0e-2);

      //Configure geometry
      if(geometry->configure(config,
			     verbose > 2 ? verbose : 1) != PEN_MESHBODY_GEO_SUCCESS){
	  printf("Unexpected Error: Unable to construct the geometry. "
		 "Please, report this error.\n");
	delete geometry;
	return -4;
      }

      //Set the geometry to the simulation context
      int err = contextSim.setGeometry(geometry);
      if(err != 0){
	printf("Unexpected Error: Unable to set geometry in the "
	       "simulation context. Please, report this error.\n");
	delete geometry;
	return -4;
      }

      //Run context configuration step with geometry
      if(contextSim.configureWithGeo(contextConf,
				     verbose > 2 ? verbose : 1) != pen_context::SUCCESS){
	printf("simAnodeDistrib: Error at simulation context "
	       "initialization with geometry. "
	       "Report this error.\n");
	delete geometry;
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
      photonVR.configure(configVRsplitting,*geometry,verbose > 2 ? verbose : 1);
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

      //Create tallies for each thread
      std::vector<measurements::measurement<double,1>> spectrums(nThreads);
      std::vector<measurements::measurement<double,2>> spatialDistribs(nThreads);

      //Init tallies
      spectrum.setDimHeader(0, "Energy (eV)");
      spectrum.setValueHeader("Value (prob)");

      spatialDistrib.setDimHeader(0, "X (cm)");
      spatialDistrib.setDimHeader(1, "Y (cm)");
      spatialDistrib.setValueHeader("Value (prob)");
      
      for(size_t i = 0; i < nThreads; ++i){
	err = spectrums[i].init(spectrum);
	if(err != 0){
	  printf("Unexpected error: Unable to copy "
		 "spectrum tally configuration. Please, report this.");
	  delete geometry;
	  return -6;	  
	}
	err = spatialDistribs[i].init(spatialDistrib);
	if(err != 0){
	  printf("Unexpected error: Unable to copy "
		 "spatial distribution tally configuration. "
		 "Please, report this.");
	  delete geometry;
	  return -6;	  
	}	
      }

      
      // ** Run simulations
      //*********************

#ifdef _PEN_USE_THREADS_
      std::vector<std::thread> threads;
      for(size_t ith = 0; ith < nThreads; ++ith){
	threads.emplace_back(runAnodeDistribSimulation,
			     nHists/nThreads,
			     eEnergy,
			     std::ref(contextSim),
			     std::ref(photonVR),
			     std::ref(spectrums[ith]),
			     std::ref(spatialDistribs[ith]),
			     std::ref(seeds1[ith]), std::ref(seeds2[ith]),
			     verbose);
      }

      //Wait until threads finish and join results
      for(size_t ith = 0; ith < nThreads; ++ith){
	threads[ith].join();
	spectrum.add(spectrums[ith]);
	spatialDistrib.add(spatialDistribs[ith]);
      }
#else
      //Simulate using a single thread
      runAnodeDistribSimulation(nHists, eEnergy, context, photonVR,
				spectrums[0], spatialDistribs[0],
				seeds1[0], seeds2[0],
				verbose);

      spectrum.add(spectrums[0]);
      spatialDistrib.add(spatialDistribs[0]);
#endif
      delete geometry;
      return 0;
    }

    void createAnode(std::ostream& out,
		     const double angle,
		     const unsigned matIndex,
		     double dx, double dy, double dz,
		     const std::string& name,
		     const std::string& parentName,
		     const bool numObjects,
		     const vector3D<double> center){
      
      //Creates a anode mesh with dimensions "dx" x "dy" x "dz" oriented
      // to the "Y" positive axis. The "dy" dimension is measured from the
      // center of the bevel to the back wall of the anode.
      //
      // Input:
      //   + angle: Anode angle in Deg
      //   + imat : Anode material index
      //   + name: Body name for the anode
      //   + parentName: Parent body name
      //   + dx,dy,dz: Anode dimensions in the X, Y and Z axis
      //   + center: Position of the center of the bevel wall
      //   + numObjects: If true, removes the geometry header
      //                 with the total number of objects
      //
      // Output:
      //   + o: Output stream where the geometry is written
      //   + transforms: Configuration with the required transforms to conform the anode
      
      //
      //
      //  Vertex are constructed as follows 
      //
      // 2*---------------------*0         n Z
      //  |_________dy________ /___center  |
      //  |                   /            |
      // 3*------------------*1            |------> Y
      //
      //                             
      // 4     0   n Y
      // *-----*   | 
      // |     |   |
      // |     |   |-----> X
      // |     |   
      // |     |   
      // *-----*   
      // 6 dx  2
      //
      //

      //Check dimensions
      if(dx <= 0.0)
	dx = 1.0;
      if(dy <= 0.0)
	dy = 1.0;
      if(dz <= 0.0)
	dz = 1.0;

      //Number of vertex
      constexpr unsigned nVertex = 8;

      //Number of triangles
      constexpr unsigned nFaces = 12;
      
      //Check if the number of objects must be printed
      if(numObjects){
	out << "# Number of objects:\n 1" << std::endl;
      }

      // ** Print object name and header
      out << "# Object: " << name << std::endl;
      out << "#MAT      #NFACES     #NVERTEX     #NAME"
	"        #PARENT NAME    #N VERTEX GROUPS"
	  << std::endl;;
      
      //Material index
      out << " " << std::to_string(matIndex);

      //Number of faces
      out << "   " << std::to_string(nFaces);
      
      //Number of vertex
      out << "   " << std::to_string(nVertex);

      //Filter name
      out << "   " << name;

      //Parent name
      out << "   " << parentName;

      //Vertex groups
      out << "   " << 5 << std::endl;

      // ** Print vertex groups
      out << "# VERTEX GROUPS\n"
      "#NAME  #NVERTEX\n"
	" back   0004\n"
	" 0002\n"
	" 0003\n"
	" 0006\n"
	" 0007\n"
	"#NAME  #NVERTEX\n"
	" front_down   0002\n"
	" 0001\n"
	" 0005\n"
	"#NAME  #NVERTEX\n"
	" front_up   0002\n"
	" 0000\n"
	" 0004\n"
	"#NAME  #NVERTEX\n"
	" top   0004\n"
	" 0000\n"
	" 0002\n"
	" 0004\n"
	" 0006\n"
	"#NAME  #NVERTEX\n"
	" bot   0004\n"
	" 0001\n"
	" 0003\n"
	" 0005\n"
	" 0007\n"
	"# VERTEX LIST\n"
	"# Index  (X Y Z)" << std::endl;

      // ** Print vertex positions
      out << "# VERTEX LIST" << std::endl;
      out << "# Index (X Y Z)" << std::endl;

      //Calculate anode transformation
      // 
      //   n  ------      --------
      //   |  |    |      |     / 
      // dz|  |    |  ->  |    /  
      //   u  ------      -----<-->
      //                        ds
      // tan = ds/dz -> ds = tan*dz      
      
      constexpr double pi = 3.141592653589793;
      constexpr double deg2rad = pi/180.0;
      const double angleRad = deg2rad*angle;
      const double ds05 = std::tan(angleRad)*dz/2.0;

      //Create vertex
      
      // +X
      out << " 0 " << center.x + dx/2.0 << " " << center.y + ds05 << " " << center.z + dz/2.0 << std::endl;
      out << " 1 " << center.x + dx/2.0 << " " << center.y - ds05 << " " << center.z - dz/2.0 << std::endl;      
      out << " 2 " << center.x + dx/2.0 << " " << center.y - dy   << " " << center.z + dz/2.0 << std::endl;      
      out << " 3 " << center.x + dx/2.0 << " " << center.y - dy   << " " << center.z - dz/2.0 << std::endl;      

      //-X
      out << " 4 " << center.x - dx/2.0 << " " << center.y + ds05 << " " << center.z + dz/2.0 << std::endl;
      out << " 5 " << center.x - dx/2.0 << " " << center.y - ds05 << " " << center.z - dz/2.0 << std::endl;      
      out << " 6 " << center.x - dx/2.0 << " " << center.y - dy   << " " << center.z + dz/2.0 << std::endl;      
      out << " 7 " << center.x - dx/2.0 << " " << center.y - dy   << " " << center.z - dz/2.0 << std::endl;      

      // Create faces
      out << "# FACES(triangles)" << std::endl;
      out << " 000 004 006" << std::endl;
      out << " 000 006 002" << std::endl;
      out << " 003 002 006" << std::endl;
      out << " 003 006 007" << std::endl;
      out << " 007 006 004" << std::endl;
      out << " 007 004 005" << std::endl;
      out << " 005 001 003" << std::endl;
      out << " 005 003 007" << std::endl;
      out << " 001 000 002" << std::endl;
      out << " 001 002 003" << std::endl;
      out << " 005 004 000" << std::endl;
      out << " 005 000 001" << std::endl;
      out << "#" << std::endl;
    }

  } // namespace xray
} // namespace penred
