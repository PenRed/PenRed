 
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
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

#include "pen_muen.hh"


int pen_muen::calculate(const double Emin,
			const double Emax,
			const unsigned nBins,
			const double tolerance,
			const double simTime,
			const char* matFilename,
			std::vector<double>& EData,
			std::vector<double>& muenData){

  // Calculates the muen values for the specified energy
  // interval in a single material
  //
  // Input:
  //
  // Emin  -> Minimum energy to be calculated
  // Emax  -> Maximum energy to be calculated
  // nBins -> Number of bins to split the interval [Emin,Emax]
  // tolerance -> Threshold for the muen calculation
  // simTime -> Time limit for the calculus of a single point
  // matFilename -> Material file path
  //
  // Output:
  //
  // EData -> Energy values for each calculated bin
  // muenData -> Muen values for each calculated bin

  if(matFilename == nullptr){
    printf("pen_muen:calculate:Error: Material filename not provided\n");
    return -1;
  }

  //Check configuration parameters
  if(tolerance <= 1.0e-6){
    printf("pen_muen:calculate:Error: Tolerance must be positive and"
	   " greater than zero: %E\n",tolerance);
    return -2;
  }
  
  //Create a dummy geometry
  pen_dummyGeo geometry;

  //Configure dummy geometry
  pen_parserSection dummySection;
  geometry.configure(dummySection,2);
  
  //Create elements data base
  pen_elementDataBase* elementsDB = new pen_elementDataBase;
  
  //Create a context
  pen_context context(*elementsDB);

  //Set context geometry
  context.setGeometry(&geometry);
  
  //Set the number of materials to context (1 per thread)
  unsigned int nCalcThreads = 1;
#ifdef _PEN_USE_THREADS_
  
    nCalcThreads = std::max(static_cast<unsigned int>(2),
			    std::thread::hardware_concurrency());
#endif
  
  int errmat = context.setMats<pen_material>(nCalcThreads);
  if(errmat != 0){
    printf("pen_muen:calculate:Error: Unable to create material context: %d.\n",errmat);
    return -3;
  }

  //Init materials parameters
  std::string PMFILEstr[constants::MAXMAT];
  for(unsigned int ith = 0; ith < nCalcThreads; ++ith){

    //Get the material
    pen_material& mat = context.getBaseMaterial(ith);
  
    //Configure the material
    mat.C1=0.2;
    mat.C2=0.2;
    mat.WCC=1.0e3;
    mat.WCR=1.0e3;

    mat.EABS[PEN_ELECTRON] = 50.0E0;
    mat.EABS[PEN_PHOTON]   = 50.0E0;
    mat.EABS[PEN_POSITRON] = 50.0E0;
    PMFILEstr[ith].assign(matFilename);
  }
  
  FILE* fcontext = nullptr;
  fcontext = fopen("muen-context.rep","w");
  if(fcontext == nullptr){
    printf("pen_muen:calculate:Error: Unable to create file 'muen-context.rep'\n");
    return -4;
  }

  double EMAX=1.0E9;
  int INFO = 0;
  int err = context.init(EMAX,fcontext,INFO,PMFILEstr);
  if(err != 0){
    printf("pen_muen:calculate:Error: Unable to configure context."
	   " Check 'muen-context.rep'.\n");
    return -5;
  }

  fclose(fcontext);

  double dE = (Emax-Emin)/static_cast<double>(nBins);

  EData.resize(nBins);
  muenData.resize(nBins);

  
#ifdef _PEN_USE_THREADS_
  
  std::vector<std::thread> calcThreads;

  std::atomic<unsigned int> atomicCount{0};
    
  for(size_t ith = 0; ith < nCalcThreads; ++ith){
    calcThreads.push_back(std::thread([&,ith](){

      //Get initial random seeds
      int seed1, seed2;
      rand0(ith,seed1,seed2);	

      //Get local material
      pen_material& mat = context.getBaseMaterial(ith);
      double EABS1 = mat.EABS[PEN_ELECTRON];
      
      unsigned int ibin = atomicCount++;
      while(ibin < nBins){

	//Get energy
	double E0 = Emin + dE*ibin;

	if(E0 < 50.0){
	  printf("pen_muen:calculate:Warning: Energy %E is too low. "
		 "Skipping.\n",E0);
	  continue;
	}

	//Ensure to include electron significant radiative yield
	for(long int i = static_cast<long int>(constants::NEGP)-1; i >= 0; --i){
	  if(exp(mat.EBRY[i])*context.grid.ET[i] < 1.0e-4*E0){
	    mat.EABS[PEN_ELECTRON] = std::max(context.grid.ET[i],EABS1);
	    break;
	  }
	}

	mat.EABS[PEN_POSITRON] =
	  std::min(std::max(1.0e-4*E0,5.0e3),mat.EABS[PEN_ELECTRON]);	

	
	double muenVal =
	  pen_muen::simulate(context,E0,simTime,tolerance,seed1,seed2,ith);

	EData[ibin] = E0;
	muenData[ibin] = muenVal;
	ibin = atomicCount++;
      }
	
    }));
  }

  for(size_t ith = 0; ith < nCalcThreads; ++ith){
    calcThreads[ith].join();
  }

#else

  //Get material
  pen_material& mat = context.getBaseMaterial(0);
  double EABS1 = mat.EABS[PEN_ELECTRON];
  
  //Get initial random seeds
  int seed1, seed2;
  rand0(5,seed1,seed2);
    
  for(unsigned ibin = 0; ibin < nBins; ++ibin){

    //Get energy
    double E0 = Emin + dE*ibin;

    if(E0 < 50.0){
      printf("pen_muen:calculate:Warning: Energy %E is too low. "
	     "Skipping.\n",E0);
      continue;
    }

    //Ensure to include electron significant radiative yield
    for(long int i = static_cast<long int>(constants::NEGP)-1; i >= 0; --i){
      if(exp(mat.EBRY[i])*context.grid.ET[i] < 1.0e-4*E0){
	mat.EABS[PEN_ELECTRON] = std::max(context.grid.ET[i],EABS1);
	break;
      }
    }

    mat.EABS[PEN_POSITRON] =
      std::min(std::max(1.0e-4*E0,5.0e3),mat.EABS[PEN_ELECTRON]);

    double muenVal =
      pen_muen::simulate(context,E0,simTime,tolerance,seed1,seed2,0);
	
    EData[ibin] = E0;
    muenData[ibin] = muenVal;
  }

    
#endif  
  
  return 0;
}

int pen_muen::calculate(const char** energySpectrums,
			const unsigned nSpectrums,
			const double tolerance,
			const double simTime,
			const char* matFilename,
			std::vector<pen_muen::muData>& muenData,
			unsigned verbose){

  // Calculates the muen values for the specified energy
  // spectrums in a single material
  //
  // Input:
  //
  // energySpectrums -> Files paths of the energy spectrums
  // nSpectrums -> Number of spectrums
  // tolerance -> Threshold for the muen calculation
  // simTime -> Time limit for the calculus of a single point
  // matFilename -> Material file path
  //
  // Output:
  //
  // muenData -> Calculated values for each spectrum

  if(matFilename == nullptr){
    printf("pen_muen:calculate:Error: Material filename not provided\n");
    return -1;
  }

  //Check configuration parameters
  if(tolerance <= 1.0e-6){
    printf("pen_muen:calculate:Error: Tolerance must be positive and"
	   " greater than zero: %E\n",tolerance);
    return -2;
  }
  
  //Create a dummy geometry
  pen_dummyGeo geometry;

  //Configure dummy geometry
  pen_parserSection dummySection;
  geometry.configure(dummySection,2);
  
  //Create elements data base
  pen_elementDataBase* elementsDB = new pen_elementDataBase;
  
  //Create a context
  pen_context context(*elementsDB);

  //Set context geometry
  context.setGeometry(&geometry);
  
  //Set the number of materials to context (1 per thread)
  unsigned int nCalcThreads = 1;
#ifdef _PEN_USE_THREADS_
  
    nCalcThreads = std::max(static_cast<unsigned int>(2),
			    std::thread::hardware_concurrency());
#endif
  
  int errmat = context.setMats<pen_material>(1);
  if(errmat != 0){
    printf("pen_muen:calculate:Error: Unable to create material context: %d.\n",errmat);
    return -3;
  }

  //Init materials parameters
  std::string PMFILEstr[constants::MAXMAT];

  //Get the material
  pen_material& mat = context.getBaseMaterial(0);
  
  //Configure the material
  mat.C1=0.2;
  mat.C2=0.2;
  mat.WCC=1.0e3;
  mat.WCR=1.0e3;

  mat.EABS[PEN_ELECTRON] = 50.0E0;
  mat.EABS[PEN_PHOTON]   = 50.0E0;
  mat.EABS[PEN_POSITRON] = 50.0E0;
  PMFILEstr[0].assign(matFilename);
  double EABS1 = mat.EABS[PEN_ELECTRON];
  
  
  FILE* fcontext = nullptr;
  fcontext = fopen("muen-context.rep","w");
  if(fcontext == nullptr){
    printf("pen_muen:calculate:Error: Unable to create file 'muen-context.rep'\n");
    return -4;
  }

  double EMAX=1.0E9;
  int INFO = 0;
  int err = context.init(EMAX,fcontext,INFO,PMFILEstr);
  if(err != 0){
    printf("pen_muen:calculate:Error: Unable to configure context."
	   " Check 'muen-context.rep'.\n");
    return -5;
  }

  fclose(fcontext);

  muenData.resize(nSpectrums);

  if(verbose > 1){
    printf("\n%s\n", muData::header(30).c_str());
  }
    
  for(size_t is = 0; is < nSpectrums; ++is){

    //Create and configure energy sampler
    fileSpectrum_energySampling spectrumSampler;
    double Emax;
    pen_parserSection config;
    config.set("filename", std::string(energySpectrums[is]));
    err = spectrumSampler.configure(Emax, config, 1);
    if(err != 0){
      printf("pen_muen:calculate:Error: Unable to calculate "
	     "energy sampler ('%s')\n", energySpectrums[is]);
      return -1;
    }

    //Get minimum energy
    double Emin = spectrumSampler.minE();

    if(Emin < 50.0){
      printf("pen_muen:calculate:Warning: Minimum spectrum energy "
	     "%E is too low.\n",Emin);
    }

    //Ensure to include electron significant radiative yield
    for(long int i = static_cast<long int>(constants::NEGP)-1; i >= 0; --i){
      if(exp(mat.EBRY[i])*context.grid.ET[i] < 1.0e-4*Emin){
	mat.EABS[PEN_ELECTRON] = std::max(context.grid.ET[i],EABS1);
	break;
      }
    }

    mat.EABS[PEN_POSITRON] =
      std::min(std::max(1.0e-4*Emin,5.0e3),mat.EABS[PEN_ELECTRON]);


    pen_muen::muenSimTally globalTally;
    globalTally.EDPT = 0.0;
    globalTally.EDPT2 = 0.0;  

    globalTally.ETRT = 0.0;
    globalTally.ETRT2 = 0.0;

    globalTally.RANGET = 0.0;
    globalTally.RANGET2 = 0.0;

    globalTally.ET = 0.0;
    globalTally.ET2 = 0.0;
    
    globalTally.nhist = 0;
    
#ifdef _PEN_USE_THREADS_
  
    std::vector<std::thread> calcThreads;
    std::vector<pen_muen::muenSimTally> tallies;
    tallies.resize(nCalcThreads);

    double localTolerance = tolerance*sqrt(static_cast<double>(nCalcThreads));
    unsigned long long minHists = 100000/static_cast<double>(nCalcThreads);
  
    for(size_t ith = 0; ith < nCalcThreads; ++ith){
      calcThreads.push_back(std::thread([&,ith](){

	//Get initial random seeds
	int seed1, seed2;
	rand0(ith,seed1,seed2);	

	//Simulate muen
	pen_muen::simulate(context,
			   spectrumSampler,
			   simTime,
			   localTolerance,
			   seed1,seed2,
			   tallies[ith],
			   minHists);
      
      }));
    }
  
    for(size_t ith = 0; ith < nCalcThreads; ++ith){
      calcThreads[ith].join();

      globalTally.EDPT  += tallies[ith].EDPT;
      globalTally.EDPT2 += tallies[ith].EDPT2;  

      globalTally.ETRT  += tallies[ith].ETRT;
      globalTally.ETRT2 += tallies[ith].ETRT2;

      globalTally.RANGET  += tallies[ith].RANGET;
      globalTally.RANGET2 += tallies[ith].RANGET2;

      globalTally.ET  += tallies[ith].ET;
      globalTally.ET2 += tallies[ith].ET2;
      
      globalTally.nhist += tallies[ith].nhist;

    }
  
#else

    //Get initial random seeds
    int seed1, seed2;
    rand0(ith,seed1,seed2);	

    //Simulate muen
    pen_muen::simulate(context,
		       spectrumSampler,
		       simTime,
		       tolerance,
		       seed1,seed2,
		       globalTally);

    
#endif

    //Calculate muen    
    double FNT  = 1.0/static_cast<double>(globalTally.nhist);
    double FG   = globalTally.EDPT*FNT;
    double EFG  = FNT*sqrt(fabs(globalTally.EDPT2-
    globalTally.EDPT*globalTally.EDPT*FNT)); // 1 sigma
    double PERR = 100.0*EFG/FG;

    double RMU0 = 1.0/(globalTally.RANGET*FNT);
    RMU0 = RMU0/mat.RHO;

    //Mean energy of the sampled particles
    double E0 = globalTally.ET*FNT;

    pen_muen::muData& mu = muenData[is];

    mu.muRho = RMU0;
    mu.muenRho = FG*RMU0/E0;
    mu.E0 = E0;
    mu.f = globalTally.ETRT*FNT;
    mu.fg = FG;
    mu.err = PERR;

    if(verbose > 1){
      printf("%-30s %s\n", energySpectrums[is], mu.stringify().c_str());
    }
  }
  
  return 0;
}

double pen_muen::simulate(const pen_context& context,
			  const double E0,
			  const double simTime,
			  const double tolerance,
			  int& seed1, int& seed2,
			  const unsigned ithread,
			  const unsigned long long minHists){

  //Calculates the muen value for the specified energy
  //
  // context -> Simulation context
  // E0 -> Energy to be calculated
  // simTime -> Maximum simulation time
  // tolerance -> Theshold to finish the simulation
  // seed1, seed2 -> Seed pair for this thread
  // ithread -> Thread identifier

  //Declare random generator
  pen_rand random;
  //Initialize random
  random.setSeeds(seed1,seed2);

  //Get material
  const pen_material& mat = context.readBaseMaterial(ithread);
  
  unsigned long long ntot = 1000000000;

  //Create particle stacks for electrons, positrons and gamma
  pen_particleStack<pen_particleState> stackE;
  pen_particleStack<pen_particleState> stackP;
  pen_particleStack<pen_state_gPol> stackG;
  
  //Create particle simulations
  pen_betaE betaE(context,stackE,stackG);
  pen_gamma gamma(context,stackE,stackP,stackG);
  pen_betaP betaP(context,stackE,stackG,stackP);
    
  double RMU0 = 1.0/context.range(E0,PEN_PHOTON,ithread);
  RMU0 = RMU0/mat.RHO;
    
  //Create a common initial history state
  pen_state_gPol genState;
  genState.E = E0;
  genState.WGHT = 1.0;
  genState.MAT = ithread+1; //Each thread uses a unique copy of the material
  genState.IBODY = 0;
  genState.ILB[0] = 1;
  genState.ILB[1] = 0;
  genState.ILB[2] = 0;
  genState.ILB[3] = 0;
  genState.ILB[4] = 0;
  genState.U = 0.0;
  genState.V = 0.0;
  genState.W = 1.0;

  double DSMAX = 1.0e20;

  double EDPT = 0.0;
  double ETRT = 0.0;
  double EDPT2 = 0.0;
  double ETRT2 = 0.0;
    
  //Create timers
  long long int mili2finish = static_cast<long long int>(simTime)*1000;
  pen_stopWatch finishWatch(mili2finish);
  finishWatch.start();

  double muenVal = 0.0;
  for(unsigned long long nhist = 1; nhist <= ntot; ++nhist){

    //Clear stacks
    stackE.cleans();
    stackP.cleans();
    stackG.cleans();

    //Copy initial state
    stateCopy(gamma.getState(),genState);
    const pen_state_gPol& state = gamma.readState();

    //Update material and body information
    gamma.updateMat();
    gamma.updateBody();
    
    gamma.START();

    //Interact until get a non delta interaction
    int icol;
    do{
      double ds, de;
      gamma.JUMP(ds,random,DSMAX);
      gamma.KNOCK(de,icol,random);
    }while(icol == GAMMA_DELTA);

    double EDP = E0 - state.E;
    double ETR = EDP;

    pen_state_gPol auxGammaState;
    bool continueSim = true;
    while(continueSim){

      continueSim = false;
	
      //Count all photons
      while(stackG.getNSec() > 0){
	continueSim = true;
	stackG.get(auxGammaState);
	EDP -= auxGammaState.E;
	if(auxGammaState.ILB[0] == 2) ETR -= auxGammaState.E;
      }

      //Count electrons
      while(stackE.getNSec() > 0){
	continueSim = true;
	stackE.get(betaE.getState());
	if(betaE.readState().E >= mat.EABS[PEN_ELECTRON]){
	  //Update material and body information
	  betaE.updateMat();
	  betaE.updateBody();
	  
	  betaE.START();
	  do{
	    double ds, de, X, Y, Z;
	    int icolDummy;
	    betaE.JUMP(ds,random,DSMAX);  //Call Jump
	    betaE.move(ds,de,X,Y,Z,random); //Move tha particle due soft interactions
	    if(betaE.readState().E < mat.EABS[PEN_ELECTRON])
	      break;
	    betaE.KNOCK(de,icolDummy,random);  //Call nock
	  }while(betaE.readState().E >= mat.EABS[PEN_ELECTRON]);
	}
      }

      //Count positrons
      while(stackP.getNSec() > 0){
	continueSim = true;
	stackP.get(betaP.getState());
	if(betaP.readState().E >= mat.EABS[PEN_POSITRON]){
	  //Update material and body information
	  betaP.updateMat();
	  betaP.updateBody();

	  betaP.START();
	  do{
	    double ds, de, X, Y, Z;
	    int icolDummy;
	    betaP.JUMP(ds,random,DSMAX);  //Call Jump
	    betaP.move(ds,de,X,Y,Z,random); //Move tha particle due soft interactions
	    if(betaP.readState().E < mat.EABS[PEN_POSITRON])
	      break;
	    betaP.KNOCK(de,icolDummy,random);  //Call nock
	  }while(betaP.readState().E >= mat.EABS[PEN_POSITRON]);
	}
      }
    }

    //Add history results to global counters
    EDPT += EDP;  EDPT2 += EDP*EDP;
    ETRT += ETR;  ETRT2 += ETR*ETR;

    //Ensure a minimum number of histories
    if(nhist < minHists)
      continue;
    
    //Calculate muen
    double FNT  = 1.0/static_cast<double>(nhist);
    double FG   = EDPT*FNT;
    double EFG  = FNT*sqrt(fabs(EDPT2-EDPT*EDPT*FNT)); // 1 sigma
    double PERR = 100.0*EFG/FG;
    muenVal = FG*RMU0/E0;

    //Check if the objective uncertainty has been reached
    if(PERR < tolerance){
      //Stop simulation
      break;
    }
    
    //Check if is time to finish
    const auto tnow = std::chrono::steady_clock::now();
    if(finishWatch.check(tnow)){
      //Time to finish
      break;
    }
    
  }

  random.getSeeds(seed1,seed2);
  return muenVal;
}

double pen_muen::simulate(const pen_context& context,
			  const fileSpectrum_energySampling& spectrum,
			  const double simTime,
			  const double tolerance,
			  int& seed1, int& seed2,
			  pen_muen::muenSimTally &tally,
			  const unsigned long long minHists){

  //Calculates the muen value for the specified energy
  //
  // context -> Simulation context
  // spectrum -> Energy sampler with the spectrum to be simulated
  // simTime -> Maximum simulation time
  // tolerance -> Theshold to finish the simulation
  // seed1, seed2 -> Seed pair for this thread
  // tally -> Structure to store tallied results
  
  //Declare random generator
  pen_rand random;
  //Initialize random
  random.setSeeds(seed1,seed2);

  //Get material
  const pen_material& mat = context.readBaseMaterial(0);
  
  unsigned long long ntot = 1000000000;

  //Create particle stacks for electrons, positrons and gamma
  pen_particleStack<pen_particleState> stackE;
  pen_particleStack<pen_particleState> stackP;
  pen_particleStack<pen_state_gPol> stackG;
  
  //Create particle simulations
  pen_betaE betaE(context,stackE,stackG);
  pen_gamma gamma(context,stackE,stackP,stackG);
  pen_betaP betaP(context,stackE,stackG,stackP);

  //double RMU0 = 1.0/context.range(E0,PEN_PHOTON,ithread);
  //RMU0 = RMU0/mat.RHO;
    
  //Create a common initial history state
  pen_state_gPol genState;
  genState.WGHT = 1.0;
  genState.MAT = 1; //All threads use the same copy of the material
  genState.IBODY = 0;
  genState.ILB[0] = 1;
  genState.ILB[1] = 0;
  genState.ILB[2] = 0;
  genState.ILB[3] = 0;
  genState.ILB[4] = 0;
  genState.U = 0.0;
  genState.V = 0.0;
  genState.W = 1.0;

  double DSMAX = 1.0e20;

  double EDPT = 0.0;
  double ETRT = 0.0;
  double EDPT2 = 0.0;
  double ETRT2 = 0.0;

  double RANGET = 0.0;
  double RANGET2 = 0.0;

  double ET = 0.0;
  double ET2 = 0.0;
  
  //Create timers
  long long int mili2finish = static_cast<long long int>(simTime)*1000;
  pen_stopWatch finishWatch(mili2finish);
  finishWatch.start();

  for(unsigned long long nhist = 1; nhist <= ntot; ++nhist){

    //Clear stacks
    stackE.cleans();
    stackP.cleans();
    stackG.cleans();

    //Sample energy
    spectrum.energySampling(genState.E,random);
    const double E0 = genState.E;
    ET  += E0;
    ET2 += E0*E0;

    //Calculate the range for this energy
    double RANGE = context.range(E0,PEN_PHOTON,0);
    RANGET  += RANGE;
    RANGET2 += RANGE*RANGE;
    
    //Copy initial state
    stateCopy(gamma.getState(),genState);
    const pen_state_gPol& state = gamma.readState();

    //Update material and body information
    gamma.updateMat();
    gamma.updateBody();
    
    gamma.START();

    //Interact until get a non delta interaction
    int icol;
    do{
      double ds, de;
      gamma.JUMP(ds,random,DSMAX);
      gamma.KNOCK(de,icol,random);
    }while(icol == GAMMA_DELTA);

    double EDP = E0 - state.E;
    double ETR = EDP;

    pen_state_gPol auxGammaState;
    bool continueSim = true;
    while(continueSim){

      continueSim = false;
	
      //Count all photons
      while(stackG.getNSec() > 0){
	continueSim = true;
	stackG.get(auxGammaState);
	EDP -= auxGammaState.E;
	if(auxGammaState.ILB[0] == 2) ETR -= auxGammaState.E;
      }

      //Count electrons
      while(stackE.getNSec() > 0){
	continueSim = true;
	stackE.get(betaE.getState());
	if(betaE.readState().E >= mat.EABS[PEN_ELECTRON]){
	  //Update material and body information
	  betaE.updateMat();
	  betaE.updateBody();
	  
	  betaE.START();
	  do{
	    double ds, de, X, Y, Z;
	    int icolDummy;
	    betaE.JUMP(ds,random,DSMAX);  //Call Jump
	    betaE.move(ds,de,X,Y,Z,random); //Move tha particle due soft interactions
	    if(betaE.readState().E < mat.EABS[PEN_ELECTRON])
	      break;
	    betaE.KNOCK(de,icolDummy,random);  //Call nock
	  }while(betaE.readState().E >= mat.EABS[PEN_ELECTRON]);
	}
      }

      //Count positrons
      while(stackP.getNSec() > 0){
	continueSim = true;
	stackP.get(betaP.getState());
	if(betaP.readState().E >= mat.EABS[PEN_POSITRON]){
	  //Update material and body information
	  betaP.updateMat();
	  betaP.updateBody();

	  betaP.START();
	  do{
	    double ds, de, X, Y, Z;
	    int icolDummy;
	    betaP.JUMP(ds,random,DSMAX);  //Call Jump
	    betaP.move(ds,de,X,Y,Z,random); //Move tha particle due soft interactions
	    if(betaP.readState().E < mat.EABS[PEN_POSITRON])
	      break;
	    betaP.KNOCK(de,icolDummy,random);  //Call nock
	  }while(betaP.readState().E >= mat.EABS[PEN_POSITRON]);
	}
      }
    }

    //Add history results to global counters
    EDPT += EDP;  EDPT2 += EDP*EDP;
    ETRT += ETR;  ETRT2 += ETR*ETR;

    //Ensure a minimum number of histories
    if(nhist < minHists)
      continue;
    
    //Calculate muen error
    double FNT  = 1.0/static_cast<double>(nhist);
    double FG   = EDPT*FNT;
    double EFG  = FNT*sqrt(fabs(EDPT2-EDPT*EDPT*FNT)); // 1 sigma
    double PERR = 100.0*EFG/FG;

    //Check if the objective uncertainty has been reached
    if(PERR < tolerance){
      //Stop simulation
      tally.nhist = nhist;
      break;
    }
    
    //Check if is time to finish
    const auto tnow = std::chrono::steady_clock::now();
    if(finishWatch.check(tnow)){
      //Time to finish
      tally.nhist = nhist;
      break;
    }
    
  }

  //Store tallied values
  tally.EDPT  = EDPT;
  tally.EDPT2 = EDPT2;

  tally.ETRT  = ETRT;
  tally.ETRT2 = ETRT2;

  tally.RANGET = RANGET;
  tally.RANGET2 = RANGET2;

  tally.ET  = ET;
  tally.ET2 = ET2;

  //Calculate muen
  double FNT  = 1.0/static_cast<double>(tally.nhist);
  double FG   = EDPT*FNT;
  
  double RMU0 = 1.0/(RANGET*FNT);
  RMU0 = RMU0/mat.RHO;

  double Emean = ET*FNT;
  
  double muenVal = FG*RMU0/Emean;
  
  random.getSeeds(seed1,seed2);
  return muenVal;
}

