
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#include <thread>
#include <limits>
#include <cctype>
#include <algorithm>
#include "pen_simulation.hh"

int createParticleGenerators(std::vector<pen_specificStateGen<pen_particleState>>& genericSources,
			     std::vector<pen_specificStateGen<pen_state_gPol>>& polarisedGammaSources,
			     const unsigned nthreads,
			     const pen_parserSection& config,
			     const unsigned verbose);
int createTallies(std::vector<pen_commonTallyCluster>& tallyGroups,
		  const size_t nthreads,
		  const pen_context& context,
		  const pen_parserSection& config,
		  const unsigned verbose);

int createGeometry(wrapper_geometry*& geometry,
		   const pen_parserSection& config,
		   const pen_context& context,
		   const unsigned verbose);

int createMaterials(pen_context& context,
		    std::string filenames[constants::MAXMAT],
		    const pen_parserSection& config,
		    const double globEmax,
		    const unsigned verbose);

int setVarianceReduction(pen_context& context,
			 const wrapper_geometry& geometry,
			 const pen_parserSection& config,
			 pen_VRCluster<pen_particleState>& genericVR,
			 pen_VRCluster<pen_state_gPol>& photonVR,
			 const unsigned verbose);

pen_context* createAuxContext(double EMAX,
			      const char* matFilename,
			      const unsigned verbose);

template<class stateType>
int configureSource(pen_specificStateGen<stateType>& source,
		    const unsigned& nthreads,
		    const pen_parserSection& config,
		    const unsigned verbose){
  
  //Check if specific sampler exists
  bool useSpecific = config.isSection("specific");
  
  //Get source kpar
  std::string auxkpar;
  int err = config.read("kpar",auxkpar);
  if(err != INTDATA_SUCCESS && !useSpecific){
    if(verbose > 0){
      printf("configureSource: Error: Particle type field ('kpar') not found and no specific sampling specified. String expected:\n");
      printf("\telectron\n");
      printf("\tgamma\n");
      printf("\tpositron\n");
    }
    return -1;
  } else if(!useSpecific && verbose > 1) {
    printf("Selected particle: %s\n",auxkpar.c_str());
  }

  //Obtain particle ID
  if(auxkpar.compare("electron") == 0)
    source.kpar = PEN_ELECTRON;
  else if(auxkpar.compare("gamma") == 0)
    source.kpar = PEN_PHOTON;
  else if(auxkpar.compare("positron") == 0)
    source.kpar = PEN_POSITRON;
  else if(!useSpecific){
    if(verbose > 0)
      printf("\nInvalid selected particle type: %s\n",auxkpar.c_str());
    return -2;
  }

  //Configure source
  source.configure(config,nthreads,verbose);
  if(source.configureStatus() != 0){
    if(verbose > 0){
      printf("\nError on source configuration\n");
    }
    return -3;
  }
  
  return 0;
}

template<class stateType>
void simulate(penred::simulation::simConfig& config,
	      const pen_context& context,
	      pen_specificStateGen<stateType>& source,
	      pen_commonTallyCluster& tallies,
	      const pen_VRCluster<pen_particleState>& genericVR,
	      const pen_VRCluster<pen_state_gPol>& photonVR){

  // Variables list
  //
  // Input:
  //  ithread     : thread number identifier
  //  pcontext    : simulation context
  //  psource     : particle source
  //  ptallies    : simulation tallies
  //  dumptime    : time between dumps
  //  dumpFilename: dump file name
  //
  // Input/output:
  //  seed1       : First random generator seed
  //  seed2       : Second random generator seed
  //  nsimulated  : total number of simulated histories (including 'initHist')

  //Create particle stacks for electrons, positrons and gamma
  pen_particleStack<pen_particleState> stackE;
  pen_particleStack<pen_particleState> stackP;
  pen_particleStack<pen_state_gPol> stackG;
  
  //Create particles to simulate
  pen_betaE betaE(context,stackE,stackG);
  pen_gamma gamma(context,stackE,stackP,stackG);
  pen_betaP betaP(context,stackE,stackG,stackP);

  //Register VR
  if(genericVR.numVR() > 0){
    betaE.registerGenericVR(genericVR);
    gamma.registerGenericVR(genericVR);
    betaP.registerGenericVR(genericVR);
  }
  if(photonVR.numVR() > 0){
    gamma.registerSpecificVR(photonVR);
  }

  //Perform the simulation
  penred::simulation::sampleAndSimulate(config, source, tallies,
					gamma, betaE, betaP);
  
}


void printUserInputOptions(const char* filename){
  printf("\nAvailable interactive options:\n"
	 "- get status    : Prints the complete status of all threads\n"
	 "- get speeds    : Prints speed in hist/s for each thread\n"
	 "- get simulated : Prints the simulated histories for each thread\n"
	 "- help          : Displays this text\n"
	 "\n"
	 "To request options, write each desired instruction on an \n"
	 "independent line in a file named '%s'.\n"
	 "\n"
	 "It is recommended to create the instruction file with a different\n"
	 "name and then rename it to '%s' to avoid processing when the file\n"
	 "while still being edited. For example, status can be obtained using: \n"
	 "\n"
	 " echo \"get status\" > userInstructions.txt\n\n", filename, filename);
}

void checkUserInput(const char* filename,
		    const std::string sourceName,
		    const std::vector<penred::simulation::simConfig>& simConfigs){

  unsigned inputsDone = 0;

  //Wait until the simulation of actual source starts for all threads
  for(;;){
    for(const penred::simulation::simConfig& simConf : simConfigs){
      if(simConf.getCurrentSourceName() != sourceName){
	std::this_thread::sleep_for(std::chrono::seconds(2));
	continue;
      }
    }
    //All threads started the source
    break;
  }

  printUserInputOptions(filename);

  //Check for user input periodically
  for(;;){

    //Check if all sources have been finished
    bool finish = true;
    for(const penred::simulation::simConfig& simConf : simConfigs){
      if(!simConf.isSourceFinished()){
	finish = false;
	break;
      }
    }

    if(finish)
      break;

    //Simulation still running, process input
    FILE* fin = nullptr;
    fin = fopen(filename, "r");
    if(fin != nullptr){
      //Input file created, read it
      char line[500];
      unsigned long nLine;
      while(pen_getLine(fin, 500, line, nLine) == 0){
	//Process input
	char word1[500];
	char word2[500];
	sscanf(line, " %s %s ", word1, word2);
	if(strcmp(word1,"get") == 0){
	  if(strcmp(word2, "status") == 0){
	    //Print status of all threads
	    for(const penred::simulation::simConfig& simConf : simConfigs){
	      std::cout << simConf.stringifyState() << std::endl;
	    }
	  }else if(strcmp(word2, "speeds") == 0 ||
		   strcmp(word2, "speed" ) == 0){
	    //Print mean speeds of all threads
	    for(const penred::simulation::simConfig& simConf : simConfigs){
	      double meanSpeed = simConf.getSpeedInSource();
	      std::string out(simConf.threadAndSourcePrefix());
	      out += std::to_string(meanSpeed) + " hist/s\n";
	      std::cout << out;
	    }
	    std::cout << std::endl;
	  }else if(strcmp(word2, "simulated") == 0){
	    //Print simulated histories in each thread
	    for(const penred::simulation::simConfig& simConf : simConfigs){
	      unsigned long long simulated  = simConf.getSimulatedInSource();
	      unsigned long long toSimulate = simConf.getToSimulateInSource();
	      std::string out(simConf.threadAndSourcePrefix());
	      out += "Simulated " + std::to_string(simulated) +
		"/" + std::to_string(toSimulate) + " histories\n";
	      std::cout << out;
	    }
	    std::cout << std::endl;	    
	  }
	}
	else if(strcmp(word1,"help") == 0){
	  printUserInputOptions(filename);
	}
      }

      //Close file and remove it
      fclose(fin);
      std::string newFilename =
	"input" +
	std::to_string(inputsDone++) +
	".txt";
      std::rename(filename, newFilename.c_str());
    }
    std::this_thread::sleep_for(std::chrono::seconds(10));
  }
}

int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s config-filename\n",argv[0]);
    return 0;
  }

  printf("***************************************************************\n");
  printf(" PenRed version: 1.9.4 (19-February-2024) \n");
  printf(" Copyright (c) 2019-2023 Universitat Politecnica de Valencia\n");
  printf(" Copyright (c) 2019-2023 Universitat de Valencia\n");
  printf(" Reference: Computer Physics Communications, 267 (2021) 108065\n"
         "            https://doi.org/10.1016/j.cpc.2021.108065\n");
  printf(" This is free software; see the source for copying conditions.\n"
	 " There is NO warranty; not even for MERCHANTABILITY or\n"
         " FITNESS FOR A PARTICULAR PURPOSE.\n");
  printf(" Please, report bugs and suggestions at our github repository\n"
	 "         https://github.com/PenRed/PenRed\n");
  printf("***************************************************************\n\n");
  
  if(strcmp(argv[1],"--version") == 0 || strcmp(argv[1],"-v") == 0){
    return 0;
  }
  
  bool addDumps = false;
  if(argc >= 3){
    if(strcmp(argv[2],"--addDumps") == 0){ //The argv[1] contains the configuration file
      printf("Adding specified dumps\n");
      if(argc < 5){
	printf("Error: At least two dumps are required\n");
	return -1;
      }
      addDumps = true;
    }
  }
  
  unsigned verbose = 2;

  //Create a timer to measure the expended time in initialization 
  pen_timer initializationTimer;
  
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  //Initialize MPI
  int rank;
  int mpiSize;
  int MThProvided;
  int MPIinitErr = MPI_Init_thread(nullptr, nullptr, MPI_THREAD_SERIALIZED,&MThProvided);
  if(MPIinitErr != MPI_SUCCESS){
    printf("Unable to initialize MPI. Error code: %d\n",MPIinitErr);
    return -1;
  }
  if(MThProvided != MPI_THREAD_SERIALIZED){
    printf("Warning: The MPI implementation used doesn't provide"
	   "support for serialized thread communication.\n"
	   "This could produce unexpected behaviours or performance issues.\n");
  }
  
  //Get current process "rank" and the total number of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

  if(verbose > 1){
    if(rank == 0){
      printf("%d MPI processes started.\n",mpiSize);
    }
  }
  
  //Redirect stdout to individual log files
  char stdoutFilename[100];
  snprintf(stdoutFilename,100,"rank-%03d.log",rank);
  printf("Rank %d: Redirect 'stdout' to file '%s'\n",rank,stdoutFilename);
  if(freopen(stdoutFilename, "w", stdout) == NULL){
    printf("Rank %d: Can't redirect stdout to '%s'\n",rank,stdoutFilename);
  }
  else if(verbose > 1){
    printf("\n**** Rank %d Log ****\n\n",rank);
  }
    
#endif
  // ***************************** MPI END ********************************** //
  
  //**************************
  // Parse configuration
  //**************************

  //Parse configuration file
  pen_parserSection config;
  std::string errorLine;
  unsigned long errorLineNum;
  int err = parseFile(argv[1],config,errorLine,errorLineNum);

  //printf("Configuration:\n");
  //printf("%s\n", config.stringify().c_str());
  
  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n",err);
    printf("Error message: %s\n",pen_parserError(err));
    printf("Error located at line %lu, at text: %s\n",
	   errorLineNum,errorLine.c_str());
    return -1;
  }

  // Get comon simulation configuration
  //*************************************
  penred::simulation::simConfig baseSimConfig;
  int errSimConfig = baseSimConfig.configure("simulation",config);
  if(errSimConfig != penred::simulation::errors::SUCCESS){
    printf("Error: Unable to parse 'simulation' section: Invalid seeds\n");
    return -1;
  }
  //Extract verbose level
  verbose = baseSimConfig.verbose;

  if(verbose > 1)
    printf("%s\n", baseSimConfig.stringifyConfig().c_str());

  // Get user interactive options
  //*************************************
  std::string userInteractiveFilename;
  bool userInteractiveEnabled;
  if(config.read("simulation/interactive",userInteractiveEnabled) != INTDATA_SUCCESS){
    userInteractiveEnabled = true;
  }
  if(config.read("simulation/interactive-file",
		 userInteractiveFilename) != INTDATA_SUCCESS){
    userInteractiveFilename = "userInstructions.txt";
  }

  if(verbose > 1){
    printf("User interactive %s ", userInteractiveEnabled ? "enabled" : "disabled");
    if(userInteractiveEnabled){
      printf("via filename '%s'\n", userInteractiveFilename.c_str());
    }else{
      printf("\n");
    }
  }

  // Get export image formats
  //***************************
  pen_parserSection formatsSection;
  std::vector<std::string> imageFormats;
  std::array<bool,pen_imageExporter::nFormats()> enabledFormats = {false};
  if(config.read("simulation/images",formatsSection) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No image format specified to export . Specify a format via the "
	     "'simulation/images' section to enable image export.\n");
    }
  }else{

    //Get keys in the format section
    formatsSection.ls(imageFormats);

    for(std::string& format : imageFormats){

      //Read if this format is enabled
      bool enabled;
      if(formatsSection.read(format.c_str(),enabled) != INTDATA_SUCCESS){
	enabled = false;
      }

      if(enabled){
	std::transform(format.begin(), format.end(), format.begin(),
		       [](unsigned char c){ return std::toupper(c); });

	//Check if this format exists
	if(!pen_imageExporter::isFormat(format.c_str())){
	  if(verbose > 0){
	    printf("Error: Unknown image format '%s'\n"
		   " Available formats:\n",format.c_str());
	    for(const char* formatName : pen_imageExporter::formatNames){
	      printf(" -%s\n",formatName);
	    }
	  }
	  return -2;	  
	}

	//Enable it
	enabledFormats[static_cast<int>(pen_imageExporter::toFormat(format.c_str()))] = true;
      }
    }

    if(verbose > 1){
      printf("Image export enabled for the following formats:\n");
      for(size_t i = 0; i < enabledFormats.size(); ++i){
	if(enabledFormats[i]){
	  printf(" -%s\n",
		 pen_imageExporter::toString(static_cast<pen_imageExporter::formatTypes>(i)));
	}
      }
    }
  }
  
  // Get number of threads to use
  //*******************************
  int auxThreads;
  if(config.read("simulation/threads",auxThreads) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("\n\nNumber of threads not specified, only one thread will be used.\n\n");
    }
    auxThreads = 1;
  }
  // ************************** MULTI-THREADING ***************************** //
#ifndef _PEN_USE_THREADS_
  else{
    printf("\n\nMulti-threading has not been activated during compilation"
	   ", only one thread will be used\n\n");
  }
  auxThreads = 1;
#else

  if(auxThreads <= 0){

    if(verbose > 1)
      printf("Automatic selection of threads enabled\n");
    
    unsigned int nConcurrency = std::thread::hardware_concurrency();

    if(nConcurrency > 0)
      auxThreads = nConcurrency;
    else{
      printf("The hardware concurrency value is not well defined or not computable."
	     " One thread will be used");
      auxThreads = 1;
    }
  }
  
#endif
  // ************************ MULTI-THREADING END *************************** //

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_

  //Check if the number of threads has been specified for this rank
  int auxThreadsMPI;
  char keyMPIthreads[60];
  if(sprintf(keyMPIthreads,"simulation/rank/%d/threads",rank) < 0){
    if(verbose > 0)
      printf("Error: Unable to create the configuration key to read the number of threads in the MPI rank %d\n",rank);
    return -2;
  }
  
  if(config.read(keyMPIthreads,auxThreadsMPI) == INTDATA_SUCCESS){
    if(verbose > 1){
      printf("Number of threads in rank %d set to %d.\n",rank,auxThreadsMPI);
    }

    if(auxThreadsMPI <= 0){
      if(verbose > 0)
	printf("Warning: Invalid number of threads, the value will not be considered\n");
    }else{
      auxThreads = auxThreadsMPI;
    }
  }
  
#endif
  // ***************************** MPI END ********************************** //
  
  unsigned nthreads = (unsigned)auxThreads;
  
  //If add dumps is enabled, set the number of
  //threads to 2 to be able to add the dump files
  if(addDumps){
    if(verbose > 1){
      printf("\n *** Number of simulating threads overwritten because 'addDumps' is enabled ***\n");
    }
    nthreads = 2;
  }

  if(verbose > 1){
    printf("\nNumber of simulating threads: %u\n",nthreads);
    if(nthreads > 1){
      printf("Initial random seeds will be selected using \"rand0\" function to ensure truly independent sequences of random numbers.\n");
    }
  }

  // Create a simulation config for each thread
  //*********************************************
  std::vector<penred::simulation::simConfig> simConfigs(nthreads);

  //Copy basic configuration
  for(unsigned i = 0; i < nthreads; i++){
    simConfigs[i].iThread = i;
    simConfigs[i].copyCommonConfig(baseSimConfig);
    //Set, by default, std::cout as output stream
    simConfigs[i].setOutstream(std::cout);
  }

  // Get initial seed pair number
  //*******************************
  int nseedPair;
  if(config.read("simulation/seedPair",nseedPair) != INTDATA_SUCCESS){
     nseedPair = -1;
  }
  else if(nseedPair < 0 || nseedPair > 1000){
    if(verbose > 0){
      printf("Invalid initial seed pair number %d\n",nseedPair);
      printf("Available seed pair range is [0,1000]\n");
    }
    return -2;
  }
  else if(verbose > 1){
    printf("Selected rand0 seed pair number: %d\n",nseedPair);
  }
  
  // Get CPU affinity option
  //*******************************
  bool CPUaffinity = false;
  if(config.read("simulation/thread-affinity",CPUaffinity) != INTDATA_SUCCESS){
    CPUaffinity = false;
  }

  if(verbose > 1){
    if(CPUaffinity){
      printf("Threads CPU affinity enabled\n");
    }
    else{
      printf("Threads CPU affinity disabled.\n");
      printf("    Enable it with: 'simulation/thread-affinity true'\n");
    }
  }

  // Get output options
  //*******************************
  bool ASCIIResults = true;
  if(config.read("simulation/ascii-results",ASCIIResults) != INTDATA_SUCCESS){
    ASCIIResults = true;
  }
  if(verbose > 1){
    if(!ASCIIResults){
      printf("ASCII results write disabled\n");
    }
  }

  bool finalDump = false;
  if(config.read("simulation/finalDump",finalDump) != INTDATA_SUCCESS){
    finalDump = false;
  }

  // Get recovery dump filename
  //***************************************
  std::string dump2read;  
  if(config.read("simulation/dump2read",dump2read) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("\n\nNo recovery dump filename specified.\n\n");
    }
    dump2read.clear();
  }

  // Check if the user only needs to read dump and write ASCII tally report
  bool dump2ascii;
  if(config.read("simulation/dump2ascii",dump2ascii) != INTDATA_SUCCESS){
    dump2ascii = false;
  }
  else if(dump2ascii && dump2read.length() == 0){
    if(verbose > 0){
      printf("\n\nError: Conversion from dump to ascii required but no dump file specified (field 'simulation/dump2read').\n\n");
    }
    return -2;
  }

  // Get context log filename
  //***************************************
  std::string contextlogfile;  
  if(config.read("simulation/contextlogfile",contextlogfile) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("\n\nNo context log filename specified.\n\n");
    }
    contextlogfile.clear();
  }

  // Get load balance configuration
  //**********************************

  std::chrono::seconds::rep balanceInterval;
  // ******************************* LB ************************************ //
  double balanceIntervald = 1.0e9;
  std::string LBhost;
  std::string LBurl;

#ifdef _PEN_USE_LB_HTTP_
  bool LBhttp = false;
#endif
  int LBport = -1;
  int LBworker = 0;
  std::string CAfile, certFile, keyFile, password, hostname;  
#ifdef _PEN_USE_LB_
  int errAuxEsp = 99;
  if((errAuxEsp = config.read("loadBalance/balance-interval",balanceIntervald)) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("\n\nTime between load balances not specified.\n");
      printf("Read of parameter 'loadBalance/balance-interval' "
	     "failed because: %s\n\n",
	     pen_parserError(errAuxEsp));
    }
  }
  else{
    if(balanceIntervald <= 10.0){
      balanceIntervald = 10.1;
      if(verbose > 1)
	printf("Balance interval must be greater than 10s.\n");
    }
    if(verbose > 1)
      printf("Balance interval set to %E.\n",balanceIntervald);
  }
  
  if(config.read("loadBalance/host",LBhost) == INTDATA_SUCCESS){
    if(config.read("loadBalance/port",LBport) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("Unable to read external balance port. Integer expected\n");
	return -2;
      }
    }
    if(config.read("loadBalance/worker",LBworker) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("Unable to read external balance worker ID. Integer expected\n");
	return -2;
      }
    }
    //Try to read CA cert, client cert, client key, key password and expected hostname
    config.read("loadBalance/CA-cert",CAfile);
    config.read("loadBalance/cert-file",certFile);
    config.read("loadBalance/key-file",keyFile);
    config.read("loadBalance/key-password",password);
    config.read("loadBalance/hostname",hostname);
    
    if(verbose > 1){
      printf("Configured external load balance server at %s:%d\n",LBhost.c_str(),LBport);
      printf(" SSL configuraiton:\n"
	     "       CA: %s\n"
	     "     cert: %s\n"
	     "      key: %s\n"
	     " password: %s\n"
	     " hostname: %s\n",
	     CAfile.c_str(),certFile.c_str(),keyFile.c_str(),
	     password.c_str(),hostname.c_str());
    }
    
  }
// ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
  else if(config.read("loadBalance/url",LBurl) == INTDATA_SUCCESS){
    if(verbose > 1)
      printf("Enabled external HTTP load balance server at url: %s\n",LBurl.c_str());
    if(config.read("loadBalance/worker",LBworker) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("Unable to read external balance worker ID. Integer expected\n");
	return -2;
      }
    }    
    LBhttp = true;
  }
#endif
// **********************   HTTP END   **********************
  else{
    if(verbose > 1){
      printf("External load balance host has not been specifeid\n");
    }
  }
  
  
  // ***************************** LB END ********************************** //
#endif
  balanceInterval = static_cast<std::chrono::seconds::rep>(balanceIntervald);
  
  //*******************************
  // Obtain initial random seeds
  //*******************************

  const size_t nRand0Seeds = 1001;

  int seedPos0 = nseedPair; //First seed pair to use
  //Set required seeds equal to the number of threads
  int nReqSeeds = nthreads; 

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_

  //Each MPI process requires a couple of initial seeds
  //equal to the number of local threads.

  nReqSeeds = 0;
  for(int irank = 0; irank < mpiSize; ++irank){

    // Send/Receive the number of threads used by the rank 'irank'
    unsigned iRankThreads = nthreads;
    MPI_Bcast(&iRankThreads,1,MPI_UNSIGNED,irank,MPI_COMM_WORLD);

    if(irank == rank){ //this rank did the last broadcast
      //Store the the first seed position to be used
      if(seedPos0 >= 0)
	seedPos0 = (seedPos0 + nReqSeeds) % nRand0Seeds;
      else
	seedPos0 = nReqSeeds % nRand0Seeds; //Avoid negative values of seedPos0
    }
    
    nReqSeeds += iRankThreads;
  }
  
#endif
  // ***************************** MPI END ********************************** //

  //First, ensure that rand0 can provide sufficient
  //initial seeds for all MPI processes and threads.
  
  if(nReqSeeds > 1001){
    printf("Error: Unsuficient initial seeds for all processes "
	   "and seeds (%d required).\n"
	   "        Please, use less threads to use, as maximum,"
	   " 1001 initial seeds.\n",nReqSeeds);
    return -3;
  }
  
  if(nthreads > 1){
    for(unsigned i = 0; i < nthreads; i++){
      int seedPos = i;
      if(seedPos0 >= 0)
	seedPos = (seedPos0+i) % nRand0Seeds;
      simConfigs[i].setSeeds(seedPos);
    }
  }
  else{    
    if(seedPos0 >= 0){
      simConfigs[0].setSeeds(seedPos0);
    }    
  }
  
  //****************************
  // Create particle generators
  //****************************
  std::vector<pen_specificStateGen<pen_particleState>> genericSources;
  std::vector<pen_specificStateGen<pen_state_gPol>> polarisedGammaSources;
      
  //Create and configure sources
  err = createParticleGenerators(genericSources,
				 polarisedGammaSources,
				 nthreads,config,verbose);

  if(err != 0){
    if(verbose > 0){
      if(err > 0){
	printf("\n");
	printf("Errors at sources creation and configuration.\n");
	printf("\n");
      }
    }
    return -3;
  }
  if(verbose > 0)
    printf("\n");

  //Set interval between checks
  for(auto& source : genericSources){
    source.setCheckTime(balanceInterval);
    source.setLBthreshold(static_cast<unsigned long long>(balanceIntervald/2));
    //Check if external balance has been configured
    // ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
    if(LBhttp){
      int ec;
      ec=source.setLBserver(LBworker,LBurl.c_str(),verbose);
      if(ec != 0){
	if(verbose > 0){
	  printf("Unable to configure load balance HTTP server\n"
		 "    Error code: %d\n",ec);
	  return -3;
	}
      }      
    }
    else
#endif
// **********************   HTTP END   **********************
    if(LBhost.length() > 0){
      int ec;
      ec=source.setLBserver(LBworker,
			    LBhost.c_str(),
			    std::to_string(LBport++).c_str(),
			    CAfile.length() == 0 ? nullptr : CAfile.c_str(),
			    certFile.length() == 0 ? nullptr : certFile.c_str(),
			    keyFile.length() == 0 ? nullptr : keyFile.c_str(),
			    password.length() == 0 ? nullptr : password.c_str(),
			    hostname.length() == 0 ? nullptr : hostname.c_str(),
			    verbose);
      if(ec != 0){
	if(verbose > 0){
	  printf("Unable to configure load balance server\n"
		 "    Error code: %d\n",ec);
	  return -3;
	}
      }
    //Notice that each source uses a unique port
    }
  }
  for(auto& source : polarisedGammaSources){
    source.setCheckTime(balanceInterval);
    source.setLBthreshold(static_cast<unsigned long long>(balanceIntervald/2));
    //Check if external balance has been configured
    // ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
    if(LBhttp){
      int ec;
      ec=source.setLBserver(LBworker,LBurl.c_str(),verbose);
      if(ec != 0){
	if(verbose > 0){
	  printf("Unable to configure load balance HTTP server\n"
		 "    Error code: %d\n",ec);
	  return -3;
	}
      }      
    }
    else
#endif
// **********************   HTTP END   **********************    
    if(LBhost.length() > 0){
      int ec;
      ec=source.setLBserver(LBworker,
			    LBhost.c_str(),
			    std::to_string(LBport++).c_str(),
			    CAfile.length() == 0 ? nullptr : CAfile.c_str(),
			    certFile.length() == 0 ? nullptr : certFile.c_str(),
			    keyFile.length() == 0 ? nullptr : keyFile.c_str(),
			    password.length() == 0 ? nullptr : password.c_str(),
			    hostname.length() == 0 ? nullptr : hostname.c_str(),
			    verbose);
      if(ec != 0){
	if(verbose > 0){
	  printf("Unable to configure load balance server\n"
		 "    Error code: %d\n",ec);
	  return -3;
	}
      }
    //Notice that each source uses a unique port
    }
  }

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  //Set MPI interval between checks
  for(auto& source : genericSources){
    source.setCheckTimeMPI(balanceInterval);
  }
  for(auto& source : polarisedGammaSources){
    source.setCheckTimeMPI(balanceInterval);
  }
#endif
  // ******************************* MPI ************************************ //

  
  //****************************
  // Context
  //****************************

  //Create elements data base
  pen_elementDataBase* elementsDB = new pen_elementDataBase;
  
  //Create simulation context
  pen_context context(*elementsDB);
  
  //Get global maximum energy
  double globEmax = -1.0;
  for(unsigned i = 0; i < genericSources.size(); i++){
    double localEmax = genericSources[i].maxEnergy();
    if(genericSources[i].kpar == PEN_POSITRON){
      //  ----  Positrons eventually give annihilation gamma rays. The maximum
      //        energy of annihilation photons is .lt. 1.21*(E0+me*c**2).
      localEmax = 1.21*(localEmax+5.12E5);
    }
    if(globEmax < localEmax)
      globEmax = localEmax;
  }
  for(unsigned i = 0; i < polarisedGammaSources.size(); i++){
    if(globEmax < polarisedGammaSources[i].maxEnergy())
      globEmax = polarisedGammaSources[i].maxEnergy();
  }

  if(verbose > 1){
    printf("Maximum global energy: %14.2E eV\n",globEmax);
  }

  //****************************
  // Materials
  //****************************

  printf("\n\n------------------------------------\n\n");
  printf(" **** Materials ****\n");
  printf(" *******************\n\n");
  
  std::string matFilenames[constants::MAXMAT];
  
  if(createMaterials(context,matFilenames,config,globEmax,verbose) != 0)
    return -4;


  //Initialize context with specified materials
  FILE* fcontext = nullptr;
  if(contextlogfile.length() > 0){
    //fcontext = fopen("context.rep","w");
    fcontext = fopen(contextlogfile.c_str(),"w");
    if(fcontext == nullptr){
      printf("Error: unable to open context log file '%s'\n",contextlogfile.c_str());
    }
  }
  if(context.init(globEmax,fcontext,verbose,matFilenames) != PEN_SUCCESS){
    if(fcontext != nullptr) fclose(fcontext);
    printf("Error at context initialization.\n");
    return -5;
  }
  if(fcontext != nullptr) fclose(fcontext);
  
  //****************************
  // Geometry parameters
  //****************************

  printf("\n\n------------------------------------\n\n");
  printf(" **** Geometry parameters\n\n");

  wrapper_geometry* geometry;
  
  //Create and configure geometry
  if(createGeometry(geometry,config,context,verbose) != 0)
    return -6;

  //Set geometry to context
  context.setGeometry(geometry);
  
  //Set source geometry
  for(unsigned i = 0; i < genericSources.size(); i++)
    genericSources[i].setGeometry(geometry);  
  for(unsigned i = 0; i < polarisedGammaSources.size(); i++)
    polarisedGammaSources[i].setGeometry(geometry);  

  //------------------------------------------------------------
  
  //Check if all required materials for the specified geometry have been created
  bool usedMaterials[constants::MAXMAT+1];
  geometry->usedMat(usedMaterials);

  unsigned missingMats = 0;
  for(unsigned i = 1; i <= constants::MAXMAT; i++){  //Avoid void material (0)
    if(matFilenames[i-1].length() == 0 && usedMaterials[i]){
      if(verbose > 0){
	printf("Error: Material %d is required by the geometry, but has not been configured.\n", i);
      }
      missingMats++;
    }
    else if(matFilenames[i-1].length() > 0 && !usedMaterials[i]){
      if(verbose > 0){
	printf("Warning: Material %d has been configured but is not used at specified geometry.\n",i);
      }
    }
  }

  if(missingMats > 0){
    printf("Error: Missing %d materials to use the specified geometry.\n",missingMats);
    return -6;
  }

  //Update absortion energies ussing geometry information
  err = context.updateEABS();
  if(err != 0){
    printf("Error calculating absorption energies\n");
    return -8;
  }
  
  //****************************
  // Create tallies
  //****************************
  std::vector<pen_commonTallyCluster> talliesVect;

  //Create one tally group per thread
  talliesVect.resize(nthreads);

  //Init tallies
  if(createTallies(talliesVect, nthreads, context, config, verbose) != 0)
    return -9;

  printf("\n%d tallies created for each thread.\n",talliesVect[0].numTallies());
  if(talliesVect[0].numTallies() == 0){
    printf("The simulation will not extract any information.\n");
    return 0;
  }
  
  // Add dumps command option
  //*******************************
  if(addDumps){

    bool firstLoad = true;
    unsigned long long totalSimHists = 0;
    for(int idump = 3; idump < argc; ++idump){

      unsigned long long simulated;
      int seed1,seed2;
      int lastSource;
      unsigned long long sourceHists;
      
      const char* filename = argv[idump];
      
      if(firstLoad){
	firstLoad = false;
	//Read the first dump file
	printf("\nLoading first dump file: '%s'\n",filename);
	int errDump = talliesVect[0].readDumpfile(filename,
						  simulated,
						  seed1,seed2,
						  lastSource,
						  sourceHists,
						  verbose);

	if(errDump != 0){
	  if(verbose > 0){
	    printf("Error loading dumped data file '%s': %d\n",filename,errDump);
	  }
	  return -10;
	}	
      }
      else{

	//Read and add the dump file
	printf("\nAdding dump file: '%s'\n",filename);
	int errDump = talliesVect[1].readDumpfile(filename,
						  simulated,
						  seed1,seed2,
						  lastSource,
						  sourceHists,
						  verbose);	
	if(errDump != 0){
	  if(verbose > 0){
	    printf("Error loading dumped data file '%s': %d\n",filename,errDump);
	  }
	  return -10;
	}

	int errSum = talliesVect[0].sum(talliesVect[1],verbose);

	if(errSum != 0){
	  if(verbose > 0){
	    printf("Error adding dumped data from file '%s': %d\n",filename,errSum);
	  }
	  return -11;
	}
      }

      totalSimHists += simulated;
      if(lastSource >= 0){
	if(verbose > 1){
	  printf("Warning: Dump file '%s' generated from a "
		 "non finished simulation\n",filename);
	}
      }
      if(verbose > 1){
	printf(" - Simulated histories: %llu\n",simulated);
	printf(" - Last seeds: %d %d\n",seed1,seed2);
      }
    }

    printf("\nTotal simulated histories: %llu\n",totalSimHists);
    
    if(ASCIIResults){
      talliesVect[0].saveData(totalSimHists);
    }
    
    talliesVect[0].dump2file("mergedDump.dump",totalSimHists,-1,-1,-1,0ull,verbose);

    
    
    return 0;
  }
  //*******************************

  // Check if we are restoring a simulation from dumpfile  
  if(dump2read.length() > 0){

    if(verbose > 1){
      printf("Load dump files '%s' for %d threads\n",dump2read.c_str(), nthreads);
    }

    std::string auxstr(dump2read); 
    std::size_t found = auxstr.find_last_of("/\\");

    //Read dump file for each thread
    for(unsigned i = 0; i < nthreads; i++){
      //Create filename

      std::string filenameDump;
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_

      if(found != std::string::npos){
	filenameDump = auxstr.substr(0,found+1) + std::string("MPI") + std::to_string(rank) +
	  std::string("-th") + std::to_string(i) + auxstr.substr(found+1);
      }else{
	filenameDump = std::string("MPI") + std::to_string(rank) +
	  std::string("-th") + std::to_string(i) + dump2read;
      }
#else
  // ***************************** MPI END ********************************** //
      if(found != std::string::npos){
	filenameDump = auxstr.substr(0,found+1) + std::string("th") + std::to_string(i) + auxstr.substr(found+1);
      } else{
	filenameDump = std::string("th") + std::to_string(i) + dump2read;
      }
#endif      

      //Read dump file
      int errDump = simConfigs[i].readTallyDump(filenameDump.c_str(),talliesVect[i]);
      
      if(errDump != 0){
	if(verbose > 0){
	  printf("Error loading dumped data for thread %d: %d\n",i,errDump);
	}
	return -10;
      }

      if(verbose > 1){
	//Report data in dump file and finish execution
	printf(" *** Dump information of thread %u:\n\n",i);
	printf("  Total simulated histories: %llu\n",
	       simConfigs[i].getSimulatedInFinished() +
	       simConfigs[i].getInitiallySimulatedInFirstSource());
	int seed1, seed2;
	simConfigs[i].getSeeds(seed1,seed2);
	printf("                 Last seeds: %9d %9d\n",seed1,seed2);
	printf("    Next source to simulate: %d\n",simConfigs[i].getFirstSourceIndex());
	printf(" Source histories simulated: %llu\n",
	       simConfigs[i].getInitiallySimulatedInFirstSource());	
      }
      
      if(dump2ascii){
	// Save the data of this thread
	talliesVect[i].saveData(simConfigs[i].getSimulatedInFinished() +
				simConfigs[i].getInitiallySimulatedInFirstSource());
      }
    }

    if(dump2ascii){
      //Finish execution
      return 0;
    }
    
    //Skip the required iterations from each source
    int iSource = 0;
    std::vector<unsigned long long>toSkip(nthreads);
    // Generic sources
    for(auto& source : genericSources){

      //Calculate the iterations per thread to skip for this source
      for(size_t ithread = 0; ithread < nthreads; ++ithread){
	if(simConfigs[ithread].getFirstSourceIndex() == iSource){
	  toSkip[ithread] = simConfigs[ithread].getInitiallySimulatedInFirstSource();
	}else{
	  toSkip[ithread] = 0;
	}
      }

      source.skip(toSkip);
      ++iSource;
    }
    // Polarised sources
    for(auto& source : polarisedGammaSources){

      //Calculate the iterations per thread to skip for this source
      for(size_t ithread = 0; ithread < nthreads; ++ithread){
	if(simConfigs[ithread].getFirstSourceIndex() == iSource){
	  toSkip[ithread] = simConfigs[ithread].getInitiallySimulatedInFirstSource();
	}else{
	  toSkip[ithread] = 0;
	}
      }
      source.skip(toSkip);
      ++iSource;
    }
    
  }

  
  
  //****************************
  // Variance Reduction 
  //****************************
  pen_VRCluster<pen_particleState> genericVR;
  pen_VRCluster<pen_state_gPol> photonVR;
  int vrRet = setVarianceReduction(context,*geometry,config,
				   genericVR,photonVR,verbose);
  if(vrRet < 0){
    if(verbose > 0){
      printf("Error on variance reduction section.\n");
      printf("             Error code: %d\n",vrRet);
    }
    return -11;
  }
  if(vrRet > 0 && verbose > 1){
    printf("No variance reduction technics enabled.\n");
  }  
    
  //****************************
  // History loop
  //****************************

  fflush(stdout);
  
  // ************************** MULTI-THREADING ***************************** //
#ifdef _PEN_USE_THREADS_
  std::vector<std::thread> simThreads;
#endif
  // ************************ MULTI-THREADING END *************************** //

  pen_timer timer;
  double time0 = CPUtime();
  
  for(unsigned i = 0; i < nthreads; i++){
    talliesVect[i].run_beginSim();
  }

  if(simConfigs[0].getFirstSourceIndex() > 0){
    if(verbose > 1)
      printf("Skip already simulated sources (%d)\n",simConfigs[0].getFirstSourceIndex());
  }

  if(verbose > 1){
    printf("Initialization processing time: %E s\n", initializationTimer.timer());
  }

  //Substract initialization time to maximum simulation time
  long long int initTime = static_cast<long long int>(initializationTimer.timer());
  for(unsigned ithread = 0; ithread < nthreads; ++ithread)
    simConfigs[ithread].maxSimTime -= initTime*1000ll;
  
  //Iterate over generic sources
  for(unsigned iSource = simConfigs[0].getFirstSourceIndex();
      iSource < genericSources.size(); ++iSource){

    //Check remaining simulation time
    if(simConfigs[0].maxSimTime <= 0)
      break; //Finish the simulation
    

  // ************************** MULTI-THREADING ***************************** //
#ifdef _PEN_USE_THREADS_
    //Run simulations for each thread
    if(nthreads > 1){
      // Multi-thread
      //****************
      
      for(unsigned ithread = 0; ithread < nthreads; ++ithread){
	// Simulate
	simThreads.push_back(std::thread(simulate<pen_particleState>,
					 std::ref(simConfigs[ithread]),
					 std::ref(context),
					 std::ref(genericSources[iSource]),
					 std::ref(talliesVect[ithread]),
					 std::ref(genericVR),
					 std::ref(photonVR)));


#ifdef _PEN_UNIX_
	//Set affinity for unix systems using pthreads
	if(CPUaffinity){
	  cpu_set_t cpuset;
	  CPU_ZERO(&cpuset);
	  CPU_SET(ithread,&cpuset);
	  int raff = pthread_setaffinity_np(simThreads.back().native_handle(),
					    sizeof(cpu_set_t), &cpuset);
	  if(raff != 0){
	    printf(" Unable to set affinity for thread %d\n",ithread);
	  }
	  else{
	    printf(" Affinity for thread %d set to CPU %d\n",ithread,ithread);
	  }
	}
#endif      
      }

      //Iterate for user instructions if required
      if(userInteractiveEnabled)
	checkUserInput(userInteractiveFilename.c_str(),
		       genericSources[iSource].name,
		       simConfigs);
    
      //Join threads
      for(unsigned ithread = 0; ithread < nthreads; ++ithread){
	simThreads[ithread].join();
	//Update remaining simulation time
	simConfigs[0].maxSimTime =
	  std::min(simConfigs[0].maxSimTime, simConfigs[ithread].maxSimTime);
      }

      //Update maximum simulation times
      for(unsigned ithread = 0; ithread < nthreads; ++ithread){
	simConfigs[ithread].maxSimTime = simConfigs[0].maxSimTime;
      }
      
      //Clear threads
      simThreads.clear();    
    }
    else{

#endif
  // ************************ MULTI-THREADING END *************************** //

      // Single thread
      //****************
      
      // Simulate      
      simulate(simConfigs[0], context,
	       genericSources[iSource],
	       talliesVect[0],
	       genericVR, photonVR);
      
#ifdef _PEN_USE_THREADS_
    }
#endif
    
  }

  //Iterate over polarized sources
  for(unsigned iSource = 0; iSource < polarisedGammaSources.size(); iSource++){

    //Check remaining simulation time
    if(simConfigs[0].maxSimTime <= 0)
      break; //Finish the simulation

    
    //Run simulations for each thread

  // ************************** MULTI-THREADING ***************************** //
#ifdef _PEN_USE_THREADS_
    if(nthreads > 1){
      // Multi-thread
      //***************
      
      for(unsigned ithread = 0; ithread < nthreads; ++ithread){
	// Simulate
	simThreads.push_back(std::thread(simulate<pen_state_gPol>,
					 std::ref(simConfigs[ithread]),
					 std::ref(context),
					 std::ref(polarisedGammaSources[iSource]),
					 std::ref(talliesVect[ithread]),
					 std::ref(genericVR),
					 std::ref(photonVR)));

#ifdef _PEN_UNIX_
	//Set affinity for unix systems using pthreads
	if(CPUaffinity){
	  cpu_set_t cpuset;
	  CPU_ZERO(&cpuset);
	  CPU_SET(ithread,&cpuset);
	  int raff = pthread_setaffinity_np(simThreads.back().native_handle(),
					    sizeof(cpu_set_t), &cpuset);
	  if(raff != 0){
	    printf(" Unable to set affinity for thread %d\n",ithread);
	  }
	  else{
	    printf(" Affinity for thread %d set to CPU %d\n",ithread,ithread);
	  }
	}
#endif
	
      }

      //Iterate for user instructions if required
      if(userInteractiveEnabled)
	checkUserInput(userInteractiveFilename.c_str(),
		       polarisedGammaSources[iSource].name,
		       simConfigs);      

      //Join threads
      for(unsigned ithread = 0; ithread < nthreads; ++ithread){
	simThreads[ithread].join();
	//Update remaining simulation time
	simConfigs[0].maxSimTime =
	  std::min(simConfigs[0].maxSimTime, simConfigs[ithread].maxSimTime);

      }

      //Update maximum simulation times
      for(unsigned ithread = 0; ithread < nthreads; ++ithread){
	simConfigs[ithread].maxSimTime = simConfigs[0].maxSimTime;
      }

      //Clear threads
      simThreads.clear();
      
    }
    else{
#endif
  // ************************ MULTI-THREADING END *************************** //

      // Single thread
      //****************
      
      // Simulate 
      simulate<pen_state_gPol>(simConfigs[0],
			       context,
			       polarisedGammaSources[iSource],
			       talliesVect[0],
			       genericVR,
			       photonVR);
      
#ifdef _PEN_USE_THREADS_
    }
#endif
    
  }

  
  //Calculate the total number of simulated hists
  unsigned long long totalHists = 0.0;
  for(unsigned ithread = 0; ithread < nthreads; ++ithread){
    totalHists += simConfigs[ithread].getSimulatedInFinished();
  }
    
  //End simulations
  for(unsigned ithread = 0; ithread < nthreads; ithread++){
    talliesVect[ithread].run_endSim(simConfigs[ithread].getSimulatedInFinished());
  }

  double simtime = timer.timer();
  double CPUendSim = CPUtime();
  double usertime = CPUendSim-time0;

  if(verbose > 1){
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
    printf("Rank %u: Simulation finished, starting results reduce step\n",rank);
#else    
  // ***************************** MPI END ********************************** //
    printf("Simulation finished, starting results reduce step\n");
#endif
  }

  //If enabled, print final dumps
  if(finalDump){
    for(unsigned ithread = 0; ithread < nthreads; ithread++){
      int seeds1, seeds2;
      simConfigs[ithread].getSeeds(seeds1, seeds2);
      talliesVect[ithread].dump2file(simConfigs[ithread].dumpFilename.c_str(),
				     simConfigs[ithread].getSimulatedInFinished(),
				     seeds1,seeds2,-1,0ull,verbose);
    }
  }
  
  //Sum tallies of all threads
  for(unsigned ithread = 1; ithread < nthreads; ithread++){
    talliesVect[0].sum(talliesVect[ithread],verbose);
  }

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  unsigned long long localHists = totalHists;
  //Reduce the results stored at the first thread of all ranks
  talliesVect[0].reduceMPI(totalHists,MPI_COMM_WORLD,verbose);

  double postProcessTime = CPUtime() - CPUendSim;
  
  //Save results stored at the first thread of rank 0 process
  if(rank == 0){
    if(ASCIIResults)
      talliesVect[0].saveData(totalHists);
    else
      talliesVect[0].dump2file("results.dump",totalHists,-1,-1,-1,0ull,verbose);
  }

  //Print local report information

  printf("\n\n");
  printf("\n*********** Rank %03u *************\n",rank);
  printf("Simulated histories: %20llu \n",localHists);
  printf("Simulation real time: %12.4E s\n",simtime);
  printf("Simulation user time: %12.4E s\n",usertime);
#if defined(_MSC_VER)
  printf("Histories per second: %12.4E\n", static_cast<double>(localHists) / usertime);
#else
  printf("Histories per second and thread: %12.4E\n",static_cast<double>(localHists)/usertime);
  printf("Histories per second: %12.4E\n",static_cast<double>(localHists)/(usertime/double(nthreads)));
#endif
  printf("Results processing time: %12.4E s\n",postProcessTime);
  printf("\n*********** END REPORT *************\n");
  fflush(stdout);
  

  //Calculate total user time and used threads
  double globalUserTime = 0.0;
  MPI_Reduce(&usertime,&globalUserTime,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  unsigned totalThreads = 0;
  MPI_Reduce(&nthreads,&totalThreads,1,MPI_UNSIGNED,MPI_SUM,0,MPI_COMM_WORLD);

  double globalSimTime = 0.0;
  //Obtain maximum simulation time
  MPI_Reduce(&simtime,&globalSimTime,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  //Close stdout
  fclose(stdout);

  //Print joined reports on stderr to print it on global log
  //Print ordered local report information
  for(int ip = 0; ip < mpiSize; ++ip){

    if(rank == ip){
      fprintf( stderr,"\n\n");
      fprintf( stderr,"\n*********** Rank %03u *************\n",rank);
      fprintf( stderr,"Simulated histories: %20llu\n",localHists);
      fprintf( stderr,"Simulation real time: %12.4E s\n",simtime);
      fprintf( stderr,"Simulation user time: %12.4E s\n",usertime);
#if defined(_MSC_VER)
      fprintf(stderr, "Histories per second: %12.4E\n", static_cast<double>(localHists) / usertime);
#else
      fprintf( stderr,"Histories per second and thread: %12.4E\n",static_cast<double>(localHists)/usertime);
      fprintf( stderr,"Histories per second: %12.4E\n",static_cast<double>(localHists)/(usertime/double(nthreads)));
#endif
      fprintf( stderr,"Results processing time: %12.4E s\n",postProcessTime);
      fprintf( stderr,"\n*********** END REPORT *************\n");
      fflush(stderr);
      std::this_thread::sleep_for (std::chrono::seconds(1));
    }
    MPI_Barrier(MPI_COMM_WORLD); //Sync all MPI processes
  }

  //Print global report
  if(rank == 0){
    fprintf( stderr,"\n\n");
    fprintf( stderr,"\n*********** Global Report *************\n");
    fprintf( stderr,"Simulated histories: %20llu\n",totalHists);
    fprintf( stderr,"Simulation real time: %12.4E s\n",globalSimTime);
    fprintf( stderr,"Simulation user time: %12.4E s\n",globalUserTime);
#if defined(_MSC_VER)
    fprintf(stderr, "Histories per second: %12.4E\n", static_cast<double>(totalHists) / globalUserTime);
#else
    fprintf( stderr,"Histories per second and thread: %12.4E\n",static_cast<double>(totalHists)/globalUserTime);
    fprintf( stderr,"Histories per second: %12.4E\n",static_cast<double>(totalHists)/(globalUserTime/double(totalThreads)));
#endif
    fprintf( stderr,"Results processing time: %12.4E s\n",postProcessTime);  
  }
  
#else
  // ***************************** MPI END ********************************** //

  double postProcessTime = CPUtime() - CPUendSim;
  
  //Save results stored at first thread
  if(ASCIIResults)
    talliesVect[0].saveData(totalHists);
  else
    talliesVect[0].dump2file("results.dump",totalHists,-1,-1,-1,0ull,verbose);

  for(size_t i = 0; i < enabledFormats.size(); ++i){
    if(enabledFormats[i]){
      talliesVect[0].exportImage(totalHists,
				 static_cast<pen_imageExporter::formatTypes>(i),false);
    }
  }
  
  printf("\n\nSimulated histories: %20llu\n",totalHists);
  printf("Simulation real time: %12.4E s\n",simtime);
  printf("Simulation user time: %12.4E s\n",usertime);
#if defined(_MSC_VER)
  printf("Histories per second: %12.4E\n", static_cast<double>(totalHists) / usertime);
#else
  printf("Histories per second and thread: %12.4E\n",static_cast<double>(totalHists)/usertime);
  printf("Histories per second: %12.4E\n",static_cast<double>(totalHists)/(usertime/double(nthreads)));
#endif
  printf("Results processing time: %12.4E s\n",postProcessTime);
  
#endif

  // ******************************* LB ************************************ //
#ifdef _PEN_USE_LB_
  
  //Print load balance reports
  if(genericSources.size() > 0){
    printf("Printing load balance reports for generic sources...");
    for(const auto& source : genericSources){
      FILE* fout = nullptr;

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
      char prefix[20];
      sprintf(prefix,"rank-%03d-",rank);
      if(rank == 0){
	std::string filenameLBMPI(prefix);
	filenameLBMPI += source.name + std::string("-MPI-generic-LBreport.rep");      
	fout = fopen(filenameLBMPI.c_str(),"w");
	source.task.printReportMPI(fout);
	fclose(fout);
      }
      std::string filenameLBTh(prefix);
      filenameLBTh += source.name + std::string("-generic-LBreport.rep");
#else
      std::string filenameLBTh(source.name);
      filenameLBTh += std::string("-generic-LBreport.rep");
#endif
  // ***************************** MPI END ********************************** //      
      
      fout = fopen(filenameLBTh.c_str(),"w");
      source.task.printReport(fout);
      fclose(fout);
    }
    printf(" Done!\n");
  }
  if(polarisedGammaSources.size() > 0){
    printf("Printing load balance reports for polarised sources...");
    for(const auto& source : polarisedGammaSources){
      FILE* fout = nullptr;

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
      char prefix[20];
      sprintf(prefix,"rank-%03d-",rank);
      std::string filenameLBMPI(prefix);
      filenameLBMPI += source.name + std::string("-MPI-polarised-LBreport.rep");      
      fout = fopen(filenameLBMPI.c_str(),"w");
      source.task.printReportMPI(fout);
      fclose(fout);

      std::string filenameLBTh(prefix);
      filenameLBTh += source.name + std::string("-polarised-LBreport.rep");
#else
      std::string filenameLBTh(source.name);
      filenameLBTh += std::string("-polarised-LBreport.rep");
#endif
  // ***************************** MPI END ********************************** //      
      
      fout = fopen(filenameLBTh.c_str(),"w");
      source.task.printReport(fout);
      fclose(fout);
    }
    printf(" Done!\n");
  }
#endif
  // ***************************** LB END ********************************** //

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  // Finalize the MPI environment.
  MPI_Finalize();
#endif
  // ***************************** MPI END ********************************** //
  
  //Free memory
  delete geometry;
  geometry = nullptr;

  delete elementsDB;
  elementsDB = nullptr;
  
  return 0;
}

int createParticleGenerators(std::vector<pen_specificStateGen<pen_particleState>>& genericSources,
			     std::vector<pen_specificStateGen<pen_state_gPol>>& polarisedGammaSources,
			     const unsigned nthreads,
			     const pen_parserSection& config,
			     const unsigned verbose){

  int errG = 0;
  int errP = 0;
  int err = 0;
  
  //Extract source sections
  pen_parserSection genericSourceSection;
  pen_parserSection polarisedSourceSection;
  errG = config.readSubsection("sources/generic",genericSourceSection);
  errP = config.readSubsection("sources/polarized",polarisedSourceSection);
  if(errG != INTDATA_SUCCESS && errP != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createParticleGenerators: Error: Fields 'sources/generic' nor 'sources/polarized' don't exists. No source specified\n");
    }
    return -2;
  }

  //Extract source names
  std::vector<std::string> genericSourceNames;
  std::vector<std::string> polarisezSourceNames;
  if(errG == INTDATA_SUCCESS)
    genericSourceSection.ls(genericSourceNames);
  if(errP == INTDATA_SUCCESS)
    polarisedSourceSection.ls(polarisezSourceNames);

  unsigned nGenericSources = genericSourceNames.size();
  unsigned nPolarizedSources = polarisezSourceNames.size();  

  if(nGenericSources < 1 && nPolarizedSources < 1){
    if(verbose > 0)
      printf("createParticleGenerators: Error: Simulation requires, at last, one particle source.\n");
    return -3;    
  }

  //Resize source vectors to store all defined sourceSection
  genericSources.clear();
  if(nGenericSources > 0){
    //Use swap because soures can't be moved
    std::vector<pen_specificStateGen<pen_particleState>>
      aux(nGenericSources);    
    genericSources.swap(aux);
  }
  
  polarisedGammaSources.clear();
  if(nPolarizedSources > 0){
    //Use swap because soures can't be moved
    std::vector<pen_specificStateGen<pen_state_gPol>>
      aux(nPolarizedSources);    
    polarisedGammaSources.swap(aux);
  }

#if (defined _PEN_USE_MPI_ && defined _PEN_USE_LB_ )
  int nextTag = 5;
#endif

  //Iterate over all generic sources
  for(unsigned i = 0; i < nGenericSources; i++){
    //Get source section
    pen_parserSection genSection;
    if(genericSourceSection.readSubsection(genericSourceNames[i],genSection) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("createParticleGenerators: Error: Unable to extract section for generic source '%s'\n",genericSourceNames[i].c_str());
      }
      err++;
      continue;
    }
    else{
      //Set source name
      genericSources[i].name.assign(genericSourceNames[i]);

      //Load history number
      double nhists;
      if(genSection.read("nhist",nhists) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("createParticleGenerators: Error: Unable to read field 'nhist' "
		 "for generic source '%s'\n",genericSourceNames[i].c_str());
	}
	err++;
	continue;
      }

      if(nhists <= 0.5){
	if(verbose > 0){
	  printf("createParticleGenerators: Error on generic source %s. "
		 "Number of histories must be greater than zero\n",genericSourceNames[i].c_str());
	}
	err++;
	continue;	
      }
      
      //Load source
      if(configureSource(genericSources[i],
			 nthreads,
			 genSection,verbose) != 0){
	if(verbose > 0)
	  printf("createParticleGenerators: Error: Can't create and "
		 "configure source '%s'.\n",genericSourceNames[i].c_str());	
	err++;
      }

      //Init the source task
      unsigned long long uihists = static_cast<unsigned long long>(nhists+0.1);
      uihists = std::max(uihists,1llu);

  // ******************************* LB ************************************ //
#if (defined _PEN_USE_MPI_ && defined _PEN_USE_LB_ )
      int errTask = genericSources[i].initTask(nthreads,uihists,MPI_COMM_WORLD,
					       nextTag,nextTag+1,verbose);
      nextTag += 2;
#else
      int errTask = genericSources[i].initTask(nthreads,uihists,verbose);      
#endif
  // ***************************** LB END ********************************** //
      if(errTask != 0){
	if(verbose > 0)
	  printf("createParticleGenerators: Error on generic source %s. "
		 "Unable to init source task\n"
		 "   Error code: %d\n",
		 genericSourceNames[i].c_str(),errTask);
	err++;
	continue;
      }

      if(verbose > 1)
	printf("Histories to simulate at source %s: %llu\n",
	       genericSources[i].name.c_str(),genericSources[i].toDo());      
    }
  }

  //Iterate over all specific sources
  for(unsigned i = 0; i < nPolarizedSources; i++){
    //Get source section
    pen_parserSection genSection;
    if(polarisedSourceSection.readSubsection(polarisezSourceNames[i],genSection) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("createParticleGenerators: Error: Unable to extract section for polarized gamma source '%s'\n",polarisezSourceNames[i].c_str());
      }
      err++;
      continue;
    }
    else{
      //Set source name
      polarisedGammaSources[i].name.assign(polarisezSourceNames[i]);

      //Load history number
      double nhists;
      if(genSection.read("nhist",nhists) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("createParticleGenerators: Error: Unable to read field 'nhist' for polarized gamma source '%s'\n",polarisezSourceNames[i].c_str());
	}
	err++;
	continue;
      }

      if(nhists <= 0.5){
	if(verbose > 0){
	  printf("createParticleGenerators: Error on polarized gamma source %s. Number of histories must be greater than zero\n",polarisezSourceNames[i].c_str());
	}
	err++;
	continue;	
      }
      
      //Load source
      if(configureSource(polarisedGammaSources[i],
			 nthreads,
			 genSection,verbose) != 0){
	if(verbose > 0)
	  printf("createParticleGenerators: Error: Can't create and configure source '%s'.\n",polarisezSourceNames[i].c_str());	
	err++;
      }

      //Init the source task
      unsigned long long uihists = static_cast<unsigned long long>(nhists);
      uihists = std::max(uihists,1llu);

      
#if defined(_PEN_USE_MPI_) && defined(_PEN_USE_LB_)
      int errTask = polarisedGammaSources[i].initTask(nthreads,uihists,MPI_COMM_WORLD,
						      nextTag,nextTag+1,verbose);
      nextTag += 2;
#else
      int errTask = polarisedGammaSources[i].initTask(nthreads,uihists,verbose);      
#endif
      if(errTask != 0){
	if(verbose > 0)
	  printf("createParticleGenerators: Error on polarised source %s. "
		 "Unable to init source task\n"
		 "   Error code: %d\n",
		 polarisedGammaSources[i].name.c_str(),errTask);
	err++;
	continue;
      }

      if(verbose > 1)
	printf("Histories to simulate at source %s: %llu\n",
	       polarisedGammaSources[i].name.c_str(),
	       polarisedGammaSources[i].toDo());
    }
  }
      
  return err;
}

int createTallies(std::vector<pen_commonTallyCluster>& tallyGroups,
		  const size_t nthreads,
		  const pen_context& context,
		  const pen_parserSection& config,
		  const unsigned verbose){

  //Extract source section
  pen_parserSection talliesSection;
  int err = config.readSubsection("tallies",talliesSection);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createTallies: Error: Configuration 'tallies' section doesn't exists.\n");
    }
    return -1;
  }

  if(nthreads < 1){
    printf("createTallies: Error: Simulation requires, at last, one thread for execution.");
    return -2;
  }

  // ************************** MULTI-THREADING ***************************** //
#ifdef _PEN_USE_THREADS_
  //Create a vector of threads
  std::vector<std::thread> threads;
#endif
  // ************************** MULTI-THREADING ***************************** //

  //Get materials array
  const abc_material* mats[constants::MAXMAT];
  context.getMatBaseArray(mats);
  
  // ************************** MULTI-THREADING ***************************** //
#ifdef _PEN_USE_THREADS_
  //Configure non zero threads with no verbose
  for(unsigned j = 1; j < nthreads; j++){
    tallyGroups[j].name.assign("common");
    //Configure tallies
    threads.push_back(tallyGroups[j].configure_async(context.readGeometry(),
						     mats,
						     j,
						     talliesSection,
						     std::min(verbose,unsigned(1))));
  }
#endif
  // ************************** MULTI-THREADING ***************************** //
  
  //Configure main thread with verbose option
  tallyGroups[0].name.assign("common");    
  //Configure tallies
  tallyGroups[0].configure(context.readGeometry(),
			   mats,
			   0,
			   talliesSection,
			   verbose);
  
  // ************************** MULTI-THREADING ***************************** //
#ifdef _PEN_USE_THREADS_
  //Sincronize all threads
  for(unsigned j = 0; j < nthreads-1; j++){
    threads[j].join();
  }
#endif
  // ************************** MULTI-THREADING ***************************** //
  
  //Check errors
  unsigned failedClusters = 0;
  for(unsigned j = 0; j < nthreads; j++){
    err = tallyGroups[j].configureStatus();
    if(err != 0){
      if(verbose > 0)
	printf("createTallies: Error on tally cluster %u "
	       "creation and configuration (err code %d).\n",j,err);
      failedClusters++;
    }
  }

  if(failedClusters > 0)
    return -3;
  
  //Share configuration from thread 0 cluster to other threads
  failedClusters = 0;
  for(unsigned j = 0; j < nthreads; j++){

    err = tallyGroups[j].shareConfig(tallyGroups[0], verbose);
    if(err != 0){
      if(verbose > 0)
	printf("createTallies: Error on tally cluster %u. "
	       "Unable to get configuration from thread 0 (err code %d).\n",j,err);
      failedClusters++;
    }
  }

  if(failedClusters > 0)
    return -4;  
  
  return 0;
}

int createGeometry(wrapper_geometry*& geometry,
		   const pen_parserSection& config,
		   const pen_context& context,
		   const unsigned verbose){
  
  //Get geometry section
  pen_parserSection geometrySection;
  int err = config.readSubsection("geometry",geometrySection);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: Configuration 'geometry' section doesn't exist.\n");
    }
    return -2;
  }

  //Append material information to geometry section
  
  for(unsigned imat = 0; imat < context.getNMats(); imat++){
    char key[400];
    sprintf(key,"materials/mat%03d/ID",imat+1);
    geometrySection.set(key,(int)imat+1);
    sprintf(key,"materials/mat%03d/density",imat+1);
    geometrySection.set(key,context.readBaseMaterial(imat).readDens());
  }
  
  //Get geometry type
  std::string geoType;
  if(geometrySection.read("type",geoType) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createGeometry: Error: field 'geometry/type' not specified. String expected.\n");
    }
    return -3;
  }

  //Create geometry
  geometry = penGeoRegister_create(geoType.c_str());
  if(geometry == nullptr){
    if(verbose > 0){
      printf("createGeometry: Error creating geometry instance of type '%s'\n", geoType.c_str());
    }
    return -4;
  }

  //Configure geometry  
  geometry->name.assign("geometry");    
  if(geometry->configure(geometrySection,verbose) != 0){
    if(verbose > 0)
      printf("createGeometry: Error: Fail on geometry configuration.\n");
    return -5;
  }
  
  //Check errors
  if(geometry->configureStatus() != 0){
    if(verbose > 0)
      printf("createGeometry: Error: Fail on geometry configuration.\n");
    return -5;
  }

  return 0;
}

int createMaterials(pen_context& context,
		    std::string filenames[constants::MAXMAT],
		    const pen_parserSection& config,
		    const double globEmax,
		    const unsigned verbose){

  //Try to read global particle ranges for all particles
  std::array<double,constants::nParTypes> maxRanges{-1.0};
  pen_parserSection rangesSection;
  if(config.readSubsection("material-ranges",rangesSection) == INTDATA_SUCCESS){

    if(verbose > 1){
      printf("Global maximum ranges specified:\n");      
    }
    
    //Read maximum ranges for each particle type
    for(unsigned ip = 0; ip < constants::nParTypes; ++ip){
      double partRange;
      if(rangesSection.read(particleName(ip),partRange) == INTDATA_SUCCESS){
	maxRanges[ip] = partRange;
	if(verbose > 1){
	  printf("  + %15s: %15.3E cm\n", particleName(ip), maxRanges[ip]);
	}
      }
    }
      
  }

  std::array<double,constants::nParTypes> defaultEabs{1.0e3};
  pen_parserSection eabsSection;
  if(config.readSubsection("material-eabs",eabsSection) == INTDATA_SUCCESS){

    if(verbose > 1){
      printf("Global absorption energies specified:\n");      
    }
    
    //Read maximum ranges for each particle type
    for(unsigned ip = 0; ip < constants::nParTypes; ++ip){
      double partEabs;
      if(eabsSection.read(particleName(ip),partEabs) == INTDATA_SUCCESS){
	defaultEabs[ip] = partEabs;
	if(verbose > 1){
	  printf("  + %15s: %15.3E cm\n", particleName(ip), defaultEabs[ip]);
	}
      }
    }
      
  }

  
  //Extract materials section
  pen_parserSection matSection;
  if(config.readSubsection("materials",matSection) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("createMaterials: Error: Configuration 'materials' section doesn't exist.\n");
    }
    return -1;
  }

  //Get "materials" subsections (each one correspond to one material)
  std::vector<std::string> matNames;
  matSection.ls(matNames);

  if(verbose > 1){
    printf("Number of materials: %u\n", unsigned(matNames.size()));
    
    for(unsigned i = 0; i < matNames.size(); i++)
      printf("\t%s\n",matNames[i].c_str());
  }
  
  if(matNames.size() == 0){
    if(verbose > 0){
      printf("Error: No material found at configuration. Simulation requires, at last, one material.\n");
    }
    return -2;
  }
  
  if(matNames.size() > constants::MAXMAT){
    if(verbose > 0){
      printf("Error: %d materials found. The maximum number of materials is %d.\n",int(matNames.size()),constants::MAXMAT);
    }
    return -2;
  }
  
  //Set materials to context
  int errmat = context.setMats<pen_material>(matNames.size());
  if(errmat != 0)
    {
      printf("Error creating materials: %d.\n",errmat);
      return -7;
    }
  
  //Iterate over each material
  int err = 0;
  for(unsigned i = 0; i < matNames.size(); i++){

    printf("\n\n------------------------------------\n\n");
    printf(" **** Material '%s'\n\n",matNames[i].c_str());

    //Extract material subsection
    pen_parserSection oneMatSec;
    if(matSection.readSubsection(matNames[i],oneMatSec) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("createMaterials: Error: Unable to read mat section 'materials/%s'.\n",matNames[i].c_str());
      }
      err++;
      continue;
    }
  
    //Get material number
    int index;
    //Get absortion energies
    if(oneMatSec.read("number",index) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("createMaterials: Error: Unable to read 'materials/%s/number'. Integer expected.\n",matNames[i].c_str());
      }
      err++;
      continue;
    }

    if(index < 1 || index > int(constants::MAXMAT)){
      if(verbose > 0){
	printf("createMaterials: Error: Material index ('materials/%s/number') out of range (1-%d).\n",matNames[i].c_str(),constants::MAXMAT);
      }
      err++;
      continue;
    }

    //Calculate material index
    int j = index-1;

    //Read material filename
    //#######################
    if(oneMatSec.read("filename",filenames[j]) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("createMaterials: Error: Unable to read 'materials/%s/filename'. String expected.\n",matNames[j].c_str());
      }
      err++;
    }
    
    //Get material
    pen_material& mat = context.getBaseMaterial(j);

    //Set default C1, C2, WCC and WCR
    mat.C1=0.05;
    mat.C2=0.05;
    mat.WCC=1.0e3;
    mat.WCR=1.0e3;    
        
    //Get absortion energies
    //#######################

    std::array<double,constants::nParTypes> localRanges = maxRanges;
    pen_context* rangeContext = nullptr;
    for(unsigned ip = 0; ip < constants::nParTypes; ++ip){
      double partRange;
      double eabs;
      std::string key = std::string("range/")+particleName(ip);
      //Try to read local range for this particle
      if(oneMatSec.read(key,partRange) == INTDATA_SUCCESS){
	localRanges[ip] = partRange;
      }

      key = std::string("eabs/")+particleName(ip);
      //Assign absortion energy by range, explicit value or default value
      if(localRanges[ip] > 0.0){

	//Create range context if has not previously created
	if(rangeContext == nullptr)
	  rangeContext = createAuxContext(globEmax,
					  filenames[j].c_str(),
					  verbose);
	if(rangeContext == nullptr){
	  if(verbose > 0){
	    printf("createMaterials: Error: Unable to create "
		   "auxiliary context to calculate range limits.\n");
	  }
	  err++;
	  break;
	}
	  
	// + Range
	double topE = globEmax;
	double lowE = 50;
	double objectiveRange = localRanges[ip];
	unsigned nTries = 0;
	do{
	  double midE = (topE+lowE)/2.0;
	  double range = rangeContext->range(midE,static_cast<pen_KPAR>(ip),0);
	  if(range == objectiveRange){
	    topE = midE;
	    lowE = midE;
	  }
	  else if(range > objectiveRange){
	    topE = midE;
	  }else{
	    lowE = midE;
	  }
	  ++nTries;
	}while(nTries < 1000000 && topE/lowE > 1.001);
	if(topE == globEmax)
	  mat.EABS[ip] = 1.0e35;
	else
	  mat.EABS[ip] = lowE;
	
      }else if(oneMatSec.read(key,eabs) == INTDATA_SUCCESS){
	// + Explicit value
	mat.EABS[ip] = eabs;
      }else{
	// + Default value
	mat.EABS[ip] = defaultEabs[ip];
      }
    }

    //Free range context if has been created
    if(rangeContext != nullptr){
      delete rangeContext;
    }
    
    // ** WARNING: DEPRECATED **
    // ** This options is mantained for retrocompatibility
    double eabs;
    if(oneMatSec.read("eabs_e-",eabs) == INTDATA_SUCCESS){
      if(verbose > 0){
	printf("createMaterials: Warning: Parameter 'materials/%s/eabs_e-' is "
	       "deprecated and will be removed. Use "
	       "'materials/%s/eabs/electron' instead.\n",
	       matNames[j].c_str(),matNames[j].c_str());
      }
      mat.EABS[PEN_ELECTRON] = eabs;
    }
    if(oneMatSec.read("eabs_e+",eabs) == INTDATA_SUCCESS){
      if(verbose > 0){
	printf("createMaterials: Warning: Parameter 'materials/%s/eabs_e+' is "
	       "deprecated and will be removed. Use "
	       "'materials/%s/eabs/positron' instead.\n",
	       matNames[j].c_str(),matNames[j].c_str());
      }
      mat.EABS[PEN_POSITRON] = eabs;      
    }
    if(oneMatSec.read("eabs_gamma",eabs) == INTDATA_SUCCESS){
      if(verbose > 0){
	printf("createMaterials: Warning: Parameter 'materials/%s/eabs_gamma' is "
	       "deprecated and will be removed. Use "
	       "'materials/%s/eabs/gamma' instead.\n",
	       matNames[j].c_str(),matNames[j].c_str());
      }
      mat.EABS[PEN_PHOTON] = eabs;     
    }
    //**************************
    
    //Get C1 and C2
    //##################
    double C1,C2;
    if(oneMatSec.read("C1",C1) == INTDATA_SUCCESS){
      mat.C1 = C1;
    }
    if(oneMatSec.read("C2",C2) == INTDATA_SUCCESS){
      mat.C2 = C2;
    }

    //Read WCC and WCR
    //##################

    //Set default values
    mat.WCC = std::min(5e3,mat.EABS[PEN_ELECTRON]/100.0);
    mat.WCR = std::min(5e3,mat.EABS[PEN_PHOTON]/100.0);

    //Read values provided by the user
    double WCC, WCR;
    if(oneMatSec.read("WCC",WCC) == INTDATA_SUCCESS){
      mat.WCC = WCC;
    }
    if(oneMatSec.read("WCR",WCR) == INTDATA_SUCCESS){
      mat.WCR = WCR;
    }
    
    if(verbose > 0){
      printf("      C1 =%11.4E       C2 =%11.4E\n",mat.C1,mat.C2);
      printf("     WCC =%11.4E eV,   WCR =%11.4E eV\n",mat.WCC,(mat.WCR > 10.0E0 ? mat.WCR : 10.0E0));
      printf("  electron EABS: %11.4E eV\n",mat.EABS[PEN_ELECTRON]);
      printf("     gamma EABS: %11.4E eV\n",mat.EABS[PEN_PHOTON]);
      printf("  positron EABS: %11.4E eV\n",mat.EABS[PEN_POSITRON]);

      printf("\n Material filename: '%s'.\n",filenames[j].c_str());
      printf("\n Material number: %d.\n",index);      
    }
    
    printf("\n\n------------------------------------\n\n");
  }

  return err;
}

int setVarianceReduction(pen_context& context,
			 const wrapper_geometry& geometry,
			 const pen_parserSection& config,
			 pen_VRCluster<pen_particleState>& genericVR,
			 pen_VRCluster<pen_state_gPol>& photonVR,
			 const unsigned verbose){

  //Extract variance reduction section
  pen_parserSection VRSection;
  if(config.readSubsection("VR",VRSection) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No variance reduction specified.\n");
    }
    return 1;
  }

  int err;

  if(verbose > 1)printf("\n");
  
  //Get materials used by the current geometry
  bool usedMat[constants::MAXMAT+1];
  geometry.usedMat(usedMat);
  
  //*********************
  // Interaction Forcing
  //*********************
  
  //Try to Read interaction forcing parameters for materials and bodies.
  //Notice that body configuration will overwrite material one.

  // Materials
  //************
  
  std::vector<std::string> IFmats;
  VRSection.ls("IForcing/materials",IFmats);

  if(IFmats.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Material interaction forcing:\n\n");
      printf("  Mat |    Particle    | Interaction |  IF factor  | Weight Range\n");      
    }
    
    for(unsigned imname = 0; imname < IFmats.size(); imname++){

      //Get ibody name section
      pen_parserSection matSection;
      std::string matSecKey = std::string("IForcing/materials/") + IFmats[imname];
      if(VRSection.readSubsection(matSecKey,matSection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("setVarianceReduction: 'VR/%s' is not a section, skip this material.\n",matSecKey.c_str());
	}
	continue;
      }

      // Material index
      //***********************

      //Read index
      int imat;
      err = matSection.read("mat-index",imat);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Material index not specified for material section '%s'. Integer expected\n",IFmats[imname].c_str());
	}
      }

      //Check if material index is in range
      if(imat < 1 || imat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Specified index (%d) for material section '%s' out of range.\n",imat,IFmats[imname].c_str());
	}
	return -1;
      }

      //Check if material is used at current geometry
      if(!usedMat[imat]){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Specified index (%d) for material section '%s' is not used at current geometry.\n",imat,IFmats[imname].c_str());
	}
	return -2;
      }

      // IF Parameters
      //****************

      std::string particle;
      unsigned kpar;
      int icol;
      double weightL,weightU,forcer;
      
      //Read particle name
      err = matSection.read("particle",particle);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("setVarianceReduction: Missing field 'particle' for interaction forcing on material section '%s'. String expected.\n",IFmats[imname].c_str());
	return -3;
      }
      //Get kpar
      kpar = particleID(particle.c_str());
      if(kpar == ALWAYS_AT_END){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Unknown particle type specified on interaction forcing for material section '%s': '%s'\n",IFmats[imname].c_str(),particle.c_str());
	}
	return -4;
      }

      //Get interaction index
      err = matSection.read("interaction",icol);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Unable to read field 'interaction' for interaction forcing on material section '%s'. Integer expected.\n",IFmats[imname].c_str());
	return -5;
      }

      //Check interaction
      if(icol < 0 || icol >= (int)constants::MAXINTERACTIONS){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Interaction index (%d) out of range at interaction forcing on material section '%s'.\n",icol,IFmats[imname].c_str());
	return -6;
      }

      //Get forcing factor
      err = matSection.read("factor",forcer);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Unable to read field 'factor' for interaction forcing on material section '%s', interaction %d. Real number expected.\n",IFmats[imname].c_str(),icol);
	return -7;
      }

      //Read weight range
      err = matSection.read("min-weight",weightL);
      if(err != INTDATA_SUCCESS){
	weightL = 0.0;
	if(verbose > 2){
	  printf("field 'min-weight' no specified, real number expected.\n");
	  printf(" minimum weight will be set to 0.0\n");
	}
      }
      err = matSection.read("max-weight",weightU);
      if(err != INTDATA_SUCCESS){
	weightU = 1.0e6;
	if(verbose > 2){
	  printf("field 'max-weight' no specified, real number expected.\n");
	  printf(" maximum weight will be set to 1.0e6\n");
	}
      }

      if(weightL < 0.0 || weightU <= weightL){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Invalid weight range for interaction forcing on material '%s', interaction %d:\n",IFmats[imname].c_str(),icol);
	  printf("               minimum weight: %12.4E\n",weightL);
	  printf("               maximum weight: %12.4E\n",weightU);
	}
	return -8;	
      }

      //Print information
      if(verbose > 1){
	printf(" %4u  %15.15s   %11d   %11.4E   %12.4E - %12.4E\n",imat,particle.c_str(),icol,forcer,weightL,weightU);
      }

      //Set parameters for bodies with this material index
      for(unsigned ibody = 0; ibody < geometry.getBodies(); ibody++){
	
	if(geometry.getMat(ibody) != (unsigned)imat) continue;

	
	if(ibody >= context.NBV){
	  if(verbose > 0){
	    printf("setVarianceReduction: Error: Maximum body index for IF reached (%u)\n",context.NBV);
	  }
	  return -8;
	}


	//Set interaction forcing in context
	context.setForcing(forcer, static_cast<pen_KPAR>(kpar),
			   icol, ibody, weightL, weightU);
      }
    }
    
  }
  else if(verbose > 1){
    printf("No materials specified to use interaction forcing.\n");
  }
  

  // Bodies
  //************
  
  std::vector<std::string> IFbodies;
  VRSection.ls("IForcing/bodies",IFbodies);

  if(IFbodies.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Body interaction forcing:\n\n");
      printf(" Body |    Particle    | Interaction |  IF factor  | Weight Range\n");      
    }

    for(unsigned ibname = 0; ibname < IFbodies.size(); ibname++){

      //Get ibody name section
      pen_parserSection bodySection;
      std::string bodySecKey = std::string("IForcing/bodies/") + IFbodies[ibname];
      if(VRSection.readSubsection(bodySecKey,bodySection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("setVarianceReduction: 'VR/%s' is not a section, skip this body.\n",bodySecKey.c_str());
	}
	continue;
      }

      //Read body alias
      std::string bodyAlias;
      unsigned ibody;
      err = bodySection.read("body",bodyAlias);
      if(err == INTDATA_SUCCESS){
	//Check if specified body exists ang get the corresponding index
	ibody = geometry.getIBody(bodyAlias.c_str());
	if(ibody >= geometry.getBodies()){
	  if(verbose > 0){
	    printf("setVarianceReduction: Body '%s' doesn't exists in loaded geometry.\n",bodyAlias.c_str());
	  }
	  return -9;
	}
      }
      else{
	//Try to read body as integer index
	int auxibody;
	err = bodySection.read("body",auxibody);
	if(err != INTDATA_SUCCESS){
	  if(verbose > 0){
	    printf("setVarianceReduction: Error: Unable to read 'VR/%s/body'. Integer or string expected.\n",bodySecKey.c_str());
	  }
	  return -9;
	}
	if(auxibody < 0 || auxibody >= (int)context.NBV){
	  if(verbose > 0)
	    printf("setVarianceReduction: Error: Specified body index (%d) out of range.\n",auxibody);
	  return -9;
	}
	ibody = (unsigned)auxibody;
      }

      if(ibody >= context.NBV){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Maximum body index for IF is (%u)\n",context.NBV);
	  printf("                  specified index: %u\n",ibody);
	}
	return -9;	  
      }
      
      //Read kpar, icol and weight limits
      std::string particle;
      unsigned kpar;
      int icol;
      double weightL,weightU,forcer;
      
      //Read particle name
      err = bodySection.read("particle",particle);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("setVarianceReduction: Missing field 'particle' for interaction forcing on body section '%s'. String expected.\n",IFbodies[ibname].c_str());
	return -10;
      }
      //Get kpar
      kpar = particleID(particle.c_str());
      if(kpar == ALWAYS_AT_END){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Unknown particle type specified on interaction forcing for body section '%s': '%s'\n",IFbodies[ibname].c_str(),particle.c_str());
	}
	return -11;
      }

      //Get interaction index
      err = bodySection.read("interaction",icol);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Unable to read field 'interaction' for interaction forcing on body section '%s'. Integer expected.\n",IFbodies[ibname].c_str());
	return -12;
      }

      //Check interaction
      if(icol < 0 || icol >= (int)constants::MAXINTERACTIONS){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Interaction index (%d) out of range at interaction forcing on body section '%s'. Integer expected.\n",icol,IFbodies[ibname].c_str());
	return -13;
      }

      //Get forcing factor
      err = bodySection.read("factor",forcer);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Unable to read field 'factor' for interaction forcing on body section '%s', interaction %d. Real number expected.\n",IFbodies[ibname].c_str(),icol);
	return -14;
      }

      //Read weight range
      err = bodySection.read("min-weight",weightL);
      if(err != INTDATA_SUCCESS){
	weightL = 0.0;
	if(verbose > 2){
	  printf("field 'min-weight' no specified, real number expected.\n");
	  printf(" minimum weight will be set to 0.0\n");
	}
      }
      err = bodySection.read("max-weight",weightU);
      if(err != INTDATA_SUCCESS){
	weightU = 1.0;
	if(verbose > 2){
	  printf("field 'max-weight' no specified, real number expected.\n");
	  printf(" maximum weight will be set to 1.0\n");
	}
      }

      if(weightL < 0.0 || weightU <= weightL){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Invalid weight range for interaction forcing on body section '%s', interaction %d:\n",IFbodies[ibname].c_str(),icol);
	  printf("               minimum weight: %12.4E\n",weightL);
	  printf("               maximum weight: %12.4E\n",weightU);
	}
	return -15;	
      }

      //Print information
      if(verbose > 1){
	printf(" %4u  %15.15s   %11d   %11.4E   %12.4E - %12.4E\n",ibody,particle.c_str(),icol,forcer,weightL,weightU);
      }

      //Set interaction forcing in context
      context.setForcing(forcer, static_cast<pen_KPAR>(kpar),
			 icol, ibody, weightL, weightU);
    }
  }
  else if(verbose > 1){
    printf("No bodies specified to use interaction forcing.\n");
  }

  //Print final interaction forcing information
  if(verbose > 1){
    printf("\n\nEnabled Interaction forcing:\n\n");
    printf(" Body |    Particle    | Interaction |  IF factor  | Weight Range\n");
    for(unsigned i = 0; i < context.NBV; i++){
      for(unsigned j = 0; j < constants::nParTypes; j++){
	if(!context.LFORCE[i][j]){continue;}
	for(unsigned k = 0; k < constants::MAXINTERACTIONS; k++){
	  if(context.FORCE[i][j][k] > 1.0){
	    printf(" %4u  %15.15s   %11d   %11.4E   %12.4E - %12.4E\n",
		   i,particleName(j),k,context.FORCE[i][j][k],
		   context.WRANGES[i][2*j],context.WRANGES[i][2*j+1]);
	  }
	}
      }
    }
    printf("\n\n");
  }


  //**************************
  // Bremsstrahlung splitting
  //**************************

  // Materials
  //************

  std::vector<std::string> bremssMat;
  VRSection.ls("bremss/materials",bremssMat);

  if(bremssMat.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Material bremsstrahlung splitting:\n\n");
      printf("  Mat  | Bremss splitting\n");      
    }
    
    for(unsigned imname = 0; imname < bremssMat.size(); imname++){

      //Get ibody name section
      pen_parserSection matSection;
      std::string matSecKey = std::string("bremss/materials/") + bremssMat[imname];
      if(VRSection.readSubsection(matSecKey,matSection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("setVarianceReduction: 'VR/%s' is not a section, skip this material.\n",matSecKey.c_str());
	}
	continue;
      }

      // Material index
      //***********************

      //Read index
      int imat;
      err = matSection.read("mat-index",imat);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Material index not specified for material '%s'. Integer expected\n",bremssMat[imname].c_str());
	}
      }

      //Check if material index is in range
      if(imat < 1 || imat > (int)constants::MAXMAT){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Specified index (%d) for material '%s' out of range.\n",imat,bremssMat[imname].c_str());
	}
	return -16;
      }

      //Check if material is used at current geometry
      if(!usedMat[imat]){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Specified index (%d) for material '%s' is not used at current geometry.\n",imat,bremssMat[imname].c_str());
	}
	return -17;
      }

      //Get splitting factor
      int splitting;
      err = matSection.read("splitting",splitting);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Unable to read field 'splitting' for bremsstrahlung splitting on material '%s'. Integer expected.\n",bremssMat[imname].c_str());
	return -18;
      }

      //Check splitting factor
      if(splitting < 1){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Invalid bremsstrahlung splitting factor (%d).\n",splitting);
	return -19;
      }

      //Set splitting factor for bodies with this material index
      for(unsigned ibody = 0; ibody < geometry.getBodies(); ibody++){
	
	if(geometry.getMat(ibody) != (unsigned)imat) continue;

	if(ibody >= context.NBV){
	  if(verbose > 0){
	    printf("setVarianceReduction: Error: Maximum body index for IF reached (%u)\n",context.NBV);
	  }
	  return -19;
	}
	
	if(context.LFORCE[ibody][PEN_ELECTRON] || context.LFORCE[ibody][PEN_POSITRON]){
	  context.IBRSPL[ibody] = (unsigned)splitting;	  
	}
      }

      //Print configuration
      if(verbose > 1){
	printf(" %5d   %5d\n", imat,splitting);
      }
      
    }
  }
  else if(verbose > 1){
    printf("No material with bremsstrahlung splitting enabled.\n");
  }
  
  // Bodies
  //************

  std::vector<std::string> bremssBodies;
  VRSection.ls("bremss/bodies",bremssBodies);

  if(bremssBodies.size() > 0){
    if(verbose > 1){
      printf("\n\n **** Body bremsstrahlung splitting:\n\n");
      printf(" Body  | Bremss splitting\n");      
    }
    
    for(unsigned ibname = 0; ibname < bremssBodies.size(); ibname++){

      //Get ibody name section
      pen_parserSection bodySection;
      std::string bodySecKey = std::string("bremss/bodies/") + bremssBodies[ibname];
      if(VRSection.readSubsection(bodySecKey,bodySection) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("setVarianceReduction: 'VR/%s' is not a section, skip this body.\n",bodySecKey.c_str());
	}
	continue;
      }

      // Body index
      //***********************

      //Check if specified body exists
      unsigned ibody = geometry.getIBody(bremssBodies[ibname].c_str());
      if(ibody >= geometry.getBodies()){
	if(verbose > 0){
	  printf("setVarianceReduction: Body '%s' doesn't exists in loaded geometry.\n",bremssBodies[ibname].c_str());
	}
	return -20;
      }
      else if(ibody >= context.NBV){
	if(verbose > 0){
	  printf("setVarianceReduction: Error: Maximum body index for IF is (%u)\n",context.NBV);
	  printf("                  specified index: %u\n",ibody);
	}
	return -20;	  
      }

      //Check if in this body IF has been enabled
      if(!context.LFORCE[ibody][PEN_ELECTRON] && !context.LFORCE[ibody][PEN_POSITRON]){
	if(verbose > 0){
	  printf("setVarianceReduction: Body '%s' (index %u) has not enabled interaction forcing.\n",bremssBodies[ibname].c_str(),ibody);
	}
	return -21;
      }

      // Splitting factor
      //***********************

      //Get splitting factor
      int splitting;
      err = bodySection.read("splitting",splitting);
      if(err != INTDATA_SUCCESS){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Unable to read field 'splitting' for bremsstrahlung splitting on body '%s'. Integer expected.\n",bremssBodies[ibname].c_str());
	return -22;
      }

      //Check splitting factor
      if(splitting < 1){
	if(verbose > 0)
	  printf("setVarianceReduction: Error: Invalid bremsstrahlung splitting factor (%d).\n",splitting);
	return -23;
      }

      //Set splitting factor for specified body
      context.IBRSPL[ibody] = (unsigned)splitting;	  

      //Print configuration
      if(verbose > 1){
	printf(" %5d   %5d\n", ibody,splitting);
      }
      
    }
  }
  else if(verbose > 1){
    printf("No bodies with bremsstrahlung splitting enabled.\n");
  }

  if(verbose > 1){
    printf("\n\nFinal bremsstrahlung splitting:\n\n");
    printf(" Body  | Bremss splitting\n");      
    for(unsigned ibody = 0; ibody < context.NBV; ibody++){
      
      if(context.IBRSPL[ibody] > 1)
	printf(" %5u   %5u\n", ibody,context.IBRSPL[ibody]);
    }
    printf("\n\n");
  }

  //**************************
  // Other VR techniques
  //**************************

  pen_parserSection VRgeneric;
  if(VRSection.readSubsection("generic",VRgeneric) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No generic variance reduction specified.\n");
    }
  }
  else{
    genericVR.configure(VRgeneric,geometry,verbose);
    genericVR.name.assign("Generic-VR");
    if(genericVR.configureStatus() != 0){
      printf("setVarianceReduction: Error: Unable to configure "
	     "generic variance reduction.");
      return -24;
    }
  }

  pen_parserSection VRphoton;
  if(VRSection.readSubsection("photon",VRphoton) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No photon based variance reduction specified.\n");
    }
  }else{
    photonVR.name.assign("Photon-VR");
    photonVR.configure(VRphoton,geometry,verbose);
    if(photonVR.configureStatus() != 0){
      printf("setVarianceReduction: Error: Unable to configure "
	     "photon specific variance reduction.");
      return -25;
    }    
  }
  //pen_VRxraysplitting dummy;
  //printf("\n*****************************");
  //printf("\n\nRegister return: %d\n\n",dummy.registerStatus());
  //printf("*****************************\n");
  
  return 0;
}

pen_context* createAuxContext(double EMAX,
			      const char* matFilename,
			      const unsigned verbose){

  //Create elements data base
  pen_elementDataBase elements;
  //Create a context
  pen_context* context = nullptr;
  context = new pen_context(elements);
  if(context == nullptr){
    if(verbose > 0){
      printf("createAuxContext: Error allocating auxilary context.\n");      
    }
    return nullptr;
  }

  //Set the number of materials to context (1)
  int errmat = context->setMats<pen_material>(1);
  if(errmat != 0){
    if(verbose > 0){
      printf("createAuxContext: Error at context "
	     "material creation: %d.\n",errmat);
    }
    delete context;
    return nullptr;
  }
  
  //Get the material
  pen_material& mat = context->getBaseMaterial(0);

  //Configure the material
  mat.C1=0.2;
  mat.C2=0.2;
  mat.WCC=1.0e3;
  mat.WCR=1.0e3;

  mat.EABS[PEN_ELECTRON] = 50.0E0;
  mat.EABS[PEN_PHOTON]   = 50.0E0;
  mat.EABS[PEN_POSITRON] = 50.0E0;

  FILE* fcontext = nullptr;
  
  if(verbose > 0){
    fcontext = fopen("rangeContext.rep","w");
    if(fcontext == nullptr){
      printf("createAuxContext: Error: unable to create "
	     "file 'rangeContext.rep'\n");
      delete context;
      return nullptr;
    }
  }
  
  int INFO = 1;
  std::string PMFILEstr[constants::MAXMAT];
  PMFILEstr[0].assign(matFilename);
  int err = context->init(EMAX,nullptr,INFO,PMFILEstr);
  if(err != 0){
    if(verbose > 0){
      printf("createAuxContext: Error: Unable to configure range context."
	     "More details can be found in 'rangeContext.rep' file.\n");
    }
    delete context;
    return nullptr;
  }

  return context;
  
}
