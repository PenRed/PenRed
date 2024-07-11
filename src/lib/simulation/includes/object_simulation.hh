
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
//        sanolgi@upvnet.upv.es
//    
//

 
#ifndef __PEN_OBJECT_SIMULATION__
#define __PEN_OBJECT_SIMULATION__

#include "base_simulation_functions.hh"
#include <ostream>
#include <iostream>
#include <memory>
#include <thread>

namespace penred{

  namespace simulation{

    template<class contextType, class stateType>
    inline void _simulate(penred::simulation::simConfig& config,
			  const contextType& context,
			  pen_specificStateGen<stateType>& source,
			  pen_commonTallyCluster& tallies,
			  const pen_VRCluster<pen_particleState>& genericVR,
			  const pen_VRCluster<pen_state_gPol>& photonVR){

      //Perform the simulation
      penred::simulation::sampleAndSimulateContext(config, context, source, tallies,
						   genericVR, photonVR);
  
    }    
  
    template <class contextType>
    class simulator{

    private:

      std::array<bool,pen_imageExporter::nFormats()> enabledFormats;
      unsigned nThreads;
      int nSeedPair;
      bool ASCIIResults, finalDump;

      pen_parserSection contextConfig;
      pen_parserSection particleSourcesConfig;
      pen_parserSection geometryConfig;
      pen_parserSection talliesConfig;
      pen_parserSection VRConfig;

      template<class stateType>
      int configureSource(pen_specificStateGen<stateType>& source,
			  const pen_parserSection& config,
			  std::ostream& out,
			  const unsigned verbose){
  
	//Check if specific sampler exists
	bool useSpecific = config.isSection("specific");
  
	//Get source kpar
	std::string auxkpar;
	int err = config.read("kpar",auxkpar);
	if(err != INTDATA_SUCCESS && !useSpecific){
	  if(verbose > 0 && out){
	    out << "configureSource: Error: Particle type field ('kpar') not found "
	      "and no specific sampling specified. String expected:"<< std::endl;
	    for(unsigned i = 0; i < constants::nParTypes; ++i){
	      out << "\t" << particleName(i) << std::endl;
	    }
	  }
	  return -1;
	} else if(!useSpecific && verbose > 1 && out) {
	  out << "Selected particle: " << auxkpar << std::endl;
	}

	//Get particle ID
	unsigned id = particleID(auxkpar.c_str());
	if(!isKpar(id)){
	  if(verbose > 0 && out)
	    out << "\nInvalid selected particle type: " << auxkpar << std::endl;
	  return -2;
	}
	source.kpar = static_cast<pen_KPAR>(id);

	//Configure source
	source.configure(config,nThreads,verbose);
	if(source.configureStatus() != 0){
	  if(verbose > 0 && out){
	    out << "\nError on source configuration" << std::endl;
	  }
	  return -3;
	}
  
	return 0;
      }

      int createParticleGenerators(std::vector<pen_specificStateGen<pen_particleState>>& genericSources,
				   std::vector<pen_specificStateGen<pen_state_gPol>>& polarisedGammaSources,
				   std::ostream& out,
				   const unsigned verbose){

	if(particleSourcesConfig.empty()){
	  if(verbose > 0 && out){
	    out << "Error: Missing particle sources configuration." << std::endl;
	  }
	  return errors::ERROR_MISSING_SOURCE_CONFIGURATION;
	}
	
	int errG = 0;
	int errP = 0;
	int err = 0;
  
	//Extract source sections
	pen_parserSection genericSourceSection;
	pen_parserSection polarisedSourceSection;
	errG = particleSourcesConfig.readSubsection("generic", genericSourceSection);
	errP = particleSourcesConfig.readSubsection("polarized", polarisedSourceSection);
	if(errG != INTDATA_SUCCESS && errP != INTDATA_SUCCESS){
	  if(verbose > 0 && out){
	    out << "createParticleGenerators: Error: Fields 'generic' nor "
	      "'polarized' don't exists. No source specified" << std::endl;;
	  }
	  return errors::ERROR_NO_SOURCE;
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
	  if(verbose > 0 && out){
	    out << "createParticleGenerators: Error: Simulation requires, "
	      "at last, one particle source." << std::endl;
	  }
	  return errors::ERROR_NO_SOURCE;
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
	    if(verbose > 0 && out){
	      out << "createParticleGenerators: Error: Unable to extract "
		"section for generic source '" << genericSourceNames[i] << "'" << std::endl;
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
	      if(verbose > 0 && out){
		out << "createParticleGenerators: Error: Unable to read field 'nhist' "
		  "for generic source '" << genericSourceNames[i] << "'" << std::endl;
	      }
	      err++;
	      continue;
	    }

	    if(nhists <= 0.5){
	      if(verbose > 0 && out){
		out << "createParticleGenerators: Error on generic source "
		    << genericSourceNames[i] << ". "
		  "Number of histories must be greater than zero" << std::endl;
	      }
	      err++;
	      continue;	
	    }
      
	    //Load source
	    if(configureSource(genericSources[i],
			       genSection,out,verbose) != 0){
	      if(verbose > 0 && out){
		out << "createParticleGenerators: Error: Can't create and "
		  "configure source '" << genericSourceNames[i] << "'." << std::endl;
	      }
	      err++;
	    }

	    //Init the source task
	    unsigned long long uihists = static_cast<unsigned long long>(nhists+0.1);
	    uihists = std::max(uihists,1llu);

	    // ******************************* LB ************************************ //
#if (defined _PEN_USE_MPI_ && defined _PEN_USE_LB_ )
	    int errTask = genericSources[i].initTask(nThreads,uihists,MPI_COMM_WORLD,
						     nextTag,nextTag+1,verbose);
	    nextTag += 2;
#else
	    int errTask = genericSources[i].initTask(nThreads,uihists,verbose);      
#endif
	    // ***************************** LB END ********************************** //
	    if(errTask != 0){
	      if(verbose > 0 && out){
		out << "createParticleGenerators: Error on generic source "
		    << genericSourceNames[i] << ". Unable to init source task\n"
		  "   Error code: " << errTask << std::endl;
	      }
	      err++;
	      continue;
	    }

	    if(verbose > 1 && out){
	      out << "Histories to simulate at source "
		  << genericSources[i].name << ": "
		  << genericSources[i].toDo() << std::endl;
	    }
	  }
	}

	//Iterate over all specific sources
	for(unsigned i = 0; i < nPolarizedSources; i++){
	  //Get source section
	  pen_parserSection genSection;
	  if(polarisedSourceSection.readSubsection(polarisezSourceNames[i],genSection) != INTDATA_SUCCESS){
	    if(verbose > 0 && out){
	      out << "createParticleGenerators: Error: Unable to extract section "
		"for polarized gamma source '"
		  << polarisezSourceNames[i] << "'" << std::endl;
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
	      if(verbose > 0 && out){
		out << "createParticleGenerators: Error: Unable to read field "
		  "'nhist' for polarized gamma source '" << polarisezSourceNames[i]
		    << "'" << std::endl;
	      }
	      err++;
	      continue;
	    }

	    if(nhists <= 0.5){
	      if(verbose > 0 && out){
		out << "createParticleGenerators: Error on polarized gamma source "
		    << polarisezSourceNames[i]
		    << ". Number of histories must be greater than zero" << std::endl;
	      }
	      err++;
	      continue;	
	    }
      
	    //Load source
	    if(configureSource(polarisedGammaSources[i],
			       genSection, out, verbose) != 0){
	      if(verbose > 0 && out){
		out << "createParticleGenerators: Error: Can't "
		  "create and configure source '"
		    << polarisezSourceNames[i] << "'." << std::endl;
	      }
	      err++;
	    }

	    //Init the source task
	    unsigned long long uihists = static_cast<unsigned long long>(nhists);
	    uihists = std::max(uihists,1llu);

      
#if defined(_PEN_USE_MPI_) && defined(_PEN_USE_LB_)
	    int errTask = polarisedGammaSources[i].initTask(nThreads,uihists,MPI_COMM_WORLD,
							    nextTag,nextTag+1,verbose);
	    nextTag += 2;
#else
	    int errTask = polarisedGammaSources[i].initTask(nThreads,uihists,verbose);      
#endif
	    if(errTask != 0){
	      if(verbose > 0 && out){
		out << "createParticleGenerators: Error on polarised source "
		    << polarisedGammaSources[i].name << ". Unable to init source task\n"
		  "   Error code: " << errTask << std::endl;
	      }
	      err++;
	      continue;
	    }

	    if(verbose > 1 && out)
	      out << "Histories to simulate at source "
		  << polarisedGammaSources[i].name << ": "
		  << polarisedGammaSources[i].toDo() << std::endl;
	  }
	}
      
	return err;
      }

      int createGeometry(std::shared_ptr<wrapper_geometry>& geometry,
			 const pen_parserSection& matInfo,
			 std::ostream& out,
			 const unsigned verbose){

	if(geometryConfig.empty()){
	  if(verbose > 0 && out){
	    out << "Error: Missing geometry configuration." << std::endl;
	  }
	  return errors::ERROR_MISSING_GEOMETRY_CONFIGURATION;
	}
	
	//Get geometry section
	pen_parserSection geometrySection = geometryConfig;

	//Append material information to geometry section
	geometrySection.addSubsection("materials", matInfo);
  
	//Get geometry type
	std::string geoType;
	if(geometrySection.read("type",geoType) != INTDATA_SUCCESS){
	  if(verbose > 0 && out){
	    out << "createGeometry: Error: field 'type' not "
	      "specified. String expected.\n" << std::endl;
	  }
	  return errors::ERROR_MISSING_TYPE;
	}

	//Create geometry
	geometry.reset(penGeoRegister_create(geoType.c_str()));
	if(geometry == nullptr){
	  if(verbose > 0 && out){
	    out << "createGeometry: Error creating geometry instance of type '"
		<< geoType <<  "'. Unknown type." << std::endl;
	  }
	  return errors::ERROR_UNKNOWN_TYPE;
	}

	//Configure geometry  
	geometry->name.assign("geometry");    
	if(geometry->configure(geometrySection,verbose) != 0){
	  if(verbose > 0 && out){
	    out << "createGeometry: Error: Geometry configuration failed." << std::endl;
	  }
	  return errors::ERROR_AT_GEOMETRY_CONFIGURATION;
	}
  
	//Check errors
	if(geometry->configureStatus() != 0){
	  if(verbose > 0 && out){
	    out << "createGeometry: Error: Geometry configuration failed." << std::endl;
	  }
	  return errors::ERROR_AT_GEOMETRY_CONFIGURATION;
	}

	return errors::SUCCESS;
      }


      int createTallies(std::vector<pen_commonTallyCluster>& tallyGroups,
			const contextType& context,
			std::ostream& out,
			const unsigned verbose){

	if(talliesConfig.empty()){
	  if(verbose > 0 && out){
	    out << "Error: Missing tally configuration." << std::endl;
	  }	  
	  return errors::ERROR_MISSING_TALLY_CONFIGURATION;
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
	for(unsigned j = 1; j < nThreads; j++){
	  tallyGroups[j].name.assign("common");
	  //Configure tallies
	  threads.push_back(tallyGroups[j].configure_async(context.readGeometry(),
							   mats,
							   j,
							   talliesConfig,
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
				 talliesConfig,
				 verbose);
  
	// ************************** MULTI-THREADING ***************************** //
#ifdef _PEN_USE_THREADS_
	//Sincronize all threads
	for(unsigned j = 0; j < nThreads-1; j++){
	  threads[j].join();
	}
#endif
	// ************************** MULTI-THREADING ***************************** //
  
	//Check errors
	unsigned failedClusters = 0;
	for(unsigned j = 0; j < nThreads; j++){
	  int err = tallyGroups[j].configureStatus();
	  if(err != 0){
	    if(verbose > 0 && out){
	      out << "createTallies: Error on tally cluster " << j <<
		" creation and configuration (err code " << err << ")." << std::endl;
	    }
	    failedClusters++;
	  }
	}

	if(failedClusters > 0)
	  return errors::ERROR_ON_TALLIES_CONFIGURATION;
  
	//Share configuration from thread 0 cluster to other threads
	failedClusters = 0;
	for(unsigned j = 0; j < nThreads; j++){

	  int err = tallyGroups[j].shareConfig(tallyGroups[0], verbose);
	  if(err != 0){
	    if(verbose > 0 && out)
	      out << "createTallies: Error on tally cluster " << j <<
		". Unable to get configuration from thread 0 (err code "
		  << err << ")." << std::endl;
	    failedClusters++;
	  }
	}

	if(failedClusters > 0)
	  return errors::ERROR_ON_TALLIES_CONFIGURATION;  
  
	return errors::SUCCESS;
      }


      int setVarianceReduction(const std::shared_ptr<wrapper_geometry>& geometry,
			       pen_VRCluster<pen_particleState>& genericVR,
			       pen_VRCluster<pen_state_gPol>& photonVR,
			       std::ostream& out,
			       const unsigned verbose){

	if(VRConfig.empty()){
	  if(verbose > 0 && out){
	    out << "No variance reduction specified." << std::endl;
	  }
	  return errors::SUCCESS;
	}

	if(verbose > 1 && out)
	  out << std::endl;
  
	pen_parserSection VRgeneric;
	if(VRConfig.readSubsection("generic",VRgeneric) != INTDATA_SUCCESS){
	  if(verbose > 1 && out){
	    out << "No generic variance reduction specified." << std::endl;
	  }
	}
	else{
	  genericVR.configure(VRgeneric,*geometry,verbose);
	  genericVR.name.assign("Generic-VR");
	  if(genericVR.configureStatus() != 0){
	    if(verbose > 0 && out){
	      out << "setVarianceReduction: Error: Unable to configure "
		"generic variance reduction." << std::endl;
	    }
	    return errors::ERROR_ON_VR_CONFIGURATION;
	  }
	}

	pen_parserSection VRphoton;
	if(VRConfig.readSubsection("photon",VRphoton) != INTDATA_SUCCESS){
	  if(verbose > 1 && out){
	    out << "No photon based variance reduction specified." << std::endl;
	  }
	}else{
	  photonVR.name.assign("Photon-VR");
	  photonVR.configure(VRphoton,*geometry,verbose);
	  if(photonVR.configureStatus() != 0){
	    if(verbose > 0 && out){
	      out << "setVarianceReduction: Error: Unable to configure "
		"photon specific variance reduction." << std::endl;
	    }
	    return errors::ERROR_ON_VR_CONFIGURATION;
	  }
	}
  
	return errors::SUCCESS;
      }
      
    public:

      simConfig baseSimConfig;

      static constexpr const char* defaultSimPath      = "simulation";
      static constexpr const char* defaultSourcePath   = "sources";
      static constexpr const char* defaultGeometryPath = "geometry";
      static constexpr const char* defaultTallyPath    = "tallies";
      static constexpr const char* defaultVRPath       = "VR";

      static std::string versionMessage(){

	return std::string("***************************************************************\n"
			   " PenRed version: 1.11.0b (11-July-2024) \n"
			   " Copyright (c) 2019-2024 Universitat Politecnica de Valencia\n"
			   " Copyright (c) 2019-2024 Universitat de Valencia\n"
			   " Reference: Computer Physics Communications, 267 (2021) 108065\n"
			   "            https://doi.org/10.1016/j.cpc.2021.108065\n"
			   " This is free software; see the source for copying conditions.\n"
			   " There is NO warranty; not even for MERCHANTABILITY or\n"
			   " FITNESS FOR A PARTICULAR PURPOSE.\n"
			   " Please, report bugs and suggestions at our github repository\n"
			   "         https://github.com/PenRed/PenRed\n"
			   "***************************************************************\n");    
      }
      
      simulator() : nThreads(1),
		    nSeedPair(-1),
		    ASCIIResults(true),
		    finalDump(false)
      {
	std::fill(enabledFormats.begin(), enabledFormats.end(), false);
      }
      

      // Set functions
      inline int enableFormat(const std::string& format){

	//Check if it is a known format
	if(!pen_imageExporter::isFormat(format.c_str()))
	  return -1;

	//Enable it
	enabledFormats[static_cast<int>(pen_imageExporter::toFormat(format.c_str()))] = true;

	return errors::SUCCESS;
      }
      inline void setThreads(const unsigned nThreadsIn){
	if(nThreadsIn == 0){
	  unsigned int nConcurrency = std::thread::hardware_concurrency();
	  if(nConcurrency > 0)
	    nThreads = nConcurrency;
	  else
	    nThreads = 1;
	}
	else{
	  nThreads = nThreadsIn;
	}
      }
      inline int setSeedPair(const int seedPairIn){
	if(seedPairIn > 1000){
	  return errors::ERROR_INVALID_SEED_PAIR;
	}
	nSeedPair = seedPairIn;
	return errors::SUCCESS;
      }
      
      inline void enableASCII(const bool enable){ ASCIIResults = enable; }
      inline void enableFinalDump(const bool enable){ finalDump = enable; }

      //Set configuration functions
      inline int setContextConfig(const pen_parserSection& config,
				  const std::string& prefix = ""){
	if(prefix.empty())
	  contextConfig = config;
	else{
	  if(config.readSubsection(prefix.c_str(), contextConfig) != INTDATA_SUCCESS)
	    return errors::ERROR_MISSING_PATH;
	}
	return errors::SUCCESS;
      }
      inline int setSourceConfig(const pen_parserSection& config,
				 const std::string& prefix = ""){
	if(prefix.empty())
	  particleSourcesConfig = config;
	else{
	  if(config.readSubsection(prefix.c_str(), particleSourcesConfig) != INTDATA_SUCCESS)
	    return errors::ERROR_MISSING_PATH;
	}
	return errors::SUCCESS;
      }
      inline int setGeometryConfig(const pen_parserSection& config,
				   const std::string& prefix = ""){
	if(prefix.empty())
	  geometryConfig = config;
	else{
	  if(config.readSubsection(prefix.c_str(), geometryConfig) != INTDATA_SUCCESS)
	    return errors::ERROR_MISSING_PATH;
	}
	return errors::SUCCESS;
      }
      inline int setTallyConfig(const pen_parserSection& config,
				const std::string& prefix = ""){
	if(prefix.empty())
	  talliesConfig = config;
	else{
	  if(config.readSubsection(prefix.c_str(), talliesConfig) != INTDATA_SUCCESS)
	    return errors::ERROR_MISSING_PATH;
	}
	return errors::SUCCESS;
      }
      inline int setVRConfig(const pen_parserSection& config,
			     const std::string& prefix = ""){
	if(prefix.empty())
	  VRConfig = config;
	else{
	  if(config.readSubsection(prefix.c_str(), VRConfig) != INTDATA_SUCCESS)
	    return errors::ERROR_MISSING_PATH;
	}
	return errors::SUCCESS;
      }      

      // Get functions
      inline simConfig& getSimConfig() { return baseSimConfig; }
      
      //Simulation configuration function
      int configureSimConfig(const pen_parserSection& config,
			     const std::string& prefix = "",
			     const unsigned verbose = 1){
	int errSimConfig = baseSimConfig.configure(prefix.c_str(),config);
	if(errSimConfig != penred::simulation::errors::SUCCESS){
	  if(verbose > 1)
	    printf("simulator: configureSimConfig: Error: Unable to parse "
		   "simulation configuration section: Invalid seeds\n");
	  return errors::ERROR_INVALID_SEEDS;
	}

	return errors::SUCCESS;
      }

      int configure(const pen_parserSection& config,
		    std::ostream& out = std::cout,
		    const std::string& prefixSimConfig = defaultSimPath,
		    const std::string& prefixSourceConfig = defaultSourcePath,
		    const std::string& prefixGeometryConfig = defaultGeometryPath,
		    const std::string& prefixTallyConfig = defaultTallyPath,
		    const std::string& prefixVRConfig = defaultVRPath){


	// ** Simulation configuration
	if(configureSimConfig(config,prefixSimConfig) != errors::SUCCESS){
	  return -1;
	}

	unsigned verbose = baseSimConfig.verbose;

	if(verbose > 1 && out){
	  out << baseSimConfig.stringifyState() << std::endl;
	}

	// ** Number of threads

	int auxThreads;
	std::string path = prefixSimConfig + "/threads";
	if(config.read(path,auxThreads) != INTDATA_SUCCESS){
	  if(verbose > 1){
	    out << "\n\nNumber of threads not specified, "
	      "only one thread will be used.\n\n" << std::endl;
	  }
	  auxThreads = 1;
	}
	// ************************** MULTI-THREADING ***************************** //
#ifndef _PEN_USE_THREADS_
	else{
	  if(verbose > 1){
	    out << "\n\nMulti-threading has not been activated during compilation"
	      ", only one thread will be used\n\n" << std::endl;
	  }
	}
	auxThreads = 1;
#else
	if(auxThreads <= 0){

	  if(verbose > 1)
	    out << "Automatic selection of threads enabled\n" << std::endl;
	  
	  unsigned int nConcurrency = std::thread::hardware_concurrency();

	  if(nConcurrency > 0)
	    auxThreads = nConcurrency;
	  else{
	    out << "The hardware concurrency value is not well defined "
	      "or not computable. One thread will be used" << std::endl;
	    auxThreads = 1;
	  }
	}
  
#endif
	// ************************ MULTI-THREADING END *************************** //

	//Set threads
	setThreads(static_cast<unsigned>(auxThreads));
	if(verbose > 1){
	  out << "\nNumber of simulating threads: " << nThreads << std::endl;
	  if(nThreads > 1){
	    out << "Initial random seeds will be selected using "
	      "\"rand0\" function to ensure truly independent sequences "
	      "of random numbers.\n" << std::endl;
	  }
	}

	// Get initial seed pair number
	//*******************************
	int nseedPairAux;
	path = prefixSimConfig + "/seedPair"; 
	if(config.read(path,nseedPairAux) == INTDATA_SUCCESS){
	  if(setSeedPair(nseedPairAux) != errors::SUCCESS){
	    if(verbose > 0){
	      out << "Invalid initial seed pair number " << nseedPairAux << std::endl;
	      out << "Available seed pair range is [0,1000]" << std::endl;
	      out << "Seed pair will be unchanged. " << std::endl;
	    }
	  }
	  else if(verbose > 1){
	    out << "Selected rand0 seed pair number: " << nSeedPair << std::endl;
	  }
	}else{
	  nSeedPair = -1;
	}

	// Get output options
	//*******************************
	bool ASCIIResultsAux = true;
	path = prefixSimConfig + "/ascii-results"; 
	if(config.read(path,ASCIIResultsAux) == INTDATA_SUCCESS){
	  enableASCII(ASCIIResultsAux);
	}

	bool finalDumpAux = false;
	path = prefixSimConfig + "/finalDump"; 
	if(config.read("simulation/finalDump",finalDumpAux) != INTDATA_SUCCESS){
	  enableFinalDump(finalDumpAux);
	}

	// Save context config
	//*******************************
	contextConfig = config;
	
	// Get source configuration
	//*******************************
	if(config.readSubsection(prefixSourceConfig.c_str(), particleSourcesConfig) != INTDATA_SUCCESS){
	  if(verbose > 0 && out){
	    out << "Error: Missing source configuration section '"
		<< prefixSourceConfig << "'" << std::endl;
	  }
	  return errors::ERROR_MISSING_SOURCE_CONFIGURATION;
	}

	// Get geometry configuration
	//*******************************
	if(config.readSubsection(prefixGeometryConfig.c_str(), geometryConfig) != INTDATA_SUCCESS){
	  if(verbose > 0 && out){
	    out << "Error: Missing geometry configuration section '"
		<< prefixGeometryConfig << "'" << std::endl;
	  }
	  return errors::ERROR_MISSING_GEOMETRY_CONFIGURATION;
	}

	// Get tally configuration
	//*******************************
	if(config.readSubsection(prefixTallyConfig.c_str(), talliesConfig) != INTDATA_SUCCESS){
	  if(verbose > 0 && out){
	    out << "Error: Missing tallies configuration section '"
		<< prefixTallyConfig << "'" << std::endl;
	  }
	  return errors::ERROR_MISSING_TALLY_CONFIGURATION;
	}

	// Get VR configuration
	//*******************************
	config.readSubsection(prefixVRConfig.c_str(), VRConfig);

	return errors::SUCCESS;
      }

      inline int configFromFile(const std::string& filename,
				std::ostream& out = std::cout,
				const std::string& prefixSimConfig = defaultSimPath,
				const std::string& prefixSourceConfig = defaultSourcePath,
				const std::string& prefixGeometryConfig = defaultGeometryPath,
				const std::string& prefixTallyConfig = defaultTallyPath,
				const std::string& prefixVRConfig = defaultVRPath){

	//Parse configuration file
	pen_parserSection config;
	std::string errorLine;
	unsigned long errorLineNum;
	int err = parseFile(filename.c_str(),config,errorLine,errorLineNum);
  
	if(err != INTDATA_SUCCESS){
	  out << "Error parsing configuration.\n"
	    "Error code: " << err <<  "\n"
	    "Error message: " << pen_parserError(err) << "\n"
	    "Error located at line " << errorLineNum
	      << ", at text: " << errorLine << std::endl;
	  return errors::ERROR_PARSING_CONFIG;
	}

	return configure(config, out,
			 prefixSimConfig,
			 prefixSourceConfig,
			 prefixGeometryConfig,
			 prefixTallyConfig,
			 prefixVRConfig);
      }
      
      int simulate(std::ostream& out = std::cout){

	//Create a timer to measure the expended time in initialization 
	pen_timer initializationTimer;

	// Get verbose level
	const unsigned verbose = baseSimConfig.verbose;
	
	// Create a simulation config for each thread
	//*********************************************
	std::vector<penred::simulation::simConfig> simConfigs(nThreads);

	// ** Copy basic configuration
	for(unsigned i = 0; i < nThreads; i++){
	  simConfigs[i].iThread = i;
	  simConfigs[i].copyCommonConfig(baseSimConfig);
	  //Set, by default, std::cout as output stream
	  simConfigs[i].setOutstream(out);
	}

	// ** Set random seeds
	constexpr size_t nRand0Seeds = 1001;

	if(nThreads > 1){
	  for(unsigned i = 0; i < nThreads; i++){
	    int seedPos = i;
	    if(nSeedPair >= 0)
	      seedPos = (nSeedPair+i) % nRand0Seeds;
	    simConfigs[i].setSeeds(seedPos);
	  }
	}
	else{    
	  if(nSeedPair >= 0){
	    simConfigs[0].setSeeds(nSeedPair);
	  }
	}

	//****************************
	// Create particle generators
	//****************************
	std::vector<pen_specificStateGen<pen_particleState>> genericSources;
	std::vector<pen_specificStateGen<pen_state_gPol>> polarisedGammaSources;

	//Create and configure sources
	int err = createParticleGenerators(genericSources,
					   polarisedGammaSources,
					   out, verbose);
	
	if(err != 0){
	  if(verbose > 0 && out){
	    if(err > 0){
	      out << "\nErrors at sources creation"
		" and configuration.\n" << std::endl;
	    }
	  }
	  return errors::ERROR_AT_SOURCE_CONFIGURATION;
	}

	if(verbose > 1 && out)
	  out << std::endl;

	std::chrono::seconds::rep balanceInterval =
	  static_cast<std::chrono::seconds::rep>(1.0e9);

	//Set interval between checks
	for(auto& source : genericSources){
	  source.setCheckTime(balanceInterval);
	  source.setLBthreshold(static_cast<unsigned long long>(balanceInterval/2));
	}
	for(auto& source : polarisedGammaSources){
	  source.setCheckTime(balanceInterval);
	  source.setLBthreshold(static_cast<unsigned long long>(balanceInterval/2));
	}

	//****************************
	// Context
	//****************************

	//Create simulation context
	contextType context;
  
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

	if(verbose > 1 && out){
	  out << "Maximum global energy: " << globEmax << " eV" << std::endl;
	}

	//****************************
	// Materials
	//****************************

	if(verbose > 1 && out){
	  out << "\n\n------------------------------------\n\n"
	    " **** Materials ****\n"
	    " *******************\n" << std::endl; 
	}

	//Initialize context with no geometry
	pen_parserSection matInfoSection;
	if(context.configure(globEmax,
			     contextConfig,
			     matInfoSection,
			     verbose) != pen_context::SUCCESS){
	  return errors::ERROR_AT_CONTEXT_CONFIGURATION;
	}

	//****************************
	// Geometry parameters
	//****************************

	if(verbose > 1 && out){
	  out << "\n\n------------------------------------\n\n"
	    " **** Geometry parameters\n" << std::endl;
	}

	std::shared_ptr<wrapper_geometry> geometry;
  
	//Create and configure geometry
	int errGeo = createGeometry(geometry, matInfoSection, out, verbose);
	if(errGeo != errors::SUCCESS)
	  return errGeo;

	//Set geometry to context
	if(context.setGeometry(geometry.get()) != 0){
	  if(verbose > 0 && out){
	    out << "Error setting geometry to context." << std::endl;
	  }
	  return errors::ERROR_SETTING_GEOMETRY_TO_CONTEXT;
	}
  
	//Set source geometry
	for(unsigned i = 0; i < genericSources.size(); i++)
	  genericSources[i].setGeometry(geometry.get());  
	for(unsigned i = 0; i < polarisedGammaSources.size(); i++)
	  polarisedGammaSources[i].setGeometry(geometry.get());  

	//------------------------------------------------------------

	//Run the second step of context configuration with the geometry
	err = context.configureWithGeo(contextConfig, verbose);
	if(err != 0){
	  return errors::ERROR_AT_CONTEXT_CONFIGURATION_WITH_GEOMETRY;
	}

	//****************************
	// Create tallies
	//****************************
	std::vector<pen_commonTallyCluster> talliesVect;

	//Create one tally group per thread
	talliesVect.resize(nThreads);

	//Init tallies
	if(createTallies(talliesVect, context, out, verbose) != 0)
	  return errors::ERROR_CREATING_TALLIES;

	out << "\n" << talliesVect[0].numTallies() << " tallies created for each thread." << std::endl;
	if(talliesVect[0].numTallies() == 0){
	  out << "The simulation will not extract any information. Abort it." << std::endl;
	  return errors::SUCCESS;
	}

	
	//****************************
	// Variance Reduction 
	//****************************
	pen_VRCluster<pen_particleState> genericVR;
	pen_VRCluster<pen_state_gPol> photonVR;
	int vrRet = setVarianceReduction(geometry,
					 genericVR, photonVR,
					 out, verbose);
	if(vrRet != errors::SUCCESS){
	  if(verbose > 0 && out){
	    out << "Error on variance reduction section.\n"
	      "             Error code: " << vrRet << std::endl;
	  }
	  return vrRet;
	}

	//****************************
	// History loop
	//****************************

	out.flush();

	// ************************** MULTI-THREADING ***************************** //
#ifdef _PEN_USE_THREADS_
	std::vector<std::thread> simThreads;
#endif
	// ************************ MULTI-THREADING END *************************** //

	pen_timer timer;
	double time0 = CPUtime();
  
	for(unsigned i = 0; i < nThreads; i++){
	  talliesVect[i].run_beginSim();
	}

	if(verbose > 1 && out){
	  out << "Initialization processing time: "
	      << initializationTimer.timer()
	      << " s" << std::endl;
	}	

	//Substract initialization time to maximum simulation time
	long long int initTime = static_cast<long long int>(initializationTimer.timer());
	for(unsigned ithread = 0; ithread < nThreads; ++ithread)
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
	  if(nThreads > 1){
	    // Multi-thread
	    //****************
      
	    for(unsigned ithread = 0; ithread < nThreads; ++ithread){
	      // Simulate
	      simThreads.push_back(std::thread(_simulate<contextType, pen_particleState>,
					       std::ref(simConfigs[ithread]),
					       std::ref(context),
					       std::ref(genericSources[iSource]),
					       std::ref(talliesVect[ithread]),
					       std::ref(genericVR),
					       std::ref(photonVR)));
	    }
    
	    //Join threads
	    for(unsigned ithread = 0; ithread < nThreads; ++ithread){
	      simThreads[ithread].join();
	      //Update remaining simulation time
	      simConfigs[0].maxSimTime =
		std::min(simConfigs[0].maxSimTime, simConfigs[ithread].maxSimTime);
	    }

	    //Update maximum simulation times
	    for(unsigned ithread = 0; ithread < nThreads; ++ithread){
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
	    penred::simulation::sampleAndSimulateContext(simConfigs[0], context,
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
	  if(nThreads > 1){
	    // Multi-thread
	    //***************
      
	    for(unsigned ithread = 0; ithread < nThreads; ++ithread){
	      // Simulate
	      simThreads.push_back(std::thread(_simulate<contextType, pen_state_gPol>,
					       std::ref(simConfigs[ithread]),
					       std::ref(context),
					       std::ref(polarisedGammaSources[iSource]),
					       std::ref(talliesVect[ithread]),
					       std::ref(genericVR),
					       std::ref(photonVR)));
	
	    }

	    //Join threads
	    for(unsigned ithread = 0; ithread < nThreads; ++ithread){
	      simThreads[ithread].join();
	      //Update remaining simulation time
	      simConfigs[0].maxSimTime =
		std::min(simConfigs[0].maxSimTime, simConfigs[ithread].maxSimTime);

	    }

	    //Update maximum simulation times
	    for(unsigned ithread = 0; ithread < nThreads; ++ithread){
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
	    penred::simulation::sampleAndSimulateContext(simConfigs[0], context,
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
	for(unsigned ithread = 0; ithread < nThreads; ++ithread){
	  totalHists += simConfigs[ithread].getSimulatedInFinished();
	}
    
	//End simulations
	for(unsigned ithread = 0; ithread < nThreads; ithread++){
	  talliesVect[ithread].run_endSim(simConfigs[ithread].getSimulatedInFinished());
	}

	double simtime = timer.timer();
	double CPUendSim = CPUtime();
	double usertime = CPUendSim-time0;

	if(verbose > 1 && out){
	  out << "SImulation finished, starting results reduce step" << std::endl;
	}

	//If enabled, print final dumps
	if(finalDump){
	  for(unsigned ithread = 0; ithread < nThreads; ithread++){
	    int seeds1, seeds2;
	    simConfigs[ithread].getSeeds(seeds1, seeds2);
	    talliesVect[ithread].dump2file(simConfigs[ithread].dumpFilename.c_str(),
					   simConfigs[ithread].getSimulatedInFinished(),
					   seeds1,seeds2,-1,0ull,verbose);
	  }
	}

	//Sum tallies of all threads
	for(unsigned ithread = 1; ithread < nThreads; ithread++){
	  talliesVect[0].sum(talliesVect[ithread],verbose);
	}

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
  
	out << "\n\nSimulated histories: " << totalHists << "\n"
	    << "Simulation real time: " << simtime << " s\n"
	    << "Simulation user time: " << usertime << " s\n"
	    << "Histories per second and thread: " << static_cast<double>(totalHists)/usertime << "\n"
	    << "Histories per second: " << static_cast<double>(totalHists)/(usertime/double(nThreads)) << "\n"
	    << "Results processing time: " << postProcessTime <<" s" << std::endl;

	// ******************************* LB ************************************ //
#ifdef _PEN_USE_LB_
  
	//Print load balance reports
	if(genericSources.size() > 0){
	  if(verbose > 1 && out)
	    out << "Printing load balance reports for generic sources...";
	  for(const auto& source : genericSources){
	    FILE* fout = nullptr;
	    std::string filenameLBTh(source.name);
	    filenameLBTh += std::string("-generic-LBreport.rep");
	    fout = fopen(filenameLBTh.c_str(),"w");
	    source.task.printReport(fout);
	    fclose(fout);	    
	  }
	  if(verbose > 1 && out)
	    out << " Done!" << std::endl;
	}

	if(polarisedGammaSources.size() > 0){
	  if(verbose > 1 && out)
	    out << "Printing load balance reports for polarised sources...";
	  for(const auto& source : polarisedGammaSources){
	    FILE* fout = nullptr;
	    std::string filenameLBTh(source.name);
	    filenameLBTh += std::string("-polarised-LBreport.rep");      
	    fout = fopen(filenameLBTh.c_str(),"w");
	    source.task.printReport(fout);
	    fclose(fout);
	  }
	  if(verbose > 1 && out)
	    out << " Done!" << std::endl;
	}
#endif

	return errors::SUCCESS;
      }
      
    };

  } // namespace simulation
  
} //namespace penred

#endif
