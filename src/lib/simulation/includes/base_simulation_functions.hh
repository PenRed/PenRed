
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

 
#ifndef __PEN_BASE_SIMULATION_FUNCTIONS__
#define __PEN_BASE_SIMULATION_FUNCTIONS__

#include <mutex>
#include <iostream>
#include <sstream>
#include <utility>

#include "PenRed.hh"
#include "pen_samplers.hh"
#include "pen_tallies.hh"
#include "pen_geometries.hh"
#include "pen_vr.hh"

namespace penred{

  namespace simulation{

    //Define tally function type
    typedef std::function<void(const pen_particleState&,  //Resulting state
      const pen_KPAR,  //State kpar
      const unsigned long long&, //History number
      const int        //Returned value by 'simulatePartCond' function
      )> tallyFuncType;
    //Define sample function type
    template<class particleState>
    using sampleFuncType =
      std::function<void(particleState&,      //Generated state
      pen_KPAR&,           //Generated kpar
      unsigned long long&, //History increment
      const unsigned,      //Thread number
      pen_rand&)>;         //Random generator instance

    namespace finishTypes{
      enum finishTypes{
	NO_CONDITION = 0,
	DETECTOR_REACHED = 1,
      };
    }

    namespace errors{
      enum errors{
	SUCCESS = 0,
	GEOMETRY_NOT_REACHED,
	KPAR_NOT_FOUND,
	ERROR_INVALID_SEEDS,
	ERROR_LOADING_DUMP,
	ERROR_SECTION_NOT_FOUND,
	ERROR_INVALID_SEED_PAIR,
	ERROR_NO_SOURCE,
	ERROR_MISSING_PATH,
	ERROR_MISSING_SOURCE_CONFIGURATION,
	ERROR_MISSING_GEOMETRY_CONFIGURATION,
	ERROR_MISSING_TALLY_CONFIGURATION,
	ERROR_AT_SOURCE_CONFIGURATION,
	ERROR_AT_CONTEXT_CONFIGURATION,
	ERROR_AT_CONTEXT_CONFIGURATION_WITH_GEOMETRY,
	ERROR_MISSING_CONTEXT_CONFIGURATION,
	ERROR_MISSING_TYPE,
	ERROR_UNKNOWN_TYPE,
	ERROR_AT_GEOMETRY_CONFIGURATION,
	ERROR_SETTING_GEOMETRY_TO_CONTEXT,
	ERROR_CREATING_TALLIES,
	ERROR_ON_TALLIES_CONFIGURATION,
	ERROR_ON_VR_CONFIGURATION,
	ERROR_PARSING_CONFIG,
	ERROR_SIMULATION_RUNNING,
      };
    }

    namespace simFlags{
      enum simFlags{
	SUCCESS = 0,
	SCHEDULED,
	NOTHING_TO_DO,
	END_OF_SOURCE,
	INITIALIZING,
	RUNNING,
	MAX_TIME_REACHED,
	ERROR_UNABLE_TO_START_WORKER,
      };
    }

    //Auxiliary functions
    template<class stateType>
    inline void updateHists(unsigned long long& nhists,
			    pen_specificStateGen<stateType>& source,
			    const unsigned ithread,
			    const unsigned verbose){
      //Update assigned histories
      unsigned long long newNhists = source.toDo(ithread);
      if(newNhists != nhists){
	if(verbose > 1)
	  penred::logs::logger::
	    printf(penred::logs::SIMULATION,
		   "Thread %u: Source '%s': Number of histories to do "
		   "updated from %llu to %llu.\n",
		   ithread,source.name.c_str(),nhists,newNhists);
	nhists = newNhists;
      }
    }

    template<class stateType>
    inline void checkpoint(pen_specificStateGen<stateType>& source,
			   const unsigned ithread,
			   const unsigned verbose){
      int errCP = source.checkPoint(verbose);
      if(errCP != 0 && verbose > 0){
	penred::logs::logger::
	  printf(penred::logs::SIMULATION,
		 "Thread %u: Source '%s': Error at load balancing "
		 "checkpoint. Error code: %d\n",
		 ithread,source.name.c_str(),errCP);
      }
      else if(errCP == 0 && verbose > 1){
	penred::logs::logger::
	  printf(penred::logs::SIMULATION,
		 "Thread %u: Source '%s': Load balance checkpoint "
		 "done.\n",
		 ithread,source.name.c_str());
      }
    }

    struct simState{

    private:
      mutable std::mutex locker;
      //Actual simulating source
      std::string actualSourceName;
      //Initial already simulated histories in this source
      unsigned long long initiallySimulatedInSource;
      //Histories already simulated in actual source
      unsigned long long simulatedInSource;
      //Histories to simulate in actual source
      unsigned long long toSimulateInSource;
      //Total number of simulated histories in previous finished sources
      unsigned long long finishedSimulated;
      //Simulation init and end time of the actual source
      std::chrono::steady_clock::time_point sourceStart;
      std::chrono::steady_clock::time_point sourceFinish;
      //Actual source simulation speed
      double sourceSpeed;
      //Estimated time to finish the simulation of actual source
      double sourceETA;
      //Flags if the simulation of the actual source has started
      bool sourceStarted;
      //Flags if the simulation of the actual source has finished
      bool sourceFinished;
      //Simulation status flag for the actual source
      std::atomic<simFlags::simFlags> sourceSimFlag;           
      
    public:

      simState() : actualSourceName("none"),
		   initiallySimulatedInSource(0),
		   simulatedInSource(0),
		   toSimulateInSource(0),
		   finishedSimulated(0),
		   sourceSpeed(1.0),
		   sourceETA(0.0),
		   sourceStarted(false),
		   sourceFinished(false),
		   sourceSimFlag(simFlags::SUCCESS)
      {}

      //Init and end source function
      inline void nextSource(const std::string& sourceName,
			     const unsigned long long& simulated,
			     const unsigned long long& toSimulate){

	std::lock_guard<std::mutex> guard(locker);

	actualSourceName = sourceName;
	initiallySimulatedInSource = simulated;
	simulatedInSource = initiallySimulatedInSource;
	toSimulateInSource = toSimulate;
	sourceSpeed = 0.0;
	sourceETA = 0.0;
	sourceFinished = false;

	//Flag the source as started and save start timestamp
	sourceStarted = true;
	sourceStart = std::chrono::steady_clock::now();
	
	sourceSimFlag = simFlags::INITIALIZING;
      }

      inline void finish() {
	std::lock_guard<std::mutex> guard(locker);
	if(sourceStarted && !sourceFinished){
	  sourceFinished = true;
	  sourceFinish = std::chrono::steady_clock::now();
	  finishedSimulated += simulatedInSource;
	}
      }

      //Elapsed time, in seconds, since last start
      inline double sinceStart() const{
	std::lock_guard<std::mutex> guard(locker);

	auto now = std::chrono::steady_clock::now();
	return static_cast<double>(std::chrono::duration_cast<
				   std::chrono::seconds>(now-sourceStart).count());
      }

      //Simulation time, in seconds, for current source
      inline double simTime() const{
	std::lock_guard<std::mutex> guard(locker);
	if(sourceFinished){
	  return static_cast<double>(std::chrono::duration_cast< std::chrono::seconds>
				     (sourceFinish-sourceStart).count());
	}else{
	  auto now = std::chrono::steady_clock::now();
	  return static_cast<double>(std::chrono::duration_cast<
				     std::chrono::seconds>(now-sourceStart).count());
	}
      }
      
      //Report functions
      inline void report(const unsigned long long& hist) {
	std::lock_guard<std::mutex> guard(locker);

	//Update source history counter
	simulatedInSource = hist-finishedSimulated;

	auto now = std::chrono::steady_clock::now();
	const double elapsed = static_cast<double>(std::chrono::duration_cast<
						   std::chrono::seconds>(now-sourceStart).count());
	sourceSpeed =
	  static_cast<double>(simulatedInSource-initiallySimulatedInSource)/elapsed;
	if(toSimulateInSource < simulatedInSource)
	  sourceETA = 0.0;
	else{
	  sourceETA =
	    (static_cast<double>(toSimulateInSource - simulatedInSource))/sourceSpeed;
	}
      }

      //Set functions
      inline void updateToSimulate(const unsigned long long& toSimulate){
	std::lock_guard<std::mutex> guard(locker);

	toSimulateInSource = toSimulate;

	if(toSimulateInSource < simulatedInSource)
	  sourceETA = 0.0;
	else{
	  sourceETA =
	    (static_cast<double>(toSimulateInSource - simulatedInSource))/sourceSpeed;
	}
      }

      //Set finished sources simulated histories
      inline void setSimulatedInFinished(const unsigned long long& nHists){
	std::lock_guard<std::mutex> guard(locker);
	finishedSimulated = nHists;
      }
      
      inline void setFlag(const simFlags::simFlags& f) { sourceSimFlag = f; }

      //Get functions
      inline unsigned long long getSimulatedInSource() const {
	std::lock_guard<std::mutex> guard(locker);
	return simulatedInSource;
      }

      inline unsigned long long getToSimulateInSource() const {
	std::lock_guard<std::mutex> guard(locker);
	return toSimulateInSource;
      }

      inline unsigned long long getSimulatedInFinished() const {
	std::lock_guard<std::mutex> guard(locker);
	return finishedSimulated;
      }

      inline unsigned long long getInitiallySimulated() const {
	std::lock_guard<std::mutex> guard(locker);
	return initiallySimulatedInSource;
      }

      inline unsigned long long getTotalSimulated() const {
	std::lock_guard<std::mutex> guard(locker);
	if(sourceFinished)
	  return finishedSimulated;
	else
	  return finishedSimulated + simulatedInSource;
      }
      
      inline double getSpeedInSource() const {
	std::lock_guard<std::mutex> guard(locker);
	return sourceSpeed;
      }
      
      inline double getSourceETA() const {
	std::lock_guard<std::mutex> guard(locker);
	return sourceETA;
      }

      inline std::string stringify() const {
	std::lock_guard<std::mutex> guard(locker);

	char aux[300];
	sprintf(aux,"Simulated %llu/%llu histories in the current source.\n"
		"   Simulation speed         : %15.5E hist/s\n"
		"   Estimated time remaining : %15.5E s\n",
		simulatedInSource, toSimulateInSource, 
		sourceSpeed, sourceETA);
	
	return std::string(aux);	
      }

      inline simFlags::simFlags readFlag() const { return sourceSimFlag; }
      
      inline bool isStarted() const { return sourceStarted; } 
      inline bool isFinished() const { return sourceFinished; }

    };    

    //Structure with configuration parameters
    struct simConfig : public penred::logs::logger{

    public:
      
      template <class T>
      friend simConfig& operator<<(simConfig& out, const T& o);

      friend simConfig& operator<<(simConfig& out, const char* o);

      static constexpr const struct Endl{
      } endl{};
      
    private:
      simState status;

      //Initial and last source seeds 
      int iSeed1, iSeed2;
      int lSeed1, lSeed2;
      
      //Recovery/dump related config
      int firstSourceIndex; 
      unsigned long long simulatedHistsInFirstSource;

      //Simulated sources list
      mutable std::mutex sourcesVecLock;      
      std::vector<std::string> simulatedSources;
      std::string actualSource;

      //Auxiliary string stream
      std::stringstream auxOut;

      static constexpr bool noFinishSim(const unsigned long long){ return true; } 
      
    public:

      //Thread number
      unsigned iThread;
      long long int dumpTime;    //In ms
      std::string dumpFilename;
      //Enable/disable write partial results
      bool writePartial;
      //Maximum simulation time
      long long int maxSimTime;  //In ms

      //Finish simulation function.
      //Returns true if the simulation must continue and false to stop it
      std::function<bool(const unsigned long long)> fSimFinish;
      
      //Verbose level
      unsigned verbose;

      simConfig();
      simConfig(const unsigned& iThreadIn,
		const long long int& dumpTimeIn,
		const std::string& dumpFilenameIn,
		const bool& writePartialIn,
		const long long int& maxSimTimeIn,
		const int& iSeed1, const int& iSeed2,
		const unsigned& verboseIn);

      //Print functions
      inline std::string threadPrefix() const {
	return std::string("Thread ") +
	  std::to_string(iThread) +
	  std::string(": ");
      }
      inline std::string sourcePrefix() const {
	return std::string("Source ") +
	  actualSource +
	  std::string(": ");
      }
      inline std::string threadAndSourcePrefix() const {
	return threadPrefix() + sourcePrefix();
      }

      inline void flush(bool appendEndl = false){
	//Flush auxiliary stream to final stream

	//Read the whole auxiliary stream
	std::string s(std::istreambuf_iterator<char>(auxOut), {});
	if(appendEndl)
	  s += '\n';
	//Write it to the output stream
	cout << s;

	//Flush output stream
	cout.flush();	
      }
      
      //Sim status functions
      inline void nextSource(const std::string& name,
			     const unsigned long long& toSimulate){
	std::lock_guard<std::mutex> guard(sourcesVecLock);
	  
	if(simulatedSources.size() == 0){
	  //Is the first source to be simulated
	  //Specify the simulation histories already simulated
	  status.nextSource(name, simulatedHistsInFirstSource, toSimulate);
	}else{
	  status.nextSource(name, 0, toSimulate);
	}

	if(verbose > 1){
	  auxOut << threadPrefix() << "Starting simulation of source "
		 << name << " with initial history "
		 << status.getTotalSimulated()
		 << " and seeds " << lSeed1 << " " << lSeed2
		 << ". Histories already simulated: "
		 << status.getSimulatedInSource()
		 << ". Remaining histories: "
		 << status.getToSimulateInSource();
	  flush(true);
	}
	

	//Update initial seeds to last ones
	iSeed1 = lSeed1;
	iSeed2 = lSeed2;

	//Add source to simulated vector
	simulatedSources.push_back(name);
	actualSource = name;
      }

      inline void finish() {
	status.finish();
      }

      inline double simTime() const{
	return status.simTime();
      }
      
      inline void report(const unsigned long long& hist){
	status.report(hist);
      }

      inline void updateToSimulate(const unsigned long long& toSimulate){
	status.updateToSimulate(toSimulate);
      }
      
      //Set functions      
      inline void setSimulatedInFinished(const unsigned long long& nHists){
	status.setSimulatedInFinished(nHists);
      }

      inline void setSeeds(const int& seed1, const int& seed2){
	iSeed1 = seed1;
	lSeed1 = seed1;

	iSeed2 = seed2;
	lSeed2 = seed2;	
      }

      inline void setSeeds(const int& seedPos){
	rand0(seedPos, iSeed1, iSeed2);
	lSeed1 = iSeed1;
	lSeed2 = iSeed2;
      }

      inline void setInitSeedsToRand(pen_rand& random){
	random.setSeeds(iSeed1, iSeed2);
      }

      inline void setFinalSeedsFromRand(const pen_rand& random){
	random.getSeeds(lSeed1,lSeed2);
      }
      
      inline void setInitialSource(const int& initSource,
				   const unsigned long long& simulated){
	firstSourceIndex = initSource;
	simulatedHistsInFirstSource = simulated;
      }

      inline size_t getNSourceNames() const {
	std::lock_guard<std::mutex> guard(sourcesVecLock);	
	return simulatedSources.size();
      }
      
      inline int getCurrentSourceIndex() const {
	std::lock_guard<std::mutex> guard(sourcesVecLock);	
	return firstSourceIndex + static_cast<int>(simulatedSources.size())-1;
      }

      inline std::string getCurrentSourceName() const {
	return actualSource;
      }

      inline int getFirstSourceIndex() const {
	return firstSourceIndex;
      }

      inline unsigned long long getInitiallySimulatedInFirstSource() const {
	return simulatedHistsInFirstSource;
      }
      
      //Status get functions
      inline simState& getStatus() { return status; }
      
      inline unsigned long long getSimulatedInSource() const {
	return status.getSimulatedInSource();
      }

      inline unsigned long long getToSimulateInSource() const {
	return status.getToSimulateInSource();
      }

      inline unsigned long long getSimulatedInFinished() const {
	return status.getSimulatedInFinished();
      }
      
      inline unsigned long long getInitiallySimulated() const {
	return status.getInitiallySimulated();
      }

      inline unsigned long long getTotalSimulated() const {
	return status.getTotalSimulated();
      }

      inline void getSeeds(int& seed1, int& seed2) const {
	seed1 = lSeed1;
	seed2 = lSeed2;	
      }
      

      inline double getSpeedInSource() const {
	return status.getSpeedInSource();
      }      

      inline double sinceSourceStart() const {
	return status.sinceStart();
      }
      
      inline double getSourceETA() const {
	return status.getSourceETA();
      }

      inline std::string stringifyState() const {
	return threadAndSourcePrefix() + status.stringify();
      }

      std::string stringifyConfig() const;

      inline bool isSourceStarted() const { return status.isStarted(); } 
      inline bool isSourceFinished() const { return status.isFinished(); }

      int configure(const pen_parserSection& config);
      inline int configure(const char* prefix, const pen_parserSection& config){

	//Parse only if the 'prefix' section is included in the 'config'
	//Otherwise, the instance remains unchanged
	pen_parserSection prefixSection;
	if(config.readSubsection(prefix, prefixSection) == INTDATA_SUCCESS){
	  return configure(prefixSection);
	}
	return errors::SUCCESS;
      }

      inline int readTallyDump(const std::string& readFilename,
			       pen_commonTallyCluster& tallies){

	//Read dump file
	unsigned long long simulated;
	int initialSource;
	unsigned long long simulatedInSource;
	int auxSeed1, auxSeed2;
	int errDump = tallies.readDumpfile(readFilename.c_str(),
					   simulated,
					   auxSeed1,auxSeed2,
					   initialSource,
					   simulatedInSource,
					   verbose);
	if(errDump != 0){
	  return errDump;
	}

	//Set seeds read from dump
	setSeeds(auxSeed1,auxSeed2);
	//Set initial source index and histories already simulated
	if(initialSource >= 0)
	  setInitialSource(initialSource,simulatedInSource);
	else
	  setInitialSource(0,0);
	//Negative source index is interpreted as completed simulation

	//Update simulated histories in finished sources
	if(simulated > simulatedInSource)
	  setSimulatedInFinished(simulated - simulatedInSource);
	else
	  setSimulatedInFinished(0);
	
	return errDump;
      }

      inline void copyCommonConfig(const simConfig& c){

	//Copy only common configuration parameters
	//to all workers in the same simulation
	iSeed1 = c.iSeed1;
	iSeed2 = c.iSeed2;

	lSeed1 = c.lSeed1;
	lSeed2 = c.lSeed2;
	
	dumpTime = c.dumpTime;
	dumpFilename = c.dumpFilename;
	writePartial = c.writePartial;
	maxSimTime = c.maxSimTime;
	verbose = c.verbose;
      }

      inline std::vector<simConfig> generateThreadConfigs(const size_t nThreads,
							  const unsigned initSeed = 0){

	std::vector<simConfig> result(nThreads);
	for(size_t i = 0; i < nThreads; ++i){
	  result[i].iThread = i;
	  result[i].copyCommonConfig(*this);
	  result[i].setSeeds(initSeed+i);
	}

	return result;
      }
      
      inline ~simConfig(){
	//Flush stream on delete
	flush();
      }

    };

    // Define extraction operators
    //*****************************
    
    template<class T>
    inline simConfig& operator<<(simConfig& sc, const T& obj) {
      // write obj to auxiliary stream
      sc.auxOut << obj;
      return sc;
    }

    //Extraction operator for strings
    template<>
    inline simConfig& operator<<(simConfig& sc, const std::string& obj) {

      //Get the last end of line character position '\n'
      std::string::size_type p = obj.rfind('\n');
      if(p == std::string::npos){
	//No end of line provided
	//Write the whole object in the stream
	sc.auxOut << obj;
      }else{
	//End of line found. Write until this point
	//and flush the stream
	sc.auxOut << obj.substr(0,p+1);
	sc.flush();

	//Write the remaining string to auxiliary buffer
	sc.auxOut << obj.substr(p+1);
      }
      
      return sc;
    }

    //Extraction operation for const char*
    inline simConfig& operator<<(simConfig& sc, const char* obj) {
      return operator<<(sc,std::string(obj));
    }

    //Extraction operation for endl
    template<>
    inline simConfig& operator<<(simConfig& sc, const simConfig::Endl& /*obj*/) {
      //Flush appending an end line character
      sc.flush(true);
      return sc;
    }
    
    //*******************************************

    //Move to geo function family:
    //
    // This function check if a sampled particle has been created in the void.
    // In that case, it tries to move the particle to reach a non void region
    //
    template <class particleType>
    bool move2geo(const unsigned long long& nhist,
		  particleType& particle,
		  pen_commonTallyCluster& tallies){

      pen_particleState& state = particle.getBaseState();
      //Check if particle has been generated at vacuum
      if(state.MAT == 0){
	//Move it to geometry system
	particle.jumpVolume();
    
	//Call tallies with move2geo collect function
	unsigned kdet = particle.getDET();
	tallies.run_move2geo(nhist,kdet,particle.kpar,
			     state,particle.DSef(),particle.DStot());

	if(state.MAT == 0){
	  tallies.run_endPart(nhist,particle.kpar,state);
	  return false;
	}
      }
      else{
	particle.updateMat();
	particle.updateBody();
      }
      return true;
    }


    template <class particleType,
	      class stateType>
    bool move2geo(const unsigned long long& nhist,
		  particleType& particle,
		  pen_commonTallyCluster& tallies,
		  pen_specificTallyCluster<stateType>& specificTallies){

      stateType& state = particle.getState();
      //Check if particle has been generated at vacuum
      if(state.MAT == 0){
	//Move it to geometry system
	particle.jumpVolume();
    
	//Call tallies with move2geo collect function
	unsigned kdet = particle.getDET();
	tallies.run_move2geo(nhist,kdet,particle.kpar,state,
			     particle.DSef(),particle.DStot());
	specificTallies.run_move2geo(nhist,kdet,particle.kpar,
				     state,particle.DSef(),particle.DStot());
    
	if(state.MAT == 0){
	  tallies.run_endPart(nhist,particle.kpar,state);
	  specificTallies.run_endPart(nhist,particle.kpar,state);
	  return false;
	}
      }
      else{
	particle.updateMat();
	particle.updateBody();
      }  
      return true;
    }

    template <class particleType>
    bool move2geo(particleType& particle){

      pen_particleState& state = particle.getBaseState();
      //Check if particle has been generated at vacuum
      if(state.MAT == 0){
	//Move it to geometry system
	particle.jumpVolume();

	if(state.MAT == 0){
	  return false;
	}
      }
      else{
	particle.updateMat();
	particle.updateBody();
      }
      return true;
    }

    //Absorb function family:
    //
    // These functions check if the particle must be absrobed
    // in its current state
    //
    template <class particleType>
    inline bool absorb(const unsigned long long& nhist,
		       particleType& particle,
		       pen_commonTallyCluster& tallies,
		       pen_rand& penRand){

      const pen_particleState& state = particle.readBaseState();

      //Check if particle must be absorbed
      if(state.E < particle.getEABS()){ 

	//Check if the particle has already been annihilated
	if(state.E > 1.0e-6){
	  // run annihilation process
	  double Eprod = particle.annihilationEDep;
	  particle.annihilate(penRand);
	  tallies.run_localEdep(nhist,particle.kpar,state,state.E+Eprod);
	}
	//Call tallies with end particle collect function
	tallies.run_endPart(nhist,particle.kpar,state);
	return true;
      }

      //Particle can't be absorbed yet
      return false;
    }

    template <class particleType,
	      class stateType>
    inline bool absorb(const unsigned long long& nhist,
		       particleType& particle,
		       pen_commonTallyCluster& tallies,
		       pen_specificTallyCluster<stateType>& specificTallies,
		       pen_rand& penRand){

      const stateType& state = particle.readState();

      //Check if particle must be absorbed
      if(state.E < particle.getEABS()){ 

	//Check if the particle has already been annihilated
	if(state.E > 1.0e-6){
	  // run annihilation process
	  double Eprod = particle.annihilationEDep;
	  particle.annihilate(penRand);
	  tallies.run_localEdep(nhist,particle.kpar,state,state.E+Eprod);
	  specificTallies.run_localEdep(nhist,particle.kpar,state,state.E+Eprod);
	}
	//Call tallies with end particle collect function
	tallies.run_endPart(nhist,particle.kpar,state);
	specificTallies.run_endPart(nhist,particle.kpar,state);      
	return true;
      }

      //Particle can't be absorbed yet
      return false;
    }

    template <class particleType>
    inline bool absorb(particleType& particle,
		       pen_rand& penRand){

      const pen_particleState& state = particle.readBaseState();

      //Check if particle must be absorbed
      if(state.E < particle.getEABS()){ 

	//Check if the particle has already been annihilated
	if(state.E > 1.0e-6){
	  // run annihilation process
	  particle.annihilate(penRand);
	}
	return true;
      }

      //Particle can't be absorbed yet
      return false;
    }

    //Simulate particle function family
    //
    // These functions simulates a particle already situated in
    // a non void region. This shoud be satisfied if the particle
    // to be simulated is a secondary produced by a previous particle
    // simulation, or the particle to be simulated has been sampled and
    // "moved" using the move2geo function. 
    //
    template <class particleType>
    void simulatePart(const unsigned long long& nhist,
		      particleType& particle,
		      pen_commonTallyCluster& tallies,
		      pen_rand& randoms){

      pen_particleState& state = particle.getBaseState();
      unsigned kdet = particle.getDET();
      const pen_KPAR kpar = particle.kpar;

      //Check if particle scaped from geometry
      if(state.MAT == 0){
	//Call tallies with end particle collect function
	tallies.run_endPart(nhist,particle.kpar,state);      
	return;
      }  
  
      //Check if the particle must be absorbed
      if(absorb(nhist,particle,tallies,randoms))
	return;

      particle.START();
  
      for(;;){

	double ds, de;
	int icol;

	//Create step data structure
	tally_StepData stepData;
    
	//Check for interaction forcing
	bool iforced = particle.readContext().isForcing(kpar,state.IBODY,state.WGHT);

	//Calculate path to next iteration
	if(iforced)
	  particle.JUMPF(ds, randoms, particle.getDSMAX());
	else
	  particle.JUMP(ds, randoms, particle.getDSMAX());

	//Call tallies with jump collect function
	tallies.run_jump(nhist,kpar,state,ds);

	//Move particle across geometry
	double lastEABS = particle.getEABS();
	particle.move(ds,stepData.softDE,
		      stepData.softX,stepData.softY,stepData.softZ,
		      randoms);    

	//Update step data
	stepData.update(particle); 
	//Call tallies with step collect function
	tallies.run_step(nhist,kpar,state,stepData);

	//Check if the particle must be absorbed at original material and body
	if(state.E < lastEABS){

	  double Eprod = particle.annihilationEDep;
	  if(state.MAT > 0){
	    // run annihilation process
	    particle.annihilate(randoms);
	  }
	  //Set to the particle state the origin material and body
	  state.IBODY = particle.lastBody();
	  state.MAT = particle.lastMat();

	  //Set position to the end of dsef
	  double XL,YL,ZL;
	  particle.lastPos(XL,YL,ZL);
	  state.X = XL+state.U*particle.DSef();
	  state.Y = YL+state.V*particle.DSef();
	  state.Z = ZL+state.W*particle.DSef();

	  tallies.run_localEdep(nhist,kpar,state,state.E+Eprod);
	  //Call tallies with end particle collect function
	  tallies.run_endPart(nhist,kpar,state);
	  return;
	}
      
	//Check if particle cross some interface
	if(particle.NCross() != 0){
      
	  //Interface crossed, Update kdet
	  kdet = particle.getDET();
      
	  //Call tallies with cross interface collect function
	  tallies.run_interfCross(nhist,kdet,kpar,state);      
	
	  //Check if material has been changed
	  if(particle.lastMat() != state.MAT){
	    //Call changed material collect function
	    tallies.run_matChange(nhist,kpar,state,particle.lastMat());
	
	    if(state.MAT == 0){
	      //New material is vacuum
	      //Call tallies with end particle collect function
	      tallies.run_endPart(nhist,kpar,state);
	      return;
	    }

	    // Particle has been stopped at interface,
	    // restart particle state
	    particle.START();

	    //Check if the particle must be absorbed
	    if(absorb(nhist,particle,tallies,randoms))
	      return;

	    //VR
	    double deVR = particle.vr_matChange(nhist,randoms,2);	
	    //VR
	    deVR += particle.vr_interfCross(nhist,randoms,2);

	    //Tally the deposited energy by VR
	    if(deVR > 0.0){
	      //We must compensate the particle weight to score correctly
	      //the energy stored in the secondary stacks by the VR call
	      double originWGHT = state.WGHT;
	      state.WGHT = 1.0;
	      tallies.run_localEdep(nhist,kpar,state,deVR);
	      state.WGHT = originWGHT;
	    }

	    //Check if the particle must be absorbed because VR
	    if(absorb(nhist,particle,tallies,randoms))
	      return;
	
	  }
	  else{
	    // Particle has been stopped at interface,
	    // restart particle state
	    particle.START();

	    //Check if the particle must be absorbed
	    if(absorb(nhist,particle,tallies,randoms))
	      return;
      
	    //VR
	    double deVR = particle.vr_interfCross(nhist,randoms,2);

	    //Tally the deposited energy by VR
	    if(deVR > 0.0){
	      //We must compensate the particle weight to score correctly
	      //the energy stored in the secondary stacks by the VR call
	      double originWGHT = state.WGHT;
	      state.WGHT = 1.0;
	      tallies.run_localEdep(nhist,kpar,state,deVR);
	      state.WGHT = originWGHT;
	    }

	    //Check if the particle must be absorbed because VR
	    if(absorb(nhist,particle,tallies,randoms))
	      return;	
	  }
	}
	else{

	  //No interface crossed, simulate interaction
	  if(iforced)
	    particle.KNOCKF(de,icol,randoms);
	  else
	    particle.KNOCK(de,icol,randoms);

	  //Call tallies with knock collect function
	  tallies.run_knock(nhist,kpar,state,icol);
	  //Call tallies with eloss collect function
	  tallies.run_localEdep(nhist,kpar,state,de);
	
	  //Check if the particle must be absorbed
	  if(absorb(nhist,particle,tallies,randoms)){
	    return;
	  }
	}
      }
    }

    template <class particleType>
    int simulatePartCond(const unsigned long long& nhist,
			 particleType& particle,
			 pen_rand& randoms,
			 const finishTypes::finishTypes& finishType,
			 const unsigned& finishValue){
  
      pen_particleState& state = particle.getBaseState();
      const pen_KPAR kpar = particle.kpar;

      //Check if particle scaped from geometry
      if(state.MAT == 0){
	return 0;
      }  
  
      //Check if the particle must be absorbed
      if(absorb(particle,randoms))
	return 0;

      particle.START();
  
      for(;;){

	double ds, de;
	int icol;

	//Create step data structure
	tally_StepData stepData;
    
	//Check for interaction forcing
	bool iforced = particle.readContext().isForcing(kpar,state.IBODY,state.WGHT);

	//Calculate path to next iteration
	if(iforced)
	  particle.JUMPF(ds, randoms, particle.getDSMAX());
	else
	  particle.JUMP(ds, randoms, particle.getDSMAX());

	//Move particle across geometry
	double lastEABS = particle.getEABS();
	particle.move(ds,stepData.softDE,
		      stepData.softX,stepData.softY,stepData.softZ,
		      randoms);    

	//Update step data
	stepData.update(particle);

	//Check if the particle must be absorbed at original material and body
	if(state.E < lastEABS){

	  if(state.MAT > 0){
	    // run annihilation process
	    particle.annihilate(randoms);
	  }
	  //Set to the particle state the origin material and body
	  state.IBODY = particle.lastBody();
	  state.MAT = particle.lastMat();

	  //Set position to the end of dsef
	  double XL,YL,ZL;
	  particle.lastPos(XL,YL,ZL);
	  state.X = XL+state.U*particle.DSef();
	  state.Y = YL+state.V*particle.DSef();
	  state.Z = ZL+state.W*particle.DSef();

	  return 0;
	}
      
	//Check if particle cross some interface
	if(particle.NCross() != 0){
      
	  if(finishType == finishTypes::DETECTOR_REACHED){      
	    //Finish on detector reached condition
	    unsigned kdet = particle.getDET();
	    if(finishValue == kdet)
	      return 1;
	  }
	
	  //Check if material has been changed
	  if(particle.lastMat() != state.MAT){
	
	    if(state.MAT == 0){
	      //New material is vacuum
	      return 0;
	    }

	    // Particle has been stopped at interface,
	    // restart particle state
	    particle.START();

	    //Check if the particle must be absorbed
	    if(absorb(particle,randoms))
	      return 0;

	    //VR
	    particle.vr_matChange(nhist,randoms,2);	
	    //VR
	    particle.vr_interfCross(nhist,randoms,2);

	    //Check if the particle must be absorbed because VR
	    if(absorb(particle,randoms))
	      return 0;
	
	  }
	  else{
	    // Particle has been stopped at interface,
	    // restart particle state
	    particle.START();

	    //Check if the particle must be absorbed
	    if(absorb(particle,randoms))
	      return 0;
      
	    //VR
	    particle.vr_interfCross(nhist,randoms,2);

	    //Check if the particle must be absorbed because VR
	    if(absorb(particle,randoms))
	      return 0;
	  }
	}
	else{

	  //No interface crossed, simulate interaction
	  if(iforced)
	    particle.KNOCKF(de,icol,randoms);
	  else
	    particle.KNOCK(de,icol,randoms);
	
	  //Check if the particle must be absorbed
	  if(absorb(particle,randoms)){
	    return 0;
	  }
	}
      }
    }

    //Simulate stack function family:
    //
    // These functions simulates "n" seconday particles stored
    // in the particle main stack belonging the particle "p" 
    //
    template<class particleType>
    inline unsigned simulateStack(const unsigned long long& hist,
				  particleType& p,
				  const unsigned& n,
				  pen_commonTallyCluster& tallies,
				  pen_rand& randoms){
      //Check the number of particles to simulate
      const unsigned toSim = std::min(n,p.nStacked());
  
      for(unsigned i = 0; i < toSim; ++i){
	//Get next state from particle stack
	p.setStateFromStack();

	//Get kdet
	unsigned kdet = p.getDET();
    
	//VR
	p.vr_particleStack(hist,randoms,2);

	//Simulate particle
	tallies.run_beginPart(hist,kdet,p.getKpar(),p.readState());
	simulatePart(hist,p,tallies,randoms);    
      }

      return toSim;
    }

    template<class particleType>
    inline unsigned simulateStackCond(const unsigned long long& hist,
				      particleType& p,
				      const unsigned& n,
				      pen_rand& randoms,
				      const finishTypes::finishTypes& finishType,
				      const unsigned& finishValue,
				      tallyFuncType& f){
  
      //Check the number of particles to simulate
      const unsigned toSim = std::min(n,p.nStacked());
  
      for(unsigned i = 0; i < toSim; ++i){
	//Get next state from particle stack
	p.setStateFromStack();

	//VR
	p.vr_particleStack(hist,randoms,2);

	//Simulate particle
	int ret = simulatePartCond(hist,p,randoms,finishType,finishValue);
	f(p.readState(), p.getKpar(), hist, ret);
      }

      return toSim;
    }

    //Simulate stacks N function family:
    //
    // These functions simulates "n" seconday particles of each
    // particle main stack belonging to each particle provided in the argument
    //
    template<class particleType>
    inline unsigned simulateStacksN(const unsigned long long& hist,
				    const unsigned& n,
				    pen_commonTallyCluster& tallies,
				    pen_rand& randoms,
				    particleType& p){
  
      //Try to simulate "n" particles in the stack
      simulateStack(hist, p, n, tallies, randoms);

      //Return the number of remaining particles in the stack
      //Notice that this is the last particle stack to process
      return p.nStacked();
    }

    template<class particleType, class... otherParticleTypes>
    inline typename std::enable_if<sizeof...(otherParticleTypes) >= 1, unsigned>::type
    simulateStacksN(const unsigned long long& hist,
		    const unsigned& n,
		    pen_commonTallyCluster& tallies,
		    pen_rand& randoms,
		    particleType& p,
		    otherParticleTypes&... others){
  
      //Try to simulate "n" particles in the stack
      simulateStack(hist, p, n, tallies, randoms);

      //Simulate other particle stacks
      unsigned remainingOthers = simulateStacksN(hist,n,tallies,randoms,others...);

      //Return the number of remaining particles in the stack
      //plus the other stacks after simulating "n" particles
      return p.nStacked() + remainingOthers;
    }

    template<class particleType>
    inline unsigned simulateStacksNCond(const unsigned long long& hist,
					const unsigned& n,
					pen_rand& randoms,
					const finishTypes::finishTypes& finishType,
					const unsigned& finishValue,
					tallyFuncType& f,
					particleType& p){
  
      //Try to simulate "n" particles in the stack
      simulateStackCond(hist, p, n, randoms, finishType, finishValue, f);

      //Return the number of remaining particles in the stack
      //Notice that this is the last particle stack to process
	return p.nStacked();      
    }

    template<class particleType, class... otherParticleTypes>
    inline typename std::enable_if<sizeof...(otherParticleTypes) >= 1, unsigned>::type
    simulateStacksNCond(const unsigned long long& hist,
			const unsigned& n,
			pen_rand& randoms,
			const finishTypes::finishTypes& finishType,
			const unsigned& finishValue,
			tallyFuncType& f,
			particleType& p,
			otherParticleTypes&... others){
  
      //Try to simulate "n" particles in the stack
      simulateStackCond(hist, p, n, randoms, finishType, finishValue, f);

      //Simulate other particle stacks
      unsigned remainingOthers = simulateStacksNCond(hist,n,randoms,
						     finishType,finishValue,
						     f,others...);
      
      //Return the number of remaining particles in the stack
      //plus the other stacks after simulating "n" particles
      return p.nStacked() + remainingOthers;
      
    }

    
    //Simulate stacks function family:
    //
    // These functions simulates all the particles in the main
    // stack belonging to each particle provided in the argument
    //
    template<class... particleTypes>
    inline void simulateStacks(const unsigned long long& hist,
			       pen_commonTallyCluster& tallies,
			       pen_rand& randoms,
			       particleTypes&... particles){

      unsigned remaining;
      do{
	//Simulate stack particles in chunks of 100 particles
	remaining = simulateStacksN(hist, 100, tallies, randoms, particles...);
      }while(remaining > 0);
  
    }

    template<class... particleTypes>
    inline void simulateStacksCond(const unsigned long long& hist,
				   pen_rand& randoms,
				   const finishTypes::finishTypes& finishType,
				   const unsigned& finishValue,
				   tallyFuncType& f,
				   particleTypes&... particles){

      unsigned remaining;
      do{
	//Simulate stack particles in chunks of 100 particles
	remaining = simulateStacksNCond(hist, 100, randoms,
					finishType, finishValue, f, particles...);
      }while(remaining > 0);
  
    }

    //Simulate shower function family
    //
    // These functions simulate a complete shower, starting from a primary
    // particle and ending with all the secondary generated particles
    //
    template<class primaryType, class... secondaryTypes>
    inline void simulateShowerKnown(const unsigned long long& hist,
				    pen_commonTallyCluster& tallies,
				    pen_rand& randoms,
				    primaryType& p,
				    secondaryTypes&... secondary){

      //Init Page
      p.page0();

      //Try to move the generated particle to geometry system
      if(move2geo(hist,p,tallies)){

	//Get detector ID
	unsigned kdet = p.getDET();	
	
	//Tally particle begins
	tallies.run_beginPart(hist,kdet,p.getKpar(),p.readState());
    
	//Simulate particle
	simulatePart(hist,p,tallies,randoms);

	//Simulate secondary particles
	simulateStacks(hist, tallies, randoms, p, secondary...);
      }

  
  
    }

    template<class primaryType, class... secondaryTypes>
    inline void simulateShowerKnownCond(const unsigned long long& hist,
					pen_rand& randoms,
					const finishTypes::finishTypes& finishType,
					const unsigned& finishValue,
					tallyFuncType& f,
					primaryType& p,
					secondaryTypes&... secondary){

      //Init Page
      p.page0();

      //Try to move the generated particle to geometry system
      if(move2geo(p)){

	//Simulate particle
	int ret = simulatePartCond(hist,p,randoms,finishType,finishValue);
	//Call tally function
	f(p.readState(), p.getKpar(), hist, ret);

	//Simulate secondary particles
	simulateStacksCond(hist, randoms, finishType, finishValue, f, p, secondary...);
      }

      
  
    }

    

    //Simulate shower by kpar family
    //
    // These functions simulate a complete shower, starting from a primary
    // particle, selected by the specified kpar, and ending with all the
    // secondary generated particles
    //

    template<class stateType, class particleType>
    inline int findAndSimulate(const unsigned long long& hist,
			       pen_commonTallyCluster& tallies,
			       pen_rand& randoms,
			       const stateType& partState,
			       const pen_KPAR& kpar,
			       particleType& p){
      
      //Check if the specified kpar match with this particle 
      if(p.getKpar() == kpar){

	//Simulate the particle

	//Set particle state
	stateCopy(p.getState(),partState);
	
	//Init Page
	p.page0();

	//Try to move the generated particle to geometry system
	if(move2geo(hist,p,tallies)){

	  //Get detector ID
	  unsigned kdet = p.getDET();	
	
	  //Tally particle begins
	  tallies.run_beginPart(hist,kdet,p.getKpar(),p.readState());
    
	  //Simulate particle
	  simulatePart(hist,p,tallies,randoms);

	  return errors::SUCCESS;
	}else{
	  return errors::GEOMETRY_NOT_REACHED;
	}
      }else{
	//No more particles available, return error
	return errors::KPAR_NOT_FOUND;
      }
      
    }

    template<class stateType, class particleType, class... otherTypes>
    inline typename std::enable_if<sizeof...(otherTypes) >= 1, int>::type
    findAndSimulate(const unsigned long long& hist,
		    pen_commonTallyCluster& tallies,
		    pen_rand& randoms,
		    const stateType& partState,
		    const pen_KPAR& kpar,
		    particleType& p,
		    otherTypes&... others){
      
      //Check if the specified kpar match with this particle 
      if(p.getKpar() == kpar){

	//Simulate the particle

	//Set particle state
	stateCopy(p.getState(),partState);

	//Init Page
	p.page0();

	//Try to move the generated particle to geometry system
	if(move2geo(hist,p,tallies)){

	  //Get detector ID
	  unsigned kdet = p.getDET();	
	
	  //Tally particle begins
	  tallies.run_beginPart(hist,kdet,p.getKpar(),p.readState());
    
	  //Simulate particle
	  simulatePart(hist,p,tallies,randoms);

	  return errors::SUCCESS;
	}else{
	  return errors::GEOMETRY_NOT_REACHED;
	}
      }else{
	//Check next particle
	return findAndSimulate(hist, tallies, randoms, partState, kpar, others...);
      }
      
    }

    template<class stateType, class particleType>
    inline int findAndSimulateCond(const unsigned long long& hist,
				   pen_rand& randoms,
				   const finishTypes::finishTypes& finishType,
				   const unsigned& finishValue,
				   tallyFuncType& f,
				   const stateType& partState,
				   const pen_KPAR& kpar,
				   particleType& p){
      
      //Check if the specified kpar match with this particle 
      if(p.getKpar() == kpar){

	//Simulate the particle

	//Set particle state
	stateCopy(p.getState(),partState);
	
	//Init Page
	p.page0();

	//Try to move the generated particle to geometry system
	if(move2geo(p)){

	  //Simulate particle
	  int ret = simulatePartCond(hist,p,randoms,finishType,finishValue);
	  f(p.readState(), p.getKpar(), hist, ret);

	  return errors::SUCCESS;
	}else{
	  return errors::GEOMETRY_NOT_REACHED;
	}
      }else{
	//No more particles available, return error
	return errors::KPAR_NOT_FOUND;
      }
      
    }

    template<class stateType, class particleType, class... otherTypes>
    inline typename std::enable_if<sizeof...(otherTypes) >= 1, int>::type
    findAndSimulateCond(const unsigned long long& hist,
			pen_rand& randoms,
			const finishTypes::finishTypes& finishType,
			const unsigned& finishValue,
			tallyFuncType& f,
			const stateType& partState,
			const pen_KPAR& kpar,
			particleType& p,
			otherTypes&... others){
      
      //Check if the specified kpar match with this particle 
      if(p.getKpar() == kpar){

	//Simulate the particle

	//Set particle state
	stateCopy(p.getState(),partState);
	
	//Init Page
	p.page0();

	//Try to move the generated particle to geometry system
	if(move2geo(p)){

	  //Simulate particle
	  int ret = simulatePartCond(hist,p,randoms,finishType,finishValue);
	  f(p.readState(), p.getKpar(), hist, ret);

	  return errors::SUCCESS;
	}else{
	  return errors::GEOMETRY_NOT_REACHED;
	}
      }else{
	//Check next particle
	return findAndSimulateCond(hist, randoms, finishType, finishValue, f,
				   partState, kpar, others...);
      }
      
    }

    template<class stateType, class... particleTypes>
    inline int simulateShowerByKpar(const unsigned long long& hist,
				    pen_commonTallyCluster& tallies,
				    pen_rand& randoms,
				    const stateType& partState,
				    const pen_KPAR& kpar,
				    particleTypes&... particles){

      int err = findAndSimulate(hist, tallies, randoms, partState, kpar, particles...);
      if(err == errors::SUCCESS){
	//Primary particle found and simulated in the geometry system.

	//Simulate secondary particles
	simulateStacks(hist, tallies, randoms, particles...);
      }
      //Return same error codes as "findAndSimulate"
      return err;      
    }


    template<class stateType, class... particleTypes>
    inline int simulateShowerByKparCond(const unsigned long long& hist,
					pen_rand& randoms,
					const finishTypes::finishTypes& finishType,
					const unsigned& finishValue,
					tallyFuncType& f,					
					const stateType& partState,
					const pen_KPAR& kpar,
					particleTypes&... particles){

      //Ensure particle state derives from base state
      static_assert(std::is_base_of<pen_particleState, stateType>::value,
		    "simulateShowerByKparCond: particle state must derive"
		    "from 'pen_particleState'");
      
      int err = findAndSimulateCond(hist, randoms, finishType, finishValue, f,
				    partState, kpar, particles...);
      if(err == errors::SUCCESS){
	//Primary particle found and simulated in the geometry system.

	//Simulate secondary particles
	simulateStacksCond(hist, randoms, finishType, finishValue, f, particles...);
      }
      //Return same error codes as "findAndSimulate"
      return err;      
    }

    

    //Sample and simulate function family
    //
    // These functions sample a initial state and simulate the complete shower
    //

    //Functions to find kpar
    template<class particleType>
    inline bool findKpar(const pen_KPAR& kpar,
			 const particleType& p){
      if(kpar == p.getKpar()){
	return true;
      }
      return false;
    }

    template<class particleType, class... otherTypes>
    inline typename std::enable_if<sizeof...(otherTypes) >= 1, bool>::type
    findKpar(const pen_KPAR& kpar,
	     const particleType& p,
	     const otherTypes&... others){
      if(kpar == p.getKpar()){
	return true;
      }
      return findKpar(kpar,others...);
    }

    //Functions to set particle stacks in the tally cluster
    template<class particleType>
    inline void setTallyStacks(pen_commonTallyCluster& tallies,
			       const particleType& p){
      //Set the particle stack in the tally cluster
      tallies.setStack(p.getKpar(), &(p.readStack()));
    }
    
    template<class particleType, class... otherTypes>
    inline typename std::enable_if<sizeof...(otherTypes) >= 1, void>::type
    setTallyStacks(pen_commonTallyCluster& tallies,
		   const particleType& p,
		   const otherTypes&... others){
      
      //Set the particle stack in the tally cluster
      tallies.setStack(p.getKpar(), &(p.readStack()));

      //Go to next particle
      setTallyStacks(tallies, others...);
    }
    
    template<class particleType, class... otherTypes>
    inline const typename particleType::typeContext&
    readContext(const particleType& p,
		const otherTypes&... /*others*/){
      return p.readContext();
    }

    //Simulation function
    template<class stateType, class... particleTypes>
    int sampleAndSimulate(simConfig& config,
			  pen_specificStateGen<stateType>& source,
			  pen_commonTallyCluster& tallies,
			  particleTypes&... particles){

      //Copy simulation thread
      const unsigned ithread = config.iThread;

      //Copy verbose level
      const unsigned verbose = config.verbose;

      //Get simulation status
      simState& status = config.getStatus();
      
      //Register the start of that thread at the source task
      //to get a history assignation
      int errStart = source.workerStart(ithread,verbose);
      if(errStart != 0){
	//Unable to start worker
	if(verbose > 0){
	  config << config.threadPrefix()
		 << "Error: Unable to start worker " << ithread
		 << " in source " << source.name << simConfig::endl;
	}

	//Init and finish the simulation
	config.nextSource(source.name, 0);
	config.finish();
	
	status.setFlag(simFlags::ERROR_UNABLE_TO_START_WORKER);
	return simFlags::ERROR_UNABLE_TO_START_WORKER;
      }

      //Get source assigned histories for this thread
      unsigned long long nhists = source.toDo(ithread);

      //Init this source simulation
      config.nextSource(source.name, nhists);
      
      //Copy total number of simulated histories
      const unsigned long long simulated = config.getTotalSimulated();      
      
      //Get last simulated history
      unsigned long long hist = simulated;
      
      //Check if there are histories to simulate
      if(nhists < 1){
	//No histories required to be simulated
	if(verbose > 1){
	  config << config.threadAndSourcePrefix()
		 << "Nothing to do" << simConfig::endl;
	}
	
	status.setFlag(simFlags::NOTHING_TO_DO);
	return simFlags::NOTHING_TO_DO;
      }
      
      //Calculate last history to simulate
      unsigned long long lastHist = simulated+nhists;

      //Read context
      const auto& context = readContext(particles...);
      
      //Create random generator
      pen_rand random;
      config.setInitSeedsToRand(random);

      //Create a stopwatch for simulation time
      pen_stopWatch simLimitWatch(config.maxSimTime);
      
      //Create a stopwatch for dumps
      pen_stopWatch dumpWatch(config.dumpTime);

      //Create a stopWatch to perform checkpoints
      pen_stopWatch checkPointWatch(source.task.getCheckTime()*1000ll);  

      //Create a stopWatch to perform speed reports
      long long int mili2report =
	static_cast<long long int>(1000*source.task.getCheckTime()/1.5);
      pen_stopWatch reportWatch(mili2report);
      
      //Start all stopWatch counts
      if(ithread == 0)
	checkPointWatch.start();
      reportWatch.start();
      dumpWatch.start();
      simLimitWatch.start();

      //Set stacks in the tally cluster
      setTallyStacks(tallies, particles...);

      //Flag simulation to running
      status.setFlag(simFlags::RUNNING);      
      
      //** Generate first state
      stateType genState;
      pen_KPAR genKpar;
      unsigned long long dhist;

      //Update tallies last hist
      tallies.run_lastHist(hist);
  
      //Sample first particle
      source.sample(genState,genKpar,dhist,ithread,random);

      //Check if sampling returns a valid kpar and if the new hist
      //is under the history limit number
      if(!findKpar(genKpar, particles...)){
	//Invalid kpar

	//Update seeds for next source
	config.setFinalSeedsFromRand(random);

	if(verbose > 0){
	  config << config.threadAndSourcePrefix()
		 << "Provided source generated a unknown kpar "
	    "on first call. Simulation finish due empty source."
		 << simConfig::endl;
	}
	
	status.setFlag(simFlags::END_OF_SOURCE);
	config.finish();	
	return simFlags::END_OF_SOURCE;
      }
      if(hist + dhist > lastHist){
	//Simulation history limit reached

	//Consider "skipped" histories as simulated
	hist = lastHist;

	//Report simulated histories
	config.report(hist);

	//Update seeds for next source
	config.setFinalSeedsFromRand(random);

	if(verbose > 0){
	  config << config.threadAndSourcePrefix()
		 << "First history increment is greater than the number of "
	    "histories to be simulated. Simulation finished."
		 << simConfig::endl;
	}
	
	//Finish the source simulation
	status.setFlag(simFlags::SUCCESS);
	config.finish();
	
	return simFlags::SUCCESS;
      }

      //Increment history counter
      hist += dhist;

      //Get detector ID
      unsigned firstKdet = context.getDET(genState.IBODY);
      //Tally first sampled particle
      tallies.run_sampledPart(hist,dhist,firstKdet,genKpar,genState);

      //Enter history loop
      unsigned simConfRepCount = 0;
      for(;;){

	int err = simulateShowerByKpar(hist, tallies, random,
				       genState, genKpar, particles...);
	if(err != errors::SUCCESS &&
	   err != errors::GEOMETRY_NOT_REACHED){
	  //Unknown particle ID signals end of source
	  status.setFlag(simFlags::END_OF_SOURCE);

	  if(verbose > 1){
	    config << config.threadAndSourcePrefix()
		   << "End of source reached." << simConfig::endl;
	  }
	  
	  //End current hist
	  tallies.run_endHist(hist);
	  break;
	}

	//Get last seeds
	int lseed1, lseed2;
	random.getSeeds(lseed1,lseed2);
    
	//Sample new particle state
	source.sample(genState,genKpar,dhist,ithread,random);
    
	//Check if is a valid particle
	if(genKpar >= ALWAYS_AT_END){

	  //End of source
	  status.setFlag(simFlags::END_OF_SOURCE);

	  if(verbose > 1){
	    config << config.threadAndSourcePrefix()
		   << "End of source reached." << simConfig::endl;
	  }
	  
	  //End last hist
	  tallies.run_endHist(hist);      
	  break;
	}


	if(dhist > 0){
	  //End previous history
	  tallies.run_endHist(hist);

	  //Get actual time
	  const auto tnow = std::chrono::steady_clock::now();
	  //Check if elapsed time reaches dump interval
	  if(dumpWatch.check(tnow)){

	    //Perform a configure report
	    config.report(hist);
	    
	    //Save dump
	    unsigned long long currentHists =
	      hist-simulated+config.getInitiallySimulated();

	    if(verbose > 1){
	      config << config.threadAndSourcePrefix()
		     << "Dump simulation at history number " << hist
		     << " with " << currentHists << "/" << nhists
		     << " simulated in actual source, with seeds "
		     << lseed1 << " " << lseed2 << simConfig::endl;
	    }
	    
	    //Dump sim
	    tallies.dump2file(config.dumpFilename.c_str(),
			      hist,
			      lseed1,lseed2,
			      config.getCurrentSourceIndex(),
			      currentHists,
			      3);
	    if(config.writePartial)
	      tallies.saveData(hist,false); //Saves data but does not repeat flush calls
	    dumpWatch.start(); //Restart watch
	  }
	  //Check if is time to make a report
	  if(reportWatch.check(tnow)){
	
	    //Report simulated histories
	    int errReport;
	    std::chrono::seconds::rep nextReport =
	      source.report(ithread,hist-simulated,&errReport,config.verbose);

	    if(nextReport < 0){
	      //Error during report
	      if(config.verbose > 0){
		config << config.threadAndSourcePrefix()
		       << "Error reporting simulated histories. "
		       << "Error code: " << errReport << simConfig::endl;
	      }
	    }
	    else{
	      //Update wait time until next report
	      if(config.verbose > 2){
		config << config.threadAndSourcePrefix()
		       << "Schedule next report "
		       << nextReport << "s" << simConfig::endl;
	      }
	      reportWatch.duration(nextReport*1000);
	    }

	    //Update assigned histories
	    updateHists(nhists,source,ithread,config.verbose);
	    lastHist = simulated+nhists;
	    reportWatch.start(); //Restart watch

	    //Update histories to simulate for this source in configuration
	    config.updateToSimulate(nhists);
	  }
	  if(ithread == 0){
	    //Check if is time to make a checkpoint
	    if(checkPointWatch.check(tnow)){
	      checkpoint(source,ithread,config.verbose);
	      checkPointWatch.start(); //Restart watch
	    }
	  }
	  //Check if the simulation must be stopped
	  if(simLimitWatch.check(tnow)){
	      
	    //Save dump
	    unsigned long long currentHists =
	      hist-simulated+config.getInitiallySimulated();

	    if(verbose > 1){
	      config << config.threadAndSourcePrefix()
		     << "Simulation time limit reached ("
		     << config.maxSimTime/1000 << "). "
		     << "Dump simulation at history number " << hist
		     << " with " << currentHists << "/" << nhists
		     << " simulated in actual source, with seeds "
		     << lseed1 << " " << lseed2 << simConfig::endl;
	    }
	    
	    //Dump sim
	    tallies.dump2file(config.dumpFilename.c_str(),
			      hist,
			      lseed1,lseed2,
			      config.getCurrentSourceIndex(),
			      currentHists,
			      3);
	    
	    if(config.writePartial)
	      tallies.saveData(hist,false); //Saves data but does not repeat flush calls
	    
	    //Finish the simulation
	    status.setFlag(simFlags::MAX_TIME_REACHED);
	    break;
	  }

	  //Check if a report is needed in the configuration
	  if(++simConfRepCount % 1000 == 0){
	    config.report(hist);
	    simConfRepCount = 0;
	  }
	    	  
	  //Increment history counter
	  hist += dhist;

	  //Check history limit
	  if(hist >= lastHist){
	    //End of simulation, ask permission to finish
	    if(source.handleFinish(ithread,hist-simulated,nhists,config.verbose)){
	      //This worker can finish, exit from the simulation loop
	      hist = lastHist;
	      status.setFlag(simFlags::SUCCESS);	      
	      break;
	    }
	    //Restart report watch
	    reportWatch.start();
	    lastHist = simulated+nhists;
	  }

	  //Check finish simulation function
	  if(!config.fSimFinish(hist)){
	    //End of simulation
	    if(verbose > 1){
	      unsigned long long currentHists =
		hist-simulated+config.getInitiallySimulated();
	      
	      config << config.threadAndSourcePrefix()
		     << "Simulation finish condition reached "
		     << "at history number " << hist
		     << " with " << currentHists << "/" << nhists
		     << " simulated in actual source, with seeds "
		     << lseed1 << " " << lseed2 << simConfig::endl;
	    }
	    
	    status.setFlag(simFlags::SUCCESS);
	    break;
	  }	  
	}

	//Get detector ID
	unsigned kdet = context.getDET(genState.IBODY);
	//Tally new sampled particle
	tallies.run_sampledPart(hist,dhist,kdet,genKpar,genState);    
      }

      //Update seeds for next source
      config.setFinalSeedsFromRand(random);

      //Perform a final report
      config.report(hist);

      //Finish source simulation
      config.finish();

      if(verbose > 1){
	int seed1, seed2;
	config.getSeeds(seed1, seed2);
	config << config.threadAndSourcePrefix()
	       << "Simulation finished at history "
	       << config.getTotalSimulated()
	       << " with " << config.getSimulatedInSource()-config.getInitiallySimulated()
	       << "/"
	       << config.getToSimulateInSource() << " histories simulated and "
	       << seed1 << " " << seed2 << " as final seeds. Total time: "
	       << config.simTime() << simConfig::endl;
      }
      
      //Update maximum simulation time
      config.maxSimTime -= static_cast<long long int>(config.simTime()*1000.0);

      return simFlags::SUCCESS;
    }

    //Simulation function
    template<class stateType, class... particleTypes>
    int sampleAndSimulateCond(simConfig& config,
			      const unsigned long long nhists,
			      const std::string& sourceName,
			      const sampleFuncType<stateType>& fsource,
			      const finishTypes::finishTypes& finishType,
			      const unsigned& finishValue,
			      tallyFuncType& ftally,
			      particleTypes&... particles){

      //Copy simulation thread
      const unsigned ithread = config.iThread;

      //Copy verbose level
      const unsigned verbose = config.verbose;

      //Get simulation status
      simState& status = config.getStatus();

      //Init this source simulation
      config.nextSource(sourceName, nhists);
      
      //Copy total number of simulated histories
      const unsigned long long simulated = config.getTotalSimulated();      
      
      //Get last simulated history
      unsigned long long hist = simulated;
      
      //Check if there are histories to simulate
      if(nhists < 1){
	//No histories required to be simulated
	if(verbose > 1){
	  config << config.threadAndSourcePrefix()
		 << "Nothing to do" << simConfig::endl;
	}
	
	status.setFlag(simFlags::NOTHING_TO_DO);
	return simFlags::NOTHING_TO_DO;
      }
      
      //Calculate last history to simulate
      unsigned long long lastHist = simulated+nhists;

      //Read context
      const auto& context = readContext(particles...);
      
      //Create random generator
      pen_rand random;
      config.setInitSeedsToRand(random);

      //Create a stopwatch for simulation time
      pen_stopWatch simLimitWatch(config.maxSimTime);
      
      //Start all stopWatch counts
      simLimitWatch.start();

      //Flag simulation to running
      status.setFlag(simFlags::RUNNING);      
      
      //** Generate first state
      stateType genState;
      pen_KPAR genKpar;
      unsigned long long dhist;
  
      //Sample first particle
      fsource(genState,genKpar,dhist,ithread,random);

      //Check if sampling returns a valid kpar and if the new hist
      //is under the history limit number
      if(!findKpar(genKpar, particles...)){
	//Invalid kpar

	//Update seeds for next source
	config.setFinalSeedsFromRand(random);

	if(verbose > 0){
	  config << config.threadAndSourcePrefix()
		 << "Provided source generated a unknown kpar "
	    "on first call. Simulation finish due empty source."
		 << simConfig::endl;
	}
	
	status.setFlag(simFlags::END_OF_SOURCE);
	config.finish();	
	return simFlags::END_OF_SOURCE;
      }

      //Locate particle in geometry
      context.readGeometry()->locate(genState);
      
      if(hist + dhist > lastHist){
	//Simulation history limit reached

	//Consider "skipped" histories as simulated
	hist = lastHist;

	//Report simulated histories
	config.report(hist);

	//Update seeds for next source
	config.setFinalSeedsFromRand(random);

	if(verbose > 0){
	  config << config.threadAndSourcePrefix()
		 << "First history increment is greater than the number of "
	    "histories to be simulated. Simulation finished."
		 << simConfig::endl;
	}
	
	//Finish the source simulation
	status.setFlag(simFlags::SUCCESS);
	config.finish();
	
	return simFlags::SUCCESS;
      }

      //Increment history counter
      hist += dhist;

      //Enter history loop
      unsigned simConfRepCount = 0;
      for(;;){

	int err = simulateShowerByKparCond(hist, random, finishType, finishValue, ftally,
					   genState, genKpar, particles...);
	if(err != errors::SUCCESS &&
	   err != errors::GEOMETRY_NOT_REACHED){
	  //Unknown particle ID signals end of source
	  status.setFlag(simFlags::END_OF_SOURCE);

	  if(verbose > 1){
	    config << config.threadAndSourcePrefix()
		   << "End of source reached." << simConfig::endl;
	  }	  
	  break;
	}

	//Get last seeds
	int lseed1, lseed2;
	random.getSeeds(lseed1,lseed2);
    
	//Sample new particle state
	fsource(genState,genKpar,dhist,ithread,random);
    
	//Check if is a valid particle
	if(genKpar >= ALWAYS_AT_END){

	  //End of source
	  status.setFlag(simFlags::END_OF_SOURCE);

	  if(verbose > 1){
	    config << config.threadAndSourcePrefix()
		   << "End of source reached." << simConfig::endl;
	  }
	  break;
	}

	//Locate particle in geometry
	context.readGeometry()->locate(genState);

	if(dhist > 0){

	  //Get actual time
	  const auto tnow = std::chrono::steady_clock::now();
	  //Check if the simulation must be stopped
	  if(simLimitWatch.check(tnow)){
	      
	    //Save dump
	    unsigned long long currentHists =
	      hist-simulated+config.getInitiallySimulated();

	    if(verbose > 1){
	      config << config.threadAndSourcePrefix()
		     << "Simulation time limit reached ("
		     << config.maxSimTime/1000 << "). "
		     << "at history number " << hist
		     << " with " << currentHists << "/" << nhists
		     << " simulated in actual source, with seeds "
		     << lseed1 << " " << lseed2 << simConfig::endl;
	    }
	    
	    //Finish the simulation
	    status.setFlag(simFlags::MAX_TIME_REACHED);
	    break;	
	  }

	  //Check if a report is needed in the configuration
	  if(++simConfRepCount % 1000 == 0){
	    config.report(hist);
	    simConfRepCount = 0;
	  }
	    	  
	  //Increment history counter
	  hist += dhist;

	  //Check history limit
	  if(hist >= lastHist){
	    //End of simulation
	    hist = lastHist;
	    status.setFlag(simFlags::SUCCESS);	      
	    break;
	  }

	  //Check finish simulation function
	  if(!config.fSimFinish(hist)){
	    //End of simulation
	    if(verbose > 1){
	      unsigned long long currentHists =
		hist-simulated+config.getInitiallySimulated();
	      
	      config << config.threadAndSourcePrefix()
		     << "Simulation finish condition reached "
		     << "at history number " << hist
		     << " with " << currentHists << "/" << nhists
		     << " simulated in actual source, with seeds "
		     << lseed1 << " " << lseed2 << simConfig::endl;
	    }
	    
	    status.setFlag(simFlags::SUCCESS);
	    break;
	  }
	}
      }

      //Update seeds for next source
      config.setFinalSeedsFromRand(random);

      //Perform a final report
      config.report(hist);

      //Finish source simulation
      config.finish();

      if(verbose > 1){
	int seed1, seed2;
	config.getSeeds(seed1, seed2);
	config << config.threadAndSourcePrefix()
	       << "Simulation finished at history "
	       << config.getTotalSimulated()
	       << " with " << config.getSimulatedInSource()-config.getInitiallySimulated()
	       << "/"
	       << config.getToSimulateInSource() << " histories simulated and "
	       << seed1 << " " << seed2 << " as final seeds. Total time: "
	       << config.simTime() << simConfig::endl;
      }
      
      //Update maximum simulation time
      config.maxSimTime -= static_cast<long long int>(config.simTime()*1000.0);

      return simFlags::SUCCESS;
    }

    //Sample and simulate function family by context
    //
    // These functions sample a initial state and simulate the complete shower
    // using directly the context informaiton. The associated particle instances
    // are created automatically.
    //
    
    template<class stateType, class contextType,
	     class... vrStateTypes,
	     size_t... I>
    inline int sampleAndSimulateContext(simConfig& config,
					const contextType& context,
					const std::index_sequence<I...>&,
					pen_specificStateGen<stateType>& source,
					pen_commonTallyCluster& tallies,
					const abc_VR<vrStateTypes>&... vr){

      //Instantiate context particles
      std::unique_ptr<penred::context::particles<contextType>> particles(new penred::context::particles<contextType>(context));

      //Register VR in particles
      particles->registerVR(vr...);

      return sampleAndSimulate(config,
			       source,
			       tallies,
			       particles->template getParticle<I>()...);      
    }

    
    template<class stateType, class contextType, class... vrStateTypes>
    inline int sampleAndSimulateContext(simConfig& config,
					const contextType& context,
					pen_specificStateGen<stateType>& source,
					pen_commonTallyCluster& tallies,
					const abc_VR<vrStateTypes>&... vr){

      return sampleAndSimulateContext(config,context,
				      std::make_index_sequence<penred::context::particles<contextType>::NTypes>{},
				      source,
				      tallies,
				      vr...);
    }

    
    template<class stateType, class contextType,
	     class... vrStateTypes,
	     size_t... I>
    inline int sampleAndSimulateCondContext(simConfig& config,
					    const contextType& context,					    
					    const std::index_sequence<I...>&,
					    const unsigned long long nhists,
					    const std::string& sourceName,
					    const sampleFuncType<stateType>& fsource,
					    const finishTypes::finishTypes& finishType,
					    const unsigned& finishValue,
					    tallyFuncType& ftally,
					    const abc_VR<vrStateTypes>&... vr){
      
      //Instantiate context particles
      std::unique_ptr<penred::context::particles<contextType>> particles(new penred::context::particles<contextType>(context));

      //Register VR in particles
      particles->registerVR(vr...);

      return sampleAndSimulateCond(config, nhists,
				   sourceName, fsource,
				   finishType, finishValue, ftally,
				   particles->template getParticle<I>()...);
    }

    template<class stateType, class contextType, class... vrStateTypes>
    inline int sampleAndSimulateCondContext(simConfig& config,
					    const contextType& context,					    
					    const unsigned long long nhists,
					    const std::string& sourceName,
					    const sampleFuncType<stateType>& fsource,
					    const finishTypes::finishTypes& finishType,
					    const unsigned& finishValue,
					    tallyFuncType& ftally,
					    const abc_VR<vrStateTypes>&... vr){

      return sampleAndSimulateCondContext(config, context,
					  std::make_index_sequence<penred::context::particles<contextType>::NTypes>{},
					  nhists,
					  sourceName,
					  fsource,
					  finishType,
					  finishValue,
					  ftally,
					  vr...);
    }
    
  } // namespace simulation
} // namespace penred
#endif
