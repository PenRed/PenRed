
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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



#ifndef __PARTICLE_GENERATION__
#define __PARTICLE_GENERATION__

#include <cmath>
#include <thread>
#include <atomic>
#include <chrono>
#include "pen_classes.hh"
#include "pen_math.hh"
#include "instantiator.hh"
#include "pen_states.hh"
#include "loadBalance.hh"

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
#include "mpi.h"
#endif
  // ***************************** MPI END ********************************** //

#define DECLARE_SAMPLER(Class) \
  public: \
  inline int registerStatus() const { return ___register_return;} \
  virtual const char* readID() const {return ___ID;}\
  private: \
  static const char* ___ID;\
  static const int ___register_return;\

#define REGISTER_SAMPLER(Class, ID) \
  const int Class::___register_return = registerSampler<Class>(static_cast<const char *>(#ID)); \
  const char* Class::___ID = static_cast<const char *>(#ID);

#define REGISTER_SPECIFIC_SAMPLER(Class, State, ID)				\
  const int Class::___register_return = registerSpecificSampler<Class,State>(static_cast<const char *>(#ID)); \
  const char* Class::___ID = static_cast<const char *>(#ID);

template <class particleState> class abc_specificSampler;
class abc_spatialSampler;
class abc_directionSampler;
class abc_energySampler;
class abc_timeSampler;
class pen_genericStateGen;
template <class particleState> class pen_specificStateGen;

enum __usedSamp{
	      USE_SPATIAL     = 1 << 0,
	      USE_DIRECTION   = 1 << 1,
	      USE_ENERGY      = 1 << 2,
	      USE_TIME        = 1 << 3,
	      USE_GENERIC     = 1 << 4,
	      USE_NONE        = 1 << 5		
};

struct pen_samplerTask{
  std::atomic<unsigned long long> iterDone;
  unsigned nworkers;
  unsigned long long toDoWorker;
  unsigned long long toDoMPI;
  unsigned long long toDo;

  pen_samplerTask() : iterDone(0), nworkers(0),
		      toDoWorker(0), toDoMPI(0), toDo(0){
  }
  
  inline int init(const size_t nw,
		  const unsigned long long nIter,
		  const char* /*logFileName*/ = nullptr,
		  const unsigned /*verbose*/ = 0){
    iterDone = 0;
    nworkers = nw;
    toDo = nIter;

    toDoMPI = toDo;
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  int mpiSize;
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  toDoMPI /= static_cast<unsigned long long>(mpiSize);
  toDoMPI += 1;
#endif
  // ******************************* MPI ************************************ //
    toDoWorker = toDoMPI/nw+1;
    return 0;
  }
  inline void skip(unsigned long long toSkip){
    unsigned long long auxIterDone = iterDone;
    if(toDo > toSkip)
      init(nworkers,toDo-toSkip);
    else
      init(nworkers,0);
    iterDone = auxIterDone;
  }
  inline unsigned long long done() const {return iterDone;}
  inline void reset(){
    iterDone = toDoWorker = toDoMPI = toDo = 0;
    nworkers = 0;
  }
  inline unsigned long long int assigned(const size_t /*iw*/){return toDoWorker;}
  inline unsigned long long int assigned(){return toDoMPI;}
  inline int workerStart(const size_t /*iw*/,
			 const unsigned /*verbose*/){return 0;}
  inline bool workerFinish(const size_t /*iw*/, int& /*why*/,
			   const unsigned /*verbose*/) const {
    return true;
  }
  inline std::chrono::seconds::rep 
  report(const size_t /*iw*/,
	 const unsigned long long nIter,
	 int* /*errRet*/,
	 const unsigned /*verbose*/) {
    iterDone += nIter;
    return std::chrono::seconds::rep(1000000000);
  }
  inline int checkPoint(const unsigned /*verbose*/) const {return 0;}
  inline std::chrono::seconds::rep
  getCheckTime() const {return std::chrono::seconds::rep(1000000000);}
  inline void setCheckTime(const std::chrono::seconds::rep) const {}
  inline void setCheckTimeMPI(const std::chrono::seconds::rep) const {}
  inline void minTime(const unsigned long long) const {}
};

template<class particleState>
class abc_specificSampler{

private:
  
  const __usedSamp usedSamplings;
  const abc_spatialSampler* pSpatial;
  const abc_directionSampler* pDirection;
  const abc_energySampler* pEnergy;
  const abc_timeSampler* pTime;
  const wrapper_geometry* geometry;

  unsigned nthread;
  
protected:
public:

  abc_specificSampler(__usedSamp usedSamplingsIn) : usedSamplings(usedSamplingsIn),
						    pSpatial(nullptr),
						    pDirection(nullptr),
						    pEnergy(nullptr),
						    pTime(nullptr)
  {}

  inline const abc_spatialSampler* spatial() const {return pSpatial;}
  inline const abc_directionSampler* direction() const {return pDirection;}
  inline const abc_energySampler* energy() const {return pEnergy;}
  inline const abc_timeSampler* time() const {return pTime;}
  inline const wrapper_geometry* geo() const {return geometry;}
  inline unsigned getThread() const {return nthread;}
  inline void setThread(const unsigned threadNum){nthread = threadNum;}
  
  inline void setSpatial(const abc_spatialSampler* newSpatial){
    pSpatial = newSpatial;
    updateSamplers(); //Handle sampler change
  }
  inline void setDirection(const abc_directionSampler* newDirection){
    pDirection = newDirection;
    updateSamplers(); //Handle sampler change
  }
  inline void setEnergy(const abc_energySampler* newEnergy){
    pEnergy = newEnergy;
    updateSamplers(); //Handle sampler change
  }
  inline void setTime(const abc_timeSampler* newTime){
    pTime = newTime;
    updateSamplers(); //Handle sampler change
  }

  inline void setGeometry(const wrapper_geometry* geometryIn){
    geometry = geometryIn;
  }

  //Virtual method to handle sapler changes
  virtual void updateSamplers(){}
  
  inline bool usesSpatial() const{
    return ((usedSamplings & USE_SPATIAL) == USE_SPATIAL);
  }
  inline bool usesDirection() const{
    return ((usedSamplings & USE_DIRECTION) == USE_DIRECTION);
  }
  inline bool usesEnergy() const{
    return ((usedSamplings & USE_ENERGY) == USE_ENERGY);
  }
  inline bool usesTime() const{
    return ((usedSamplings & USE_TIME) == USE_TIME);
  }
  inline bool usesGeneric() const{
    return ((usedSamplings & USE_GENERIC) == USE_GENERIC);
  }
  
  virtual void sample(particleState& state,
		      pen_KPAR& genKpar,
		      unsigned long long& dhist,
		      pen_rand& random) = 0;  
  
  virtual void skip(const unsigned long long /*dhists*/){}
  
  virtual int configure(double&,
			const abc_spatialSampler*,
			const abc_directionSampler*,
			const abc_energySampler*,
			const abc_timeSampler*,
			const pen_parserSection&,
			const unsigned = 0) = 0;

  inline static const char* type() {return "SPECIFIC";}
  virtual const char* readID() const = 0;
  
  virtual ~abc_specificSampler(){}
};

class abc_spatialSampler{

private:
  
protected:
  double rotation[9];
  bool rotate;
  virtual void geoSampling(double pos[3], pen_rand& random) const = 0;
public:
  double translation[3];

  abc_spatialSampler();
    
  void sample(pen_particleState& state, pen_rand& random) const;

  virtual int configure(const pen_parserSection&, const unsigned = 0) = 0;
  inline void setRotationZYZ(const double omega, const double theta, const double phi){
    createRotationZYZ(omega,theta,phi,rotation);
  }

  inline static const char* type() {return "SPATIAL";}
  virtual const char* readID() const = 0;
  
  virtual ~abc_spatialSampler(){}
};

class abc_directionSampler{

private:
  
protected:

  virtual void directionSampling(double dir[3], pen_rand& random) const = 0;
public:

  abc_directionSampler();
  
  void sample(pen_particleState& state, pen_rand& random) const;

  virtual int configure(const pen_parserSection&, const unsigned = 0) = 0;

  inline static const char* type() {return "DIRECTION";}
  virtual const char* readID() const = 0;
  
  virtual ~abc_directionSampler(){}
};

class abc_energySampler{

private:
  
protected:

  virtual void energySampling(double& Energy, pen_rand& random) const = 0;
public:

  double maxEnergy;      //Maximum energy allowed by physic models used (eV)
  double minimumEnergy;  //Minimum energy allowed by physic models used (eV)

  abc_energySampler();
  
  void sample(pen_particleState& state, pen_rand& random) const;

  virtual int configure(double&, const pen_parserSection&, const unsigned = 0) = 0;

  inline static const char* type() {return "ENERGY";}
  virtual const char* readID() const = 0;
  
  virtual ~abc_energySampler(){}
};

class abc_timeSampler{

private:
  
protected:

  virtual void timeSampling(double& time, pen_rand& random) const = 0;
public:

  abc_timeSampler();	
  
  void sample(pen_particleState& state, pen_rand& random) const;

  virtual int configure(const pen_parserSection&, const unsigned = 0) = 0;

  inline static const char* type() {return "TIME";}
  virtual const char* readID() const = 0;
  
  virtual ~abc_timeSampler(){}
};

class pen_genericStateGen{

  template<class stateType> friend class pen_specificStateGen; 
  
private:

  //Use the "Construct Members On First Use Idiom" to avoid the static
  //initialization order fiasco
  static instantiator<abc_spatialSampler>& spatialSamplers();
  static instantiator<abc_directionSampler>& directionSamplers();
  static instantiator<abc_energySampler>& energySamplers();
  static instantiator<abc_timeSampler>& timeSamplers();

  abc_spatialSampler* spatialSampler;
  abc_directionSampler* directionSampler;
  abc_energySampler* energySampler;
  abc_timeSampler* timeSampler;
  
  const wrapper_geometry* geometry;
  
  int configStatus;
  double Emax; //Maximum possible value for sampled energies
  
public:

  int sourceBody; //Body index source
  unsigned sourceMat;  //Body material source
  
  std::string name; //Source name
  bool LAGE; //Age recording status
  pen_KPAR kpar; //Generated particles kpar

// **************************** LB ********************************* //
#ifdef _PEN_USE_LB_
  typedef LB::task taskType;
#else
  typedef pen_samplerTask taskType;
#endif
// **************************** LB ********************************* //

  taskType task;

  pen_genericStateGen();

  inline int initTask(const size_t nw,
		      const unsigned long long nIter,
		      const unsigned verbose){
    return task.init(nw,nIter,name.c_str(),verbose);
  }

  
  // ******************************* MPI ************************************ //
#if defined(_PEN_USE_MPI_) && defined(_PEN_USE_LB_)
  inline int initTask(const size_t nw,
		      const unsigned long long nIter,
		      const MPI_Comm commIn,
		      const int tagReqests,
		      const int tagProcess,
		      const unsigned verbose){
    return task.init(nw,nIter,commIn,tagReqests,tagProcess,
		     name.c_str(),verbose);
  }
#endif
  // ***************************** MPI END ********************************** //

  inline void skip(const unsigned long long nIter){task.skip(nIter);}

  inline unsigned long long toDo(unsigned ithread){
    return task.assigned(ithread);
  }
  inline unsigned long long toDo(){
    return task.assigned();
  }

  inline std::chrono::seconds::rep 
  report(const size_t iw,
	 const unsigned long long nIter,
	 int* errRet,
	 const unsigned verbose){
    return task.report(iw,nIter,errRet,verbose);
  }

  inline int checkPoint(const unsigned verbose){
    return task.checkPoint(verbose);
  }
  inline int workerStart(const size_t iw,
			 const unsigned verbose){
    return task.workerStart(iw,verbose);
  }
  inline bool workerFinish(const size_t iw, int& why, const unsigned verbose){
    return task.workerFinish(iw,why,verbose);
  }

  bool handleFinish(const unsigned iw,
		    const unsigned long long nDone,
		    unsigned long long& assigned,
		    const unsigned verbose);  
  inline void setCheckTime(const std::chrono::seconds::rep t){
    task.setCheckTime(t);
  }
  inline void setLBthreshold(const unsigned long long t) {
    task.minTime(t);
  }
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  inline void setCheckTimeMPI(const std::chrono::seconds::rep t){
    task.setCheckTimeMPI(t);
  }
#endif
  // ***************************** MPI END ********************************** //
  
  inline const abc_spatialSampler* spatial() const {
    return spatialSampler;}
  inline const abc_directionSampler* direction() const {
    return directionSampler;}
  inline const abc_energySampler* energy() const {
    return energySampler;}
  inline const abc_timeSampler* time() const {
    return timeSampler;}
  
  static std::string samplersList();
  static std::string samplersList(std::vector<std::string>& spatial,
				  std::vector<std::string>& direction,
				  std::vector<std::string>& energy,
				  std::vector<std::string>& time);
  int setGeometry(const wrapper_geometry* geometryIn);

  inline int configureStatus(){return configStatus;}
  inline double maxEnergy(){return Emax;}
  
  template <class subclass>
  static int addSpatialSampler(const char* typeID){
    return spatialSamplers().addSubType<subclass>(typeID);
  }

  template <class subclass>
  static int addDirectionSampler(const char* typeID){
    return directionSamplers().addSubType<subclass>(typeID);
  }
  
  template <class subclass>
  static int addEnergySampler(const char* typeID){
    return energySamplers().addSubType<subclass>(typeID);
  }
  
  template <class subclass>
  static int addTimeSampler(const char* typeID){
    return timeSamplers().addSubType<subclass>(typeID);
  }
  
  int selectSpatialSampler(const char* ID,
			   const pen_parserSection& config,
			   const unsigned verbose = 0);

  int selectDirectionSampler(const char* ID,
			     const pen_parserSection& config,
			     const unsigned verbose = 0);
  
  int selectEnergySampler(const char* ID,
			  const pen_parserSection& config,
			  const unsigned verbose = 0);

  int selectTimeSampler(const char* ID,
			const pen_parserSection& config,
			const unsigned verbose = 0);    

  void sample(pen_particleState& state,
	      pen_rand& random) const;

  void configure(const pen_parserSection& config,
		 const unsigned verbose = 0);
  
  void clear();

  inline const char* spatialID(){
    if(spatialSampler != nullptr)
      return spatialSampler->readID();
    return nullptr;
  }

  inline const char* directionID(){
    if(directionSampler != nullptr)
      return directionSampler->readID();
    return nullptr;
  }

  inline const char* energyID(){
    if(energySampler != nullptr)
      return energySampler->readID();
    return nullptr;
  }

  inline const char* timeID(){
    if(timeSampler != nullptr)
      return timeSampler->readID();
    return nullptr;
  }
  
  virtual ~pen_genericStateGen(){clear();}
};

template <class particleState>
class pen_specificStateGen{
  
private:

  pen_genericStateGen genericGen;
  static instantiator<abc_specificSampler<particleState>>& specificSamplers(){
    static instantiator<abc_specificSampler<particleState>>* ans =
      new instantiator<abc_specificSampler<particleState>>;
    return *ans;
  }

  std::vector<abc_specificSampler<particleState>*> specificSamplerVect;

  void clearSpecificSamplers(){
    for(unsigned i = 0; i < specificSamplerVect.size(); i++){
      if(specificSamplerVect[i] != nullptr)
	delete specificSamplerVect[i];
    }

    specificSamplerVect.clear();
    useSpecific = false;
  }

  
  int configureGeneric(const pen_parserSection& config,
		       const bool specific,
		       const unsigned verbose){
    
  int err;

  //Clear generator
  clear();
  
  if(specific && verbose > 1)
    printf("Specific sampler expected.\n");
  
  //*******************************
  // Check if age must be recorded
  //*******************************

  err = config.read("record-time",LAGE);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("genericStateGen: configure: Field ('record-time') not found. Age recording disabled.\n");
    }
    LAGE = false;
  }

  if(verbose > 1){
    printf("Age recording: %s\n", LAGE ? "enabled" : "disabled");
  }
  
  //*******************
  // Samplers sections
  //*******************
    
  //Extract 'spatial' sampler field
  pen_parserSection spatialSection;
  bool useSpatial = true;
  err = config.readSubsection("spatial",spatialSection);
  if(err != INTDATA_SUCCESS){
    if(specific){
      useSpatial = false;
    }
    else{
      if(verbose > 0){
	printf("genericStateGen: configure: Error: '%s': Unable to read field 'spatial'.\n",name.c_str());
      }
      return -1;
    }
  }

  //Extract 'direction' sampler field
  pen_parserSection directionSection;
  bool useDirection = true;
  err = config.readSubsection("direction",directionSection);
  if(err != INTDATA_SUCCESS){
    if(specific){
      useDirection = false;
    }
    else{
      if(verbose > 0){
	printf("genericStateGen: configure: Error: '%s': Unable to read field 'direction'.\n",name.c_str());
      }
      return -2;
    }
  }

  //Extract 'energy' sampler field
  pen_parserSection energySection;
  bool useEnergy = true;
  err = config.readSubsection("energy",energySection);
  if(err != INTDATA_SUCCESS){
    if(specific){
      useEnergy = false;
    }
    else{
      if(verbose > 0){    
	printf("genericStateGen: configure: Error: '%s': Unable to read field 'energy'.\n",name.c_str());
      }
      return -3;
    }
  }
  
  //Extract 'time' sampler field
  pen_parserSection timeSection;
  bool UseTimeSampling = true;
  err = config.readSubsection("time",timeSection);
  if(err != INTDATA_SUCCESS){
    //No time sampler
    UseTimeSampling = false;
  }


  //**************
  // Samplers IDs
  //**************    
    
  //Get spatial sampler ID
  std::string spatialType;
  if(useSpatial){
    err = spatialSection.read("type",spatialType);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("genericStateGen:configure: Error: '%s': Unable to read field 'spatial/type'. String expected.\n",name.c_str());
      }
      return -1;
    }
  }

  //Get direction sampler ID
  std::string directionType;
  if(useDirection){
    err = directionSection.read("type",directionType);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("genericStateGen:configure: Error: '%s': Unable to read field 'direction/type'. String expected.\n",name.c_str());
      }
      return -2;
    }
  }

  //Get energy sampler ID
  std::string energyType;
  if(useEnergy){
    err = energySection.read("type",energyType);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("genericStateGen:configure: Error: '%s': Unable to read field 'energy/type'. String expected.\n",name.c_str());
      }
      return -3;
    }
  }
    
  std::string timeType;
  if(UseTimeSampling){
    //Get time sampler ID
    err = timeSection.read("type",timeType);
    if(err != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("genericStateGen:configure: Error: '%s': Unable to read field 'time/type'. String expected.\n",name.c_str());
      }
      return -4;
    }
  }

  
  //********************************
  // Create and configure samplers
  //********************************

  if(verbose > 1){
    printf("\n------------------------------------\n");
    printf("\n **** Source '%s'\n",name.c_str());
    printf("\n------------------------------------\n\n");
  }

  if(useSpatial){
    if(verbose > 1){printf("\nSpatial sampler '%s':\n\n",spatialType.c_str());}
    err = selectSpatialSampler(spatialType.c_str(),spatialSection,verbose);
    if(verbose > 1){printf("\n------------------------------------\n\n");}     
    if(err != INTDATA_SUCCESS){
      clear();
      return -1;
    }
  }

  if(useDirection){
    if(verbose > 1){printf("\nDirection sampler '%s':\n\n",directionType.c_str());}
    err = selectDirectionSampler(directionType.c_str(),directionSection,verbose);
    if(verbose > 1){printf("\n------------------------------------\n\n");} 
    if(err != INTDATA_SUCCESS){
      clear();
      return -2;
    }
  }

  if(useEnergy){
    if(verbose > 1){printf("\nEnergy sampler '%s':\n\n",energyType.c_str());}    
    err = selectEnergySampler(energyType.c_str(),energySection,verbose);
    if(verbose > 1){printf("\n------------------------------------\n\n");}
    if(err != INTDATA_SUCCESS){
      clear();
      return -3;
    }
  }

  if(UseTimeSampling){
    if(verbose > 1){printf("\nTime sampler '%s':\n\n",timeType.c_str());}
    err = selectTimeSampler(timeType.c_str(),timeSection,verbose);
    if(verbose > 1){printf("\n------------------------------------\n\n");}
    if(err != INTDATA_SUCCESS){
      clear();
      return -4;
    }
  }

  //********************************
  // Optional parameters
  //********************************  

  //Get source body
  int sbody;
  sourceBody = -1;
  err = config.read("source-body",sbody);
  if(err == INTDATA_SUCCESS && sbody >= 0){
    if(sbody >= int(pen_geoconst::NB)){
      if(verbose > 0){
	printf("genericStateGen: configure: Warning: '%s': Selected source body (%d) out of range [0-%u)\n",name.c_str(),sbody,pen_geoconst::NB);
      }
    }
    else{
      sourceBody = sbody;
      if(verbose > 1){
	printf("Selected source body: %d\n",sourceBody);
      }
    }
  } else if(verbose > 1){
    printf("No source body selected\n");
  }

  //Get source material
  int smat;
  sourceMat = 0;
  err = config.read("source-material",smat);
  if(err == INTDATA_SUCCESS){
    if(smat < 1 || smat > int(constants::MAXMAT)){
      if(verbose > 0){
	printf("genericStateGen: configure: Warning: '%s': Selected source material (%d) out of range (0-%u]\n",name.c_str(),smat,constants::MAXMAT);	
      }
    }
    else{
      sourceMat = unsigned(smat);
      if(verbose > 1){
	printf("Selected source material: %u\n",sourceMat);
      }
    }
  } else if(verbose > 1){
    printf("No source material selected\n");
  }
  
  if(verbose > 1){

    printf("\n");
    printf("Generator '%s' configured with generic samplers:\n",name.c_str());
    printf("Spatial   -> %s\n", spatialID());
    printf("Direction -> %s\n", directionID());
    printf("Energy    -> %s\n", energyID());
    printf("Time      -> %s\n", timeID());
  }

  return 0;
}
  

  bool useGeneric;
  bool useSpecific;

  const wrapper_geometry*& geometry;
  
  int& configStatus;
  double& Emax;
  
public:

  typedef pen_genericStateGen::taskType taskType;
  
  int& sourceBody; //Body index source
  unsigned& sourceMat;  //Body material source
  
  std::string& name; //Source name
  bool& LAGE; //Age recording status
  pen_KPAR& kpar; //Generated particles kpar

  taskType& task;
  
  pen_specificStateGen() : useGeneric(false),
			   useSpecific(false),
			   geometry(genericGen.geometry),
			   configStatus(genericGen.configStatus),
			   Emax(genericGen.Emax),
			   sourceBody(genericGen.sourceBody),
			   sourceMat(genericGen.sourceMat),
			   name(genericGen.name),
			   LAGE(genericGen.LAGE),
  			   kpar(genericGen.kpar),
			   task(genericGen.task)
{}

  inline int initTask(const size_t nw,
		      const unsigned long long nIter,
		      const unsigned verbose){
    return genericGen.initTask(nw,nIter,verbose);
  }

  
  // ******************************* MPI ************************************ //
#if defined(_PEN_USE_MPI_) && defined(_PEN_USE_LB_)
  inline int initTask(const size_t nw,
		      const unsigned long long nIter,
		      const MPI_Comm commIn,
		      const int tagReqests,
		      const int tagProcess,
		      const unsigned verbose){
    return genericGen.initTask(nw,nIter,commIn,tagReqests,tagProcess,verbose);
  }
#endif
  // ***************************** MPI END ********************************** //

  inline void skip(const unsigned long long nIter){
    genericGen.skip(nIter);
  }

  inline unsigned long long toDo(unsigned ithread){
    return genericGen.toDo(ithread);
  }
  inline unsigned long long toDo(){
    return genericGen.toDo();
  }

  inline std::chrono::seconds::rep 
  report(const size_t iw,
	 const unsigned long long nIter,
	 int* errRet,
	 const unsigned verbose){
    return genericGen.report(iw,nIter,errRet,verbose);
  }

  inline int checkPoint(const unsigned verbose){
    return genericGen.checkPoint(verbose);
  }
  inline int workerStart(const size_t iw,
			 const unsigned verbose){
    return genericGen.workerStart(iw,verbose);
  }  
  inline bool workerFinish(const size_t iw, int& why, const unsigned verbose){
    return genericGen.workerFinish(iw,why,verbose);
  }

  inline bool handleFinish(const unsigned iw,
			   const unsigned long long nDone,
			   unsigned long long& assigned,
			   const unsigned verbose){
    return genericGen.handleFinish(iw,nDone,assigned,verbose);
  }
  inline void setCheckTime(const std::chrono::seconds::rep t){
    genericGen.setCheckTime(t);
  }
  inline void setLBthreshold(const unsigned long long t) {
    genericGen.setLBthreshold(t);
  }
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  inline void setCheckTimeMPI(const std::chrono::seconds::rep t){
    genericGen.setCheckTimeMPI(t);
  }
#endif
  // ***************************** MPI END ********************************** //
  
  inline double maxEnergy(){return Emax;}
  inline const abc_spatialSampler* spatial() const {
    return genericGen.spatialSampler;}
  inline const abc_directionSampler* direction() const {
    return genericGen.directionSampler;}
  inline const abc_energySampler* energy() const {
    return genericGen.energySampler;}
  inline const abc_timeSampler* time() const {
    return genericGen.timeSampler;}
  inline const abc_specificSampler<particleState>* specific(unsigned ithread) const {
    if(ithread > specificSamplerVect.size())
      throw std::out_of_range("abc_specificSampler: specific: Specified thread out of range");
    return specificSamplerVect[ithread];}
  
  static std::string samplersList(){
    std::string aux(pen_genericStateGen::samplersList());
    aux += " --- Specific:\n";
    aux += specificSamplers().typesList();
    return aux;
  }
  static std::string samplersList(std::vector<std::string>& spatial,
				  std::vector<std::string>& direction,
				  std::vector<std::string>& energy,
				  std::vector<std::string>& time,
				  std::vector<std::string>& specific){
    std::string aux(pen_genericStateGen::samplersList(spatial,direction,energy,time));
    aux += "Specific:\n";
    aux += specificSamplers().typesList(specific);
    return aux;
  }
  
  inline int setGeometry(const wrapper_geometry* geometryIn){
    int ret = genericGen.setGeometry(geometryIn);
    if(ret == 0){
      for(unsigned ithread = 0; ithread < specificSamplerVect.size(); ithread++){
	specificSamplerVect[ithread]->setGeometry(geometry);
      }
    }
    return ret;
  }

  inline bool usesGeneric() const {return useGeneric;}
  inline bool usesSpecific() const {return useSpecific;}
  
  
  template <class subclass>
  static int addSpecificSampler(const char* typeID){
    return specificSamplers().template addSubType<subclass>(typeID);
  }
  
  int selectSpecificSampler(const char* ID,
			    const unsigned nthreads,
			    const pen_parserSection& config,
			    const unsigned verbose = 0){

    //Clear previous specific samplers
    clearSpecificSamplers();

    if(nthreads < 1){
      if(verbose > 0){
	printf("SelectSpecificSampler: Error: Number of threads must be greater than zero.\n");
      }
      return -1;
    }

    if(verbose > 1){
      printf(" Configure specific sampler with %u threads.\n",nthreads);
    }

    //Resize samplers vector
    specificSamplerVect.resize(nthreads,nullptr);
    
    //Create sampler instances
    //*************************
    for(unsigned i = 0; i < nthreads; i++){
      specificSamplerVect[i] = specificSamplers().createInstance(ID);
      if(specificSamplerVect[i] == nullptr){
	if(verbose > 0){
	  printf("SelectSpecificSampler: Error: Unable to create specific sampler '%s' for thread %d\n",ID,i);
	}
	clearSpecificSamplers();
	return -2;
      }
      //Set thread to sampler
      specificSamplerVect[i]->setThread(i);
    }

    //Check if all required samplers has been provided
    //***************************************************

    // *** Spatial
    
    if(specificSamplerVect[0]->usesSpatial() && spatial() == nullptr){
      if(verbose > 0){
	printf("SelectSpecificSampler: Error: Specific sampler '%s' requires a spatial sampler.\n",ID);
      }
      clearSpecificSamplers();
      return -3;
    }
    else if(spatial() != nullptr &&
	    !specificSamplerVect[0]->usesSpatial() &&
	    !specificSamplerVect[0]->usesGeneric()
	    ){
      if(verbose > 2){
	printf("SelectSpecificSampler: Warning: Spatial sampler selected but not required by specific sampler. Could be ignored.\n");
      }
    }

    // *** Direction

    if(specificSamplerVect[0]->usesDirection() && direction() == nullptr){
      if(verbose > 0){
	printf("SelectSpecificSampler: Error: Specific sampler '%s' requires a direction sampler.\n",ID);
      }
      clearSpecificSamplers();
      return -4;
    }
    else if(direction() != nullptr &&
	    !specificSamplerVect[0]->usesDirection() &&
	    !specificSamplerVect[0]->usesGeneric()
	    ){
      if(verbose > 2){
	printf("SelectSpecificSampler: Warning: Direction sampler selected but not required by specific sampler. Could be ignored.\n");
      }
    }    
    
    // *** Energy
    
    if(specificSamplerVect[0]->usesEnergy() && energy() == nullptr){
      if(verbose > 0){
	printf("SelectSpecificSampler: Error: Specific sampler '%s' requires an energy sampler.\n",ID);
      }
      clearSpecificSamplers();
      return -5;
    }
    else if(energy() != nullptr &&
	    !specificSamplerVect[0]->usesEnergy() &&
	    !specificSamplerVect[0]->usesGeneric()
	    ){
      if(verbose > 2){
	printf("SelectSpecificSampler: Warning: Energy sampler selected but not required by specific sampler. Could be ignored.\n");
      }
    }    

    // *** Time
    
    if(specificSamplerVect[0]->usesTime() && time() == nullptr){
      if(verbose > 0){
	printf("SelectSpecificSampler: Error: Specific sampler '%s' requires a time sampler.\n",ID);
      }
      clearSpecificSamplers();
      return -6;
    }
    else if(time() != nullptr &&
	    !specificSamplerVect[0]->usesTime() &&
	    !specificSamplerVect[0]->usesGeneric()
	    ){
      if(verbose > 2){
	printf("SelectSpecificSampler: Warning: Time sampler selected but not required by specific sampler. Could be ignored.\n");
      }
    }
    
    // *** Generic

    if(specificSamplerVect[0]->usesGeneric()){
      if(spatial() == nullptr ||
	 direction() == nullptr ||
	 energy() == nullptr){
	
	if(verbose > 0){
	  printf("SelectSpecificSampler: Error: Specific sampler '%s' requires generic sampling. So, spatial, direction and energy generic samplers must be provided.\n",ID);
	}
	clearSpecificSamplers();
	return -7;
      }
      useGeneric = true;
    }

    //Configure specific samplers
    //******************************
    unsigned auxVerbose = verbose;
    for(unsigned i = 0; i < nthreads; i++){
      int errConfig = specificSamplerVect[i]->configure(Emax,
							spatial(),
							direction(),
							energy(),
							time(),
							config,
							auxVerbose);
      if(errConfig != 0){
	if(verbose > 0){
	  printf("SelectSpecificSampler: Error: Unable to configure specific sampler '%s' for thread %u\n",ID,i);
	}
	clearSpecificSamplers();
	return -8;
      }
      if(verbose >= 1)
	auxVerbose = 1;
    }
    
    useSpecific = true;
    return 0;
  }
  
  inline int selectSpatialSampler(const char* ID,
				  const pen_parserSection& config,
				  const unsigned verbose = 0){
    int ret = genericGen.selectSpatialSampler(ID,config,verbose);
    if(ret == 0){
      for(unsigned ithread = 0; ithread < specificSamplerVect.size(); ithread++){
	specificSamplerVect[ithread]->setSpatial(spatial());
      }
    }
    return ret;
  }

  inline int selectDirectionSampler(const char* ID,
				    const pen_parserSection& config,
				    const unsigned verbose = 0){
    int ret = genericGen.selectDirectionSampler(ID,config,verbose);
    if(ret == 0){
      for(unsigned ithread = 0; ithread < specificSamplerVect.size(); ithread++){
	specificSamplerVect[ithread]->setDirection(direction());
      }
    }    
    return ret;
  }
  
  inline int selectEnergySampler(const char* ID,
				 const pen_parserSection& config,
				 const unsigned verbose = 0){
    int ret = genericGen.selectEnergySampler(ID,config,verbose);
    if(ret == 0){
      for(unsigned ithread = 0; ithread < specificSamplerVect.size(); ithread++){
	specificSamplerVect[ithread]->setEnergy(energy());
      }
    }
    return ret;
  }

  inline int selectTimeSampler(const char* ID,
			       const pen_parserSection& config,
			       const unsigned verbose = 0){
    int ret = genericGen.selectTimeSampler(ID,config,verbose);
    if(ret == 0){
      for(unsigned ithread = 0; ithread < specificSamplerVect.size(); ithread++){
	specificSamplerVect[ithread]->setTime(time());
      }
    }
    return ret;
  }

  void sample(particleState& state,
	      pen_KPAR& genKpar,
	      unsigned long long& dhist,
	      const unsigned thread, 
	      pen_rand& random){

    genKpar = kpar;
    
    //Check if generic sampling is required
    if(useGeneric){
      //Perform generic sampling
      dhist = 1;
      genericGen.sample(state,random);

      //If specified, perform specific sampling
      if(useSpecific){
	specificSamplerVect[thread]->sample(state,genKpar,dhist,random);
      }
    }
    else{ //No generic sampling
      //If specified, perform specific sampling
      if(useSpecific){
	specificSamplerVect[thread]->sample(state,genKpar,dhist,random);

	//Locate particle in geometry
	geometry->locate(state);
      }
      else{
	//No sampler specified!
	printf("pen_specificStateGen:sample: Error: No generic nor specific sampler specified!\n");
	state.reset();
	genKpar = ALWAYS_AT_END;
	dhist = 0;	
      }
    }    
  }

  void skip(const unsigned long long dhists, const unsigned thread){
    //Only specific samplers can skip histories.
    //Notice that the history skip will not change the random seeds
    if(useSpecific){
      specificSamplerVect[thread]->skip(dhists);
    }
    skip(dhists);
  }

  void skip(const std::vector<unsigned long long> dhists){

    const size_t nthreads = specificSamplerVect.size();
    if(dhists.size() < nthreads)
      return;
    unsigned long long totalHists = 0.0;
    for(size_t ithread = 0; ithread < nthreads; ++ithread){
      totalHists += dhists[ithread];
      if(useSpecific){
	specificSamplerVect[ithread]->skip(dhists[ithread]);
      }
    }
    skip(totalHists);
  }
  
  void clear(){
      genericGen.clear();
      clearSpecificSamplers();
      useGeneric = false;
  }

  inline int configureStatus(){return configStatus;}

  void configure(const pen_parserSection& config,
		 const unsigned nthreads,
		 const unsigned verbose = 0){

    //Check if specific field exists
    bool specificSpecified = config.isSection("specific");
    
    //Configure generic part
    int err = configureGeneric(config,specificSpecified,verbose);
    if(err != 0){
      if(verbose > 0){
	printf("specificStateGen: configure: Error: Unable to configure generic state generation part.\n");
      }
      clear();
      configStatus = -1;
      return;
    }

    //Check if specific sapling has been enabled
    if(!specificSpecified){
      if(verbose > 1){
	printf("Section 'specific' not found. No specific sampler will be used\n");
	useGeneric = true;
      }
      configStatus = 0;
      return;
    }

    //Extract 'specific' sampler field
    pen_parserSection specificSection;
    if(config.readSubsection("specific",specificSection) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("specificStateGen: configure: Error: Unable to read 'specific' field. Section expected\n");
      }
      clear();
      configStatus = -2;
      return;
    }

    //Get specific sampler ID
    std::string specificType;
    if(specificSection.read("type",specificType) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("specificStateGen: configure: Error: Unable to read field 'specific/type'. String expected.\n");
      }
      clear();      
      configStatus = -5;
      return;
    }

    if(verbose > 1){printf("\nSpecific sampler '%s':\n\n",specificType.c_str());}
    if(selectSpecificSampler(specificType.c_str(),
			     nthreads,
			     specificSection,
			     verbose) != 0){
      clear();
      configStatus = -6;
      return;
    }
    
    if(verbose > 1){
      printf("specificStateGen: configure: Info: Generator '%s' configured with %u threads and specific sampler:\n",name.c_str(),nthreads);
      printf("Specific   -> %s\n", specificID());
    }
    
    configStatus = 0;
    return;
  }

  inline const char* specificID(){
    if(specificSamplerVect.size() > 0)
      return specificSamplerVect[0]->readID();
    return nullptr;
  }

  inline const char* spatialID(){
    return genericGen.spatialID();
  }

  inline const char* directionID(){
    return genericGen.directionID();
  }

  inline const char* energyID(){
    return genericGen.energyID();
  }

  inline const char* timeID(){
    return genericGen.timeID();
  }
  
  virtual ~pen_specificStateGen(){clear();}
  
};

template <class subclass>
static int registerSampler(const char* typeID){
  //Check what kind of sampler is "subclass"
  if(strcmp(subclass::type(), "SPATIAL") == 0){
    return pen_genericStateGen::addSpatialSampler<subclass>(typeID);
  }
  else if(strcmp(subclass::type(), "DIRECTION") == 0){
    return pen_genericStateGen::addDirectionSampler<subclass>(typeID);
  }
  else if(strcmp(subclass::type(), "ENERGY") == 0){
    return pen_genericStateGen::addEnergySampler<subclass>(typeID);
  }
  else if(strcmp(subclass::type(), "TIME") == 0){
    return pen_genericStateGen::addTimeSampler<subclass>(typeID);
  }
  else{
    //Unknown sampler type
    return -99;
  }
}


template <class subclass,class particleState>
static int registerSpecificSampler(const char* typeID){
  //Check what kind of sampler is "subclass"
  if(strcmp(subclass::type(), "SPECIFIC") == 0){
    return pen_specificStateGen<particleState>::template addSpecificSampler<subclass>(typeID);
  }
  else{
    //Unknown sampler type
    return -99;
  }
}

//Include defined samplers
#include "specificSamplers.hh"
#include "spatialSamplers.hh"
#include "directionSamplers.hh"
#include "energySamplers.hh"
#include "timeSamplers.hh"
  
#endif
