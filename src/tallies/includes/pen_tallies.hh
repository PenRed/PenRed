
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
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#ifndef __PEN_TALLIES__
#define __PEN_TALLIES__

#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <thread>
#include <stdexcept>
#include <functional>

#include "pen_auxiliar.hh"

// ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
#include "mpi.h"
#endif
// ***************************** MPI END ********************************** //


#include "pen_classes.hh"
#include "pen_math.hh"
#include "instantiator.hh"
#include "pen_data.hh"
#include "pen_dump.hh"
#include "pen_image.hh"
#include "pen_geometries.hh"
#include "pen_states.hh"

// Auxiliar structures
//////////////////////////

struct tally_StepData{
  double dsef;
  double dstot;
  double softDE;
  double softX;
  double softY;
  double softZ;
  unsigned originIBODY;
  unsigned originMAT;

  template <class particleType>
  inline void update(const particleType& particle){
    dsef = particle.DSef();
    dstot = particle.DStot();
    originMAT = particle.lastMat();
    originIBODY = particle.lastBody();
  }  
  
};


// Register macros
///////////////////////

#define DECLARE_TALLY(Class,State)			   \
  public:						   \
  static int registerStatus();				   \
  const char* readID() const;				   \
  static const char* ___ID;\
  static volatile const int ___register_return;\
  inline int sum(const pen_genericTally<State>& sumtally){\
  const Class& derived = dynamic_cast<const Class&>(sumtally);\
  return sumTally(derived);\
  }\
  inline int shareConfig(const pen_genericTally<State>& sharingTally){ \
  const Class& derived = dynamic_cast<const Class&>(sharingTally);\
  return sharedConfig(derived);\
  }\
  private:

#define REGISTER_COMMON_TALLY(Class, ID) \
  volatile const int Class::___register_return = pen_commonTallyCluster::addTally<Class>(static_cast<const char *>(#ID)); \
  const char* Class::___ID = static_cast<const char *>(#ID);		\
  int Class::registerStatus() { return ___register_return;}		\
  const char* Class::readID() const { return ___ID;}

#define REGISTER_SPECIFIC_TALLY(Class, ID) \
  volatile const int Class::___register_return = pen_specificTallyCluster::addTally<Class>(static_cast<const char *>(#ID)); \
  const char* Class::___ID = static_cast<const char *>(#ID);		\
  int Class::registerStatus() { return ___register_return;}		\
  const char* Class::readID() const { return ___ID;}

// Enumeration flags
//////////////////////////

enum __usedFunc{
	      USE_BEGINSIM     = 1 << 0,
	      USE_ENDSIM       = 1 << 1,
	      USE_SAMPLEDPART  = 1 << 2,
	      USE_ENDHIST      = 1 << 3,
	      USE_MOVE2GEO     = 1 << 4,
	      USE_BEGINPART    = 1 << 5,
	      USE_ENDPART      = 1 << 6,
	      USE_JUMP         = 1 << 7,
	      USE_STEP         = 1 << 8,
	      USE_INTERFCROSS  = 1 << 9,
	      USE_MATCHANGE    = 1 << 10,
	      USE_KNOCK        = 1 << 11,
	      USE_LOCALEDEP    = 1 << 12,
	      USE_LASTHIST     = 1 << 13,
	      INVALID_FUNCTION = 1 << 14
		
};


inline __usedFunc operator|(__usedFunc a, __usedFunc b)
{return static_cast<__usedFunc>(static_cast<int>(a) | static_cast<int>(b));}

// Tally class
//////////////////////////

template <class stateType>
class pen_genericTally : public penred::logs::logger{

private:

  unsigned nthread;
  std::string name;
  std::string OutputDirPath;
  const __usedFunc usedFunctions;

  //Create a vector of image exporters
  std::vector<pen_imageExporter> imageExporters;

  //Create an array to store particle stacks
  std::array<const abc_particleStack*, constants::nParTypes> stacks;
protected:

  //Create dump instance
  pen_dump dump;

  //Include sub dump
  inline void addSubDump(pen_genericTally<stateType>& subTally){
    dump.toDump(subTally.dump);
  }

  //Include measurement in dump
  template<class valType, size_t dim>
  inline void toDump(penred::measurements::measurement<valType, dim>& m){
    dump.toDump(m.getData().data(), m.getNBins());
    dump.toDump(m.getData2().data(), m.getNBins());
  }

  //Create filenames with output dir path, thread and MPI number
  std::string createFilename(const char* filename) const{

    // ******************************* MPI ********************************** //
#ifdef _PEN_USE_MPI_
    int rank;
    //Add the MPI rank number to filename
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    return std::string(OutputDirPath + name + std::string("-MPI") +
		       std::to_string(rank) + std::string("-th") +
		       std::to_string(nthread) + std::string("-") +
		       std::string(filename));

    // ***************************** MPI END ******************************** //
#else
    
    return std::string(OutputDirPath + name + std::string("-th") +
		       std::to_string(nthread) + std::string("-") +
		       std::string(filename));

#endif    
  }  
  
  //Create private fopen method to append tally name and thread
  //identifier to generated and read files
  inline FILE* fopen(const char* filename, const char* mode) const{

    if(filename == nullptr || mode == nullptr)
      return nullptr;

    std::string finalFilename;
    FILE* fout = nullptr;
    
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
    int rank;
    //Add the MPI rank number to filename
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    finalFilename = OutputDirPath + name + std::string("-MPI") +
      std::to_string(rank) + std::string("-th") + std::to_string(nthread) +
      std::string("-") + std::string(filename);

  // ***************************** MPI END ********************************** //
#else
    
    finalFilename = OutputDirPath + name + std::string("-th") +
      std::to_string(nthread) + std::string("-") + std::string(filename);

#endif


    fout = ::fopen(finalFilename.c_str(),mode);
    //Check if the file can be opened
    if(fout == nullptr){
      //Try to open the same file name but with a 'not-found' prefix
      printf("pen_genericTally: fopen: Error: Unable to open "
	     "file ('%s'): \n%s\n",
	     mode, finalFilename.c_str());
      
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_

      finalFilename = std::string("not-found-") + name + std::string("-MPI") +
	std::to_string(rank) + std::string("-th") + std::to_string(nthread) +
	std::string("-") + std::string(filename);
      
  // ***************************** MPI END ********************************** //
#else
      
      finalFilename = std::string("not-found-") + name + std::string("-th") +
	std::to_string(nthread) + std::string("-") + std::string(filename);
      
#endif

      printf("Trying to redirect to: '%s'\n",finalFilename.c_str());
      fflush(stdout);

      fout = ::fopen(finalFilename.c_str(),mode);
    }
    //Return resulting file pointer
    return fout;
  }

public:

  pen_genericTally(__usedFunc functionsToUse) : nthread(0),
						name(""),
						usedFunctions(functionsToUse)
  {
    std::fill(stacks.begin(), stacks.end(), nullptr);
  }

  virtual const char* readID() const = 0;
  inline const std::string& readName() const { return name;}
  inline unsigned getThread() const {return nthread;}
  inline const std::string& readOutputDirPath() const { return OutputDirPath;}
  inline void setName(const char* newName) { name.assign(newName);}
  inline void setName(const std::string& newName) { name.assign(newName);}
  inline void setOutputDirPath(const char* newName) { OutputDirPath.assign(newName);}
  inline void setOutputDirPath(const std::string& newName) { OutputDirPath.assign(newName);}
  inline void setThread(const unsigned threadNum){nthread = threadNum;}
  inline void setStack(const pen_KPAR kpar, const abc_particleStack* stack){stacks[kpar] = stack;}
  inline const abc_particleStack* readStack(const pen_KPAR kpar) const {return stacks[kpar];}

  template<class T>
  inline void addImage(const char* filename,
		       const unsigned nDim,
		       const unsigned* elements,
		       const float* delements,
		       const double* origin,
		       std::function<T(unsigned long long, size_t)> fin){
    imageExporters.emplace(imageExporters.end(), fin);
    imageExporters.back().setDimensions(nDim,elements,delements);
    imageExporters.back().setOrigin(origin);
    imageExporters.back().baseName = createFilename(filename);
  }

  template<class T>
  inline void addImage(const char* filename,
		       const unsigned nDim,
		       const unsigned* elements,
		       const float* delements,
		       const double* origin,
		       std::function<T(unsigned long long, size_t, T&)> fin){
    imageExporters.emplace(imageExporters.end(), fin);
    imageExporters.back().setDimensions(nDim,elements,delements);
    imageExporters.back().setOrigin(origin);
    imageExporters.back().baseName = createFilename(filename);
  }  

  inline void exportImage(const unsigned long long nhists,
			  const pen_imageExporter::formatTypes format) const {
    for(auto& exporter : imageExporters)
      exporter.exportImage(nhists,format);
  }
  
  void getUsedMethods(std::vector<std::string>& funcNames) const{
    if(has_beginSim())
      funcNames.push_back(std::string("beginSim"));
    if(has_endSim())
      funcNames.push_back(std::string("endSim"));  
    if(has_sampledPart())
      funcNames.push_back(std::string("sampledPart"));  
    if(has_endHist())
      funcNames.push_back(std::string("endHist"));
    if(has_move2geo())
      funcNames.push_back(std::string("move2geo"));
    if(has_beginPart())
      funcNames.push_back(std::string("beginPart"));  
    if(has_endPart())
      funcNames.push_back(std::string("endPart"));  
    if(has_jump())
      funcNames.push_back(std::string("jump"));      
    if(has_step())
      funcNames.push_back(std::string("step"));
    if(has_interfCross())
      funcNames.push_back(std::string("interfCross"));      
    if(has_matChange())
      funcNames.push_back(std::string("matChange"));      
    if(has_knock())
      funcNames.push_back(std::string("knock"));      
    if(has_localEdep())
      funcNames.push_back(std::string("localEdep"));
    if(has_lastHist())
      funcNames.push_back(std::string("lastHist"));    
  }
  
  virtual void tally_beginSim(){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_beginSim' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  virtual void tally_endSim(const unsigned long long /*nhist*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_endSim' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_sampledPart(const unsigned long long /*nhist*/,
				 const unsigned long long /*dhist*/,
				 const unsigned /*kdet*/,
				 const pen_KPAR /*kpar*/,
				 const stateType& /*state*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_sampledPart' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_endHist(const unsigned long long /*nhist*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_endHist' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }

  virtual void tally_move2geo(const unsigned long long /*nhist*/,
			      const unsigned /*kdet*/,
			      const pen_KPAR /*kpar*/,
			      const stateType& /*state*/,
			      const double /*dsef*/,
			      const double /*dstot*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_move2geo' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_beginPart(const unsigned long long /*nhist*/,
			       const unsigned /*kdet*/,
			       const pen_KPAR /*kpar*/,
			       const stateType& /*state*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_beginPart' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_endPart(const unsigned long long /*nhist*/,
			     const pen_KPAR /*kpar*/,
			     const stateType& /*state*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_endPart' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_localEdep(const unsigned long long /*nhist*/,
			       const pen_KPAR /*kpar*/,
			       const stateType& /*state*/,
			       const double /*dE*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_localEdep' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_step(const unsigned long long /*nhist*/,
			  const pen_KPAR /*kpar*/,
			  const stateType& /*state*/,
			  const tally_StepData& /*stepData*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_step' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_interfCross(const unsigned long long /*nhist*/,
				 const unsigned /*kdet*/,
				 const pen_KPAR /*kpar*/,
				 const stateType& /*state*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_interfCross' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_matChange(const unsigned long long /*nhist*/,
			       const pen_KPAR /*kpar*/,
			       const stateType& /*state*/,
			       const unsigned /*prevMat*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_matChange' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_jump(const unsigned long long /*nhist*/,
			  const pen_KPAR /*kpar*/,
			  const stateType& /*state*/,
			  const double /*ds*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_jump' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual void tally_knock(const unsigned long long /*nhist*/,
			   const pen_KPAR /*kpar*/,
			   const stateType& /*state*/,
			   const int /*icol*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_knock' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }

  virtual void tally_lastHist(const unsigned long long /*lasthist*/){
    char error[1000];
    sprintf(error,"Error on tally %s of type %s: \n Trying to call 'tally_lastHist' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }
  
  virtual int configure(const wrapper_geometry& /*geometry*/,
			const abc_material* const /*materials*/[constants::MAXMAT],
			const pen_parserSection& /*config*/,
			const unsigned /*verbose*/) = 0;
  
  virtual void saveData(const unsigned long long /*nhist*/) const = 0;
  virtual int sum(const pen_genericTally<stateType>&) = 0;
  virtual int shareConfig(const pen_genericTally<stateType>&) = 0;
  int sharedConfig(const pen_genericTally<stateType>&){return 0;}
  virtual void flush() = 0;
  
  inline __usedFunc readFlags() const {return usedFunctions;}
  
  inline bool has_beginSim() const{
    return ((usedFunctions & USE_BEGINSIM) == USE_BEGINSIM);
  }
  inline bool has_endSim() const{
    return ((usedFunctions & USE_ENDSIM) == USE_ENDSIM);
  }
  inline bool has_sampledPart() const{
    return ((usedFunctions & USE_SAMPLEDPART) == USE_SAMPLEDPART);
  }
  inline bool has_endHist() const{
    return ((usedFunctions & USE_ENDHIST) == USE_ENDHIST);
  }
  inline bool has_move2geo() const{
    return ((usedFunctions & USE_MOVE2GEO) == USE_MOVE2GEO);
  }
  inline bool has_beginPart() const{
    return ((usedFunctions & USE_BEGINPART) == USE_BEGINPART);
  }
  inline bool has_endPart() const{
    return ((usedFunctions & USE_ENDPART) == USE_ENDPART);
  }
  inline bool has_jump() const{
    return ((usedFunctions & USE_JUMP) == USE_JUMP);
  }
  inline bool has_step() const{
    return ((usedFunctions & USE_STEP) == USE_STEP);
  }
  inline bool has_interfCross() const{
    return ((usedFunctions & USE_INTERFCROSS) == USE_INTERFCROSS);
  }
  inline bool has_matChange() const{
    return ((usedFunctions & USE_MATCHANGE) == USE_MATCHANGE);
  }
  inline bool has_knock() const{
    return ((usedFunctions & USE_KNOCK) == USE_KNOCK);
  }
  inline bool has_localEdep() const{
    return ((usedFunctions & USE_LOCALEDEP) == USE_LOCALEDEP);
  }
  inline bool has_lastHist() const{
    return ((usedFunctions & USE_LASTHIST) == USE_LASTHIST);
  }
  

  inline int writeDump(unsigned char*& pout, size_t& pos, const unsigned verbose = 0) {
    //Flush tally
    flush();
    //Dump it
    return dump.dump(pout,pos,0,verbose);
  }

  inline int readDump(const unsigned char* const pin, size_t& pos, const unsigned verbose = 0){
    return dump.read(pin,pos,verbose);
  }
  
  virtual ~pen_genericTally(){}
  
};

// Auxiliar functions
///////////////////////

bool __tallySort (pen_genericTally<pen_particleState>* i,
		  pen_genericTally<pen_particleState>* j);

namespace penred{
  namespace tally{
    
    //Check registered types
    template <class stateType>
    bool checkRegistered(const unsigned verbose);
  }
}

// Tally cluster classes
//////////////////////////

class pen_commonTallyCluster : public penred::logs::logger{

private:
  unsigned nthread;
protected:

  typedef pen_genericTally<pen_particleState>* ptallyType;
  typedef std::vector<ptallyType>::iterator tallyIterator;
  
  //Tallies instantiators
  static instantiator<pen_genericTally<pen_particleState>>& genericTallies();
  
  //Registered tallies vector
  std::vector<pen_genericTally<pen_particleState>*> tallies;

  //Vectors containing tally pointers with specific collection functions
  std::vector<pen_genericTally<pen_particleState>*> tallies_beginSim;
  std::vector<pen_genericTally<pen_particleState>*> tallies_endSim;
  std::vector<pen_genericTally<pen_particleState>*> tallies_sampledPart;
  std::vector<pen_genericTally<pen_particleState>*> tallies_endHist;
  std::vector<pen_genericTally<pen_particleState>*> tallies_move2geo;
  std::vector<pen_genericTally<pen_particleState>*> tallies_beginPart;
  std::vector<pen_genericTally<pen_particleState>*> tallies_endPart;
  std::vector<pen_genericTally<pen_particleState>*> tallies_localEdep;
  std::vector<pen_genericTally<pen_particleState>*> tallies_step;
  std::vector<pen_genericTally<pen_particleState>*> tallies_interfCross;
  std::vector<pen_genericTally<pen_particleState>*> tallies_matChange;
  std::vector<pen_genericTally<pen_particleState>*> tallies_jump;
  std::vector<pen_genericTally<pen_particleState>*> tallies_knock;
  std::vector<pen_genericTally<pen_particleState>*> tallies_lastHist;

  int configStatus;

  pen_commonTallyCluster* mpiBuffer;

  int createTally(const char* OutDir, const char* ID,
		  const char* name,
		  const wrapper_geometry& geometry,
		  const abc_material* const materials[constants::MAXMAT], 
		  const pen_parserSection& config,
		  const unsigned verbose);  
public:

  std::string name;

  pen_commonTallyCluster() : nthread(0),
			     configStatus(0),
			     mpiBuffer(nullptr),
			     name("unnamed")
  {}

  inline int configureStatus(){return configStatus;}
  
  inline unsigned numTallies(){return tallies.size();}
  
  static std::string typesList(){
    return genericTallies().typesList();
  }

  static std::string typesList(std::vector<std::string>& list){
    return genericTallies().typesList(list);
  }

  static pen_genericTally<pen_particleState>* createInstance(const char* type){
    return genericTallies().createInstance(type);
  }  
  static pen_genericTally<pen_particleState>* createInstance(const std::string type){
    return genericTallies().createInstance(type.c_str());
  }

  inline unsigned getThread() const {return nthread;}
  
  template <class subclass>
  static int addTally(const char* typeID){
    return genericTallies().addSubType<subclass>(typeID);
  }

#ifdef _PEN_USE_THREADS_
  std::thread configure_async(const wrapper_geometry* geometry,
			      const abc_material* const materials[constants::MAXMAT],
			      const unsigned threadNum,
			      const pen_parserSection& config,
			      const unsigned verbose = 0);
#endif
  
  void configure(const wrapper_geometry* geometry,
		 const abc_material* const materials[constants::MAXMAT],
		 const unsigned threadNum,
		 const pen_parserSection config,
		 const unsigned verbose = 0);

  inline void setStack(const pen_KPAR kpar, const abc_particleStack* stack){
    for(tallyIterator i = tallies.begin();
	i != tallies.end(); ++i)
      (*i)->setStack(kpar,stack);
  }
  
  inline void run_beginSim(){
    for(tallyIterator i = tallies_beginSim.begin();
	i != tallies_beginSim.end(); ++i)
      (*i)->tally_beginSim();
  }

  inline void run_endSim(const unsigned long long nhist){
    for(tallyIterator i = tallies_endSim.begin();
	i != tallies_endSim.end(); ++i)
      (*i)->tally_endSim(nhist);
  }
  
  inline void run_sampledPart(const unsigned long long nhist,
			      const unsigned long long dhist,
			      const unsigned kdet,
			      const pen_KPAR kpar,
			      const pen_particleState& state){

    for(tallyIterator i = tallies_sampledPart.begin();
	i != tallies_sampledPart.end(); ++i)
      (*i)->tally_sampledPart(nhist,dhist,kdet,kpar,state);
  }
  
  inline void run_endHist(const unsigned long long nhist){

    for(tallyIterator i = tallies_endHist.begin();
	i != tallies_endHist.end(); ++i)
      (*i)->tally_endHist(nhist);
  }

  inline void run_move2geo(const unsigned long long nhist,
			   const unsigned kdet,
			   const pen_KPAR kpar,
			   const pen_particleState& state,
			   const double dsef,
			   const double dstot){
    
    for(tallyIterator i = tallies_move2geo.begin();
	i != tallies_move2geo.end(); ++i)
      (*i)->tally_move2geo(nhist,kdet,kpar,state,dsef,dstot);
  }
  
  inline void run_beginPart(const unsigned long long nhist,
			    const unsigned kdet,
			    const pen_KPAR kpar,
			    const pen_particleState& state){

    for(tallyIterator i = tallies_beginPart.begin();
	i != tallies_beginPart.end(); ++i)
      (*i)->tally_beginPart(nhist,kdet,kpar,state);
  }
  
  inline void run_endPart(const unsigned long long nhist,
			  const pen_KPAR kpar,
			  const pen_particleState& state){

    for(tallyIterator i = tallies_endPart.begin();
	i != tallies_endPart.end(); ++i)
      (*i)->tally_endPart(nhist,kpar,state);
  }
  
  inline void run_localEdep(const unsigned long long nhist,
			const pen_KPAR kpar,
			const pen_particleState& state,
			const double dE){

    for(tallyIterator i = tallies_localEdep.begin();
	i != tallies_localEdep.end(); ++i)
      (*i)->tally_localEdep(nhist,kpar,state,dE);
  }
  
  inline void run_step(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const tally_StepData& stepData){

    for(tallyIterator i = tallies_step.begin();
	i != tallies_step.end(); ++i)
      (*i)->tally_step(nhist,kpar,state,stepData);    
  }
  
  inline void run_interfCross(const unsigned long long nhist,
			      const unsigned kdet,
			      const pen_KPAR kpar,
			      const pen_particleState& state){

    for(tallyIterator i = tallies_interfCross.begin();
	i != tallies_interfCross.end(); ++i)
      (*i)->tally_interfCross(nhist,kdet,kpar,state);
  }
  
  inline void run_matChange(const unsigned long long nhist,
			    const pen_KPAR kpar,
			    const pen_particleState& state,
			    const unsigned prevMat){

    for(tallyIterator i = tallies_matChange.begin();
	i != tallies_matChange.end(); ++i)
      (*i)->tally_matChange(nhist,kpar,state,prevMat);
  }
  
  inline void run_jump(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const double ds){

    for(tallyIterator i = tallies_jump.begin();
	i != tallies_jump.end(); ++i)
      (*i)->tally_jump(nhist,kpar,state,ds);
  }
  
  inline void run_knock(const unsigned long long nhist,
			const pen_KPAR kpar,
			const pen_particleState& state,
			const int icol){

    for(tallyIterator i = tallies_knock.begin();
	i != tallies_knock.end(); ++i)
      (*i)->tally_knock(nhist,kpar,state,icol);
  }

  inline void run_lastHist(const unsigned long long lastHist){

    for(tallyIterator i = tallies_lastHist.begin();
	i != tallies_lastHist.end(); ++i)
      (*i)->tally_lastHist(lastHist);
  }
  
  inline void saveData(const unsigned long long nhist, const bool doflush = true){

    for(tallyIterator i = tallies.begin();
	i != tallies.end(); ++i){
      if(doflush)
	(*i)->flush();
      (*i)->saveData(nhist);
    }
  }

  inline void exportImage(const unsigned long long nhist,
			  const pen_imageExporter::formatTypes format,
			  const bool doflush = true){

    for(tallyIterator i = tallies.begin();
	i != tallies.end(); ++i){
      if(doflush)
	(*i)->flush();
      (*i)->exportImage(nhist,format);
    }
  }  
  
  int writeDump(unsigned char*& pdump,
		size_t& dim,
		const unsigned long long nhist,
		const int seed1, const int seed2,
		const int lastSource,
		const unsigned long long sourceHists,
		const unsigned verbose = 0);
  int readDump(const unsigned char* const pdump,
	       size_t& pos,
	       unsigned long long& nhist,
	       int& seed1, int& seed2,
	       int& lastSource,
	       unsigned long long& sourceHists,
	       const unsigned verbose = 0);

  int dump2file(const char* filename,
		const unsigned long long nhist,
		const int seed1, const int seed2,
		const int lastSource,
		const unsigned long long sourceHists,
		const unsigned verbose = 0);

  int readDumpfile(const char* filename,
		   unsigned long long& nhist,
		   int& seed1, int& seed2,
		   int& lastSource,
		   unsigned long long& sourceHists,
		   const unsigned verbose = 0);
  
  int sum(pen_commonTallyCluster& cluster, const unsigned verbose = 0);
  int shareConfig(const pen_commonTallyCluster& cluster, const unsigned verbose = 0);
  
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  int reduceMPI(unsigned long long& simulatedHists,
		const MPI_Comm comm,
		const unsigned verbose);
#endif
  // ***************************** MPI END ********************************** //
  
  void clear();

  ~pen_commonTallyCluster();
  
};

template <class stateType>
class pen_specificTallyCluster : public penred::logs::logger{

private:
  unsigned nthread;
protected:

  typedef pen_genericTally<stateType>* ptallyType;
  typedef typename std::vector<ptallyType>::iterator tallyIterator;
  
  //Tallies instantiators
  static instantiator<pen_genericTally<stateType>>& specificTallies(){
    static instantiator<pen_genericTally<stateType>>* ans =
      new instantiator<pen_genericTally<stateType>>;
    return *ans;  
  }
  
  //Registered tallies vector
  std::vector<pen_genericTally<stateType>*> tallies;

  //Vectors containing tally pointers with specific collection functions
  std::vector<pen_genericTally<stateType>*> tallies_beginSim;
  std::vector<pen_genericTally<stateType>*> tallies_endSim;
  std::vector<pen_genericTally<stateType>*> tallies_sampledPart;
  std::vector<pen_genericTally<stateType>*> tallies_endHist;
  std::vector<pen_genericTally<stateType>*> tallies_move2geo;
  std::vector<pen_genericTally<stateType>*> tallies_beginPart;
  std::vector<pen_genericTally<stateType>*> tallies_endPart;
  std::vector<pen_genericTally<stateType>*> tallies_localEdep;
  std::vector<pen_genericTally<stateType>*> tallies_step;
  std::vector<pen_genericTally<stateType>*> tallies_interfCross;
  std::vector<pen_genericTally<stateType>*> tallies_matChange;
  std::vector<pen_genericTally<stateType>*> tallies_jump;
  std::vector<pen_genericTally<stateType>*> tallies_knock;
  std::vector<pen_genericTally<stateType>*> tallies_lastHist;
  int configStatus;

  pen_commonTallyCluster* mpiBuffer;
  
  int createTally(const char* ID,
		  const char* tallyname,
		  const wrapper_geometry& geometry,
		  const abc_material* const materials[constants::MAXMAT],
		  const pen_parserSection& config,
		  const unsigned verbose){

  
    //Check if tally name already exists
    std::string strName;
    if(tallyname == nullptr || strName.assign(tallyname).length() == 0){
      if(verbose > 0){
	printf("specificTallyCluster: createTally: Error: empty tally name.\n");
      }
      return -3;
    }
    for(unsigned i = 0; i < tallies.size(); i++){
      if(strName.compare(tallies[i]->readName()) == 0){
	if(verbose > 0){
	  printf("specificTallyCluster: createTally: Error: Tally name '%s' already used.\n", tallyname);
	}
	return -4;
      }
    }
  
    //Create specified tally by ID
    pen_genericTally<stateType>* ptally = nullptr;  
    ptally = specificTallies().createInstance(ID);
    if(ptally == nullptr){
      if(verbose > 0){
	printf("specificTallyCluster: createTally: Error: unable to create a tally of type '%s'.\n",ID);
      }
      return -1;
    }

    //Set tally name
    ptally->setName(tallyname);
    
    //Set thread
    ptally->setThread(nthread);

    //Configure tally
    int errConfig = ptally->configure(geometry,materials,config,verbose);
    if(errConfig != 0){
      delete ptally;
      if(verbose > 0){
	printf("specificTallyCluster: createTally: Error: tally '%s' of type '%s' failed on configuration step.\n",tallyname,ID);
      }
      return -2;
    }
  
    //Append tally pointer to necessary vectors
    tallies.push_back(ptally);
    
    if(ptally->has_beginSim()){ tallies_beginSim.push_back(ptally);}
    if(ptally->has_endSim()){ tallies_endSim.push_back(ptally);}
    if(ptally->has_sampledPart()){ tallies_sampledPart.push_back(ptally);}
    if(ptally->has_endHist()){ tallies_endHist.push_back(ptally);}
    if(ptally->has_move2geo()){ tallies_move2geo.push_back(ptally);}
    if(ptally->has_beginPart()){ tallies_beginPart.push_back(ptally);}
    if(ptally->has_endPart()){ tallies_endPart.push_back(ptally);}
    if(ptally->has_localEdep()){ tallies_localEdep.push_back(ptally);}
    if(ptally->has_step()){ tallies_step.push_back(ptally);}
    if(ptally->has_interfCross()){ tallies_interfCross.push_back(ptally);}
    if(ptally->has_matChange()){ tallies_matChange.push_back(ptally);}
    if(ptally->has_jump()){ tallies_jump.push_back(ptally);}
    if(ptally->has_knock()){ tallies_knock.push_back(ptally);}
    if(ptally->has_lastHist()){ tallies_lastHist.push_back(ptally);}
    
    //Sort vector tallies by name
    std::sort(tallies.begin(), tallies.end(), __tallySort);
    
    return 0;
  
  }
  
public:

  std::string name;

  pen_specificTallyCluster() : nthread(0),
			       configStatus(0),
			       mpiBuffer(nullptr),
			       name("unnamed")
  {}

  inline unsigned numTallies(){return tallies.size();}
  inline int configureStatus(){return configStatus;}

  inline unsigned getThread() const {return nthread;}
  
  template <class subclass>
  static int addTally(const char* typeID){
    return specificTallies().template addSubType<subclass>(typeID);
  }

  int sum(pen_specificTallyCluster<stateType>& cluster, const unsigned verbose = 0){
    //Check if both clusters have the same ammount of tallies
    unsigned size1 = tallies.size();
    unsigned size2 = cluster.tallies.size();
    if(size1 != size2){
      if(verbose > 0){
	printf("pen_specificTallyCluster: sum: Number of tallies doesn't match.\n");
      }      
      return -1;
    }

    //Check if both clusters have the same tally types and names
    for(unsigned i = 0; i < size1; i++){
      if(tallies[i]->readName().compare(cluster.tallies[i]->readName()) != 0){
	if(verbose > 0){
	  printf("pen_specificTallyCluster: sum: Tallies names doesn't match:\n");
	  printf("                          %s\n",tallies[i]->readName().c_str());
	  printf("                          %s\n",cluster.tallies[i]->readName().c_str());	  
	}      	
	return i+1;
      }
      if(strcmp(tallies[i]->readID(), cluster.tallies[i]->readID()) != 0){
	if(verbose > 0){
	  printf("pen_specificTallyCluster: sum: Tallies IDs doesn't match:\n");
	  printf("                          %s\n",tallies[i]->readID());
	  printf("                          %s\n",cluster.tallies[i]->readID());
	}
	return i+1;
      }
    }

    //Sum all tallies
    int err = 0;
    for(unsigned i = 0; i < size1; ++i){

      //Flush both tallies
      tallies[i]->flush();
      cluster.tallies[i]->flush();

      //Sum tallies
      int err2 = tallies[i]->sum(*(cluster.tallies[i]));
      if(err2 != 0){ 
	if(verbose > 0){
	  printf("pen_specificTallyCluster: sum: Error adding up tallies %s (position %d)\n",tallies[i]->readName().c_str(),i);
	  printf("                               Error code: %d\n", err2);
	}
	++err;
      }
    }
    return err;
  }
  
  void clear(){
    for(std::size_t i = 0; i < tallies.size(); i++){
      delete tallies[i];
    }
    tallies.clear();
    tallies_beginSim.clear();
    tallies_endSim.clear();
    tallies_sampledPart.clear();
    tallies_endHist.clear();
    tallies_move2geo.clear();
    tallies_beginPart.clear();
    tallies_endPart.clear();
    tallies_localEdep.clear();
    tallies_step.clear();
    tallies_interfCross.clear();
    tallies_matChange.clear();
    tallies_jump.clear();
    tallies_knock.clear();
    tallies_lastHist.clear();

    if(mpiBuffer != nullptr)
      delete mpiBuffer;
    mpiBuffer = nullptr;
    
    nthread = 0;
    configStatus = 0;
  }  
  
  void configure(const wrapper_geometry* geometry,
		 const abc_material* const materials[constants::MAXMAT],
		 const unsigned threadNum,
		 const pen_parserSection config,
		 const unsigned verbose = 0){

    int err = 0;
    
    //Clear previous configuration
    clear();

    if(verbose > 1){
      printf("\n **** Tally group '%s'\n",name.c_str());
    }
    
    //Set thread
    nthread = threadNum;
  
    //Extract tally names
    std::vector<std::string> tallyNames;
    config.ls(tallyNames);
  
    //Iterate for each tally name
    for(unsigned i = 0; i < tallyNames.size(); i++){

      if(verbose > 1){
	printf("\n------------------------------------\n\n");
	printf("\nTally '%s':\n\n",tallyNames[i].c_str());
      }
    
      //Get subsection for tally 'i'
      pen_parserSection tallySec;
      if(config.readSubsection(tallyNames[i],tallySec) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("specificTallyCluster: configure: Error: unable to read section '%s' to configure tally.\n",tallyNames[i].c_str());
	}
	err++;
	continue;
      }

      //Try to read tally type ID
      std::string tallyID;
      if(tallySec.read("type",tallyID) != INTDATA_SUCCESS){
	if(verbose > 0){
	  printf("specificTallyCluster: configure: Error: unable to read field %s/type. String expected\n",tallyNames[i].c_str());
	  err++;
	  continue;
	}
      }

      if(verbose > 1){
	printf("Tally type: '%s'\n\n",tallyID.c_str());
      }
      
      //Try to create and configure tally 'i'
      if(createTally(tallyID.c_str(),tallyNames[i].c_str(),*geometry,materials,tallySec,verbose) != 0){
	if(verbose > 0){
	  printf("specificTallyCluster: configure: Error: Unable to create and configure tally '%s' of type '%s'.\n",tallyNames[i].c_str(),tallyID.c_str());
	  err++;
	  continue;
	}
      }    
    }
    if(verbose > 1){printf("\n------------------------------------\n\n");}

    if(err > 0){
      if(verbose > 0){
	printf("specificTallyCluster: configure: Error: %d tallies creation failed.\n",err);
      }
    }

    //Print created tallies
    if(verbose > 1){
      printf("\nCreated specific tallies (name and type):\n\n");
      for(unsigned i = 0; i < tallies.size(); i++){
	printf("  %20s -> %20s\n",
	       tallies[i]->readName().c_str(),
	       tallies[i]->readID());
      }
    }

    configStatus = err;

  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
  //If this tally cluster runs on thread 0, we need a buffer cluster
  //to reduce the results of all MPI processes
  if(nthread == 0){
    mpiBuffer = new pen_commonTallyCluster;
    mpiBuffer->name = std::string("MPI_buffer_") + name;
    //Configure tally cluster with no verbose and a non zero thread
    //to avoid recursive creation of mpi buffers
    mpiBuffer->configure(geometry,
			 materials,
			 1,
			 config,
			 0);
    if(mpiBuffer->configureStatus() > 0){
      if(verbose > 0){
	printf("commonTallyCluster: configure: Error creating %d tallies on MPI buffer cluster.\n",err);
      }
    }
  }
#endif
  // ***************************** MPI END ********************************** //
    
  }

  std::thread configure_async(const wrapper_geometry* geometry,
			      const abc_material* const materials[constants::MAXMAT],
			      const unsigned threadNum,
			      const pen_parserSection& config,
			      const unsigned verbose = 0){
    return std::thread(configure(&pen_specificTallyCluster::configure,this,
				 geometry,materials,threadNum,config,verbose));
  }

  inline void run_beginSim(){

    //Default all tally logs to simulation log
    for(auto t : tallies)
      t->setDefaultLog(penred::logs::SIMULATION);
    
    for(tallyIterator i = tallies_beginSim.begin();
	i != tallies_beginSim.end(); ++i)
      (*i)->tally_beginSim();
  }

  inline void run_endSim(const unsigned long long nhist){

    //Default all tally logs to configuration log
    for(auto t : tallies)
      t->setDefaultLog(penred::logs::CONFIGURATION);
    
    for(tallyIterator i = tallies_endSim.begin();
	i != tallies_endSim.end(); ++i)
      (*i)->tally_endSim(nhist);
  }
  
  inline void run_sampledPart(const unsigned long long nhist,
			      const unsigned long long dhist,
			      const unsigned kdet,
			      const pen_KPAR kpar,
			      const stateType& state){

    for(tallyIterator i = tallies_sampledPart.begin();
	i != tallies_sampledPart.end(); ++i)
      (*i)->tally_sampledPart(nhist,dhist,kdet,kpar,state);
  }
  
  
  inline void run_endHist(const unsigned long long nhist){

    for(tallyIterator i = tallies_endHist.begin();
	i != tallies_endHist.end(); ++i)
      (*i)->tally_endHist(nhist);
  }

  inline void run_move2geo(const unsigned long long nhist,
			   const unsigned kdet,
			   const pen_KPAR kpar,
			   const stateType& state,
			   const double dsef,
			   const double dstot){
    
    for(tallyIterator i = tallies_move2geo.begin();
	i != tallies_move2geo.end(); ++i)
      (*i)->tally_move2geo(nhist,kdet,kpar,state,dsef,dstot);
  }
  
  inline void run_beginPart(const unsigned long long nhist,
			    const unsigned kdet,
			    const pen_KPAR kpar,
			    const stateType& state){

    for(tallyIterator i = tallies_beginPart.begin();
	i != tallies_beginPart.end(); ++i)
      (*i)->tally_beginPart(nhist,kdet,kpar,state);
  }
  
  inline void run_endPart(const unsigned long long nhist,
			  const pen_KPAR kpar,
			  const stateType& state){

    for(tallyIterator i = tallies_endPart.begin();
	i != tallies_endPart.end(); ++i)
      (*i)->tally_endPart(nhist,kpar,state);
  }
  
  inline void run_localEdep(const unsigned long long nhist,
			const pen_KPAR kpar,
			const stateType& state,
			const double dE){

    for(tallyIterator i = tallies_localEdep.begin();
	i != tallies_localEdep.end(); ++i)
      (*i)->tally_localEdep(nhist,kpar,state,dE);
  }
  
  inline void run_step(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const stateType& state,
		       const tally_StepData& stepData){

    for(tallyIterator i = tallies_step.begin();
	i != tallies_step.end(); ++i)
      (*i)->tally_step(nhist,kpar,state,stepData);    
  }
  
  inline void run_interfCross(const unsigned long long nhist,
			      const unsigned kdet,
			      const pen_KPAR kpar,
			      const stateType& state){

    for(tallyIterator i = tallies_interfCross.begin();
	i != tallies_interfCross.end(); ++i)
      (*i)->tally_interfCross(nhist,kdet,kpar,state);
  }
  
  inline void run_matChange(const unsigned long long nhist,
			    const pen_KPAR kpar,
			    const stateType& state,
			    const unsigned prevMat){

    for(tallyIterator i = tallies_matChange.begin();
	i != tallies_matChange.end(); ++i)
      (*i)->tally_matChange(nhist,kpar,state,prevMat);
  }
  
  inline void run_jump(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const stateType& state,
		       const double ds){

    for(tallyIterator i = tallies_jump.begin();
	i != tallies_jump.end(); ++i)
      (*i)->tally_jump(nhist,kpar,state,ds);
  }
  
  inline void run_knock(const unsigned long long nhist,
			const pen_KPAR kpar,
			const stateType& state,
			const int icol){

    for(tallyIterator i = tallies_knock.begin();
	i != tallies_knock.end(); ++i)
      (*i)->tally_knock(nhist,kpar,state,icol);
  }

  inline void run_lastHist(const unsigned long long lastHist){

    for(tallyIterator i = tallies_lastHist.begin();
	i != tallies_lastHist.end(); ++i)
      (*i)->tally_lastHist(lastHist);
  }
  
  inline void saveData(const unsigned long long nhist, const bool doflush = true){

    for(tallyIterator i = tallies.begin();
	i != tallies.end(); ++i){
      if(doflush)
	(*i)->flush();
      (*i)->saveData(nhist);
    }
  }

  inline void exportImage(const unsigned long long nhist,
			  const pen_imageExporter::formatTypes format,
			  const bool doflush = true){

    for(tallyIterator i = tallies.begin();
	i != tallies.end(); ++i){
      if(doflush)
	(*i)->flush();
      (*i)->exportImage(nhist,format);
    }
  }  

  int writeDump(unsigned char*& pdump,
		size_t& dim,
		const pen_rand& random ,
		const unsigned verbose = 0) {

    if(tallies.size() < 1){
      if(verbose > 0){
	printf("commontallyCluster: writeDump: Error: Any tally has been registered.\n");
      }
      return -1;
    }
    
    //Create a pointers to store dumped data of each tally
    std::vector<unsigned char*> vdumps;
    std::vector<size_t> dumpsDim;
    vdumps.resize(tallies.size());
    dumpsDim.resize(tallies.size());

    size_t dumpSize = 0;
    dumpSize += 2*sizeof(int32_t); //To store last random seeds
    //Get seeds
    int32_t seeds[2]; 
    random.getSeeds(seeds[0],seeds[1]);
      
      
    //Write all dumps and get its size, name and ID size
    for(unsigned i = 0; i < tallies.size(); i++){
      int err;
      err = tallies[i]->writeDump(vdumps[i],dumpsDim[i],verbose);
      if(err != PEN_DUMP_SUCCESS){
	//Free allocated data
	for(unsigned j = 0; j < i; j++){
	  free(vdumps[i]);
	  vdumps[i] = nullptr;
	}
	if(verbose > 0){
	  printf("Error dumping data of tally %s with name %s (position %d)\n",
		 tallies[i]->readID(),tallies[i]->readName().c_str(), i);
	}
	return -2;
      }
      dumpSize += dumpsDim[i];
      dumpSize += 2*sizeof(uint32_t); //To store name and ID length
      dumpSize += tallies[i]->readName().length(); //To store name
      dumpSize += strlen(tallies[i]->readID()); //To store ID      
    }
    
    //Allocate memory
    if(dumpSize < 1){
      if(verbose > 0){
	printf("commontallyCluster: writeDump: Error: No data to dump.\n");
      }      
      return -3;
    }

    pdump = nullptr;
    pdump = (unsigned char*) malloc(dumpSize);
    if(pdump == nullptr){
      if(verbose > 0){
	printf("commontallyCluster: writeDump: Error: Allocation fail.\n");
      }      
      return -4;
    }

    //Write all data to dump buffer
    size_t pos = 0;

    //  last random seeds
    memcpy(&pdump[pos],seeds,2*sizeof(int32_t));
    pos += 2*sizeof(int32_t);
    
    for(unsigned i = 0; i < tallies.size(); i++){
      uint32_t nameSize = (uint32_t)tallies[i]->readName().length();
      uint32_t IDsize = (uint32_t)strlen(tallies[i]->readID());

      //Save name size
      memcpy(&pdump[pos],&nameSize,sizeof(uint32_t));
      pos += sizeof(uint32_t);
      //Save name
      if(nameSize > 0){
	memcpy(&pdump[pos],tallies[i]->readName().c_str(),nameSize);
	pos += nameSize;
      }

      //Save ID size
      memcpy(&pdump[pos],&IDsize,sizeof(uint32_t));
      pos += sizeof(uint32_t);
      //Save ID
      if(IDsize > 0){
	memcpy(&pdump[pos],tallies[i]->readID(),IDsize);
	pos += IDsize;
      }

      //Store dump
      memcpy(&pdump[pos],vdumps[i],dumpsDim[i]);
      pos += dumpsDim[i];
    }

    //Check dumped data
    if(pos != dumpSize){
      if(verbose > 0){
	printf("commontallyCluster: writeDump: Error: Dumped data size doesn't match with allocated.\n");
	printf("                                dumped: %lu\n",pos);
	printf("                                dumped: %lu\n",dumpSize);
      }
      free(pdump);
      pdump = nullptr;
      return -5;      
    }

    //Free memory
    for(unsigned i = 0; i < tallies.size(); i++)
      free(vdumps[i]);

    //Store final dimension
    dim = dumpSize;
    
    return 0;
  }
  int readDump(const unsigned char* const pdump,
	       size_t& pos,
	       int& seed1, int& seed2,
	       const unsigned verbose = 0){

    if(tallies.size() < 1){
      if(verbose > 0){
	printf("commontallyCluster: readDump: Error: Any tally has been registered.\n");
      }
      return -1;
    }

    //Read seeds
    int32_t seeds[2];
    memcpy(seeds,&pdump[pos],2*sizeof(int32_t));
    pos += 2*sizeof(int32_t);

    seed1 = seeds[0];
    seed2 = seeds[1];

    //Read tally dumped data
    for(unsigned i = 0; i < tallies.size(); i++){
      uint32_t nameSize, IDsize;

      //Read name size
      memcpy(&nameSize,&pdump[pos],sizeof(uint32_t));
      pos += sizeof(uint32_t);      
      if(nameSize < 1){
	if(verbose > 0){
	  printf("commontallyCluster: readDump: Error: Dumped tally %d has null name size.\n",i);
	}
	return -2;
      }
      
      //Compare names
      if(memcmp(&pdump[pos], tallies[i]->readName().c_str(),nameSize) != 0){
	if(verbose > 0){
	  printf("commontallyCluster: readDump: Error: Tally names doesn't match at read tally number %d.\n",i);
	  printf("                             Dumped name: %*s\n",nameSize,&pdump[pos]);
	  printf("                           Expected name: %s\n",tallies[i]->readName().c_str());
	}
	return -2;	
      }
      pos += nameSize;
      

      //Read ID size
      memcpy(&IDsize,&pdump[pos],sizeof(uint32_t));
      pos += sizeof(uint32_t);      
      if(IDsize < 1){
	if(verbose > 0){
	  printf("commontallyCluster: readDump: Error: Dumped tally %d has null ID size.\n",i);
	}
	return -3;
      }

      //Compare IDs
      if(memcmp(&pdump[pos], tallies[i]->readID(),IDsize) != 0){
	if(verbose > 0){
	  printf("commontallyCluster: readDump: Error: Tally IDs doesn't match at read tally number %d.\n",i);
	  printf("                             Dumped ID: %*s\n",IDsize,&pdump[pos]);
	  printf("                           Expected ID: %s\n",tallies[i]->readID());
	}
	return -3;	
      }
      pos += IDsize;      

      int err = tallies[i]->readDump(pdump,pos,verbose);
      if(err != PEN_DUMP_SUCCESS){
	printf("commontallyCluster: readDump: Error reading dumped data of tally %d.\n",i);
	return -4;
      }
    }

    return 0;
  }

  int dump2file(const char* filename,
		const pen_rand& random,
		const unsigned verbose = 0){

    unsigned char* pdump = nullptr;
    size_t dumpSize;

    int err;
    
    //Create dump
    err = writeDump(pdump,dumpSize,random,verbose);
    if(err != 0){
      if(verbose > 0){
	printf("commontallyCluster: dump2file: Error creating dump.\n");
      }
      return -1;
    }

    std::string auxstr(filename);
    std::size_t found = auxstr.find_last_of("/\\");

    std::string dumpFilename;

    
    if(found != std::string::npos){
      dumpFilename = auxstr.substr(0,found+1) + std::string("th") + std::to_string(nthread) + auxstr.substr(found+1);
    }else{
      dumpFilename = std::string("th") + std::to_string(nthread) + filename;
    }

    FILE* fdump = nullptr;
    fdump = fopen(dumpFilename.c_str(), "wb");
    if(fdump == nullptr){
      if(verbose > 0){
	printf("commontallyCluster: dump2file: Error creating dump file %s.\n",dumpFilename.c_str());
      }
      free(pdump);
      return -2;
    }

    //Write data
    size_t nwrite = fwrite(pdump,1,dumpSize,fdump);
    if(nwrite != dumpSize){
      if(verbose > 0){
	printf("commontallyCluster: dump2file: Failed to write all data.\n");
      }
      free(pdump);
      return -3;
    }

    //Close file
    fclose(fdump);
    
    //Free memory
    free(pdump);
    
    return 0;
  }

  int readDumpfile(const char* filename,
		   int& seed1, int& seed2,
		   const unsigned verbose = 0){

    unsigned char* pdump = nullptr;
    size_t dumpSize;

    int err;

    //Read data
    FILE* fdump = nullptr;
    fdump = fopen(filename, "rb");
    if(fdump == nullptr){
      if(verbose > 0){
	printf("commontallyCluster: readDumpfile: Error opening dump file %s.\n",filename);
      }
      return -2;
    }

    // Get file size
    //***************

    unsigned char buffer[5000];
    dumpSize = 0;
    size_t nread = 0;
    while((nread = fread(buffer,sizeof(unsigned char),5000,fdump)) > 0){
      dumpSize += nread;
    }

    //Return to file beginning
    rewind(fdump);

    // Allocate memory to store the entire file
    //*******************************************

    pdump = (unsigned char*) malloc(dumpSize);
    if(pdump == nullptr){
      if(verbose > 0){
	printf("commontallyCluster: readDumpfile: Bad alloc.\n");
      }
      return -4;
    }
    
    // Read data
    //*************
    nread = fread(pdump,1,dumpSize,fdump);
    if(nread != dumpSize){
      if(verbose > 0){
	printf("commontallyCluster: readDumpfile: Failed to read all data file.\n");
      }
      free(pdump);
      return -5;
    }

    //Close file
    fclose(fdump);


    // Extract data
    //***************
    
    nread = 0;
    err = readDump(pdump,nread,seed1,seed2,verbose);
    if(err != 0){
      if(verbose > 0){
	printf("commontallyCluster: readDumpfile: Failed to extract data.\n");
      }
      free(pdump);
      return -6;
    }

    //Free memory
    free(pdump);
    
    return 0;
  }
  
  ~pen_specificTallyCluster(){clear();}
  
};

//Include defined tallies
#include "genericTallies.hh"
#include "specificTallies.hh"

#endif
