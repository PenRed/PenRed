
//
//
//    Copyright (C) 2020-2024 Universitat de València - UV
//    Copyright (C) 2020-2024 Universitat Politècnica de València - UPV
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

#ifndef __PEN_VR__
#define __PEN_VR__

#include <cstdio>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

#include "pen_classes.hh"
#include "pen_auxiliar.hh"
#include "pen_math.hh"
#include "pen_data.hh"
#include "instantiator.hh"
#include "pen_states.hh"

// Register macros
///////////////////////

#define DECLARE_VR(Class)			\
  public:					\
  static int registerStatus();			\
  const char* readID() const;			\
  static const char* ___ID;			\
  static volatile const int ___register_return;	\
  private:

#define REGISTER_VR(Class,stateType, ID)				\
  volatile const int Class::___register_return = registerGenericVR<Class,stateType>(static_cast<const char*>(#ID)); \
  const char* Class::___ID = static_cast<const char *>(#ID);		\
  int Class::registerStatus() { return ___register_return;}		\
  const char* Class::readID() const { return ___ID;}

// Enumeration flags
//////////////////////////

enum __usedVRFunc{
  VR_USE_PARTICLESTACK  = 1 << 0,
  VR_USE_MATCHANGE      = 1 << 1,
  VR_USE_INTERFCROSS    = 1 << 2,
  VR_INVALID_FUNCTION   = 1 << 3
};

template <class stateType>
class pen_genericVR : public penred::logs::logger{

private:
  std::string name;
  const __usedVRFunc usedFunctions;

public:
  pen_genericVR(__usedVRFunc functionsToUse) : name(""),
					       usedFunctions(functionsToUse)
  {}

  virtual const char* readID() const = 0;
  inline const std::string& readName() const {return name;}
  inline void setName(const char* newName) { name.assign(newName);}
  inline void setName(const std::string& newName) { name.assign(newName);}

  inline __usedVRFunc readFlags() const {return usedFunctions;}
  
  inline bool has_particleStack() const{
    return ((usedFunctions & VR_USE_PARTICLESTACK) == VR_USE_PARTICLESTACK);
  }
  inline bool has_matChange() const{
    return ((usedFunctions & VR_USE_MATCHANGE) == VR_USE_MATCHANGE);
  }
  inline bool has_interfCross() const{
    return ((usedFunctions & VR_USE_INTERFCROSS) == VR_USE_INTERFCROSS);
  }

  inline void getUsedMethods(std::vector<std::string>& funcNames) const{
    if(has_particleStack())
      funcNames.push_back(std::string("particleStack"));      
    if(has_matChange())
      funcNames.push_back(std::string("matChange"));
    if(has_interfCross())
      funcNames.push_back(std::string("interfCross"));
  }

  virtual void vr_particleStack(const unsigned long long /*nhist*/,
				const pen_KPAR /*kpar*/,
				const unsigned /*kdet*/,
				stateType& /*state*/,
				std::vector<stateType>& /*stack*/,
				unsigned& /*created*/,
				const unsigned /*available*/,
				pen_rand& /*random*/) const{
    char error[1000];
    sprintf(error,"Error on VR %s of type %s: \n Trying to call 'vr_particleStack' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();    
  }

  virtual void vr_matChange(const unsigned long long /*nhist*/,
			    const pen_KPAR /*kpar*/,
			    const unsigned /*prevMat*/,			    
			    stateType& /*state*/,
			    std::vector<stateType>& /*stack*/,
			    unsigned& /*created*/,
			    const unsigned /*available*/,
			    pen_rand& /*random*/) const{
    char error[1000];
    sprintf(error,"Error on VR %s of type %s: \n Trying to call 'vr_matChange' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();    
  }

  virtual void vr_interfCross(const unsigned long long /*nhist*/,
			      const pen_KPAR /*kpar*/,
			      const unsigned /*kdet*/,
			      stateType& /*state*/,
			      std::vector<stateType>& /*stack*/,
			      unsigned& /*created*/,
			      const unsigned /*available*/,
			      pen_rand& /*random*/) const{
    char error[1000];
    sprintf(error,"Error on VR %s of type %s: \n Trying to call 'vr_interfCross' but has not been implemented.\n",readName().c_str(),readID());
    printf("%s\n",error);
    fflush(stdout);
    throw std::bad_function_call();
  }


  virtual int configure(const pen_parserSection& /*config*/,
			const wrapper_geometry& /*geometry*/,
			const unsigned /*verbose*/) = 0;

  virtual ~pen_genericVR(){};
};

namespace penred{
  namespace vr{
    
    //Check registered types
    template <class stateType>
    bool checkRegistered(const unsigned verbose);
  }
}

// Auxiliar functions
///////////////////////

template<class stateType>
bool __VRSort (pen_genericVR<stateType>* i,
	       pen_genericVR<stateType>* j) {
  return (i->readName().compare(j->readName()) < 0);
}

// VR cluster classes
//////////////////////////
template<class stateType>
class pen_VRCluster : public abc_VR<stateType>{
private:
protected:
  typedef pen_genericVR<stateType>* pVRtype;
  typedef typename std::vector<pVRtype>::const_iterator VRIterator;

  //VR instantiators
  static instantiator<pen_genericVR<stateType>>& genericVRs();

  //Registered VR vector
  std::vector<pVRtype> VRs;

  //Vectors containing VR pointers with specific execution functions
  std::vector<pVRtype> VRs_particleStack;
  std::vector<pVRtype> VRs_matChange;
  std::vector<pVRtype> VRs_interfCross;

  int configStatus;

  int createVR(const char* ID,
	       const char* VRname,
	       const wrapper_geometry& geometry,
	       const pen_parserSection& config,
	       const unsigned verbose);
public:
  std::string name;

  pen_VRCluster() : configStatus(0),
		    name("unnamed")
  {}
  
  inline int configureStatus() const {return configStatus;}
  
  inline unsigned numVR() const {return VRs.size();}

  static std::string typesList(){
    return genericVRs().typesList();
  }

  static std::string typesList(std::vector<std::string>& list){
    return genericVRs().typesList(list);
  }

  template <class subclass>
  static int addVR(const char* typeID){
    return genericVRs().template addSubType<subclass>(typeID);
  }

  void configure(const pen_parserSection& config,
		 const wrapper_geometry& geometry,
		 const unsigned verbose);

  inline void run_particleStack(const unsigned long long nhist,
				const pen_KPAR kpar,
				const unsigned kdet,
				stateType& state,
				std::vector<stateType>& stack,
				unsigned& created,
				const unsigned available,
				pen_rand& random) const{
    for(VRIterator i = VRs_particleStack.begin();
	i != VRs_particleStack.end(); ++i)
      (*i)->vr_particleStack(nhist,kpar,kdet,state,stack,
			     created,available,random);
  }

  inline void run_matChange(const unsigned long long nhist,
			    const pen_KPAR kpar,
			    const unsigned prevMat,			    
			    stateType& state,
			    std::vector<stateType>& stack,
			    unsigned& created,
			    const unsigned available,
			    pen_rand& random) const{
    for(VRIterator i = VRs_matChange.begin();
	i != VRs_matChange.end(); ++i)
      (*i)->vr_matChange(nhist,kpar,prevMat,state,stack,
			 created,available,random);
  }

  inline void run_interfCross(const unsigned long long nhist,
			      const pen_KPAR kpar,
			      const unsigned kdet,
			      stateType& state,
			      std::vector<stateType>& stack,
			      unsigned& created,
			      const unsigned available,
			      pen_rand& random) const{
    for(VRIterator i = VRs_interfCross.begin();
	i != VRs_interfCross.end(); ++i)
      (*i)->vr_interfCross(nhist,kpar,kdet,state,stack,
			   created,available,random);
  }

  void clear();
  
  ~pen_VRCluster();
};

template<class stateType>
instantiator<pen_genericVR<stateType>>& pen_VRCluster<stateType>::genericVRs(){
  static instantiator<pen_genericVR<stateType>>* ans =
    new instantiator<pen_genericVR<stateType>>;
  return *ans;    
}

template<class stateType>
int pen_VRCluster<stateType>::createVR(const char* ID,
				       const char* VRname,
				       const wrapper_geometry& geometry,
				       const pen_parserSection& config,
				       const unsigned verbose){

  //Check if tally name already exists
  std::string strName;
  if(VRname == nullptr || strName.assign(VRname).length() == 0){
    if(verbose > 0){
      penred::logs::logger::printf("VRCluster: createVR: Error: empty VR name.\n");
    }
    return -3;
  }
  for(unsigned i = 0; i < VRs.size(); i++){
    if(strName.compare(VRs[i]->readName()) == 0){
      if(verbose > 0){
	penred::logs::logger::printf("VRCluster: createVR: Error: VR name '%s' "
				     "already used.\n", VRname);
      }
      return -4;
    }
  }

  //Create specified VR by ID
  pen_genericVR<stateType>* pVR = nullptr;  
  pVR = genericVRs().createInstance(ID);
  if(pVR == nullptr){
    if(verbose > 0){
      penred::logs::logger::printf("VRCluster: createVR: Error: unable to create a "
				   "VR of type '%s'.\n",ID);
    }
    return -1;
  }

  //Set VR name
  pVR->setName(VRname);

  //Configure VR
  int errConfig = pVR->configure(config,geometry,verbose);
  if(errConfig != 0){
    delete pVR;
    if(verbose > 0){
      penred::logs::logger::printf("VRCluster: createVR: Error: VR "
				   "'%s' of type '%s' "
				   "failed on configuration step.\n",
				   VRname,ID);
    }
    return -2;
  }

  //Append VR pointer to necessary vectors
  VRs.push_back(pVR);

  if(pVR->has_particleStack()){VRs_particleStack.push_back(pVR);}
  if(pVR->has_matChange()){VRs_matChange.push_back(pVR);}
  if(pVR->has_interfCross()){VRs_interfCross.push_back(pVR);}

  //Sort vector VR by name
  std::sort(VRs.begin(), VRs.end(), __VRSort<stateType>);

  return 0;
  
}

template<class stateType>
void pen_VRCluster<stateType>::configure(const pen_parserSection& config,
					 const wrapper_geometry& geometry,
					 const unsigned verbose){
  int err = 0;

  //Clear previous configuration
  clear();

  if(verbose > 1){
    penred::logs::logger::printf("\n------------------------------------\n");
    penred::logs::logger::printf("\n **** VR group '%s'\n",name.c_str());
  }

  //Check registered types to ensure static library linking of the register variable
  if(!penred::vr::checkRegistered<stateType>(verbose)){
    if(verbose > 0){
      printf("Warning: Some VR types are not properly registered\n");
    }
  }

  //Extract VR names
  std::vector<std::string> vrNames;
  config.ls(vrNames);

  //Iterate for each tally name
  for(unsigned i = 0; i < vrNames.size(); ++i){

    if(verbose > 1){
      penred::logs::logger::printf("\n------------------------------------\n\n");
      penred::logs::logger::printf("\nVR '%s':\n\n",vrNames[i].c_str());
    }

    //Get subsection for VR 'i'
    pen_parserSection vrSec;
    if(config.readSubsection(vrNames[i],vrSec) != INTDATA_SUCCESS){
      if(verbose > 0){
	penred::logs::logger::printf("VRCluster: configure: Error: unable to read "
				     "section '%s' to configure VR.\n",
				     vrNames[i].c_str());
      }
      err++;
      continue;
    }

    //Try to read VR type ID
    std::string vrID;
    if(vrSec.read("type",vrID) != INTDATA_SUCCESS){
      if(verbose > 0){
	penred::logs::logger::printf("VRCluster: configure: Error: unable to read "
				     "field %s/type. String expected\n",
				     vrNames[i].c_str());
	err++;
	continue;
      }
    }

    if(verbose > 1){
      penred::logs::logger::printf("VR type: '%s'\n\n",vrID.c_str());
    }

    //Try to create and configure VR 'i'
    if(createVR(vrID.c_str(),vrNames[i].c_str(),
		geometry,vrSec,verbose) != 0){
      if(verbose > 0){
	penred::logs::logger::printf("VRCluster: configure: Error: Unable to "
				     "create and configure VR '%s' of type '%s'.\n",
				     vrNames[i].c_str(),vrID.c_str());
	err++;
	continue;
      }
    }
  }
  if(verbose > 1){
    penred::logs::logger::printf("\n------------------------------------\n\n");
  }

  if(err > 0){
    if(verbose > 0){
      penred::logs::logger::printf("VRCluster: configure: Error: %d VR "
				   "creation failed.\n",err);
    }
  }

  //Print created tallies
  if(verbose > 1){
    penred::logs::logger::printf("\nCreated common VR (name and type):\n\n");
    for(unsigned i = 0; i < VRs.size(); i++){
      penred::logs::logger::printf("  %20s -> %20s\n",
				   VRs[i]->readName().c_str(),
				   VRs[i]->readID());
    }
  }

  configStatus = err;  
}

template<class stateType>
void pen_VRCluster<stateType>::clear(){
  for(std::size_t i = 0; i < VRs.size(); ++i)
    delete VRs[i];
  VRs.clear();
  VRs_particleStack.clear();
  VRs_matChange.clear();
  VRs_interfCross.clear();

  configStatus = 0;  
}

template<class stateType>
pen_VRCluster<stateType>::~pen_VRCluster(){clear();}

template <class subclass, class stateType>
static int registerGenericVR(const char* typeID){
  return pen_VRCluster<stateType>::template addVR<subclass>(typeID);
}

//Include defined VR
#include "../generic/includes/genericVR.hh"
#include "../specific/includes/specificVR.hh"

#endif
