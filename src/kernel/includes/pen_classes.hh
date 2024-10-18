
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
//        vicente.gimenez@uv.es
//    
//


#ifndef __PEN_CLASSES__
#define __PEN_CLASSES__

#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <type_traits>
#include <array>
#include <functional>
#include <algorithm>
#include <numeric>
#include <memory>

#include "../states/pen_baseState.hh"
#include "pen_constants.hh"
#include "../errors/pen_errors.hh"
#include "../rands/includes/pen_randoms.hh"
#include "../parsers/internalData/includes/pen_data.hh"

#include "math_classes.hh"

template<class particleType, class contextType, class materialType> class abc_interaction;
template<class stateType, class contextType> class pen_particle;
template<class stateType> class pen_particleStack;
class wrapper_geometry;

//-------------------
// Materials
//-------------------

class abc_material : public penred::logs::logger{

private:
  
protected:
  
  bool initialized;

  //Stores the required identifier to get the data to construct this material.
  //For example, a filename or an ID in a data base
  std::string dataPath;
public:

  //  ----  Material densities and its reciprocals.
  double DEN;    //Density
  double RDEN;   //1/Density
  
  //  ----  Absorption energies, EABS(KPAR,MAT).
  double EABS[constants::nParTypes];  
  
  abc_material() : initialized(false), DEN(1.0), RDEN(1.0){
    for(unsigned i = 0; i < constants::nParTypes; i++){
      EABS[i] = 50.0; //eV
    }
  }
  inline void setDataPath(const std::string& newDataPath) { dataPath = newDataPath; }
  inline std::string readDataPath() const { return dataPath; }
  inline bool initDone() const { return initialized;}
  inline double getEABS(const unsigned kpar) const {return EABS[kpar];}
  inline void setEABS(const unsigned kpar, const double eabs){
    if(eabs > 0.0){
      EABS[kpar] = eabs;
    }
  }
  inline void setDens(const double dens){
    if(dens > 1.0e-17){
      DEN = dens;
      RDEN = 1.0/dens;
    }
  }

  inline double readDens() const {return DEN;}
  inline double readIDens() const {return RDEN;}

  // ** Functions to be implemented in each class type
  virtual double meanZ() const = 0;
  virtual void composition(std::vector<std::pair<unsigned, double>>&) const = 0;
  
  virtual ~abc_material(){};
  
};

//-------------------
// Geometry
//-------------------

class wrapper_geometry : public penred::logs::logger{

protected:
  int configStatus;

public:

  std::string name;
    
  wrapper_geometry() : configStatus(0),
		       name("unnamed") {}

  inline int configureStatus() const {return configStatus;}
  inline virtual const char* getType() const {return "UNKNOWN";}
  
  virtual void locate(pen_particleState& state) const = 0;
  virtual void step(pen_particleState& state, double DS, double &DSEF, double &DSTOT, int &NCROSS) const = 0;
  virtual int configure(const pen_parserSection& config, const unsigned verbose) = 0;
  virtual void usedMat(bool[constants::MAXMAT+1]) const = 0;
  virtual double getEabs(const unsigned ibody, const unsigned kpar) const = 0;
  virtual double getDSMAX(const unsigned ibody) const = 0;
  virtual unsigned getDET(const unsigned ibody) const = 0;
  virtual unsigned getMat(const unsigned ibody) const = 0;
  virtual unsigned long getElements() const = 0;
  virtual unsigned getBodies() const = 0;
  virtual unsigned getIBody(const char* elementName) const = 0;
  virtual std::string getBodyName(const unsigned ibody) const = 0;
  virtual void getOffset(double* offset) const { offset[0] = 0.0; offset[1] = 0.0; offset[2] = 0.0; }

  //Method to try to convert the wrapper to a specific geometry type
  template<class geoType>
  const geoType* convertTo() const {
    return dynamic_cast<const geoType*>(this);
  }


  //Methods for nested geometries
  virtual size_t nInternalGeometries() const { return 0; }
  virtual const wrapper_geometry* getInternalGeo(const size_t /*index*/) const { return nullptr; }
  
  template<class geoType>
  const geoType* getInternalGeoType (size_t& geoPos, const size_t firstPos = 0) const {
    const size_t nInternal = nInternalGeometries();
    if(firstPos >= nInternal)
      return nullptr;
    for(size_t i = firstPos; i < nInternal; ++i){
      //Get internal geometry pointer
      const wrapper_geometry* internalGeo = getInternalGeo(i);
      //Try to convert it to the desired class
      const geoType* pGeoType = internalGeo->convertTo<geoType>();
      if(pGeoType != nullptr){
	geoPos = i;
	return pGeoType;
      }
    }
    geoPos = nInternal;
    return nullptr;
  }
  
  virtual ~wrapper_geometry(){}

  
};

//-------------------
// Variance reduction
//-------------------

template<class stateType>
class abc_VR : public penred::logs::logger{

public:
  virtual void run_particleStack(const unsigned long long nhist,
				 const pen_KPAR kpar,
				 const unsigned kdet,
				 stateType& state,
				 std::vector<stateType>& stack,
				 unsigned& created,
				 const unsigned available,
				 pen_rand& random) const = 0;

  virtual void run_matChange(const unsigned long long nhist,
			     const pen_KPAR kpar,
			     const unsigned prevMat,			    
			     stateType& state,
			     std::vector<stateType>& stack,
			     unsigned& created,
			     const unsigned available,
			     pen_rand& random) const = 0;

  virtual void run_interfCross(const unsigned long long nhist,
			       const pen_KPAR kpar,
			       const unsigned kdet,
			       stateType& state,
			       std::vector<stateType>& stack,
			       unsigned& created,
			       const unsigned available,
			       pen_rand& random) const = 0;

  virtual ~abc_VR(){}
};

//-------------------
// Interactions
//-------------------

template<class particleType, class contextType, class materialType>
class abc_interaction{

 protected:
  const int ID;

 public:

  abc_interaction() : ID(-1) {}
  abc_interaction(const int id) : ID(id) {}
  virtual void init(const contextType&){};
  virtual double iMeanFreePath(const materialType&, const particleType&) const = 0;
  virtual int interact(const contextType&, const materialType&, particleType&, double&, pen_rand&) const = 0;
  virtual int interactF(const contextType&, const materialType&, particleType&, double&, pen_rand&) const = 0;  
  inline int getID() const {return ID;}
  virtual ~abc_interaction(){}
};

//-------------------
// Particles
//-------------------

template<class stateType, class contextType, class materialType>
class abc_particle{

private:
  std::vector<pen_particleState> genericStates;
  std::vector<stateType> specificStates;
protected:

  //Particle energy at previous JUMP
  double ELAST1;

  //  **** Particle interaction constants
  
  double P[constants::MAXINTERACTIONS], ST;

  // ****  Energy grid variables for the next knock call.
  unsigned int KE;
  double XEL, XE, XEK;

  // **** Variance reduction variables

  double P0[constants::MAXINTERACTIONS];
  bool   LFORC[constants::MAXINTERACTIONS];

  //     KSOFTI ... take value 1 (0) if soft energy loss is (not) required
  int KSOFTI;

  //Particle state (This variable will be used to store particles in the stack)
  stateType state;

  //Particle current material pointer
  const materialType* pmat;
  
  //Particle last movement variables
  double dsef;   //Distance traveled on same material (Non void zones)
  double dstot;  //Total traveled distance (Included void zones)
  int ncross; //Is non zero if some interface has been crossed

  //Particle previous position variables
  unsigned MATL;
  unsigned IBODYL;
  double XL;
  double YL;
  double ZL;

  //Particle local parameters
  double EABS;
  double DSMAXbody;
  unsigned KDET;

  //Generic VR cluster
  const abc_VR<pen_particleState>* genericVR;
  //Specific VR cluster
  const abc_VR<stateType>* specificVR;
  //Reference to particle stack
  pen_particleStack<stateType>& stack;
  
public:

  typedef stateType typeState;
  typedef contextType typeContext;
  
  const contextType& context;
  const pen_KPAR kpar;
  const unsigned int interactions;

  //Energy deposited at material on particle annihilation at rest (eV)
  const double annihilationEDep;
  //  ---- 

  abc_particle(const contextType& contextIn,
	       const pen_KPAR KPAR,
	       const unsigned int nInt,
	       const double annihilationEDepIn,
	       pen_particleStack<stateType>& stackIn) : KSOFTI(0),
							pmat(nullptr),
							dsef(0.0),
							dstot(0.0),
							ncross(0),
							MATL(0),
							IBODYL(0),
							XL(0.0),
							YL(0.0),
							ZL(0.0),
							EABS(0.0),
							DSMAXbody(1.0e35),
							KDET(0),
							genericVR(nullptr),
							specificVR(nullptr),
							stack(stackIn),
							context(contextIn),
							kpar(KPAR),
							interactions(nInt),
							annihilationEDep(annihilationEDepIn)
  {
    //Check if kpar is in range [0,nParTypes)
    if(KPAR < 0 || KPAR >= constants::nParTypes)
      {
	char error[300];
	sprintf(error,"Particle kpar (%d) out of range [0,%d).\n  Define an identifier for this particle at file 'kernel/particles/includes/pen_particles_ID.hh'",KPAR,constants::nParTypes);
	throw std::out_of_range(error);      
      }

    if(interactions > constants::MAXINTERACTIONS){
	char error[300];
	sprintf(error,"Number of interactions of particle with kpar (%d) out of range [0,%d)\n  Increment constant 'MAXINTERACTIONS' in constants namespace",KPAR,constants::MAXINTERACTIONS);
	throw std::out_of_range(error);      
    }

    //Check if context material is compatible
    if(!contextIn.template compatible<materialType>()){
      char error[300];
      sprintf(error,"Incompatible material type at particle instantation (kpar = %d).",KPAR);
      throw std::invalid_argument(error);
    }

    //Reserve space in the auxiliary stacks
    genericStates.resize(constants::NMS);
    specificStates.resize(constants::NMS);
  }

  virtual void START() = 0;
  virtual void JUMP(double &DS, pen_rand& penRand, const double DSMAX) = 0;
  virtual void JUMPF(double &DS, pen_rand& penRand, const double DSMAX) = 0;
  virtual void KNOCK(double &DE, int &ICOL, pen_rand& penRand) = 0;
  virtual void KNOCKF(double &DE, int &ICOL, pen_rand& penRand) = 0;

  virtual void softEloss(double& X,
			 double& Y,
			 double& Z,
			 double& DE,
			 pen_rand& penRand) = 0;
  
  virtual void dpage() = 0;
  virtual void page0() = 0;

  inline unsigned getDET() const {return KDET;}
  inline double getDSMAX() const {return DSMAXbody;}
  
  inline double getEABS() const {return EABS;}
  
  inline double DSef() const{return dsef;}
  inline double DStot() const{return dstot;}
  inline int NCross() const{return ncross;}

  inline void updateMat(){
    if(state.MAT > 0)
      pmat = context.template readMaterialPointer<materialType>(state.MAT-1);
    else
      pmat = nullptr;
  }
  inline void setMat(const unsigned imat){    
    state.MAT = imat;
    updateMat();
  }

  inline void updateBody(){
    if(state.IBODY < context.readGeometry()->getBodies()){
      EABS      = context.getEABS(state.IBODY,kpar);
      DSMAXbody = context.readGeometry()->getDSMAX(state.IBODY);
      KDET      = context.readGeometry()->getDET(state.IBODY);
    }
    else{
      EABS = 0.0;
      DSMAXbody = 1.0e35;
      KDET = 0;
    }
  }
  
  
  inline void setStep(const double dsefIn,
		      const double dstotIn,
		      const int ncrossIn){
    dsef = dsefIn;
    dstot = dstotIn;
    ncross = ncrossIn;
  }
  
  inline unsigned lastMat() const {return MATL;}
  inline unsigned lastBody() const {return IBODYL;}
  inline void lastPos(double& X, double& Y, double& Z){
    X = XL; Y = YL; Z = ZL;
  }

  inline void setLastMat(const unsigned lastMatIn) {MATL = lastMatIn;}
  inline void setLastBody(const unsigned lastBodyIn) {IBODYL = lastBodyIn;}
  inline void setLastPos(const double X, const double Y, const double Z) {
    XL = X; YL = Y; ZL = Z;
  }
  inline void saveLastPos(){
    XL = state.X; YL = state.Y; ZL = state.Z;
  }
  
  inline void jumpVolume(){

    //Store last position information
    IBODYL = state.IBODY;
    MATL = state.MAT;
    
    XL = state.X;
    YL = state.Y;
    ZL = state.Z;

    //Jump until particle crosses some interface or scapes from the geometry
    context.readGeometry()->step(state,1.0e30,dsef,dstot,ncross);
    //Calculate new particle age 
    if(state.LAGE){dpage();}

    if(MATL != state.MAT)
      updateMat();
    if(IBODYL != state.IBODY)
      updateBody();
  }
  
  inline void move(const double ds,
		   double& de,
		   double& softX,
		   double& softY,
		   double& softZ,
		   pen_rand& penRand){

    //Save actual position
    IBODYL = state.IBODY;
    MATL = state.MAT;
    
    XL = state.X;
    YL = state.Y;
    ZL = state.Z;

    //Move the particle
    context.readGeometry()->step(state, ds, dsef, dstot, ncross);
    //Calculate new particle age 
    if(state.LAGE){dpage();}

    //Check if material needs to be updated
    if(MATL != state.MAT)
      updateMat();
    if(IBODYL != state.IBODY)
      updateBody();
    
    //Check if soft energy deposition is required
    if(reqSoftELoss() == 1){
      //Calculate soft energy loss
      softEloss(softX,softY,softZ,de,penRand);
    }else{
      de = 0.0;
    }
  }
  
  virtual void annihilate(pen_rand&){}
  

  inline double xel() const {return XEL;}
  inline double xe() const {return XE;}
  inline double xek() const {return XEK;}
  inline unsigned ke() const {return KE;}

  inline void getGrid(unsigned& ke,
		      double& xel,
		      double& xe,
		      double& xek) const {
    ke = KE;
    xel = XEL; xe = XE; xek = XEK;
  }
  
  inline int reqSoftELoss() const {return KSOFTI;}
  
  inline pen_KPAR getKpar() const {
    return kpar;
  }

  inline void setBaseState(const pen_particleState& newState){
    state = newState;
  }
  
  inline void setState(const stateType& newState){
    state = newState;
  }

  inline pen_particleState& getBaseState(){
    return state;
  }
  
  inline stateType& getState(){
    return state;
  }

  inline const pen_particleState& readBaseState() const {
    return state;
  }
  
  inline const stateType& readState() const{
    return state;
  }

  inline unsigned nStacked() const{
    return stack.getNSec();
  }

  inline void setStateFromStack(){
    
    stack.get(state);

    //Update body and material
    updateBody();
    updateMat();
    
  }

  inline const pen_particleStack<stateType>& readStack() const { return stack; }
  
  inline const contextType& readContext() const {return context;}

  inline void registerGenericVR(const abc_VR<pen_particleState>& vrIn){
    genericVR = &vrIn;
  }

  inline void registerSpecificVR(const abc_VR<stateType>& vrIn){
    specificVR = &vrIn;
  }

  void vr_particleStack(const unsigned long long nhist,
			pen_rand& random,
			const unsigned verbose);  

  double vr_matChange(const unsigned long long nhist,
		      pen_rand& random,
		      const unsigned verbose);

  double vr_interfCross(const unsigned long long nhist,
			pen_rand& random,
			const unsigned verbose);
  
  void baseClear(){
    KSOFTI = 0;
    pmat = nullptr;
    dsef = 0.0;
    dstot = 0.0;
    ncross = 0;
    MATL = 0;
    IBODYL = 0;
    XL = 0.0;
    YL = 0.0;
    ZL = 0.0;
    state.reset();
  }
  
  virtual ~abc_particle(){}
  
  //CJUMP0
  //CJUMP1
  //CEGRID (only XEL, XE, XEK, and KE)  
};

//-------------------
// Particle stacks
//-------------------

class abc_particleStack{

 protected:
  unsigned int NSEC;
  
 public:
  constexpr abc_particleStack() : NSEC(0) {}

  inline unsigned int getNSec() const {return NSEC;}
  inline void cleans(){NSEC = 0;}

  virtual pen_particleState readBaseState(const unsigned i) const = 0;
};

template<class stateType> class pen_particleStack : public abc_particleStack{

 protected:
  std::vector<stateType> states;
 public:
  
  pen_particleStack() : abc_particleStack() {
    states.resize(constants::NMS);
  }

  void store(const stateType& state)
  {    
    if(NSEC < constants::NMS)
      {
	states[NSEC] = state;
	NSEC++;
      }
    else
      {
	penred::logs::logger::printf(penred::logs::SIMULATION,
				     "pen_particleStack:store:Warning: Stack full\n");
	//Stack is full remove particle with less energy
	unsigned int lessEpos = 0;
	double minE = 1.0e35;
	for(unsigned int i = 0; i < NSEC; i++)
	  {
	    if(states[i].E < minE)
	      {
		lessEpos = i;
		minE = states[i].E;
	      }
	  }
	//Store new particle in less energy position
	if(minE < state.E)
	  {
	    states[lessEpos] = state;
	  }
	//Create warning
	penError(ERR_store_FULL_STACK);
      }
  }
  
  inline unsigned get(stateType& state){
    if(NSEC > 0){
      NSEC--;
      state = states[NSEC];
      return NSEC;
    }
    return 0;
  }

  inline stateType readState(const unsigned i) const {
    return states[i];
  }

  inline pen_particleState readBaseState(const unsigned i) const{
    return pen_particleState(readState(i));
  }
  
  //CERSEC
  //SECST
  
};

//-------------------
// Context
//-------------------

class wrapper_context : public penred::logs::logger{

public:
  
  // Configuration functions
  //-------------------------

  //Configuration with no geometry set
  virtual int configure(const double EMAX,
			const pen_parserSection& config,
			pen_parserSection& matInfo,
			const unsigned verbose) = 0;

  //Configuration with geometry set
  virtual int configureWithGeo(const pen_parserSection& config,
			       const unsigned verbose) = 0;
  //-----------------------------

  // Functions to get particle transport characteristics 
  //----------------------------------------------------------
  
  //Function 'range' is intended to return the range in the specified
  //material for a particle with the specified energy and type
  virtual double range(const double E, const pen_KPAR kpar, const unsigned M) const = 0;

  inline void findRange(double& lowE, double& topE, const double objectiveRange,
			const pen_KPAR kpar, const unsigned M,
			const unsigned long maxTries = 1000000, const double tol = 0.001){
    
    //Function 'findRange' returns the energy interval corresponding
    //to the specified range. This one is calculated for one particle
    //type and material
    //
    // * Input: 
    // lowE : minimum inital energy
    // topE : maximum initial energy
    // objectiveRange: Range to be found
    // kpar : particle type
    // M    : Material index in context [0,nMaterials)
    // maxTries : Maximum iterations to found the objective range
    // tol  : Maximum relative difference between the final lowE and topE.
    //
    // * Output:
    // lowE : final minimum energy of the resulting interval
    // topE : final maximum energy of the resulting interval

    
    unsigned nTries = 0;
    const double maxRatio = 1.0 + tol;
    do{
      double midE = (topE+lowE)/2.0;
      double r = range(midE,kpar,M);
      if(r == objectiveRange){
	topE = midE;
	lowE = midE;
      }
      else if(r > objectiveRange){
	topE = midE;
      }else{
	lowE = midE;
      }
      ++nTries;
    }while(nTries < maxTries && topE/lowE > maxRatio);    
  }

  virtual double avncol(const double E,
			const pen_KPAR kpar,
			const int icol,
			const unsigned imat) const = 0;

  virtual double avninter(const double E,
			  const pen_KPAR kpar,
			  const int icol,
			  const unsigned imat,
			  const bool calc_piecewise) const = 0;

  //Interaction forcing set function  
  virtual void setForcing(const double forcerIn,
			  const pen_KPAR kpar,
			  const int icol,
			  const unsigned ibody,
			  const double weightL,
			  const double weightU) = 0;    

  //----------------------------------------------------------
  
};

template <class baseMat>
class abc_context : public wrapper_context{
public:

  typedef baseMat matType;
  
private:

  //Array with cutoff energies
  double* maxEABS;

  unsigned geoBodies; // Controls geometry number of bodies

  // Geometry
  const wrapper_geometry* geometry;

  bool matsSet; //Controls if materials has been set

  // Materials
  baseMat* materials[constants::MAXMAT];
  unsigned nMats;

  void clearMats(){
    for(unsigned i = 0; i < constants::MAXMAT; i++){
      if(materials[i] != nullptr){
	delete materials[i];
	materials[i] = nullptr;
      }
    }
    matsSet = false;
    nMats = 0;
  }
  void clearGeo(){
    if(maxEABS != nullptr){
      delete [] maxEABS;
      maxEABS = nullptr;
    }
    geoBodies = 0;
    geometry = nullptr;
  }
  void clear(){
    clearMats();
    clearGeo();
  }
  
protected:

  
public:

  const static unsigned int NBV = pen_geoconst::NB;

  abc_context() : maxEABS(nullptr),
		  geoBodies(0),
		  geometry(nullptr),
		  matsSet(false),
		  nMats(0){
    
    //Set material pointers to null
    for(unsigned int i = 0; i < constants::MAXMAT; i++)
      materials[i] = nullptr;
  }
  
  inline void getMatBaseArray(const abc_material* mats[constants::MAXMAT]) const {

    for(unsigned i = 0; i < constants::MAXMAT; i++){
      mats[i] = materials[i];
    }
  }

  inline double getMatEABS(const unsigned imat, const unsigned kpar){
    if(imat >= nMats){
      char error[300];
      sprintf(error,"getMatEABS: %d exceeds number of available materials (%d).",imat,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return materials[imat]->getEABS(kpar);
    }
  }
  inline double getEABS(const unsigned ibody, const unsigned kpar) const{
      unsigned index = kpar+ibody*constants::nParTypes;
      return maxEABS[index];
  }
  inline unsigned getDET(const unsigned ibody) const {
    return geometry->getDET(ibody);
  }
  inline double getDSMAX(const unsigned ibody) const {
    return geometry->getDSMAX(ibody);
  }
  
  int setGeometry(const wrapper_geometry* geoIn, const unsigned verbose = 2){
    
    if(geoIn == nullptr){
      if(verbose > 0){
	printf("setGeometry: Error: Null geometry\n");
      }
      return -1;
    }

    //Clear previous geometry
    clearGeo();
    
    //Save pointer
    geometry = geoIn;

    //Get number of bodies in geometry
    geoBodies = geometry->getBodies();
    
    //Allocate memory for absorption energies of each body
    maxEABS = new double[(geoBodies+1)*constants::nParTypes];
    if(maxEABS == nullptr){
      if(verbose > 0){
	printf("setGeometry: Error: Unable to allocate materials");
      }
      return -2;
    }

    //Get materials used by the current geometry
    bool usedMat[constants::MAXMAT+1];
    geometry->usedMat(usedMat);

    //Ensure all required materials have been defined
    for(unsigned i = getNMats()+1; i < constants::MAXMAT+1; ++i){
      if(usedMat[i]){
	if(verbose > 0){
	  printf("setGeometry: Error: Geometry uses material %u, which has "
		 "not been defined in the context configuration.\n",i);
	}
	return -3;
      }
    }

    //Update absorption energies ussing geometry information
    int err = updateEABS();
    if(err != 0){
      if(verbose > 0){
	printf("setGeometry: Error updating absorption energies\n");
      }
      return -4;
    }
      
    return 0;
  }
  inline const wrapper_geometry* readGeometry() const {return geometry;}
  
  int updateEABS(const unsigned verbose = 2){
    if(geometry == nullptr){
      if(verbose > 0){
	printf("updateEABS: Error: Geometry not set, null pointer stored.\n");
      }
      return -1;
    }

    //Check if the number of elements has been changed
    if(geoBodies != geometry->getBodies()){
      //Geometry has been changed, run set geometry again,
      //where "updateEABS" will be called again
      if(verbose > 1){
	printf("updateEABS: Number of bodies in the geometry has changed. "
	       "Reset it via 'setGeometry'.\n");
      }
      return setGeometry(geometry);
    }

    //Fill absorption energies
    for(unsigned i = 0; i < geoBodies; i++){
      unsigned index0 = i*constants::nParTypes;

      //Get material index for this body
      int mat = geometry->getMat(i);
      if(mat < 0 || mat > (int)nMats){ //Check material bounds
	//Index out of range
	if(verbose > 0){
	  printf("updateEABS: Error: Geometry material %d is out of context "
		 "material range [0,%u]. Check the configuration to set and "
		 "configure the apropiate number of materials for this geometry.\n",
		 mat,nMats);
	}
	return -3;
      }

      //Iterate over particle types
      for(unsigned j = 0; j < constants::nParTypes; j++){
	double Egeo = geometry->getEabs(i,j);
	if(mat < 1){ // Void
	  maxEABS[index0+j] = Egeo;
	}
	else{
	  double Emat = materials[mat-1]->getEABS(j);
	  if(Egeo > Emat){ maxEABS[index0+j] = Egeo;}
	  else{ maxEABS[index0+j] = Emat;}
	}
      }
    }
    

    return 0;
  }
  
  template<class derivedMat>
  int setMats(const unsigned M){

    //Check if materials has been already created
    if(matsSet)
      return -1;
    
    //Check index
    if(M > constants::MAXMAT || M < 1)
      return -2;

    //Clear materials (just in case)
    clearMats();
    
    //Create materials
    for(unsigned i = 0; i < M; i++){
      materials[i] = new derivedMat();
      if(materials[i] == nullptr){
	clearMats();
	return -3;
      }
    }
    nMats = M;
    matsSet = true;
    return 0;
  }

  inline unsigned getNMats() const {return nMats;}
  
  inline int setMatEABS(const unsigned M,
			const unsigned kpar,
			const double eabs){

    if(M >= constants::MAXMAT)
      return -1;

    if(materials[M] == nullptr)
      return -2;

    materials[M]->setEABS(kpar,eabs);
    return 0;
  }

  //Ensure that read and get material uses only convertible classes to baseMat
  template <class derivedMat>
  inline const typename std::enable_if<std::is_convertible<derivedMat,baseMat>::value, derivedMat&>::type
  readMaterial(const unsigned M) const{
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return *(static_cast<derivedMat*>(materials[M]));
    }
  }

  template <class derivedMat>
  inline const typename std::enable_if<std::is_convertible<derivedMat,baseMat>::value, derivedMat*>::type
  readMaterialPointer(const unsigned M) const{
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return static_cast<derivedMat*>(materials[M]);
    }
  }

  
  template <class derivedMat>
  inline typename std::enable_if<std::is_convertible<derivedMat,baseMat>::value, derivedMat&>::type
  getMaterial(const unsigned int M){
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return *(static_cast<derivedMat*>(materials[M]));
    }
  }

  inline const baseMat& readBaseMaterial(const unsigned M) const{
    
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
    
      return *materials[M];
    }
  }

  inline baseMat& getBaseMaterial(const unsigned M){
    if(M >= nMats){
      char error[300];
      sprintf(error,"%d exceeds number of available materials (%d).",M,nMats);
      throw std::out_of_range(error);            
    }
    else{
      return *materials[M];
    }
  }
  
  //Check material compatibility function
  template <class derivedMat>
  bool compatible() const{
    if(!matsSet)
      return false;

    derivedMat* pderived = nullptr;
    pderived = dynamic_cast<derivedMat*>(materials[0]);

    if(pderived == nullptr){
      return false;
    }
    return true;
  }
  
  virtual ~abc_context(){
    clear();
  };
};

template<class contextType>
inline std::shared_ptr<contextType> createContext(){
  static_assert(std::is_base_of<wrapper_context, contextType>::value,
		"Error: 'createContext' cannot create contexts not"
		"derived from 'wrapper_context'");
  static_assert(std::is_base_of<abc_context<typename contextType::matType>, contextType>::value,
		"Error: Context type provided to 'createContext' is not "
		"derived from 'abc_context'");  
  return std::make_shared<contextType>();
}

//-------------------
// Energy grid
//-------------------

class abc_grid : public penred::logs::logger{

public:
  // ****  Energy grid and interpolation constants. The means of "Raw" calificative
  //       is that the variable has not been modifiqued by any transformation. For
  //       example, in a logarithmic scale, if the mean energy is 50eV, EL must be
  //       set to 50eV, not to log(50).
  
  // EMIN  : Minimum grid value
  // EL    : Raw lowest grid value (typically 0.99999*EMIN)
  // EU    : Raw last grid value
  // ET    : Raw Value of each bin (with no transformation)
  // DLEMP : Transformed ET (in a logarithmic scale, the "I" component must be log(ET[i]))
  // DLFC  : Inverse of distance between transofrmed bins
  // DLEMP1: Transformed EL (for example, in a logarithmic scale, log(EL))
  // 
  double EMIN, EL, EU, ET[constants::NEGP], DLEMP[constants::NEGP], DLEMP1, DLFC;

  bool initialized;

  abc_grid() : initialized(false){}
  
  virtual int init(double EMINu, double EMAXu) = 0;
  virtual void getInterval(const double E, int& KE, double& XEL, double& XE, double& XEK) const = 0;
  
  virtual ~abc_grid(){};
};

template<size_t dim>
class abc_genericGrid : public penred::logs::logger{
public:
  static const size_t size = dim;
  double EMIN, EL, EU, ET[size], DLEMP[size], DLEMP1, DLFC;

  bool initialized;

  abc_genericGrid() : initialized(false){}
  
  virtual int init(double EMINu, double EMAXu) = 0;
  virtual void getInterval(const double E, long int& KE,
			   double& XEL, double& XE, double& XEK) const = 0;  
  virtual ~abc_genericGrid(){};
};


//  Implementations
//-------------------

template<class stateType, class contextType, class materialType>
void abc_particle<stateType,contextType,materialType>::
vr_particleStack(const unsigned long long nhist,
		 pen_rand& random, const unsigned verbose){

  if(genericVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    genericVR->run_particleStack(nhist,kpar,KDET,
				 state,genericStates,
				 created,spaceAvailable,random);
      
    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Generic VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }
    stateType defaultState;
    for(unsigned i = 0; i < created; ++i){
      defaultState.copyBase(genericStates[i]);
      stack.store(defaultState);
    }
  }
  if(specificVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    specificVR->run_particleStack(nhist,kpar,KDET,
				  state,specificStates,
				  created,spaceAvailable,random);

    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Specific VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }

    for(unsigned i = 0; i < created; ++i){
      stack.store(specificStates[i]);
    }
  }
}

template<class stateType, class contextType, class materialType>
double abc_particle<stateType,contextType,materialType>::
vr_matChange(const unsigned long long nhist,
	     pen_rand& random, const unsigned verbose){

  double de = 0.0;
  if(genericVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    genericVR->run_matChange(nhist,kpar,MATL,
			     state,genericStates,
			     created,spaceAvailable,random);
      
    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Generic VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }
    stateType defaultState;
    for(unsigned i = 0; i < created; ++i){
      defaultState.copyBase(genericStates[i]);
      stack.store(defaultState);
      de += defaultState.E*defaultState.WGHT;
    }
  }
  if(specificVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    specificVR->run_matChange(nhist,kpar,MATL,
			      state,specificStates,
			      created,spaceAvailable,random);

    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Specific VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }

    for(unsigned i = 0; i < created; ++i){
      stack.store(specificStates[i]);
      de += specificStates[i].E*specificStates[i].WGHT;      
    }
  }
  return de;
}

template<class stateType, class contextType, class materialType>
double abc_particle<stateType,contextType,materialType>::
vr_interfCross(const unsigned long long nhist,
	       pen_rand& random, const unsigned verbose){

  double de = 0.0;
  if(genericVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    genericVR->run_interfCross(nhist,kpar,KDET,
			       state,genericStates,
			       created,spaceAvailable,random);

    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Generic VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }
    stateType defaultState;
    for(unsigned i = 0; i < created; ++i){
      defaultState.copyBase(genericStates[i]);
      stack.store(defaultState);
      de += defaultState.E*defaultState.WGHT;      
    }
  }
  if(specificVR != nullptr){
    unsigned spaceAvailable = constants::NMS-stack.getNSec();
    unsigned created = 0;
    specificVR->run_interfCross(nhist,kpar,KDET,
				state,specificStates,
				created,spaceAvailable,random);

    if(created > spaceAvailable){
      created = spaceAvailable;
      if(verbose > 1)
	printf("abc_particle: Warning: Specific VR creates "
	       "more particles (%u) than available space (%u)\n"
	       "Extra particles will be ignored\n",created,spaceAvailable);
    }

    for(unsigned i = 0; i < created; ++i){
      stack.store(specificStates[i]);
      de += specificStates[i].E*specificStates[i].WGHT;      
    }
  }
  return de;
}



#endif
