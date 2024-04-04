
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

 
#ifndef __PEN_CONTEXT_
#define __PEN_CONTEXT_

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <memory>

#ifdef _PEN_EMBEDDED_DATA_BASE_
#include "materialCreator.hh"
#endif

//Interactions enumerations

enum pen_betaE_interact{
  BETAe_HARD_ELASTIC = 0,
  BETAe_HARD_INELASTIC,
  BETAe_HARD_BREMSSTRAHLUNG,
  BETAe_HARD_INNER_SHELL,
  BETAe_DELTA,
  BETAe_SOFT_INTERACTION,
  BETAe_HARD_TOTAL
};

enum pen_gamma_interact{
  GAMMA_RAYLEIGH = 0,
  GAMMA_COMPTON,
  GAMMA_PHOTOELECTRIC,
  GAMMA_PAIR_PRODUCTION,
  GAMMA_DELTA  
};

enum pen_betaP_interact{
  BETAp_HARD_ELASTIC = 0,
  BETAp_HARD_INELASTIC,
  BETAp_HARD_BREMSSTRAHLUNG,
  BETAp_HARD_INNER_SHELL,
  BETAp_ANNIHILATION,
  BETAp_DELTA,
  BETAp_SOFT_INTERACTION,
  BETAp_HARD_TOTAL
};

//-------------------
// Context
//-------------------

class pen_context;

//Reader for material configuration
class pen_contextReaderMat : public pen_configReader<pen_contextReaderMat>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
    UNKNOWN_PARTICLE = 2,
    INVALID_ATOMIC_NUMBER = 3,
    INVALID_NAME = 4,
  };
  
  int family; //0 -> material-ranges
              //1 -> material-eabs
              //2 -> materials
              //3 -> materials/${subsection}/elements
              //4 -> materials/${subsection}/range
              //5 -> materials/${subsection}/eabs
  
  unsigned actualPartID;
  std::array<double,constants::nParTypes> maxRanges;
  std::array<double,constants::nParTypes> defaultEabs;

  std::string contextlogfile;
  
  struct materialData{
    std::string name;
    unsigned index;
    std::string filename;
    bool forceCreation;
    double density;
    double C1, C2, WCC, WCR;

    std::array<double,constants::nParTypes> maxRanges;
    std::array<double,constants::nParTypes> eabs;
    
    std::vector<penred::penMaterialCreator::massFraction> composition;
    
    inline materialData(const std::string& nameIn) : name(nameIn){
      filename = nameIn + ".mat";
      std::fill(maxRanges.begin(), maxRanges.end(), -1.0);
      std::fill(eabs.begin(), eabs.end(), -1.0e3);
    }
  };

  std::vector<materialData> mats;
  

  pen_contextReaderMat() : family(-1),
			   actualPartID(ALWAYS_AT_END){
    std::fill(maxRanges.begin(), maxRanges.end(), -1.0);
    std::fill(defaultEabs.begin(), defaultEabs.end(), 1.0e3);
  }

  int beginSectionFamily(const std::string& pathInSection,
			 const size_t size,
			 const unsigned verbose);

  int endSectionFamily(const unsigned verbose);
  int beginSection(const std::string& name,
		   const unsigned verbose);
  int endSection(const unsigned verbose);
  
  int storeElement(const std::string& pathInSection,
		   const pen_parserData& element,
		   const unsigned verbose);
  int storeString(const std::string& pathInSection,
		  const std::string& element,
		  const unsigned verbose);

};

//Reader for VR configuration
class pen_contextReaderVR : public pen_configReader<pen_contextReaderVR>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
    UNKNOWN_PARTICLE = 2,
    INVALID_MATERIAL_INDEX = 3,
    INVALID_INTERACTION_INDEX = 4,
  };
  
  int family; //0 -> VR/IForcing
              //1 -> VR/IForcing/${subsection}/bodies
              //2 -> VR/IForcing/${subsection}/materials
              //3 -> VR/bremss
              //4 -> VR/bremss/${subsection}/bodies
              //5 -> VR/bremss/${subsection}/materials
  
  struct IFdata{
    std::string name;
    std::vector<std::string> enabledBodies;
    std::vector<unsigned> enabledMats;
    pen_KPAR particleType;
    unsigned interaction;

    double minW, maxW, factor;

    inline IFdata(const std::string& nameIn) : name(nameIn){}
  };

  std::vector<IFdata> iforcing;
  
  struct bremssData{
    std::string name;
    std::vector<std::string> enabledBodies;
    std::vector<unsigned> enabledMats;

    unsigned factor;

    inline bremssData(const std::string& nameIn) : name(nameIn){}
  };

  std::vector<bremssData> bremss;

  inline pen_contextReaderVR() : family(-1){}
  
  int beginSectionFamily(const std::string& pathInSection,
			 const size_t size,
			 const unsigned verbose);

  int endSectionFamily(const unsigned verbose);
  int beginSection(const std::string& name,
		   const unsigned verbose);
  int endSection(const unsigned verbose);  
  int storeElement(const std::string& pathInSection,
		   const pen_parserData& element,
		   const unsigned verbose);
  int storeString(const std::string& pathInSection,
		  const std::string& element,
		  const unsigned verbose);

};

class pen_context : public abc_context<pen_material>{
  
public:

  enum configErrors{
    SUCCESS = 0,
    NO_MATERIAL_PROVIDED,
    TOO_MUCH_MATERIALS,
    UNABLE_TO_CREATE_MATERIALS,
    MISSING_MATERIAL_COMPOSITION,
    MATERIAL_CREATION_FAILED,
    RANGE_CONTEXT_CREATION_FAILED,
    UNUSED_MATERIAL,
    MAXIMUM_BODY_IF_REACHED,
    UNKNOWN_BODY,
    GEOMETRY_NOT_SET,
    MATERIAL_OUT_OF_RANGE,
    MATERIAL_DEFINED_TWICE,
    GEOMETRY_USES_NON_DEFINED_MATERIAL,
    UNABLE_TO_UPDATE_EABS_FROM_GEOMETRY,
    ERROR_AT_CONTEXT_INIT,
    UNABLE_TO_ACCESS_MATERIAL_FILE,
  };
  
  // ****  Energy grid and interpolation constants.
  pen_logGrid grid;

  // Element data
  std::unique_ptr<pen_elementDataBase> pelements;
  pen_elementDataBase& elements;

  // RITA random sampling variables
  CRNDG3 rndg3;

  //  ****  Variance-reduction parameters.

  //  ----  Parameter values defined for each body. NBV is the maximum
  //        number of bodies. When using PENGEOM, NBV must have the same
  //        value as the parameter NB in module PENGEOM_mod.
  const static unsigned int NBV = pen_geoconst::NB;
  //  ----  Forcing factors, FORCE(IBODY,KPAR,ICOL).
  double FORCE[NBV][constants::nParTypes][constants::MAXINTERACTIONS];
  //  ----  Interaction forcing bodies state 
  bool  LFORCE[NBV][constants::nParTypes];
  // WRANGES store interaction forcing weight ranges (wlow,wup)
  // for each particle in each body
  double WRANGES[NBV][2*constants::nParTypes];
  //  ----  Bremsstrahlung splitting numbers, IBRSPL(IBODY).
  unsigned int IBRSPL[NBV];
  
  inline pen_context() : pelements(new pen_elementDataBase()),
			 elements(*pelements.get()){
    
    //  ****  Variance-reduction parameters.
    
    // Forcing
    for(unsigned int i = 0; i < NBV; i++)
      for(unsigned int j = 0; j < constants::nParTypes; j++){
	WRANGES[i][2*j]   = 0.0;
	WRANGES[i][2*j+1] = 1.0e6;
	LFORCE[i][j] = false;
	for(unsigned int k = 0; k < constants::MAXINTERACTIONS; k++)
	  FORCE[i][j][k] = 1.0;
      }
    
    // Bremsstrahlung splitting
    for(unsigned int i = 0; i < NBV; i++)
      IBRSPL[i] = 1;

  }
  
  int init(const double EMAX, FILE *IWR, int INFO, std::string PMFILE[constants::MAXMAT]);
  int configure(const double EMAX,
		const pen_parserSection& config,
		pen_parserSection& matInfo,
		const unsigned verbose = 2) override;
  int configureWithGeo(const pen_parserSection& config,
		       const unsigned verbose = 2) override;

  inline bool isForcing(const unsigned kpar,
			const unsigned ibody,
			const double wght) const {
    unsigned kpar2 = 2*kpar;
    if(LFORCE[ibody][kpar] &&
       wght >= WRANGES[ibody][kpar2] &&
       wght <= WRANGES[ibody][kpar2+1]){
      return true;
    }
    return false;
  }

  double range(const double E,
	       const pen_KPAR kpar,
	       const unsigned M) const override;

  double avncol(const double E,
		const pen_KPAR kpar,
		const int icol,
		const unsigned imat) const override;

  double IMFP(const double E,
	      const pen_KPAR kpar,
	      const int icol,
	      const unsigned M) const;
  
  double avninter(const double E,
		  const pen_KPAR kpar,
		  const int icol,
		  const unsigned imat,
		  const bool calc_piecewise) const override;

  double getIF(const double forcerIn,
	       const pen_KPAR kpar,
	       const int icol,
	       const unsigned imat,
	       const bool calc_piecewise = true) const;

  inline void setForcing(const double forcerIn,
			 const pen_KPAR kpar,
			 const int icol,
			 const unsigned ibody,
			 const double weightL,
			 const double weightU) override{

    //Read body material
    unsigned imat = readGeometry()->getMat(ibody);

    //Avoid void regions
    if(imat == 0)
      return;

    //Handle negative forcer values
    double forcer = getIF(forcerIn, kpar, icol, imat-1);
    
    //Set interaction forcing for the specified parameters
    FORCE[ibody][kpar][icol] = forcer;

    if(WRANGES[ibody][2*kpar] < weightL)
      WRANGES[ibody][2*kpar] = weightL;
    
    if(WRANGES[ibody][2*kpar+1] > weightU)
      WRANGES[ibody][2*kpar+1] = weightU;

    LFORCE[ibody][kpar] = true;    
  }

  static pen_context* createAuxContext(double EMAX,
				       const char* matFilename,
				       const unsigned verbose);

  //Function to access element DB. Use it only for debug and testing
  inline pen_elementDataBase& getElementDB() { return elements; }

  
};

//Define the configure formats for this context

template<>
struct pen_format<pen_contextReaderMat>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

#Context
context-log/reader-description "Context report file"
context-log/reader-value "context.rep"
context-log/reader-required/type "optional"

#Materials
material-ranges/${subsection}/reader-description "Maximum range, in cm, for the specified particle applied to all materials.\nA negative value disables this value"
material-ranges/${subsection}/reader-value -1.0
material-ranges/${subsection}/reader-required/type "optional"

material-eabs/${subsection}/reader-description "Absorption energy, in eV, applied in all materials for the specified particle"
material-eabs/${subsection}/reader-value 1.0e3
material-eabs/${subsection}/reader-conditions/min/type "greater"
material-eabs/${subsection}/reader-conditions/min/value 50.0
material-eabs/${subsection}/reader-required/type "optional"

materials/${subsection}/reader-description "Sections with materials definitions"
materials/${subsection}/number/reader-value 1 
materials/${subsection}/number/reader-description "Material index in the geometry system. The index 0 is reserved to void regions"
materials/${subsection}/number/reader-conditions/noVoid/type "greater"
materials/${subsection}/number/reader-conditions/noVoid/value 0

materials/${subsection}/filename/reader-description "Material filename to read/write the material information"
materials/${subsection}/filename/reader-value "-"
materials/${subsection}/filename/reader-required/type "optional_if_exist"
materials/${subsection}/filename/reader-required/value "elements"

materials/${subsection}/force-creation/reader-description "If this value is true, and the elements section is provided,\n the material will be created although the material file exists.\nOtherwise, the material is not created unless the provided file does not exist."
materials/${subsection}/force-creation/reader-value false
materials/${subsection}/force-creation/reader-required/type "optional"

materials/${subsection}/elements/${subsection}/reader-description "The element ${subsection} must be the element atomic number (Z), and the corresponding value its fraction by weight in the created material"
materials/${subsection}/elements/${subsection}/reader-value 0.1
materials/${subsection}/elements/${subsection}/reader-required/type "optional"
materials/${subsection}/elements/${subsection}/reader-conditions/greater/type "greater"
materials/${subsection}/elements/${subsection}/reader-conditions/greater/value 0.0

materials/${subsection}/density/reader-description "Material density in g/cm**3. Only required if the material is defined in the configuration"
materials/${subsection}/density/reader-value 1.0
materials/${subsection}/density/reader-required/type "required_if_exist"
materials/${subsection}/density/reader-required/value "elements"
materials/${subsection}/density/reader-conditions/positive/type "greater"
materials/${subsection}/density/reader-conditions/positive/value 0.0

materials/${subsection}/range/${subsection}/reader-description "Range, in cm, for the specified particle in this material. Can be disabled setting a negative value"
materials/${subsection}/range/${subsection}/reader-value -1.0
materials/${subsection}/range/${subsection}/reader-required/type "optional"

materials/${subsection}/eabs/${subsection}/reader-description "Absortpion energy, in eV, for the specified particle in this material.\nIf a negative value is set, it will be replaced by the default material value or specified range."
materials/${subsection}/eabs/${subsection}/reader-value -1.0e3
materials/${subsection}/eabs/${subsection}/reader-required/type "optional"

materials/${subsection}/C1/reader-description "C1 parameter for class II transport"
materials/${subsection}/C1/reader-value 0.05
materials/${subsection}/C1/reader-required/type "optional"
materials/${subsection}/C1/reader-conditions/positive/type "positive"

materials/${subsection}/C2/reader-description "C2 parameter for class II transport"
materials/${subsection}/C2/reader-value 0.05
materials/${subsection}/C2/reader-required/type "optional"
materials/${subsection}/C2/reader-conditions/positive/type "positive"

materials/${subsection}/WCC/reader-description "WCC parameter for class II transport. If a negative value is provided,\n WCC will be set automatically by the program."
materials/${subsection}/WCC/reader-value -1.0
materials/${subsection}/WCC/reader-required/type "optional"

materials/${subsection}/WCR/reader-description "WCR parameter for class II transport. If a negative value is provided,\n WCR will be set automatically by the program."
materials/${subsection}/WCR/reader-value -1.0
materials/${subsection}/WCR/reader-required/type "optional"

)===";  
};

template<>
struct pen_format<pen_contextReaderVR>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Variance reduction

## Interaction forcing

VR/IForcing/${subsection}/reader-description "Interaction forcing configuration"
VR/IForcing/${subsection}/reader-required/type "optional"

VR/IForcing/${subsection}/particle/reader-description "Particle type to apply the interaction forcing"
VR/IForcing/${subsection}/particle/reader-value "electron"

VR/IForcing/${subsection}/interaction/reader-description "Particle specific interaction index to apply the interaction forcing.\nAll the interaction index for each particle can be found in the documentation."
VR/IForcing/${subsection}/interaction/reader-value 1
VR/IForcing/${subsection}/interaction/reader-conditions/positive/type "positive"

VR/IForcing/${subsection}/factor/reader-description "Specify how many times the interaction is more probable.\nA negative value is interpreted as the desired mean number of interactions in a mean free path."
VR/IForcing/${subsection}/factor/reader-value 10.0

VR/IForcing/${subsection}/bodies/${subsection}/reader-description "Enable/disable bodies, which alias are specified by the subsections names, where the interaction forcing is applied."
VR/IForcing/${subsection}/bodies/${subsection}/reader-required/type "optional_if_exist"
VR/IForcing/${subsection}/bodies/${subsection}/reader-required/value "materials"
VR/IForcing/${subsection}/bodies/${subsection}/reader-value true

VR/IForcing/${subsection}/materials/${subsection}/reader-description "Enable/disable materials, which alias are specified by the subsections names, where the interaction forcing is applied."
VR/IForcing/${subsection}/materials/${subsection}/reader-required/type "optional_if_exist"
VR/IForcing/${subsection}/materials/${subsection}/reader-required/value "bodies"
VR/IForcing/${subsection}/materials/${subsection}/reader-value true


VR/IForcing/${subsection}/min-weight/reader-description "The particles with a weight below this threshold will not be forced"
VR/IForcing/${subsection}/min-weight/reader-value 0.0
VR/IForcing/${subsection}/min-weight/reader-required/type "optional"
VR/IForcing/${subsection}/min-weight/reader-conditions/positive/type "positive"

VR/IForcing/${subsection}/max-weight/reader-description "The particles with a weight above this threshold will not be forced"
VR/IForcing/${subsection}/max-weight/reader-value 1.0e6
VR/IForcing/${subsection}/max-weight/reader-required/type "optional"
VR/IForcing/${subsection}/max-weight/reader-conditions/greaterThanMin/type "greater"
VR/IForcing/${subsection}/max-weight/reader-conditions/greaterThanMin/value "min-weight"

## Bremsstrahlung splitting

VR/bremss/${subsection}/reader-description "Bremsstrahlung splitting configuration"
VR/bremss/${subsection}/reader-required/type "optional"

VR/bremss/${subsection}/splitting/reader-description "Specify the final number of duplicated particles"
VR/bremss/${subsection}/splitting/reader-value 1
VR/bremss/${subsection}/splitting/reader-conditions/greaterThan0/type "greater"
VR/bremss/${subsection}/splitting/reader-conditions/greaterThan0/value 0

VR/bremss/${subsection}/bodies/${subsection}/reader-description "Enable/disable bodies, which alias are specified by the subsections names, where the Bremsstrahlung splitting is applied."
VR/bremss/${subsection}/bodies/${subsection}/reader-required/type "optional_if_exist"
VR/bremss/${subsection}/bodies/${subsection}/reader-required/value "materials"
VR/bremss/${subsection}/bodies/${subsection}/reader-value true

VR/bremss/${subsection}/materials/${subsection}/reader-description "Enable/disable materials, which alias are specified by the subsections names, where the Bremsstrahlung splitting is applied."
VR/bremss/${subsection}/materials/${subsection}/reader-required/type "optional_if_exist"
VR/bremss/${subsection}/materials/${subsection}/reader-required/value "bodies"
VR/bremss/${subsection}/materials/${subsection}/reader-value true

)===";  
};

//-----------------------------------------------
// PENELOPE common  functions
//-----------------------------------------------

void DIRECT(const double CDT,
	    const double DF,
	    double &U,
	    double &V,
	    double &W);


void DIRPOL(const double CDT,
	    double &DF,
	    const double CONS,
	    double &SP1,
	    double &SP2,
	    double &SP3,
	    double &U,
	    double &V,
	    double &W,
	    pen_rand& random);

void RELAX(const pen_elementDataBase& elements,
	   const pen_material& mat,
	   pen_particleState& state,
	   const int ICOL,
	   const int MODER,
	   const int IZ,
	   const int IS,
	   int &KS,
	   pen_particleStack<pen_particleState>& stackE,
	   pen_particleStack<pen_state_gPol>& stackG,
	   pen_rand& random);

void PANaR(const pen_particleState& betaState,
	   const double ECUT,
	   pen_particleStack<pen_state_gPol>& stackG,
	   pen_rand& penRand);

void GPHaT(const double E,
	   double &XS,
	   const pen_material& mat,
	   const pen_elementDataBase& elemDB);

#endif
