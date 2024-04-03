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

#ifndef _PEN_X_RAY_UTILITIES_
#define _PEN_X_RAY_UTILITIES_

#include "x-ray-common.hh"
#include "anode.hh"
#include "collimator.hh"
#include "filter.hh"

namespace penred{

  namespace xray{

    int filterWithCollimation(std::vector<detectedPart>& particlesIn,
			      const double filterZorigin,
			      const double filter2col,
			      const std::vector<std::pair<unsigned, double>>& filters,
			      const double coldz,
			      const double coldx1, const double coldy1,
			      const double coldx2, const double coldy2,
			      const double emin,
			      const unsigned nthreadsIn = 0,
			      const unsigned verbose = 2,
			      const std::string& geoFilename = "");

    //Reader for material configuration
    class readerHVL : public pen_configReader<readerHVL>{

    public:

      enum errors{
	SUCCESS = 0,
	UNHANDLED = 1,
      };
      
      int family; // -1 root section

      bool simAnode;
      unsigned long long nhists;
      double beamEnergy;
      double minE;
      double focalSpot;
      double anodeAngle;
      double sourcePositionZ;
      double source2det;
      double source2filter;
      double source2register;
      double detectorRad;

      unsigned anodeZ;
      double anodeDensity;
      std::string anodeMatFilename;
      bool printGeometries;

      std::vector<std::pair<unsigned, double>> filters;

      unsigned nThreads;

      readerHVL() : family(-1){
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

    int HVL(const pen_parserSection config,
      const std::vector<detectedPart>& particlesIn,
      double& hvl,
      const unsigned verbose);
    
  };
};

template<>
struct pen_format<penred::xray::readerHVL>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

#Anode configuration

## force anode simulation
anode/simulation/enabled/reader-description "Enables/disables splicit anode simulation. If disabled, particles must be provided as argument."
anode/simulation/enabled/reader-value false
anode/simulation/enabled/reader-required/type "optional"

## Anode simulation histories
anode/simulation/nhists/reader-description "Number of histories to be simulated"
anode/simulation/nhists/reader-value 1.0e6
anode/simulation/nhists/reader-required/type "required_if"
anode/simulation/nhists/reader-required/value "anode/simulation/enabled"
anode/simulation/nhists/reader-conditions/1orMore/type "greater"
anode/simulation/nhists/reader-conditions/1orMore/value 1

## Anode initial electron energy
anode/simulation/energy/reader-description "Energy, in eV, of the initial electron beam"
anode/simulation/energy/reader-value 100.0e3
anode/simulation/energy/reader-required/type "required_if"
anode/simulation/energy/reader-required/value "anode/simulation/enabled"
anode/simulation/energy/reader-conditions/not2low/type "greater"
anode/simulation/energy/reader-conditions/not2low/value 1e3

## Anode minimum tallied energy
anode/simulation/focalSpot/reader-description "Beam focal spot, in cm"
anode/simulation/focalSpot/reader-value 0.1
anode/simulation/focalSpot/reader-required/type "required_if"
anode/simulation/focalSpot/reader-required/value "anode/simulation/enabled"
anode/simulation/focalSpot/reader-conditions/greaterThan0/type "greater"
anode/simulation/focalSpot/reader-conditions/greaterThan0/value 0.0

## Anode angle
anode/simulation/angle/reader-description "Anode angle, in DEG"
anode/simulation/angle/reader-value 5.0
anode/simulation/angle/reader-required/type "required_if"
anode/simulation/angle/reader-required/value "anode/simulation/enabled"
anode/simulation/angle/reader-conditions/greaterThan0/type "greater"
anode/simulation/angle/reader-conditions/greaterThan0/value 0.0
anode/simulation/angle/reader-conditions/lesserThan90/type "lesser"
anode/simulation/angle/reader-conditions/lesserThan90/value 90.0

## Anode materials

anode/material/z/reader-description "Anode material atomic number. Set it to 0 to disable its construction"
anode/material/z/reader-value 0
anode/material/z/reader-required/type "optional_if_exist"
anode/material/z/reader-required/value "anode/material/filename"
anode/material/z/reader-conditions/gt0/type "greater"
anode/material/z/reader-conditions/gt0/value 0

anode/material/density/reader-description "Anode material density. Only used if a file is not provided."
anode/material/density/reader-value -1.0
anode/material/density/reader-required/type "optional"
anode/material/density/reader-required/value "anode/material/filename"
anode/material/density/reader-conditions/gt0/type "greater"
anode/material/density/reader-conditions/gt0/value 0.0

anode/material/filename/reader-description "Anode material filename"
anode/material/filename/reader-value "-"
anode/material/filename/reader-required/type "optional_if_exist"
anode/material/filename/reader-required/value "anode/material/z"

## Filters

filters/${subsection}/reader-description "Filters to be used"

# width
filters/${subsection}/width/reader-description "Filter width in cm"
filters/${subsection}/width/reader-value 0.1
filters/${subsection}/width/reader-conditions/gt0/type "greater"
filters/${subsection}/width/reader-conditions/gt0/value 0.0

# material
filters/${subsection}/z/reader-description "Filter atomic number Z"
filters/${subsection}/z/reader-value 13
filters/${subsection}/z/reader-conditions/gt0/type "greater"
filters/${subsection}/z/reader-conditions/gt0/value 0

# Source variables

# Z position
source/zpos/reader-description "Z position of the source. This value is ignored if the anode is simulated"
source/zpos/reader-value 0.0
source/zpos/reader-required/type "optional_if"
source/zpos/reader-required/value "anode/simulation/enabled"

# Distance variables

## Distance source to detector
distance/source-detector/reader-description "Distance source to detector, in cm"
distance/source-detector/reader-value 1.0
distance/source-detector/reader-conditions/non0/type "greater"
distance/source-detector/reader-conditions/non0/value 0.0

## Distance from source to registered anode particles
distance/source-register/reader-description "Distance from source to registered particles, in cm"
distance/source-register/reader-value 1.0
distance/source-register/reader-required/type "optional_if"
distance/source-register/reader-required/value "anode/simulation/enabled"
distance/source-register/reader-conditions/non0/type "greater"
distance/source-register/reader-conditions/non0/value 0.0

## Distance source to filter
distance/source-filter/reader-description "Distance source to filter, in cm"
distance/source-filter/reader-required/type "optional"
distance/source-filter/reader-value 1.0
distance/source-filter/reader-conditions/non0/type "greater"
distance/source-filter/reader-conditions/non0/value 0.0

# Detector variables
detector/radius/reader-description "Distance source to filter, in cm"
detector/radius/reader-required/type "optional"
detector/radius/reader-value 1.0
detector/radius/reader-conditions/non0/type "greater"
detector/radius/reader-conditions/non0/value 0.0

# Minimum tallied energy
minE/reader-description "Minimum particle energy to be tallied"
minE/reader-value 20.0e3
minE/reader-required/type "optional"
minE/reader-conditions/lesserThanE/type "lesser"
minE/reader-conditions/lesserThanE/value "anode/simulation/energy"

# Number of threads
nThreads/reader-description "Number of threads to be used"
nThreads/reader-value 0
nThreads/reader-required/type "optional"
nThreads/reader-conditions/positive/type "positive"

# Print geometries
print-geometries/reader-description "Enable/disable geometries print"
print-geometries/reader-value false
print-geometries/reader-required/type "optional"

)===";
};


#endif
