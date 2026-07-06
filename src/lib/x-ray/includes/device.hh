//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2025-2026 Vicent Giménez Alventosa
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

#ifndef _PEN_X_RAY_DEVICE_
#define _PEN_X_RAY_DEVICE_

#include <cstdlib>
#include <algorithm>

#include "x-ray-common.hh"
#include "anode.hh"
#include "collimator.hh"
#include "filter.hh"
#include "phantom.hh"
#include "pen_muen.hh"
#include "utilities.hh"

#include <thread>

namespace penred{

  namespace xray{

    class readerXRayDeviceSimulate;

    int constructDevice(std::ostream& out,
                        const double focalSpot,
                        const double source2det,
                        const double source2filter,
                        const double detectorDx,
                        const double detectorDy,
                        const double detectorDz,
                        const double inherentFilterSize,
                        const std::vector<double>& filters,
                        unsigned& initMat,
                        const vector3D<double> detectorPos =
                        vector3D<double>(0.0,0.0,0.0),
                        const bool constructAnode = false,
                        const double anodeAngle = 5.0,
                        const unsigned verbose = 1,
                        const bool PSFFilter = false);


    int constructSimDevice(const pen_parserSection& config,
                           penred::simulation::simulator<pen_context>& simula,
                           const unsigned verbose);

    int simDevice(const pen_parserSection& config,
		  const unsigned verbose);


    int simDevice(const readerXRayDeviceSimulate& reader,
		  const double maxE,
		  const simulation::sampleFuncType<pen_particleState>& fsample,
		  const penred::simulation::simConfig& baseSimConfig,
		  measurements::measurement<double, 2>& detFluence,
		  measurements::measurement<double, 2>& detEdep,
		  measurements::measurement<double, 1>& detSpec,
		  unsigned long long& simHistsOut,
		  const unsigned verbose);
    
    //Reader for device construction configuration
    class readerXRayDeviceCreate : public pen_configReader<readerXRayDeviceCreate>{

    public:

      enum errors{
	SUCCESS = 0,
	UNHANDLED = 1,
	BAD_DIMENSIONS,
      };
      
      int family; // -1 root section, 0 filter section

      bool createAnode;
      double focalSpot;
      double anodeAngle;
      vector3D<double> detectorPosition;
      double source2det;
      double source2filter;
      double detectorDx;
      double detectorDy;
      double detectorDz;
      double inherentFilterWidth;
  
      std::vector<double> filters;      

      readerXRayDeviceCreate() : family(-1){ }

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

      int beginArray(const std::string& pathInSection,
		     const size_t size,
		     const unsigned verbose);
      
      int endArray(const unsigned verbose);

      int storeArrayElement(const std::string& pathInSection,
			    const pen_parserData& element,
			    const size_t pos,
			    const unsigned verbose);
    };

    //Reader for device construction configuration
    class readerXRayDeviceSimulate : public pen_configReader<readerXRayDeviceSimulate>{

    public:

      enum errors{
	SUCCESS = 0,
	UNHANDLED = 1,
	BAD_DIMENSIONS,
	INVALID_ATOMIC_NUMBER,
      };
      
      int family; // -1 root section
                  //  0 filter section
                  //  1 materials section
                  //  2 material elements section

      bool simAnode;
      unsigned long long nHists;
      double maxTime;
      double minEnergy;
      unsigned nThreads;
      unsigned seedPair;
      
      double focalSpot;
      
      vector3D<double> detectorPosition;
      std::string PSFFile;
      vector3D<double> PSFTrans;
      vector3D<double> PSFRotation;      

      double source2det;
      double detectorDx;
      double detectorDy;
      double detectorDz;
      bool detectorIdeal;
      std::string detectorMatFile;

      double inherentFilterWidth;
      std::vector<double> filtersWidth;
      std::vector<unsigned> filtersZ;
      std::vector<std::string> filtersMatFile;
      double source2filter;

      bool storeDetectedPSF;
      bool storeFilteredPSF;
      
      double anodeAngle;
      unsigned anodeZ;
      double kvp;

      unsigned long detBinsX, detBinsY;
      unsigned long eBins;

      std::string outputPrefix;

      struct materialData{
	std::string name;
	int index;
	double density;
    
	std::vector<penred::massFraction> composition;
    
	inline materialData(const std::string& nameIn) : name(nameIn){
	}
      };

      std::vector<materialData> addedGeoMats;
      std::string addedGeoType;
      pen_parserSection addedGeoConf;
      
      readerXRayDeviceSimulate() : family(-1){ }

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

      int beginArray(const std::string& pathInSection,
		     const size_t size,
		     const unsigned verbose);
      
      int endArray(const unsigned verbose);

      int storeArrayElement(const std::string& pathInSection,
			    const pen_parserData& element,
			    const size_t pos,
			    const unsigned verbose);
    };
    

    inline int constructDevice(std::ostream& out,
			       const pen_parserSection& config,
			       const unsigned verbose = 1){

      //Read configuration
      readerXRayDeviceCreate reader;
      int err = reader.read(config, 2);
      if(err != penred::xray::readerXRayDeviceCreate::SUCCESS){
        penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                     "Error: Bad configuration values\n");
        return -2;
      }

      unsigned initMat = 1;
      return constructDevice(out,
                             reader.focalSpot,
                             reader.source2det,
                             reader.source2filter,
                             reader.detectorDx,
                             reader.detectorDy,
                             reader.detectorDz,
                             reader.inherentFilterWidth,
                             reader.filters,
                             initMat,
                             reader.detectorPosition,
                             reader.createAnode,
                             reader.anodeAngle,
                             verbose);
    }

  } // namespace xray
} // namespace penred

template<>
struct pen_format<penred::xray::readerXRayDeviceCreate>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

#Anode configuration

## Anode enabled/disabled
anode/create/reader-description "Enable/disable the anode construction"
anode/create/reader-value true

#Anode angle
anode/angle/reader-description "Sets the anode angle in DEG"
anode/angle/reader-value 5.0
anode/angle/reader-required/type "required_if"
anode/angle/reader-required/value "anode/create"
anode/angle/reader-conditions/gt0/type "greater"
anode/angle/reader-conditions/gt0/value 0.0
anode/angle/reader-conditions/lesserThan90/type "lesser"
anode/angle/reader-conditions/lesserThan90/value 90.0

#Focal spot
focalSpot/reader-description "X-ray focal spot in cm"
focalSpot/reader-value 0.1
focalSpot/reader-conditions/positive/type "positive"

#Distance source to detector
distance/detector/reader-description "Distance, in cm, from source to detector"
distance/detector/reader-value 30.0
distance/detector/reader-conditions/gt0/type "greater"
distance/detector/reader-conditions/gt0/value 0.0

#Detector size
detector/dx/reader-description "Detector size, in cm, in the X axis"
detector/dx/reader-value 50.0
detector/dx/reader-conditions/gt0/type "greater"
detector/dx/reader-conditions/gt0/value 0.0

detector/dy/reader-description "Detector size, in cm, in the Y axis"
detector/dy/reader-value 20.0
detector/dy/reader-conditions/gt0/type "greater"
detector/dy/reader-conditions/gt0/value 0.0

detector/dz/reader-description "Detector depth, in cm, in the Z axis"
detector/dz/reader-value 1.0
detector/dz/reader-required/type "optional"
detector/dz/reader-conditions/gt0/type "greater"
detector/dz/reader-conditions/gt0/value 0.0

# Detector position
detector/pos/x/reader-description "Detector position in X axis (cm)"
detector/pos/x/reader-value 0.0
detector/pos/x/reader-required/type "optional"

detector/pos/y/reader-description "Detector position in Y axis (cm)"
detector/pos/y/reader-value 0.0
detector/pos/y/reader-required/type "optional"

detector/pos/z/reader-description "Detector position in Z axis (cm)"
detector/pos/z/reader-value 0.0
detector/pos/z/reader-required/type "optional"


#Inherent filter
inherent-filter/width/reader-description "Inherent filter size, in cm"
inherent-filter/width/reader-value 0.1

#Distance source to filter
distance/filter/reader-description "Distance, in cm, from source to first filter"
distance/filter/reader-value 8.0
distance/filter/reader-conditions/gt0/type "greater"
distance/filter/reader-conditions/gt0/value 0.0
distance/filter/reader-conditions/lesserThanDet/type "lesser"
distance/filter/reader-conditions/lesserThanDet/value "distance/detector"

## Filters

filters/${subsection}/reader-description "Additional filters to be used"
filters/${subsection}/reader-required/type "optional"

# Filters width
filters/${subsection}/width/reader-description "Filter width in cm"
filters/${subsection}/width/reader-value 0.1
filters/${subsection}/width/reader-conditions/gt0/type "greater"
filters/${subsection}/width/reader-conditions/gt0/value 0.0

)===";
};

template<>
struct pen_format<penred::xray::readerXRayDeviceSimulate>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Simulation generic parameters
simulation/sim-anode/reader-description "Enable/disable anode simulation"
simulation/sim-anode/reader-value true
simulation/sim-anode/reader-required/type "optional"

simulation/histories/reader-description "Maximum number of histories to simulate"
simulation/histories/reader-value 1.0e11
simulation/histories/reader-conditions/gt0/type "greater"
simulation/histories/reader-conditions/gt0/value 0.0
simulation/histories/reader-required/type "optional_if_exist"
simulation/histories/reader-required/value "simulation/max-time"

simulation/max-time/reader-description "Maximum time, in seconds, to perform the simulation"
simulation/max-time/reader-value 1.0e35
simulation/max-time/reader-conditions/gt0/type "greater"
simulation/max-time/reader-conditions/gt0/value 0.0
simulation/max-time/reader-required/type "optional_if_exist"
simulation/max-time/reader-required/value "simulation/histories"

simulation/min-energy/reader-description "Minimum energy to be simulated and tallied, in eV"
simulation/min-energy/reader-value 1.0e3
simulation/min-energy/reader-conditions/gt/type "greater"
simulation/min-energy/reader-conditions/gt/value 50.0
simulation/min-energy/reader-required/type "optional"

simulation/nthreads/reader-description "Number of threads to be used. Set it to 0 to get all available threads."
simulation/nthreads/reader-value 0
simulation/nthreads/reader-conditions/positive/type "positive"
simulation/nthreads/reader-required/type "optional"

simulation/seedPair/reader-description "Initial seed pair to be used"
simulation/seedPair/reader-value 0
simulation/seedPair/reader-conditions/positive/type "positive"
simulation/seedPair/reader-conditions/lesser/type "lesser"
simulation/seedPair/reader-conditions/lesser/value 1001
simulation/seedPair/reader-required/type "optional"

simulation/detBins/nx/reader-description "Number of spatial detector bins in X axis"
simulation/detBins/nx/reader-value 100
simulation/detBins/nx/reader-conditions/enought/type "greater"
simulation/detBins/nx/reader-conditions/enought/value 0
simulation/detBins/nx/reader-required/type "optional"

simulation/detBins/ny/reader-description "Number of spatial detector bins in Y axis"
simulation/detBins/ny/reader-value 100
simulation/detBins/ny/reader-conditions/enought/type "greater"
simulation/detBins/ny/reader-conditions/enought/value 0
simulation/detBins/ny/reader-required/type "optional"

simulation/eBins/reader-description "Number of energetic bins to tally spectrums"
simulation/eBins/reader-value 100
simulation/eBins/reader-conditions/enought/type "greater"
simulation/eBins/reader-conditions/enought/value 3
simulation/eBins/reader-required/type "optional"

simulation/output-prefix/reader-description "Output prefix path"
simulation/output-prefix/reader-value ""
simulation/output-prefix/reader-required/type "optional"

#X-ray device common characteristics
x-ray/focal-spot/reader-description "X-ray focal spot in cm"
x-ray/focal-spot/reader-value 0.1
x-ray/focal-spot/reader-conditions/positive/type "positive"

#Detector position
x-ray/detector/position/reader-description "Position (x,y,z) of the detector's top face center, in cm."
x-ray/detector/position/reader-value [0.0,0.0,0.0]
x-ray/detector/position/reader-required/type "optional"

#Source PSF
x-ray/source/psf/path/reader-description "Path to the source phase space file"
x-ray/source/psf/path/reader-value "path/to/data.psf"
x-ray/source/psf/path/reader-required/type "optional_if"
x-ray/source/psf/path/reader-required/value "simulation/sim-anode"

x-ray/source/psf/translation/reader-description "Translation applied to each particle stored in the PSF."
x-ray/source/psf/translation/reader-value [0.0,0.0,0.0]
x-ray/source/psf/translation/reader-required/type "optional"

x-ray/source/psf/rotation/reader-description "Euler angles to rotate the PSF's particles. Rotation is performed around the Z,Y,Z axis, in that order. Notice the detector position is below the source position (-Z direction)."
x-ray/source/psf/rotation/reader-value [0.0,0.0,0.0]
x-ray/source/psf/rotation/reader-required/type "optional"

#Distance source to detector
x-ray/distance/detector/reader-description "Distance, in cm, from anode impact point to the detector."
x-ray/distance/detector/reader-value 30.0
x-ray/distance/detector/reader-conditions/gt0/type "greater"
x-ray/distance/detector/reader-conditions/gt0/value 0.0

#Detector size
x-ray/detector/dx/reader-description "Detector size, in cm, in the X axis"
x-ray/detector/dx/reader-value 50.0
x-ray/detector/dx/reader-conditions/gt0/type "greater"
x-ray/detector/dx/reader-conditions/gt0/value 0.0

x-ray/detector/dy/reader-description "Detector size, in cm, in the Y axis"
x-ray/detector/dy/reader-value 20.0
x-ray/detector/dy/reader-conditions/gt0/type "greater"
x-ray/detector/dy/reader-conditions/gt0/value 0.0

x-ray/detector/dz/reader-description "Detector depth, in cm, in the Z axis"
x-ray/detector/dz/reader-value 1.0
x-ray/detector/dz/reader-required/type "optional"
x-ray/detector/dz/reader-conditions/gt0/type "greater"
x-ray/detector/dz/reader-conditions/gt0/value 0.0

# Detector absorption
x-ray/detector/ideal/reader-description "Enable/disable ideal detection, i.e. perfect absorber"
x-ray/detector/ideal/reader-value true
x-ray/detector/ideal/reader-required/type "optional"

x-ray/detector/material/reader-description "Sets the detector material. Only for non ideal detectors"
x-ray/detector/material/reader-value "detector.mat"
x-ray/detector/material/reader-required/type "optional_if"
x-ray/detector/material/reader-required/value "x-ray/detector/ideal"

#Inherent filter
x-ray/inherent-filter/width/reader-description "Inherent filter width, in cm"
x-ray/inherent-filter/width/reader-value 0.1

#Distance source to filter
x-ray/distance/filter/reader-description "Distance, in cm, from anode impact point to first filter."
x-ray/distance/filter/reader-value 8.0
x-ray/distance/filter/reader-conditions/gt0/type "greater"
x-ray/distance/filter/reader-conditions/gt0/value 0.0
x-ray/distance/filter/reader-conditions/lesserThanDet/type "lesser"
x-ray/distance/filter/reader-conditions/lesserThanDet/value "x-ray/distance/detector"

## Filters

x-ray/filters/${subsection}/reader-description "Additional filters to be used"
x-ray/filters/${subsection}/reader-required/type "optional"

# width
x-ray/filters/${subsection}/width/reader-description "Filter width in cm"
x-ray/filters/${subsection}/width/reader-value 0.1
x-ray/filters/${subsection}/width/reader-conditions/gt0/type "greater"
x-ray/filters/${subsection}/width/reader-conditions/gt0/value 0.0

# material
x-ray/filters/${subsection}/z/reader-description "Filter atomic number Z"
x-ray/filters/${subsection}/z/reader-value 13
x-ray/filters/${subsection}/z/reader-conditions/gt0/type "greater"
x-ray/filters/${subsection}/z/reader-conditions/gt0/value 0
x-ray/filters/${subsection}/z/reader-required/type "optional_if_exist"
x-ray/filters/${subsection}/z/reader-required/value "mat-file"

x-ray/filters/${subsection}/mat-file/reader-description "Filter material file path"
x-ray/filters/${subsection}/mat-file/reader-value "-"
x-ray/filters/${subsection}/mat-file/reader-required/type "optional_if_exist"
x-ray/filters/${subsection}/mat-file/reader-required/value "z"

#Anode configuration

## Anode

#Angle
x-ray/anode/angle/reader-description "Sets the anode angle in DEG"
x-ray/anode/angle/reader-value 5.0
x-ray/anode/angle/reader-required/type "required"
x-ray/anode/angle/reader-conditions/gt0/type "greater"
x-ray/anode/angle/reader-conditions/gt0/value 0.0
x-ray/anode/angle/reader-conditions/lesserThan90/type "lesser"
x-ray/anode/angle/reader-conditions/lesserThan90/value 90.0

#Z
x-ray/anode/z/reader-description "Sets the anode atomic number"
x-ray/anode/z/reader-value 74
x-ray/anode/z/reader-required/type "optional"
x-ray/anode/z/reader-conditions/gt0/type "greater"
x-ray/anode/z/reader-conditions/gt0/value 0

#Beam energy
x-ray/kvp/reader-description "X-ray KVP value"
x-ray/kvp/reader-value 120
x-ray/kvp/reader-required/type "required_if"
x-ray/kvp/reader-required/value "simulation/sim-anode"
x-ray/kvp/reader-conditions/gt1/type "greater"
x-ray/kvp/reader-conditions/gt1/value 1.0
x-ray/kvp/reader-conditions/lesserThan1MeV/type "lesser"
x-ray/kvp/reader-conditions/lesserThan1MeV/value 1000.0

## Added geometry

# Geometry configuration section
geometry/config/${subsection}/reader-description "Added geometry configuration section"
geometry/config/${subsection}/reader-required/type "required_if_exist"
geometry/config/${subsection}/reader-required/value "geometry/materials"

# Materials
geometry/materials/${subsection}/reader-description "Materials of the added geometry"
geometry/materials/${subsection}/reader-required/type "required_if_exist"
geometry/materials/${subsection}/reader-required/value "geometry/config"

geometry/materials/${subsection}/number/reader-description "Material number in the added geometry"
geometry/materials/${subsection}/number/reader-value 1
geometry/materials/${subsection}/number/reader-conditions/gt0/type "greater"
geometry/materials/${subsection}/number/reader-conditions/gt0/value 0


geometry/materials/${subsection}/density/reader-description "Material density in g/cm**3"
geometry/materials/${subsection}/density/reader-value 1.0
geometry/materials/${subsection}/density/reader-conditions/gt0/type "greater"
geometry/materials/${subsection}/density/reader-conditions/gt0/value 0.0

geometry/materials/${subsection}/elements/reader-description "The value of ${subsection} is expected to be the element atomic number (Z), and the corresponding value its fraction by weight in the created material"
geometry/materials/${subsection}/elements/${subsection}/reader-value 1.0
geometry/materials/${subsection}/elements/${subsection}/reader-conditions/gt0/type "greater"
geometry/materials/${subsection}/elements/${subsection}/reader-conditions/gt0/value 0.0

## Added tallies

# Geometry configuration section
tallies/${subsection}/reader-description "User defined tallies. Each tally name will be prefixed with 'added_'"
tallies/${subsection}/reader-required/type "optional"


## Create PSFs
psf/detected/reader-description "If enabled, a PSF is created at the detector"
psf/detected/reader-value false
psf/detected/reader-required/type "optional"

psf/filtered/reader-description "If enabled, a PSF is created after the last filter"
psf/filtered/reader-value false
psf/filtered/reader-required/type "optional"

)===";
};


#endif
