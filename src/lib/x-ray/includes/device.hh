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
			const double source2bowtie,
			const double detectorDx,
			const double detectorDy,
			const double inherentFilterSize,
			const std::vector<double>& filters,
			std::vector<double> bowtieDz,
			const vector3D<double> center =
			vector3D<double>(0.0,0.0,0.0),
			const bool constructAnode = false,
			const double anodeAngle = 5.0,
			const unsigned verbose = 1);


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
      vector3D<double> sourcePosition;
      double source2det;
      double source2filter;
      double detectorDx;
      double detectorDy;
      double inherentFilterWidth;
  
      std::vector<double> filters;

      double source2bowtie;
      std::vector<double> bowtieDz;
      

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
      
      vector3D<double> sourcePosition;
      std::string spatialDistribFile;
      std::string energyDistribFile;
      double distrib2source;

      double source2det;
      double detectorDx;
      double detectorDy;

      double inherentFilterWidth;
      std::vector<double> filtersWidth;
      std::vector<unsigned> filtersZ;
      std::vector<std::string> filtersMatFile;
      double source2filter;

      double source2bowtie;
      std::vector<double> bowtieDz;
      unsigned bowtieZ;
      std::string bowtieMatFile;
      bool bowtieAutoDesign;
      unsigned bowtieDesignBins;
      
      double anodeAngle;
      unsigned anodeZ;
      double kvp;
  
      std::vector<double> filters;

      bool printGeo;

      unsigned detBins;
      unsigned eBins;

      double tolerance;

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
	printf("Error: Bad configuration values\n");
	return -2;
      }

      return constructDevice(out,
			     reader.focalSpot,
			     reader.source2det,
			     reader.source2filter,
			     reader.source2bowtie,
			     reader.detectorDx,
			     reader.detectorDy,
			     reader.inherentFilterWidth,
			     reader.filters,
			     reader.bowtieDz,
			     reader.sourcePosition,
			     reader.createAnode,
			     reader.anodeAngle,
			     verbose);
    }

    inline int checkSimDevice(const pen_parserSection& config,
			      unsigned& nMats,
			      const unsigned verbose){

      // ** Parse configuration
      
      //Read information from config section
      readerXRayDeviceSimulate reader;
      int err = reader.read(config,verbose);
      if(err != readerXRayDeviceSimulate::SUCCESS){
	return err;
      }

      nMats = 2; //Collimators and detector
      if(reader.simAnode)
	nMats += 1; //Anode material

      //Add inherent filter
      if(reader.inherentFilterWidth > 0.0)
	nMats += 1;

      //Add filters materials
      nMats += reader.filtersZ.size();

      //Add bowtie material
      if(reader.source2bowtie > 0.0 && reader.bowtieDz.size() > 0)
	nMats += 1;

      return err;
    }
  };
};

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
source2det/reader-description "Distance, in cm, from source to detector"
source2det/reader-value 30.0
source2det/reader-conditions/gt0/type "greater"
source2det/reader-conditions/gt0/value 0.0

#Detector size
detector/dx/reader-description "Detector size, in cm, in the X axis"
detector/dx/reader-value 50.0
detector/dx/reader-conditions/gt0/type "greater"
detector/dx/reader-conditions/gt0/value 0.0

detector/dy/reader-description "Detector size, in cm, in the Y axis"
detector/dy/reader-value 20.0
detector/dy/reader-conditions/gt0/type "greater"
detector/dy/reader-conditions/gt0/value 0.0

#Inherent filter
inherentFilter/width/reader-description "Inherent filter size, in cm"
inherentFilter/width/reader-value 0.1

#Distance source to filter
source2filter/reader-description "Distance, in cm, from source to first filter"
source2filter/reader-value 8.0
source2filter/reader-conditions/gt0/type "greater"
source2filter/reader-conditions/gt0/value 0.0
source2filter/reader-conditions/lesserThanDet/type "lesser"
source2filter/reader-conditions/lesserThanDet/value "source2det"

## Filters

filters/${subsection}/reader-description "Additional filters to be used"
filters/${subsection}/reader-required/type "optional"

# Filters width
filters/${subsection}/width/reader-description "Filter width in cm"
filters/${subsection}/width/reader-value 0.1
filters/${subsection}/width/reader-conditions/gt0/type "greater"
filters/${subsection}/width/reader-conditions/gt0/value 0.0

## Bowtie

# Distance source to bowtie
source2bowtie/reader-description "Distance, in cm, from source to bowtie"
source2bowtie/reader-value -1.0
source2bowtie/reader-conditions/greaterThanFilters/type "greater"
source2bowtie/reader-conditions/greaterThanFilters/value "source2filter"
source2bowtie/reader-required/type "required_if_exist"
source2bowtie/reader-required/value "bowtie/dz"

bowtie/dz/reader-description "Bowtie heights"
bowtie/dz/reader-value [0.5,0.5,0.4,0.3,0.2,0.3,0.4,0.5,0.5]
bowtie/dz/reader-required/type "required_if_exist"
bowtie/dz/reader-required/value "source2bowtie"

# Source position
source/pos/x/reader-description "Source position in X axis (cm)"
source/pos/x/reader-value 0.0
source/pos/x/reader-required/type "optional"

# Source position
source/pos/y/reader-description "Source position in Y axis (cm)"
source/pos/y/reader-value 0.0
source/pos/y/reader-required/type "optional"

# Source position
source/pos/z/reader-description "Source position in Z axis (cm)"
source/pos/z/reader-value 0.0
source/pos/z/reader-required/type "optional"

)===";
};

template<>
struct pen_format<penred::xray::readerXRayDeviceSimulate>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Simulation generic parameters
simulation/sim-anode/reader-description "Enable/disable anode simulation"
simulation/sim-anode/reader-value false
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

simulation/tolerance/reader-description "Tolerance to finish the simulation"
simulation/tolerance/reader-value 0.0
simulation/tolerance/reader-conditions/positive/type "positive"
simulation/tolerance/reader-required/type "optional"

simulation/print-geometry/reader-description "Enable/disable geometry print"
simulation/print-geometry/reader-value false
simulation/print-geometry/reader-required/type "optional"

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

simulation/detBins/reader-description "Number of spatial detector bins"
simulation/detBins/reader-value 100
simulation/detBins/reader-conditions/enought/type "greater"
simulation/detBins/reader-conditions/enought/value 3
simulation/detBins/reader-required/type "optional"

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

#Source position
x-ray/source/position/reader-description "Position (x,y,z) of the source, in cm. If the anode is simulated, this position corresponds to the center of the anode face where the electron beam collides.\n If a spatial distribution is used instead, the position specify where the center of the distribution is located."
x-ray/source/position/reader-value [0.0,0.0,0.0]

#Source distributions
x-ray/source/spatial-distrib-file/reader-description "Path to the spatial source distribution file"
x-ray/source/spatial-distrib-file/reader-value "path/to/distrib.dat"
x-ray/source/spatial-distrib-file/reader-required/type "optional_if"
x-ray/source/spatial-distrib-file/reader-required/value "simulation/sim-anode"

x-ray/source/energy-distrib-file/reader-description "Path to the energy source distribution file"
x-ray/source/energy-distrib-file/reader-value "path/to/distrib.dat"
x-ray/source/energy-distrib-file/reader-required/type "optional_if"
x-ray/source/energy-distrib-file/reader-required/value "simulation/sim-anode"

x-ray/source/distrib2source/reader-description "Distance between the source spatial distribution and the 'real' source in cm."
x-ray/source/distrib2source/reader-value 0.4
x-ray/source/distrib2source/reader-required/type "optional_if"
x-ray/source/distrib2source/reader-required/value "simulation/sim-anode"
x-ray/source/distrib2source/reader-conditions/positive/type "positive"

#Distance source to detector
x-ray/source2det/reader-description "Distance, in cm, from source to detector"
x-ray/source2det/reader-value 30.0
x-ray/source2det/reader-conditions/gt0/type "greater"
x-ray/source2det/reader-conditions/gt0/value 0.0

#Detector size
x-ray/detector/dx/reader-description "Detector size, in cm, in the X axis"
x-ray/detector/dx/reader-value 50.0
x-ray/detector/dx/reader-conditions/gt0/type "greater"
x-ray/detector/dx/reader-conditions/gt0/value 0.0

x-ray/detector/dy/reader-description "Detector size, in cm, in the Y axis"
x-ray/detector/dy/reader-value 20.0
x-ray/detector/dy/reader-conditions/gt0/type "greater"
x-ray/detector/dy/reader-conditions/gt0/value 0.0

#Inherent filter
x-ray/inherentFilter/width/reader-description "Inherent filter width, in cm"
x-ray/inherentFilter/width/reader-value 0.1

#Distance source to filter
x-ray/source2filter/reader-description "Distance, in cm, from source to first filter"
x-ray/source2filter/reader-value 8.0
x-ray/source2filter/reader-conditions/gt0/type "greater"
x-ray/source2filter/reader-conditions/gt0/value 0.0
x-ray/source2filter/reader-conditions/lesserThanDet/type "lesser"
x-ray/source2filter/reader-conditions/lesserThanDet/value "x-ray/source2det"

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

## Bowtie

# Distance source to bowtie
x-ray/source2bowtie/reader-description "Distance, in cm, from source to bowtie"
x-ray/source2bowtie/reader-value -1.0
x-ray/source2bowtie/reader-conditions/greaterThanFilters/type "greater"
x-ray/source2bowtie/reader-conditions/greaterThanFilters/value "x-ray/source2filter"
x-ray/source2bowtie/reader-conditions/lesserThanDet/type "lesser"
x-ray/source2bowtie/reader-conditions/lesserThanDet/value "x-ray/source2det"
x-ray/source2bowtie/reader-required/type "required_if_exist"
x-ray/source2bowtie/reader-required/value "x-ray/bowtie/dz"

x-ray/bowtie/dz/reader-description "Bowtie top face displacements"
x-ray/bowtie/dz/reader-value [0.5,0.5,0.4,0.3,0.2,0.3,0.4,0.5,0.5]
x-ray/bowtie/dz/reader-required/type "required_if_exist"
x-ray/bowtie/dz/reader-required/value "x-ray/source2bowtie"

x-ray/bowtie/z/reader-description "Bowtie material atomic number Z"
x-ray/bowtie/z/reader-value 13
x-ray/bowtie/z/reader-conditions/gt0/type "greater"
x-ray/bowtie/z/reader-conditions/gt0/value 0
x-ray/bowtie/z/reader-required/type "optional_if_exist"
x-ray/bowtie/z/reader-required/value "x-ray/bowtie/mat-file"

x-ray/bowtie/mat-file/reader-description "Bowtie filter material file path"
x-ray/bowtie/mat-file/reader-value "-"
x-ray/bowtie/mat-file/reader-required/type "optional_if_exist"
x-ray/bowtie/mat-file/reader-required/value "x-ray/bowtie/z"

x-ray/bowtie/auto-design/reader-description "Eanble/disable bowtie automatic design"
x-ray/bowtie/auto-design/reader-value false
x-ray/bowtie/auto-design/reader-required/type "optional"

x-ray/bowtie/design-bins/reader-description "Number of bins to design bowtie"
x-ray/bowtie/design-bins/reader-value 0
x-ray/bowtie/design-bins/reader-conditions/gt0/type "greater"
x-ray/bowtie/design-bins/reader-conditions/gt0/value 0
x-ray/bowtie/design-bins/reader-required/type "optional"

#Anode configuration

## Anode

#Angle
anode/angle/reader-description "Sets the anode angle in DEG"
anode/angle/reader-value 5.0
anode/angle/reader-required/type "required"
anode/angle/reader-conditions/gt0/type "greater"
anode/angle/reader-conditions/gt0/value 0.0
anode/angle/reader-conditions/lesserThan90/type "lesser"
anode/angle/reader-conditions/lesserThan90/value 90.0

#Z
anode/z/reader-description "Sets the anode atomic number"
anode/z/reader-value 74
anode/z/reader-required/type "optional"
anode/z/reader-conditions/gt0/type "greater"
anode/z/reader-conditions/gt0/value 0

#Beam energy
x-ray/kvp/reader-description "X-ray KVP value"
x-ray/kvp/reader-value 120.0e3
x-ray/kvp/reader-required/type "required_if"
x-ray/kvp/reader-required/value "simulation/sim-anode"
x-ray/kvp/reader-conditions/gt0/type "greater"
x-ray/kvp/reader-conditions/gt0/value 1.0
x-ray/kvp/reader-conditions/lesserThan90/type "lesser"
x-ray/kvp/reader-conditions/lesserThan90/value 1.0e6

## Added geometry

# Type
geometry/type/reader-description "Added geometry type"
geometry/type/reader-value "-"
geometry/type/reader-required/type "required_if_exist"
geometry/type/reader-required/value "geometry/config"

# Geometry configuration section
geometry/config/${subsection}/reader-description "Added geometry configuration section"
geometry/config/${subsection}/reader-required/type "required_if_exist"
geometry/config/${subsection}/reader-required/value "geometry/type"

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

)===";
};


#endif
