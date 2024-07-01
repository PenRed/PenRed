
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

 
#ifndef __PEN_SPATIAL_DETECTION_DISTRIB_TALLY__
#define __PEN_SPATIAL_DETECTION_DISTRIB_TALLY__

#include "pen_constants.hh"

class pen_DetectionSpatialDistrib: public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_DetectionSpatialDistrib,pen_particleState)
    
  private:

  penred::measurements::measurement<double, 4> results;
  
  unsigned idet;
  unsigned iPar;

  bool printBins, printCoord, reduceDimensions;
    
public:
    
  pen_DetectionSpatialDistrib() : pen_genericTally(USE_INTERFCROSS |
						   USE_MOVE2GEO)    
  {}
    
  void tally_interfCross(const unsigned long long nhist,
			 const unsigned kdet,
			 const pen_KPAR kpar,
			 const  pen_particleState& state);

  void tally_move2geo(const unsigned long long /*nhist*/,
		      const unsigned /*kdet*/,
		      const pen_KPAR /*kpar*/,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);  
    
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  inline void flush(){ results.flush(); }
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_DetectionSpatialDistrib& tally);
};

//Tally configuration reader
class tallyReader_DetectionSpatialDistrib : public pen_configReader<tallyReader_DetectionSpatialDistrib>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
    UNKNOWN_PARTICLE = 2,
  };

  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;

  unsigned long nx, ny, nz;

  double emin, emax;
  unsigned long nEBins;

  unsigned kdet;

  unsigned ipar;

  bool printBins, printCoord;

  int storeElement(const std::string& pathInSection,
		   const pen_parserData& element,
		   const unsigned verbose);
  
  int storeString(const std::string& pathInSection,
		  const std::string& element,
		  const unsigned verbose);  
};

  
template<>
struct pen_format<tallyReader_DetectionSpatialDistrib>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Tally "DetectionSpatialDistrib" reader configuration

reader-description "Tally to register the spatial distribution of the number of particles reaching a detector"

# Spatial dimensions

## X
spatial/xmin/reader-description "Minimum X value to tally, in cm"
spatial/xmin/reader-value -1.0e35
spatial/xmin/reader-required/type "required_if_exist"
spatial/xmin/reader-required/value "spatial/nx"

spatial/xmax/reader-description "Maximum X value to tally, in cm"
spatial/xmax/reader-value 1.0e35
spatial/xmax/reader-conditions/gt/type "greater"
spatial/xmax/reader-conditions/gt/value "spatial/xmin"
spatial/xmax/reader-required/type "required_if_exist"
spatial/xmax/reader-required/value "spatial/nx"

spatial/nx/reader-description "Number of bins in X axis"
spatial/nx/reader-value 1
spatial/nx/reader-conditions/gt0/type "greater"
spatial/nx/reader-conditions/gt0/value 0
spatial/nx/reader-required/type "optional"

## Y
spatial/ymin/reader-description "Minimum Y value to tally, in cm"
spatial/ymin/reader-value -1.0e35
spatial/ymin/reader-required/type "required_if_exist"
spatial/ymin/reader-required/value "spatial/ny"

spatial/ymax/reader-description "Maximum Y value to tally, in cm"
spatial/ymax/reader-value 1.0e35
spatial/ymax/reader-conditions/gt/type "greater"
spatial/ymax/reader-conditions/gt/value "spatial/ymin"
spatial/ymax/reader-required/type "required_if_exist"
spatial/ymax/reader-required/value "spatial/ny"

spatial/ny/reader-description "Number of bins in Y axis"
spatial/ny/reader-value 1
spatial/ny/reader-conditions/gt0/type "greater"
spatial/ny/reader-conditions/gt0/value 0
spatial/ny/reader-required/type "optional"

## Z
spatial/zmin/reader-description "Minimum Z value to tally, in cm"
spatial/zmin/reader-value -1.0e35
spatial/zmin/reader-required/type "required_if_exist"
spatial/zmin/reader-required/value "spatial/nz"

spatial/zmax/reader-description "Maximum Z value to tally, in cm"
spatial/zmax/reader-value 1.0e35
spatial/zmax/reader-conditions/gt/type "greater"
spatial/zmax/reader-conditions/gt/value "spatial/zmin"
spatial/zmax/reader-required/type "required_if_exist"
spatial/zmax/reader-required/value "spatial/nz"

spatial/nz/reader-description "Number of bins in Z axis"
spatial/nz/reader-value 1
spatial/nz/reader-conditions/gt0/type "greater"
spatial/nz/reader-conditions/gt0/value 0
spatial/nz/reader-required/type "optional"

## Detector

detector/reader-description "Detector index to tally at"
detector/reader-value 1
detector/reader-conditions/gt0/type "greater"
detector/reader-conditions/gt0/value 0

## Energy

energy/nbins/reader-description "Number of energy bins to tally spectrums"
energy/nbins/reader-value 1
energy/nbins/reader-required/type "optional"

energy/emin/reader-description "Minimum energy to be tallied"
energy/emin/reader-value 0.0
energy/emin/reader-required/type "required_if_exist"
energy/emin/reader-required/value "energy/nbins"

energy/emax/reader-description "Maximum energy to be tallied"
energy/emax/reader-value 1.0e30
energy/emax/reader-required/type "required_if_exist"
energy/emax/reader-required/value "energy/nbins"
energy/emax/reader-conditions/gt/type "greater"
energy/emax/reader-conditions/gt/value "energy/emin"

## Particle to register
particle/reader-description "Particle to be registered"
particle/reader-value "gamma"

## Print options
printBins/reader-description "Enable/disable printing bin numbers in results report"
printBins/reader-value false
printBins/reader-required/type "optional"

printCoord/reader-description "Enable/disable printing bin coordinates in results report"
printCoord/reader-value true
printCoord/reader-required/type "optional"

)===";
};


#endif
