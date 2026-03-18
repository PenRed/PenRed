
//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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

 
#ifndef __PEN_DETECTION_EDEP_TALLY__
#define __PEN_DETECTION_EDEP_TALLY__

//#include "pen_constants.hh"

class pen_DetectionEDep: public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_DetectionEDep,pen_particleState,DETECTION_EDEP,
		std::pair<double, penred::tally::Dim<4>>
		)
    
  private:

  penred::measurements::measurement<double, 4> measure;

  unsigned iBodyBind;
  unsigned idet;
  unsigned iPar;
  bool detect;

  bool printBins, printCoord, reduceDimensions;

  const wrapper_geometry* geo;
  
public:
    
  pen_DetectionEDep() : pen_genericTally(USE_STEP |
					 USE_LOCALEDEP |
					 USE_INTERFCROSS |
					 USE_MOVE2GEO |
					 USE_BEGINPART |
					 USE_SAMPLEDPART),
			geo(nullptr)
  {
    //Register results functions
    setResultsGenerator<0>
      ([this](const unsigned long long nhists) -> penred::measurements::results<double, 4>{
	penred::measurements::results<double, 4> generated;
	measure.results(nhists, generated);
	return generated;
      });
  }
    
  void tally_step(const unsigned long long nhist,
		  const pen_KPAR /*kpar*/,
		  const pen_particleState& state,
		  const tally_StepData& stepData);

  void tally_localEdep(const unsigned long long nhist,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state,
		       const double dE);

  void tally_interfCross(const unsigned long long /*nhist*/,
			 const unsigned kdet,
			 const pen_KPAR kpar,
			 const pen_particleState& /*state*/);

  void tally_move2geo(const unsigned long long nhist,
		      const unsigned kdet,
		      const pen_KPAR kpar,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);

  void tally_beginPart(const unsigned long long nhist,
		       const unsigned kdet,
		       const pen_KPAR kpar,
		       const pen_particleState& state);

  void tally_sampledPart(const unsigned long long nhist,
			 const unsigned long long /*dhist*/,
			 const unsigned kdet,
			 const pen_KPAR kpar,
			 const pen_particleState& state);  
    
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  inline void flush(){ measure.flush(); }
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_DetectionEDep& tally);

  inline void transformAndAdd(vector3D<double>& pos, double t, double val, unsigned long long nhist){
    if(geo != nullptr){
      geo->composeInvTransform(iBodyBind, pos, t);
    }
    measure.add({pos.x, pos.y, pos.z, t}, val, nhist);
  }
};

//Tally configuration reader
class tallyReader_DetectionEDep : public pen_configReader<tallyReader_DetectionEDep>{

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

  double tmin, tmax;
  unsigned long nt;

  unsigned kdet;
  std::string bodyBind;

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
struct pen_format<tallyReader_DetectionEDep>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Tally "DetectionEDep" reader configuration

reader-description "Tally to register the energy deposition distribution in the selected detector"

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

## Time

time/nbins/reader-description "Number of temporal bins"
time/nbins/reader-value 1
time/nbins/reader-required/type "optional"

time/min/reader-description "Minimum energy to be tallied"
time/min/reader-value 0.0
time/min/reader-required/type "required_if_exist"
time/min/reader-required/value "time/nbins"

time/max/reader-description "Maximum energy to be tallied"
time/max/reader-value 1.0e30
time/max/reader-required/type "required_if_exist"
time/max/reader-required/value "time/nbins"
time/max/reader-conditions/gt/type "greater"
time/max/reader-conditions/gt/value "time/min"

## Particle to register
particle/reader-description "Particle to be registered. Set it to 'all' to register energy depositions from all particle types"
particle/reader-value "all"
particle/reader-required/type "optional"

## Print options
printBins/reader-description "Enable/disable printing bin numbers in results report"
printBins/reader-value false
printBins/reader-required/type "optional"

printCoord/reader-description "Enable/disable printing bin coordinates in results report"
printCoord/reader-value true
printCoord/reader-required/type "optional"

## Object binding
binding/reader-description "Binds the tally to the specified object alias, applying its transforms"
binding/reader-value ""
binding/reader-required/type "optional"

)===";
};


#endif
