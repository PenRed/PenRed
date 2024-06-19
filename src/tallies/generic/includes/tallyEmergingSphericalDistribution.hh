
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


#ifndef __PEN_EMERGIN_SPHERICAL_DISTRIB_TALLY__
#define __PEN_EMERGIN_SPHERICAL_DISTRIB_TALLY__

#include "pen_constants.hh" 

class pen_EmergingSphericalDistrib : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_EmergingSphericalDistrib,pen_particleState)
  
private:

  //Results bined by energy, theta and phi
  std::array<penred::measurements::measurement<double, 3>, constants::nParTypes> results;

  std::array<bool, constants::nParTypes> enabled;

  bool printBins, printCoord, reduceDimensions;    
    
public:
    
  pen_EmergingSphericalDistrib() : pen_genericTally( USE_MOVE2GEO |
						     USE_MATCHANGE |
						     USE_ENDSIM){
    std::fill(enabled.begin(), enabled.end(), false);    
  }
  
  void tally_move2geo(const unsigned long long nhist,
		      const unsigned /*kdet*/,
		      const pen_KPAR kpar,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);
  
  void tally_matChange(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const unsigned /*prevMat*/);
  
  void tally_endSim(const unsigned long long /*nhist*/);
  
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],     
		const pen_parserSection& config, const unsigned verbose);
  void flush();
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_EmergingSphericalDistrib& tally);
};


//Tally configuration reader
class tallyReader_EmergingSphericalDistrib : public pen_configReader<tallyReader_EmergingSphericalDistrib>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
    UNKNOWN_PARTICLE = 2,
  };

  int family;

  double emin, emax;
  double tmin, tmax;
  double pmin, pmax;

  unsigned long ne, nt, np;

  unsigned actualPartID;
  std::array<bool, constants::nParTypes> enabledPart;

  bool printBins, printCoord;

  tallyReader_EmergingSphericalDistrib() : family(-1){
    std::fill(enabledPart.begin(), enabledPart.end(), false);
  }
    
  int storeElement(const std::string& pathInSection,
		   const pen_parserData& element,
		   const unsigned verbose);
  
  int beginSectionFamily(const std::string& pathInSection,
			 const size_t /*size*/,
			 const unsigned /*verbose*/);

  int endSectionFamily(const unsigned /*verbose*/);

  int beginSection(const std::string& name,
		   const unsigned verbose);

  int endSection(const unsigned);
};

  
template<>
struct pen_format<tallyReader_EmergingSphericalDistrib>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Tally "EmergingSphericalDistrib" reader configuration

reader-description "Tally to register the spherical distribution of emerging particles, i.e. escpaing from the geometry system"

# Spatial dimensions

## Energy

energy/nbins/reader-description "Number of energy bins to tally spectrums"
energy/nbins/reader-value 1
energy/nbins/reader-required/type "optional"

energy/min/reader-description "Minimum energy to be tallied"
energy/min/reader-value 0.0
energy/min/reader-required/type "required_if_exist"
energy/min/reader-required/value "energy/nbins"

energy/max/reader-description "Maximum energy to be tallied"
energy/max/reader-value 1.0e30
energy/max/reader-required/type "required_if_exist"
energy/max/reader-required/value "energy/nbins"
energy/max/reader-conditions/gt/type "greater"
energy/max/reader-conditions/gt/value "energy/min"

## Theta

theta/nbins/reader-description "Number of polar bins to tally"
theta/nbins/reader-value 1
theta/nbins/reader-required/type "optional"

theta/min/reader-description "Minimum polar angle to be tallied, in DEG"
theta/min/reader-value 0.0
theta/min/reader-required/type "optional"
theta/min/reader-conditions/gt/type "greater_equal"
theta/min/reader-conditions/gt/value 0.0

theta/max/reader-description "Maximum azimuthal angle to be tallied, in DEG"
theta/max/reader-value 180.0
theta/max/reader-required/type "optional"
theta/max/reader-conditions/gt/type "greater"
theta/max/reader-conditions/gt/value "theta/min"
theta/max/reader-conditions/lt/type "lesser_equal"
theta/max/reader-conditions/lt/value 180.0

## Phi

phi/nbins/reader-description "Number of azimuthal bins to tally"
phi/nbins/reader-value 1
phi/nbins/reader-required/type "optional"

phi/min/reader-description "Minimum azimuthal angle to be tallied, in DEG"
phi/min/reader-value 0.0
phi/min/reader-required/type "optional"
phi/min/reader-conditions/gt/type "greater_equal"
phi/min/reader-conditions/gt/value 0.0

phi/max/reader-description "Maximum azimuthal angle to be tallied, in DEG"
phi/max/reader-value 360.0
phi/max/reader-required/type "optional"
phi/max/reader-conditions/gt/type "greater"
phi/max/reader-conditions/gt/value "phi/min"
phi/max/reader-conditions/lt/type "lesser_equal"
phi/max/reader-conditions/lt/value 360.0

## Particle to register
particle/${subsection}/reader-description "Particle to enable/disable register"
particle/${subsection}/reader-value true

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
