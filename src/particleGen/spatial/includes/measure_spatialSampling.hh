
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

 
#ifndef __MEASURE_SPATIAL_SAMPLING__
#define __MEASURE_SPATIAL_SAMPLING__

#include <fstream>

// ** 3D SAMPLING

class measure3D_spatialSampling : public abc_spatialSampler {

  DECLARE_SAMPLER(measure3D_spatialSampling)

private:

  penred::sampling::aliasing<3> sampler;  
  
public:

  void geoSampling(double pos[3], pen_rand& random) const;

  int configure(const pen_parserSection& config, const unsigned verbose = 0);
  
};

//Tally configuration reader
class samplerReader_spatial_measure3D : public pen_configReader<samplerReader_spatial_measure3D>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
  };

  std::string filename;

  double dx,dy,dz;

  int storeElement(const std::string& pathInSection,
		   const pen_parserData& element,
		   const unsigned);  
  int storeString(const std::string& pathInSection,
		  const std::string& element,
		  const unsigned verbose);  
};

  
template<>
struct pen_format<samplerReader_spatial_measure3D>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Sampler "3D_MEASURE" reader configuration

reader-description "Spatial sampler to generate particle positions according to a 3D measruement. For example, generated by the 'DETECTION_SPATIAL_DISTRIB' tally"

filename/reader-description "File with the spatial distribution"
filename/reader-value "distrib.dat"
filename/reader-required/type "required"

dx/reader-description "Displacement in the X axis"
dx/reader-value 0.0
dx/reader-required/type "optional"

dy/reader-description "Displacement in the Y axis"
dy/reader-value 0.0
dy/reader-required/type "optional"

dz/reader-description "Displacement in the Z axis"
dz/reader-value 0.0
dz/reader-required/type "optional"

)===";
};

// ** 2D SAMPLING

class measure2D_spatialSampling : public abc_spatialSampler {

  DECLARE_SAMPLER(measure2D_spatialSampling)

private:

  penred::sampling::aliasing<2> sampler;

  unsigned xindex, yindex, zindex;

  double constantCoord;

public:

  void geoSampling(double pos[3], pen_rand& random) const;

  int configure(const pen_parserSection& config, const unsigned verbose = 0);
  
};

//Tally configuration reader
class samplerReader_spatial_measure2D : public pen_configReader<samplerReader_spatial_measure2D>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
    INVALID_VALUE = 2,
  };

  std::string filename;

  double dx,dy,dz;

  unsigned xindex, yindex, zindex;

  double constantCoord;

  int storeElement(const std::string& pathInSection,
		   const pen_parserData& element,
		   const unsigned);  
  int storeString(const std::string& pathInSection,
		  const std::string& element,
		  const unsigned verbose);  
};

  
template<>
struct pen_format<samplerReader_spatial_measure2D>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Sampler "2D_MEASURE" reader configuration

reader-description "Spatial sampler to generate particle positions according to a 2D measruement. For example, generated by the 'DETECTION_SPATIAL_DISTRIB' tally"

filename/reader-description "File with the spatial distribution"
filename/reader-value "distrib.dat"
filename/reader-required/type "required"

dx/reader-description "Displacement in the X axis"
dx/reader-value 0.0
dx/reader-required/type "optional"

dy/reader-description "Displacement in the Y axis"
dy/reader-value 0.0
dy/reader-required/type "optional"

dz/reader-description "Displacement in the Z axis"
dz/reader-value 0.0
dz/reader-required/type "optional"

plane/reader-description "Plane to be sampled (XY,XZ,YZ)"
plane/reader-value "xy"
plane/reader-required/type "required"

constant-coordinate/reader-description "Sets the constant coordinate value"
constant-coordinate/value 0.0

)===";
};

// ** 1D SAMPLING

class measure1D_spatialSampling : public abc_spatialSampler {

  DECLARE_SAMPLER(measure1D_spatialSampling)

private:

  penred::sampling::aliasing<1> sampler;

  unsigned sampleIndex;

  double constX;
  double constY;
  double constZ;

public:

  void geoSampling(double pos[3], pen_rand& random) const;

  int configure(const pen_parserSection& config, const unsigned verbose = 0);
  
};

//Tally configuration reader
class samplerReader_spatial_measure1D : public pen_configReader<samplerReader_spatial_measure1D>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
    INVALID_VALUE = 2,
  };

  std::string filename;

  double dx,dy,dz;

  unsigned sampleIndex;

  double constX;
  double constY;
  double constZ;

  int storeElement(const std::string& pathInSection,
		   const pen_parserData& element,
		   const unsigned);  

  int storeString(const std::string& pathInSection,
		  const std::string& element,
		  const unsigned verbose);  
};

  
template<>
struct pen_format<samplerReader_spatial_measure1D>{
  //By default, no format is defined
  static constexpr const char* format = R"===(

# Sampler "1D_MEASURE" reader configuration

reader-description "Spatial sampler to generate particle positions according to a 1D measruement. For example, generated by the 'DETECTION_SPATIAL_DISTRIB' tally"

filename/reader-description "File with the spatial distribution"
filename/reader-value "distrib.dat"
filename/reader-required/type "required"

dx/reader-description "Displacement in the X axis"
dx/reader-value 0.0
dx/reader-required/type "optional"

dy/reader-description "Displacement in the Y axis"
dy/reader-value 0.0
dy/reader-required/type "optional"

dz/reader-description "Displacement in the Z axis"
dz/reader-value 0.0
dz/reader-required/type "optional"

axis/reader-description "Coordinate to be sampled (x,y,z)"
axis/reader-value "x"
axis/reader-required/type "required"

x/reader-description "Sets the constant value for x coordinate"
x/reader-value 0.0
x/reader-required/type "optional"

y/reader-description "Sets the constant value for y coordinate"
y/reader-value 0.0
y/reader-required/type "optional"

z/reader-description "Sets the constant value for z coordinate"
z/reader-value 0.0
z/reader-required/type "optional"

)===";
};


#endif
