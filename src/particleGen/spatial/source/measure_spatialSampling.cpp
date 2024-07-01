
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


#include "measure_spatialSampling.hh"

// ** 3D SAMPLING

// * Reader functions

int samplerReader_spatial_measure3D::storeElement(const std::string& pathInSection,
						  const pen_parserData& element,
						  const unsigned){

  if(pathInSection.compare("dx") == 0){
    dx = element;
  }
  else if(pathInSection.compare("dy") == 0){
    dy = element;
  }
  else if(pathInSection.compare("dz") == 0){
    dz = element;
  }
  else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
}

int samplerReader_spatial_measure3D::storeString(const std::string& pathInSection,
						 const std::string& element,
						 const unsigned){

  if(pathInSection.compare("filename") == 0){
    filename = element;
  }
  else{
    return errors::UNHANDLED;
  }

  return errors::SUCCESS;
  
}

void measure3D_spatialSampling::geoSampling(double pos[3], pen_rand& random) const{

  std::array<double, 3> sampled = sampler.samplePositions(random);
  
  pos[0] = sampled[0];
  pos[1] = sampled[1];
  pos[2] = sampled[2];
  
}

int measure3D_spatialSampling::configure(const pen_parserSection& config, const unsigned verbose){


  //Read configuration
  samplerReader_spatial_measure3D reader;
  int err = reader.read(config,verbose);
  if(err != samplerReader_spatial_measure3D::SUCCESS){
    return err;
  }
  
  //Read measure file
  penred::measurements::results<double, 3> distrib;
  std::ifstream fin(reader.filename, std::ifstream::in);
  if(!fin){
    if(verbose > 0){
      printf("measure3D_spatialSampling: Error: Unable to "
	     "open distribution file '%s'",
	     reader.filename.c_str());
    }
    return -1;
  }
  
  err = distrib.read(fin);
  if(err != 0){
    if(verbose > 0){
      printf("measure3D_spatialSampling: Error: Unable to "
	     "read distribution file '%s'.\n  Error code: %d\n",
	     reader.filename.c_str(), err);
    }
    return -2;
  }

  //Init sampling
  err = sampler.init(distrib.readData(),
		     distrib.readDimBins(),
		     distrib.readLimits());

  if(err != 0){
    if(verbose > 0){
      printf("measure3D_spatialSampling: Error: Unable to "
	     "init spatial distribution with walker's aliasing.\n  Error code: %d\n",
	     err);
    }
    return -3;
  }

  //Save translation
  translation[0] = reader.dx;
  translation[1] = reader.dy;
  translation[2] = reader.dz;

  return 0;
}

REGISTER_SAMPLER(measure3D_spatialSampling,3D_MEASURE)


// ** 2D SAMPLING

// * Reader functions

int samplerReader_spatial_measure2D::storeElement(const std::string& pathInSection,
						  const pen_parserData& element,
						  const unsigned){
  
  if(pathInSection.compare("dx") == 0){
    dx = element;
  }
  else if(pathInSection.compare("dy") == 0){
    dy = element;
  }
  else if(pathInSection.compare("dz") == 0){
    dz = element;
  }
  else if(pathInSection.compare("constant-coordinate") == 0){
    constantCoord = element;
  }
  else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
}

int samplerReader_spatial_measure2D::storeString(const std::string& pathInSection,
						 const std::string& element,
						 const unsigned verbose){

  if(pathInSection.compare("filename") == 0){
    filename = element;
  }
  else if(pathInSection.compare("plane") == 0){
    std::string plane;
    plane = element;
    std::transform(plane.cbegin(), plane.cend(), plane.begin(),
		   [](unsigned char c){ return std::tolower(c); });

    if(plane.compare("xy") == 0){
      xindex = 0;
      yindex = 1;
      zindex = 2;
    }
    else if(plane.compare("yx") == 0){
      xindex = 1;
      yindex = 0;
      zindex = 2;
    }
    else if(plane.compare("xz") == 0){
      xindex = 0;
      yindex = 2;
      zindex = 1;
    }
    else if(plane.compare("zx") == 0){
      xindex = 1;
      yindex = 2;
      zindex = 0;
    }
    else if(plane.compare("yz") == 0){
      xindex = 2;
      yindex = 0;
      zindex = 1;
    }    
    else if(plane.compare("zy") == 0){
      xindex = 2;
      yindex = 1;
      zindex = 0;
    }
    else{
      if(verbose > 0){
	printf("Error: Invalid 'plane' value (%s)\n"
	       "       Valid values are: xy,yx,xz,zx,yz,zy\n",
	       plane.c_str());
      }
      return errors::INVALID_VALUE;
    }
  }
  else{
    return errors::UNHANDLED;
  }

  return errors::SUCCESS;
  
}

void measure2D_spatialSampling::geoSampling(double pos[3], pen_rand& random) const{

  std::array<double, 2> sampled = sampler.samplePositions(random);
  std::array<double, 3> aux = {sampled[0], sampled[1], constantCoord};
  
  pos[0] = aux[xindex];
  pos[1] = aux[yindex];
  pos[2] = aux[zindex];
  
}

int measure2D_spatialSampling::configure(const pen_parserSection& config, const unsigned verbose){


  //Read configuration
  samplerReader_spatial_measure2D reader;
  int err = reader.read(config,verbose);
  if(err != samplerReader_spatial_measure2D::SUCCESS){
    return err;
  }
  
  //Read measure file
  penred::measurements::results<double, 2> distrib;
  std::ifstream fin(reader.filename, std::ifstream::in);
  if(!fin){
    if(verbose > 0){
      printf("measure2D_spatialSampling: Error: Unable to "
	     "open distribution file '%s'",
	     reader.filename.c_str());
    }
    return -1;
  }
  
  err = distrib.read(fin);
  if(err != 0){
    if(verbose > 0){
      printf("measure2D_spatialSampling: Error: Unable to "
	     "read distribution file '%s'.\n  Error code: %d\n",
	     reader.filename.c_str(), err);
    }
    return -2;
  }

  //Init sampling
  err = sampler.init(distrib.readData(),
		     distrib.readDimBins(),
		     distrib.readLimits());

  if(err != 0){
    if(verbose > 0){
      printf("measure2D_spatialSampling: Error: Unable to "
	     "init spatial distribution with walker's aliasing.\n  Error code: %d\n",
	     err);
    }
    return -3;
  }

  //Save coordinates indexes and constant value
  constantCoord = reader.constantCoord;
  xindex = reader.xindex;
  yindex = reader.yindex;
  zindex = reader.zindex;

  //Save translation
  translation[0] = reader.dx;
  translation[1] = reader.dy;
  translation[2] = reader.dz;

  return 0;
}

REGISTER_SAMPLER(measure2D_spatialSampling,2D_MEASURE)


// ** 1D SAMPLING

// * Reader functions

int samplerReader_spatial_measure1D::storeElement(const std::string& pathInSection,
						  const pen_parserData& element,
						  const unsigned){
  
  if(pathInSection.compare("dx") == 0){
    dx = element;
  }
  else if(pathInSection.compare("dy") == 0){
    dy = element;
  }
  else if(pathInSection.compare("dz") == 0){
    dz = element;
  }
  else if(pathInSection.compare("x") == 0){
    constX = element;
  }
  else if(pathInSection.compare("y") == 0){
    constY = element;
  }
  else if(pathInSection.compare("z") == 0){
    constZ = element;
  }  
  else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
}

int samplerReader_spatial_measure1D::storeString(const std::string& pathInSection,
						 const std::string& element,
						 const unsigned verbose){

  if(pathInSection.compare("filename") == 0){
    filename = element;
  }
  else if(pathInSection.compare("axis") == 0){
    std::string aux;
    aux = element;
    std::transform(aux.cbegin(), aux.cend(), aux.begin(),
		   [](unsigned char c){ return std::tolower(c); });

    if(aux.compare("x") == 0){
      sampleIndex = 0;
    }
    else if(aux.compare("y") == 0){
      sampleIndex = 1;
    }
    else if(aux.compare("z") == 0){
      sampleIndex = 2;
    }
    else{
      if(verbose > 0){
	printf("Error: Invalid 'axis' value (%s)\n"
	       "       Valid values are: x,y,z\n",
	       aux.c_str());
      }
      return errors::INVALID_VALUE;
    }
  }
  else{
    return errors::UNHANDLED;
  }

  return errors::SUCCESS;
  
}

void measure1D_spatialSampling::geoSampling(double pos[3], pen_rand& random) const{

  std::array<double, 1> sampled = sampler.samplePositions(random);
  
  pos[0] = constX;
  pos[1] = constY;
  pos[2] = constZ;

  pos[sampleIndex] = sampled[0];
  
}

int measure1D_spatialSampling::configure(const pen_parserSection& config, const unsigned verbose){


  //Read configuration
  samplerReader_spatial_measure1D reader;
  int err = reader.read(config,verbose);
  if(err != samplerReader_spatial_measure1D::SUCCESS){
    return err;
  }
  
  //Read measure file
  penred::measurements::results<double, 1> distrib;
  std::ifstream fin(reader.filename, std::ifstream::in);
  if(!fin){
    if(verbose > 0){
      printf("measure1D_spatialSampling: Error: Unable to "
	     "open distribution file '%s'",
	     reader.filename.c_str());
    }
    return -1;
  }
  
  err = distrib.read(fin);
  if(err != 0){
    if(verbose > 0){
      printf("measure1D_spatialSampling: Error: Unable to "
	     "read distribution file '%s'.\n  Error code: %d\n",
	     reader.filename.c_str(), err);
    }
    return -2;
  }

  //Init sampling
  err = sampler.init(distrib.readData(),
		     distrib.readDimBins(),
		     distrib.readLimits());

  if(err != 0){
    if(verbose > 0){
      printf("measure1D_spatialSampling: Error: Unable to "
	     "init spatial distribution with walker's aliasing.\n  Error code: %d\n",
	     err);
    }
    return -3;
  }

  sampleIndex = reader.sampleIndex;
  
  constX = reader.constX;
  constY = reader.constY;
  constZ = reader.constZ;

  //Save translation
  translation[0] = reader.dx;
  translation[1] = reader.dy;
  translation[2] = reader.dz;

  return 0;
}

REGISTER_SAMPLER(measure1D_spatialSampling,1D_MEASURE)
