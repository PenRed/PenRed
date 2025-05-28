
//
//
//    Copyright (C) 2021-2023 Universitat de València - UV
//    Copyright (C) 2021-2023 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#ifndef __BRACHY_SPECIFIC_SAMPLER__
#define __BRACHY_SPECIFIC_SAMPLER__

#ifdef _PEN_USE_DICOM_

#include "PSFsource.hh"
#include "pen_geometries.hh"

class brachy_specificSampler : public abc_specificSampler<pen_particleState>{
  DECLARE_SPECIFIC_SAMPLER(brachy_specificSampler, pen_particleState)
  
  psf_specificSampler psf;
  
  std::vector<vector3D<double>> positions;
  std::vector<std::array<double,9>> rotations;
  std::vector<double> weights;
  
  //Parameter to select kind of seed to be read from DICOM file 
  int seedID;

  //Enable/disable seed rotation
  bool seedRotation;
  
  public:

  brachy_specificSampler() : abc_specificSampler<pen_particleState>(USE_NONE),
			     seedID(-1),
			     seedRotation(true){}


  void skip(const unsigned long long dhists);

  int configure(double& Emax,
		const pen_parserSection& config,
		const unsigned nthreads,
		const unsigned verbose);  

  void sample(pen_particleState& state,
	      pen_KPAR& genKpar,
	      unsigned long long& dhist,
	      pen_rand& random);

  void updateGeometry(const wrapper_geometry* geometryIn);

  int sharedConfig(const brachy_specificSampler& o);
};

#endif

#endif
