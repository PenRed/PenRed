
//
//
//    Copyright (C) 2020 Universitat de València - UV
//    Copyright (C) 2020 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

 
#ifndef __IMAGE_SPATIAL_SAMPLING__
#define __IMAGE_SPATIAL_SAMPLING__

#ifdef _PEN_USE_DICOM_

#include "pen_dicom.hh"

class image_spatialSampling : public abc_spatialSampler {

  DECLARE_SAMPLER(image_spatialSampling)

private:
  long int nx,ny,nz;
  long int nxy;
  long int nvox;
  double dx, dy, dz;
  double imageCx,imageCy,imageCz;
  double Ox,Oy,Oz;
  double* F;
  long int* K;
  double isocenter[3];
public:

  image_spatialSampling() : nx(0), ny(0), nz(0),
			    dx(0.0), dy(0.0), dz(0.0),
			    F(nullptr), K(nullptr)
  {}

  void geoSampling(double pos[3], pen_rand& random) const;

  int configure(const pen_parserSection& config, const unsigned verbose = 0);

  ~image_spatialSampling();
};

#endif
#endif
