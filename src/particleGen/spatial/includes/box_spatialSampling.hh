
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
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
//        vicente.gimenez@uv.es
//    
//

 
#ifndef __BOX_SPATIAL_SAMPLING__
#define __BOX_SPATIAL_SAMPLING__

class box_spatialSampling : public abc_spatialSampler {

  DECLARE_SAMPLER(box_spatialSampling)

private:

  double dx, dy, dz;
  double dx05, dy05, dz05;
  
public:

  box_spatialSampling() : dx(1.0),
			  dy(1.0),
			  dz(1.0),
			  dx05(0.5),
			  dy05(0.5),
			  dz05(0.5)
  {}

  void geoSampling(double pos[3], pen_rand& random) const;

  int configure(const pen_parserSection& config, const unsigned verbose = 0);
  
};

#endif
