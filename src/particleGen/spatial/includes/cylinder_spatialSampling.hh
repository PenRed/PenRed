
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es               (Sandra Oliver Gil)
//    
//

 
#ifndef __CYLINDER_SPATIAL_SAMPLING__
#define __CYLINDER_SPATIAL_SAMPLING__

class cylinder_spatialSampling : public abc_spatialSampler {

  DECLARE_SAMPLER(cylinder_spatialSampling)

private:

  double rmin, rmax, dzmin, dzmax;
  double rmin2, dr2, dzmin05, dzmax05, coverHeight;
  double p1, p2, p3;
  
public:

  static const double PI;
  static const double PI2;

  cylinder_spatialSampling() : rmin(0.0),
			       rmax(1.0),
			       dzmin(1.0),
			       dzmax(1.0),
			       rmin2(0.0),
			       dr2(1.0),
			       dzmin05(0.5),
			       dzmax05(0.5)
  {}

  void geoSampling(double pos[3], pen_rand& random) const;

  int configure(const pen_parserSection& config, const unsigned verbose = 0);
  
};

#endif
