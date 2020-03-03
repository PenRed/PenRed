
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

 
#ifndef __CONE_DIRECTION_SAMPLING__
#define __CONE_DIRECTION_SAMPLING__

class cone_directionSampling : public abc_directionSampler{
  DECLARE_SAMPLER(cone_directionSampling)
private:

  static const double deg2rad;
  double CPCT,CPST,SPCT,SPST,SPHI,CPHI,STHE,CTHE,CAPER;
  
public:
  void directionSampling(double dir[3], pen_rand& random) const;
  int configure(const pen_parserSection& config, const unsigned verbose = 0);
  
};

#endif
