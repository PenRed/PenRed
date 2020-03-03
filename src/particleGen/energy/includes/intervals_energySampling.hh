
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


#ifndef __INTERVALS_ENERGY_SAMPLING__
#define __INTERVALS_ENERGY_SAMPLING__

#include "pen_auxiliar.hh"

class intervals_energySampling : public abc_energySampler{
  DECLARE_SAMPLER(intervals_energySampling)
private:

  static const unsigned maxIntervals = 30;
  double energies[maxIntervals];
  double dE[maxIntervals];
  double cummulative[maxIntervals+1];
  unsigned nIntervals;
  
  
public:
  void energySampling(double& energy, pen_rand& random) const;
  int configure(double& Emax, const pen_parserSection& config, const unsigned verbose);
  
};

#endif
