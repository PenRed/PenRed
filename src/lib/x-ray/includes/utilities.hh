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
//        sanolgi@upvnet.upv.es
//    
//

#ifndef _PEN_X_RAY_UTILITIES_
#define _PEN_X_RAY_UTILITIES_

#include <cstdlib>
#include <algorithm>

#include "x-ray-common.hh"
#include "anode.hh"
#include "collimator.hh"
#include "filter.hh"
#include "phantom.hh"
#include "pen_muen.hh"

namespace penred{

  namespace xray{

    int filterWithCollimation(std::vector<detectedPart>& particlesIn,
			      const double filterZorigin,
			      const double filter2col,
			      const std::vector<std::pair<unsigned, double>>& filters,
			      const double coldz,
			      const double coldx1, const double coldy1,
			      const double coldx2, const double coldy2,
			      const double emin,
			      const unsigned nthreadsIn = 0,
			      const unsigned verbose = 2,
			      const std::string& geoFilename = "");

    int bowtieWithCollimation(std::vector<detectedPart>& particlesIn,
			      const double filterZorigin,
			      const double filter2col,
			      const unsigned matZ,
			      const std::vector<double>& bowtieDs,
			      const double coldz,
			      const double coldx1, const double coldy1,
			      const double coldx2, const double coldy2,
			      const double emin,
			      const unsigned nthreadsIn = 0,
			      const unsigned verbose = 2,
			      const std::string& geoFilename = "");
    
  };
};


#endif
