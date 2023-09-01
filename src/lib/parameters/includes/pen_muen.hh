 
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

#ifndef __PEN_MUEN__
#define __PEN_MUEN__

#include <thread>
#include <atomic>
#include "PenRed.hh"
#include "pen_geometries.hh"

namespace pen_muen{

  int calculate(const double Emin,
		const double Emax,
		const unsigned nBins,
		const double tolerance,
		const double simTime,
		const char* matFilename,
		std::vector<double>& EData,
		std::vector<double>& muenData);

  double simulate(const pen_context& context,
		  const double E0,
		  const double simTime,
		  const double tolerance,
		  int& seed1, int& seed2,
		  const unsigned ithread);
  
};

#endif
