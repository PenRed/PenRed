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

#ifndef __PEN_X_RAY_ANODE__
#define __PEN_X_RAY_ANODE__

#include "x-ray-common.hh"
#include <thread>
#include <functional>

namespace penred{

  namespace xray{

    void runAnodeSimulation(const unsigned long long nHists,
			    const double Einit,
			    const double beamRad,
			    const pen_context& context,
			    const pen_VRCluster<pen_state_gPol>& photonVR,
			    std::vector<detectedPart>& results,
			    int& seed1, int& seed2,
			    const bool onlyPhotons);

    //Function to simulate a monoenergetic electron beam aiming an anode
    int simAnode(const char* matFilename,
		 const double eEnergy,
		 const double eMin,
		 const double focalSpot,
		 const double angle,
		 const unsigned long long nHists,
		 double& dReg,
		 std::vector<detectedPart>& results,
		 const bool onlyPhotons = true,
		 const unsigned verbose = 1,
		 const unsigned threads2Use = 0);

    

  };

};

#endif
