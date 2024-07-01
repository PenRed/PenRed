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

#include "pen_randoms.hh"

int main(){

  //Create distribution
  constexpr unsigned nx = 53;
  constexpr unsigned ny = 87;
  constexpr unsigned nz = 23;
  constexpr unsigned nE = 100;

  constexpr unsigned nxy = nx*ny;
  constexpr unsigned nxyz = nxy*nz;  
  
  std::vector<double> dis(nx*ny*nz*nE);

  for(size_t e = 0; e < nE; ++e){
    for(size_t k = 0; k < nz; ++k){
      for(size_t j = 0; j < ny; ++j){
	for(size_t i = 0; i < nx; ++i){
	  dis[e*nxyz + k*nxy + j*nx + i] =
	    (2.0+sin(static_cast<double>(j)+0.5)*
	     cos(static_cast<double>(i)+0.5))/((e+k)/(k+1));
	}
      }
    }
  }

  //Init sampler
  penred::sampling::aliasing<4> sampler; 
  sampler.init(dis, {nx,ny,nz,nE},
	       {std::pair<double, double>(-102.54,10.06),
		std::pair<double, double>(302.5,543.87),
		std::pair<double, double>(-23,64.87),
		std::pair<double, double>(-65.2,45.76)});

  //Create random number generator
  pen_rand random;
  
  constexpr unsigned long nIter = 1000000000;
  for(size_t i = 0; i < nIter; ++i){

    //Get seeds
    int seed1, seed2;
    random.getSeeds(seed1,seed2);

    //Sample global index
    unsigned long globIndex = sampler.sample(random);

    //Reset seeds
    random.setSeeds(seed1,seed2);

    //Sample local indexes
    std::array<unsigned long, 4> localIndexes = sampler.sampleByDim(random);

    //Compare both
    if(globIndex !=
       localIndexes[3]*nxyz +
       localIndexes[2]*nxy +
       localIndexes[1]*nx +
       localIndexes[0]){
      printf("Error: Indexes mismatch!\n"
	     "       Global index : %lu\n"
	     "       Local indexes: %lu %lu %lu %lu\n",
	     globIndex,
	     localIndexes[0],
	     localIndexes[1],
	     localIndexes[2],
	     localIndexes[3]);
    }
  }

  return 0;
}
