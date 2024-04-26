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
  constexpr unsigned xbins = 100;
  constexpr unsigned ybins = 100;
  
  std::vector<double> dis(xbins*ybins);

  for(size_t j = 0; j < ybins; ++j){
    for(size_t i = 0; i < xbins; ++i){
      dis[j*xbins + i] =
	2.0+sin(static_cast<double>(j)+0.5)*cos(static_cast<double>(i)+0.5);
    }
  }

  //Init sampler
  penred::sampling::aliasing<2> sampler;
  sampler.init(dis, {100,100},
	       {std::pair<double, double>(0.0,10.0),
		std::pair<double, double>(0.0,10.0)});

  //Create random number generator
  pen_rand random;

  //Create results vector
  std::vector<double> res(xbins*ybins,0.0);
  
  constexpr unsigned long nIter = 10000000000;
  for(size_t i = 0; i < nIter; ++i){
    ++res[sampler.sample(random)];
  }

  //Get normalization factor of both distributions
  double normDis = std::accumulate(dis.cbegin(), dis.cend(), 0.0);
  double normRes = std::accumulate(res.cbegin(), res.cend(), 0.0);
  
  //Print results
  FILE* fout = fopen("results.dat", "w");
  for(size_t j = 0; j < ybins; ++j){
    for(size_t i = 0; i < xbins; ++i){
      fprintf(fout,"%E %E\n",
	      dis[j*xbins + i]/normDis,
	      res[j*xbins + i]/normRes);
    }
    fprintf(fout,"\n");
  }
  fclose(fout);
  
  //Compare both values  
  for(size_t i = 0; i < xbins*ybins; ++i){
    if(std::fabs((dis[i]/normDis - res[i]/normRes)/(dis[i]/normDis)) > 1.0e-2){
      printf("Results mismatch!\n"
	     " distribution : %15.5E\n"
	     " generated    : %15.5E\n",
	     dis[i]/normDis, res[i]/normRes);
      return 2;
    }
  }


  // ** Test sampling indexes
  //Create results vector
  std::vector<double> res2(xbins*ybins,0.0);
  
  for(size_t i = 0; i < nIter; ++i){
    const std::array<unsigned long, 2> indexes = sampler.sampleByDim(random);
    ++res2[indexes[1]*xbins + indexes[0]];
  }

  //Get normalization factor of both distributions
  double normRes2 = std::accumulate(res2.cbegin(), res2.cend(), 0.0);
  
  //Print results
  fout = fopen("results2.dat", "w");
  for(size_t j = 0; j < ybins; ++j){
    for(size_t i = 0; i < xbins; ++i){
      fprintf(fout,"%E %E\n",
	      dis[j*xbins + i]/normDis,
	      res2[j*xbins + i]/normRes2);
    }
    fprintf(fout,"\n");
  }
  fclose(fout);
  
  //Compare both values  
  for(size_t i = 0; i < xbins*ybins; ++i){
    if(std::fabs((dis[i]/normDis - res2[i]/normRes2)/(dis[i]/normDis)) > 1.0e-2){
      printf("Results 2 mismatch!\n"
	     " distribution : %15.5E\n"
	     " generated    : %15.5E\n",
	     dis[i]/normDis, res2[i]/normRes2);
      return 2;
    }
  }  

  return 0;
}
