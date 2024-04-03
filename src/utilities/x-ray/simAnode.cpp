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

#include "x-ray.hh"

int main(int argc, char** argv){

  if(argc < 6){
    printf("usage: mat-filename beamEnergy minE focalSpot anode-angle histories\n");
    return 1;
  }

  //Extract data
  const char* matFilename = argv[1];
  const double beamE = std::atof(argv[2]);
  const double eMin = std::atof(argv[3]);
  const double fs = std::atof(argv[4]); 
  const double angle = std::atof(argv[5]);
  const double nHistAux = std::atof(argv[6]);

  if(nHistAux <= 0.0){
    printf("Number of histories must be greater than 0.\n");
    return -1;
  }

  const unsigned long long nHist = static_cast<unsigned long long>(nHistAux);

  //Simulate the anode
  std::vector<penred::xray::detectedPart> results;
  double dReg;
  int err = penred::xray::simAnode(matFilename, beamE, eMin, fs, angle,
				   nHist, 100000000, dReg, results, true, 2);
  if(err != 0){
    printf("Error on anode simulation.\n");
    return -2;
  }

  //Construct a spectrum
  constexpr double dx = 200.0;
  std::vector<double> spec(10+static_cast<unsigned>(beamE/dx),0.0);
  const double fact = 1.0/(dx*static_cast<double>(nHist));
  for(const penred::xray::detectedPart& p : results){
    spec[static_cast<unsigned>(p.state.E/dx)] += p.state.WGHT*fact;
  }

  //Print the spectrum
  FILE* fout = fopen("spectrum.spc", "w");
  for(size_t i = 0; i < spec.size(); ++i){
    fprintf(fout,"%15.5E %15.5E\n", dx*static_cast<double>(i), spec[i]);
  }
  fclose(fout);
  
  return 0;  
}
