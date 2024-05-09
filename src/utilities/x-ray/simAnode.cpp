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

  if(argc < 8){
    printf("usage: %s mat-filename beamEnergy minE anode-angle "
	   "bins pixel-size histories max-angle\n", argv[0]);
    return 1;
  }

  //Extract data
  const char* matFilename = argv[1];
  const double beamE = std::atof(argv[2]);
  const double eMin = std::atof(argv[3]);
  const double angle = std::atof(argv[4]);
  const double nBinsAux = std::atof(argv[5]);
  const double pixelSize = std::atof(argv[6]);
  const double nHistAux = std::atof(argv[7]);

  double maxAngle = 0.0;
  if(argc > 8)
    maxAngle = std::atof(argv[8]);

  if(nBinsAux <= 0.0){
    printf("Number of bins must be greater than 0.\n");
    return -1;
  }
  
  if(nHistAux <= 0.0){
    printf("Number of histories must be greater than 0.\n");
    return -2;
  }

  if(maxAngle <= 0.0)
    maxAngle = 80.0;

  const unsigned long nBins = static_cast<unsigned long>(nBinsAux);
  const unsigned long long nHist = static_cast<unsigned long long>(nHistAux);

  //Create tallies
  penred::measurements::measurement<double,1> spectrum;
  penred::measurements::measurement<double,2> spatialDistrib;

  spectrum.init({nBins}, {std::pair<double,double>(eMin,beamE)});
  
  //Simulate the anode
  double dReg;
  int err = penred::xray::simAnodeDistrib(matFilename, beamE, eMin, pixelSize, angle,
					  nHist, dReg, spectrum, spatialDistrib,
					  maxAngle, 2);
  if(err != 0){
    printf("Error on anode simulation.\n");
    return -2;
  }

  //Print the spectrum
  FILE* fout = fopen("spectrum.dat", "w");
  spectrum.print(fout, nHist, 2, true, false);
  fclose(fout);

  //Print the spatial distribution
  fout = fopen("spatialDistrib.dat", "w");
  spatialDistrib.print(fout, nHist, 2, true, true);
  fclose(fout);

  printf("Data recorded at %f cm from the anode face center.\n", dReg);
  
  return 0;  
}
