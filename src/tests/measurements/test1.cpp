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

#include "math_classes.hh"
#include <random>
#include <fstream>

int main(){

  //Create a 2D measurement
  penred::measurements::measurement<double, 3> genData;
  genData.init({220,150,1},
	       {std::pair<double,double>(-1.0,20.2),
		std::pair<double,double>(-20.0,35.6),
		std::pair<double,double>(-1000.0,1000.0)});

  genData.description = "Test measurement generated data\nfrom C++ random number generation.";
  genData.setDimHeader(0, "Position x(cm) with long| text ");
  genData.setDimHeader(2, "z (cm)");
  genData.setValueHeader("Value (prob)");

  std::random_device rd;
  std::mt19937 gen(rd()); 
  std::uniform_real_distribution<> disX(-1.0, 20.2);  
  std::uniform_real_distribution<> disY(-20.0, 35.6);
  std::uniform_int_distribution<> partPerHist(10, 100);  

  constexpr unsigned long long nhists = 1000000;
  for(unsigned long long i = 0; i < nhists; ++i){

    //Get generated particles
    int genPart = partPerHist(gen);
    for(int j = 0; j < genPart; ++j){
      genData.add({disX(gen), disY(gen), 0.0},
		  1.0/static_cast<double>(genPart),
		  i);
    }
  }

  //Print resulting data
  FILE* fout = nullptr;
  fout = fopen("testResult.dat","w");
  if(fout == nullptr){
    printf("Unable to open file\n");
    return 1;
  }
  genData.print(fout, nhists, 2, true, true);
  fclose(fout);
  fout = nullptr;

  //Print resulting data
  fout = fopen("testResult2.dat","w");
  if(fout == nullptr){
    printf("Unable to open file\n");
    return 1;
  }
  genData.print(fout, nhists, 2, false, true);
  fclose(fout);  
  fout = nullptr;

  //Print resulting data
  fout = fopen("testResult3.dat","w");
  if(fout == nullptr){
    printf("Unable to open file\n");
    return 1;
  }
  genData.print(fout, nhists, 2, true, false);
  fclose(fout);  
  fout = nullptr;

  //Print resulting data
  fout = fopen("testResult4.dat","w");
  if(fout == nullptr){
    printf("Unable to open file\n");
    return 1;
  }
  genData.print(fout, nhists, 2, false, false);
  fclose(fout);  
  fout = nullptr;
  
  //Read resulting data
  penred::measurements::results<double, 3> readResults;

  std::ifstream fin("testResult.dat", std::ifstream::in);
  if(!fin){
    printf("Unable to open file\n");
    return 1;
  }
  int err = readResults.read(fin);
  if(err != 0){
    printf("Error reading results.\n"
	   "  Error code: %d\n", err);
    return 1;
  }
  fin.close();

  //Get results with no read
  penred::measurements::results<double, 3> genResults;
  genData.results(nhists, genResults);
  
  //Compare both values  
  for(size_t i = 0; i < readResults.getNBins(); ++i){
    if(std::fabs((genResults.data[i] - readResults.data[i])/genResults.data[i]) > 1.0e-4){
      printf("Results mismatch!\n"
	     " read result     : %15.5E\n"
	     " generated result: %15.5E\n",
	     readResults.data[i], genResults.data[i]);
      return 2;
    }

    if(std::fabs((genResults.sigma[i] - readResults.sigma[i])/genResults.sigma[i]) > 1.0e-2){
      printf("Sigma mismatch!\n"
	     " read sigma     : %15.5E\n"
	     " generated sigma: %15.5E\n",
	     readResults.sigma[i], genResults.sigma[i]);
      return 2;
    }
  }

  //Read resulting data with one dimension less
  penred::measurements::results<double, 2> readResultsLess;

  std::ifstream fin2("testResult2.dat", std::ifstream::in);
  if(!fin2){
    printf("Unable to open file\n");
    return 1;
  }
  err = readResultsLess.read(fin2);
  if(err != 0){
    printf("Error reading results.\n"
	   "  Error code: %d\n", err);
    return 1;
  }
  fin2.close();
  
  //Compare both values  
  for(size_t i = 0; i < readResultsLess.getNBins(); ++i){
    if(std::fabs((genResults.data[i] - readResultsLess.data[i])/genResults.data[i]) > 1.0e-4){
      printf("Results mismatch!\n"
	     " read result     : %15.5E\n"
	     " generated result: %15.5E\n",
	     readResultsLess.data[i], genResults.data[i]);
      return 2;
    }

    if(std::fabs((genResults.sigma[i] - readResultsLess.sigma[i])/genResults.sigma[i]) > 1.0e-2){
      printf("Sigma mismatch!\n"
	     " read sigma     : %15.5E\n"
	     " generated sigma: %15.5E\n",
	     readResultsLess.sigma[i], genResults.sigma[i]);
      return 2;
    }
  }

  //Print both results
  fout = nullptr;
  fout = fopen("results1.dat", "w");
  readResults.print(fout, 2, true, true);
  fclose(fout);

  fout = nullptr;
  fout = fopen("results2.dat", "w");
  readResultsLess.print(fout, 2, true, true);
  fclose(fout);
  
  return 0;
}
