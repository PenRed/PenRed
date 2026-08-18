
//
//
//    Copyright (C) 2021 Universitat de València - UV
//    Copyright (C) 2021 Universitat Politècnica de València - UPV
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


#include "fileSpectrum_energySampling.hh"


void fileSpectrum_energySampling::energySampling(double& energy,
						 pen_rand& random) const{

  // Get random number
  double rand = random.rand();

  // Get interval
  unsigned interval = seeki(&cumulative.front(),rand,nEBins);

  // Calculate sampled energy
  energy = energies[interval]+(rand-cumulative[interval])*dE[interval];
}

int fileSpectrum_energySampling::configure(double& Emax,
					   const pen_parserSection& config,
					   const unsigned verbose){

  int err;
  std::string filename;
  err = config.read("filename",filename);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("fileSpectrum_energySampling:configure:unable to read "
	     "'filename' in configuration. String expected\n");
    }
    return -1;
  }

  //Open the file
  FILE* fin = nullptr;
  fin = fopen(filename.c_str(),"r");
  if(fin == nullptr){
    if(verbose > 0){
      printf("fileSpectrum_energySampling:configure: Error: unable "
	     "to open spectrum file '%s'\n",filename.c_str());
    }
    return -2;
  }

  if(verbose > 1){
    printf("Spectrum file: %s\n",filename.c_str());
  }

  //Clear vectors
  energies.clear();
  dE.clear();
  cumulative.clear();
  nEBins = 0;

  //Set initial cumulative value (0.0)
  cumulative.push_back(0.0);
  
  //Read the data
  char line[1000];
  long unsigned nRead = 0;
  long unsigned auxRead;
  double prevE = -1;
  while(pen_getLine(fin,1000,line,auxRead) == 0){
    nRead += auxRead;
    //Read energy and probability
    double E,p;
    int nscan = sscanf(line,"%lE %lE",&E,&p);
    if(nscan != 2){
      if(verbose > 0)
	printf("fileSpectrum_energySampling:configure: Error in spectrum file."
	       " Unexpected format in line %lu: \n    %s\n",nRead,line);
      return -3;
    }

    if(E <= 0.0){
      if(verbose > 0)
	printf("fileSpectrum_energySampling:configure: Error in spectrum file."
	       " Invalid energy in line %lu: \n    %s\n",nRead,line);
      return -4;
    }

    //Store energy and energy interval
    energies.push_back(E);
    if(!std::signbit(prevE)){
      if(E < prevE){
	if(verbose > 0)
	  printf("fileSpectrum_energySampling:configure: Error in spectrum file."
		 " Ascending energy order is not accomplished in line %lu:\n"
		 "  Previous energy: %.4E\n"
		 "      Next energy: %.4E\n",nRead,prevE,E);
	return -5;
      }
      //Store the interval size
      dE.push_back(E-prevE);
    }
    //Update previous energy
    prevE = E;

    //Check if this is the last value to store
    if(std::signbit(p)){
      //Finish the read
      break;
    }
    cumulative.push_back(p);
    
  }

  if(cumulative.size() != energies.size()){
    printf("fileSpectrum_energySampling:configure: Error: "
	   " Unexpected error, please, report it. CODE: 1\n");
    return -999;    
  }
  
  //Obtain the cumulative probabilities
  nEBins = cumulative.size();
  for(unsigned i = 2; i < nEBins; ++i)
    cumulative[i] += cumulative[i-1];
  
  
  //Normalize the probabilities
  double totalProb = cumulative.back();
  if(totalProb <= 1.0e-15){
    if(verbose > 0)
      printf("fileSpectrum_energySampling:configure: Error: "
	     " Null probability\n");
    return -6;
  }
  for(double& p : cumulative)
    p /= totalProb;

  //Precalculate weighted bin widths
  for(unsigned i = 0; i < nEBins-1; ++i){
    double diffcumul = cumulative[i+1]-cumulative[i];
    if(diffcumul < 1.0e-15)
      dE[i] = 0.0;
    else
      dE[i] /= diffcumul;
  }

  //Get maximum energy
  Emax = energies.back();

  if(verbose > 1){
    printf("      Energy          Etop       LowCumul"
	   "     topCumul     prob\n");
    for(unsigned j = 0; j < nEBins-1; j++)
      printf(" %12.4E  %12.4E  %12.4E  %12.4E  %12.4E\n",
	     energies[j],energies[j+1],
	     cumulative[j],cumulative[j+1],
	     cumulative[j+1]-cumulative[j]);

    printf("Maximum possible energy: %12.4E eV\n",Emax);
  }  
  
  return 0;
}

REGISTER_SAMPLER(fileSpectrum_energySampling, FILE_SPECTRUM)
