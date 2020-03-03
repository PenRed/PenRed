
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


#include <cmath>
#include "pen_phaseSpaceFile.hh"

int main(int argc, char** argv){

 if(argc < 5){
    printf("usage: %s filename emin emax nbins\n",argv[0]);
    return 1;
  }

  
  //Create tally variables
  const int maxBins = 20000;
  double spectre[maxBins][constants::nParTypes+1];

  for(unsigned i = 0; i < maxBins; i++){
    for(unsigned j = 0; j <= constants::nParTypes; j++)
      spectre[i][j] = 0.0e0;
  }
  
  //Get emin, emax and nbins
  double emin = atof(argv[2]);
  double emax = atof(argv[3]);
  double de,ide;
  int nbins = atoi(argv[4]);

  
  //Check values
  if(emin >= emax){
    printf("Error: 'emin' must be greater than 'emax'\n");
    return -1;
  }

  if(nbins <= 0 || nbins > maxBins){
    printf("Error: Invalid number of bins (%d)\n",nbins);
    return -2;
  }

  //Calculate bin width
  de = (emax - emin) / (double)nbins;
  ide = 1.0/de;
  
  //Create a phase space file
  pen_psfreader psf;

  printf("Expected particles per chunk: %lu (%lu B)\n",psf.bufferSize(),psf.memory());
  
  //Open specified input file
  FILE* fin = nullptr;
  fin = fopen(argv[1],"rb");
  if(fin == nullptr){
    printf("Error: unable to open file '%s'.\n",argv[1]);
    return -3;
  }

  //Open output file
  FILE* fout = nullptr;
  fout = fopen("spectre.psf","w");
  if(fout == nullptr){
    printf("Error: unable to open file '%s'.\n","spectre.psf");
    return -4;
  }

  //Read input file until the end
  unsigned nchunks = 0;
  long long unsigned nhists = 0;
  while(psf.read(fin,1) == PEN_PSF_SUCCESS){
    nchunks++;
    //Iterate over read states
    pen_particleState state;
    unsigned long dhist;
    unsigned kpar;
    while(psf.get(dhist,kpar,state) > 0){

      //printf("%s\n",state.stringify().c_str());
      //Calculate bin index
      int ibin = (state.E-emin)*ide;

      nhists += dhist;
      //Check bin index
      if(ibin >= 0 && ibin < nbins){
	spectre[ibin][kpar+1] += state.WGHT;
	spectre[ibin][0] += state.WGHT;
      }
    }
  }

  printf("Read psf chunks: %u (%lu B)\n",nchunks,size_t(nchunks)*psf.memory());
  printf("Number of histories: %llu\n",nhists);
  
  //Write spectrum to output file
  fprintf(fout,"# '%s' psf spectrum\n",argv[1]);
  fprintf(fout,"# Emin, Emax, nbins:\n");
  fprintf(fout,"# %12.4E %12.4E %d:\n",emin,emax,nbins);
  fprintf(fout,"# \n");

  fprintf(fout,"#  bin  |     E(eV)    | ");

  fprintf(fout,"%15.15s | ","Total");
  for(unsigned i = 0; i < constants::nParTypes; i++){
    fprintf(fout,"%15.15s | ",particleName(i));
  }
  fprintf(fout,"\n");

  for(int i = 0; i < nbins; i++){
    fprintf(fout," %6d   %12.4E   ",i,double(i)*de+emin);
      for(unsigned j = 0; j <= constants::nParTypes; j++){
	fprintf(fout,"%15.4E   ",spectre[i][j]);
    }
      fprintf(fout,"\n");
  }

  //Close files
  fclose(fin);
  fclose(fout);
  
  return 0;
}
