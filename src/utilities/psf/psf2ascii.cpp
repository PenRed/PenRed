
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


#include "pen_phaseSpaceFile.hh"

int main(int argc, char** argv){

 if(argc < 3){
    printf("usage: %s input output\n",argv[0]);
    return 1;
  }
    
  //Create a phase space file
  pen_psfreader psf;

  //Read phase space file and write it to ASCII output
  FILE* fin = nullptr;
  fin = fopen(argv[1],"rb");
  if(fin == nullptr){
    printf("Error: unable to open input file '%s'.\n",argv[1]);
    return -1;
  }

  //Open output file
  FILE* fout = nullptr;
  fout = fopen(argv[2],"w");
  if(fout == nullptr){
    printf("Error: unable to open output file '%s'.\n",argv[2]);
    return -2;
  }

  //Print header
  fprintf(fout,"# DHIST  KPAR %s\n",baseStateHeaderNoGeo());

  //Read input file until the end
  unsigned nchunks = 0;
  while(psf.read(fin,1) == PEN_PSF_SUCCESS){
    nchunks++;
    //Iterate over read states
    pen_particleState state;
    unsigned long dhist;
    unsigned kpar;
    while(psf.get(dhist,kpar,state) > 0){
      //Print base state without body and material information
      fprintf(fout,"  %5lu  %4u %s\n",dhist,kpar,state.stringifyBaseNoGeo().c_str());
    }
  }

  //Close files
  fclose(fin);
  fclose(fout);  

  printf("Processed psf chunks: %u (%lu B)\n",nchunks, static_cast<unsigned long>(size_t(nchunks)*psf.memory()));
  
  return 0;
}
