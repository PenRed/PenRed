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
//        vicent.gimenez.alventosa@gmail.com
//    
//


#include "pen_muen.hh"

int main(int argc, const char** argv){

  //Check the arguments
  if(argc < 5){
    printf("usage: %s material-file tolerance sim-time "
	   "energy-spectrum-file-1 energy-spectrum-file-2 ... \n",argv[0]);
    return 1;
  }

  //Get configuration parameters

  double tolerance = atof(argv[2]);
  double simTime = atof(argv[3]);

  if(tolerance <= 1.0e-10){
    printf("Error: Tolerance must be positive and"
	   " greater than zero: %E\n",tolerance);
    return -1;
  }

  //Get number of provided spectrums
  const unsigned nSpectrums = static_cast<unsigned>(argc - 4);

  //Open output file
  FILE* fout =  nullptr;
  fout =  fopen("mutren.dat","w");
  if(fout == nullptr){
    printf("Error: Unable to create file 'mutren.dat'");
    return -4;
  }

  //Print header information
  fprintf(fout,"#  Photon mass energy-absorption coefficients\n");
  fprintf(fout,"#  Material file: %s\n",argv[1]);
  fprintf(fout,"#  \n");
  fprintf(fout,"#  Tolerance: %E %%\n",tolerance);
  fprintf(fout,"#  Allotted simulation time: %E s\n",simTime);
  fprintf(fout,"#  \n");
  fprintf(fout,"#%30s    Mean E       mu/rho      mu_en/rho     "
	  "(1-g)*f        f            1-g      uncert.      \n",
	  "  Spectrum filename ");
  fprintf(fout,"#%30s     (eV)       (cm^2/g)     (cm^2/g)%45.45s(%%)       \n",
	  " "," ");
  fprintf(fout,"# \n");
  fflush(fout);

  //Calculate data for all spectrums
  std::vector<pen_muen::muData> results;
  int err = pen_muen::calculate(&argv[4], nSpectrums,
				tolerance, simTime,
				argv[1],
				results, 2);
  
  for(unsigned is = 0; is < nSpectrums; ++is){

    pen_muen::muData& result = results[is];
    fprintf(fout,"%-30s   %.5E  %.5E  %.5E  %.5E  %.5E  %.5E  %.1E\n",
	    argv[4+is],
	    result.E0,result.muRho,result.muen(),
	    result.fg/result.E0,result.f/result.E0,
	    result.fg/result.f,result.err);
  }

  if(err != 0){
    printf("Error on calculation.\n");
  }

  fclose(fout);

  return 0;
}
