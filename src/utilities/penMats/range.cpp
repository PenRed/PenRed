
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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



#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include "PenRed.hh"

int main(int argc, char** argv)
{

  //Check argument number
  if(argc < 3)
    {
      printf("Usage: %s material-file E1 E2 E3 ...\n",argv[0]);
      return 1;
    }

  if(strcmp(argv[1],"--help") == 0 || strcmp(argv[1],"-help") == 0 || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--h") == 0)
    {
      printf("Usage: %s material-file E1 E2 E3 ...\n",argv[0]);
      return 0;
    }

  //Store energies
  std::vector<double> energies;
  for(int i = 2; i < argc; ++i)
    energies.push_back(atof(argv[i]));

  //Sort energies
  std::sort(energies.begin(), energies.end());
  
  for(double E : energies){
    if(E < 50.0 || E > 1.0e9){
      printf("Invalid energy (%12.4E)\n",E);
      printf("Energy must be in the range (50,1e9) eV\n");
      return -2;
    }
  }

  // Configure context
  //*********************

  //Create elements data base
  pen_elementDataBase elements;
  //Create a context
  pen_context context(elements);

  //Set the number of materials to context (1)
  int errmat = context.setMats<pen_material>(1);
  if(errmat != 0){
    printf("Error at context material creation: %d.\n",errmat);
    return -1;
  }  
  
  //Get the material
  pen_material& mat = context.getBaseMaterial(0);

  //Configure the material
  mat.C1=0.2;
  mat.C2=0.2;
  mat.WCC=1.0e3;
  mat.WCR=1.0e3;

  mat.EABS[PEN_ELECTRON] = 50.0E0;
  mat.EABS[PEN_PHOTON]   = 50.0E0;
  mat.EABS[PEN_POSITRON] = 50.0E0;

  //Configure context
  printf("  \n");
  printf("Configuring context with material '%s'. Please, wait...\n",argv[1]);

  FILE* fcontext = nullptr;
  fcontext = fopen("context.rep","w");
  if(fcontext == nullptr){
    printf("Error: unable to create file 'context.rep'\n");
    return -2;
  }
  
  double EMAX=1.0E9;
  int INFO = 2;
  std::string PMFILEstr[constants::MAXMAT];
  PMFILEstr[0].assign(argv[1]);
  int err = context.init(EMAX,fcontext,INFO,PMFILEstr);
  if(err != 0){
    printf("Error: Unable to configure context. Check 'context.rep'.\n");
    return -3;
  }
  fclose(fcontext);

  //Print ranges for for specified energies
  for(double E : energies){

    //Get energy interval
    int KE;
    double XEL, XE, XEK;
    context.grid.getInterval(E,KE,XEL,XE,XEK);

    for(size_t i = 0; i < constants::nParTypes; ++i)
      printf("# range for particle %16s with %12.4E eV: %12.4E cm\n", particleName(i), E, context.range(E,static_cast<pen_KPAR>(i),0));
  }
  
  return 0;
}
