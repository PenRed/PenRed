
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



#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include "pen_materials.hh"


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

  // Load material
  //****************
  
  //Open material file
  FILE* fmat = 0;
  fmat = fopen(argv[1],"r");
  //Check if material has been opened
  if(fmat == nullptr)
    {
      printf("Error: Can't open material file '%s'\n",argv[1]);
      return -1;
    }

  //Load material
  initStructs initStore;
  pen_elementDataBase elements;
  pen_logGrid grid;
  
  //Initialize grid energy limits
  grid.init(50.0, 1.0e9);

  pen_material* mat = new pen_material;
  mat->load(fmat, stdout, initStore,elements,grid, 0);
  if(penGetError() != PEN_SUCCESS)
    {
      printf("%s",penErrorString(penGetError()));
      return -3;
    }

  //Print ranges for for specified energies
  for(double E : energies){

    //Get energy interval
    int KE;
    double XEL, XE, XEK;
    grid.getInterval(E,KE,XEL,XE,XEK);

    for(size_t i = 0; i < constants::nParTypes; ++i)
      printf("# range for particle %16s with %12.4E eV: %12.4E cm\n", particleName(i), E, exp(mat->RANGEL[i][KE] + (mat->RANGEL[i][KE+1]-mat->RANGEL[i][KE])*XEK));
  }
  
  delete mat;
  return 0;
}
