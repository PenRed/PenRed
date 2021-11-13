
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
#include "pen_materials.hh"


int main(int argc, char** argv)
{

  //Check argument number
  if(argc < 4)
    {
      printf("Usage: %s ELOW EMAX material-file1 material-file2 ...\n",argv[0]);
      return 1;
    }
    
  //******************************************************
  //Fix the minimum number of exponent digits in MVS to 2 
#ifdef _MSC_VER
  unsigned int prev_exponent_format =
      _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
  //******************************************************      

  if(strcmp(argv[1],"--help") == 0 || strcmp(argv[1],"-help") == 0 || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--h") == 0)
    {
      printf("Usage: %s material-file1 material-file2 ...\n",argv[0]);
      return 0;
    }

  //Get Elow and Emax
  double Elow = atof(argv[1]);
  double Emax = atof(argv[2]);

  if(Elow >= Emax)
    {
      printf("Error: Emax must be greater than Elow.\n");
      printf("       Elow: %e\n",Elow);
      printf("       Emax: %e\n",Emax);
      return 2;
    }
      
  //Load materials
  for(int i = 3; i < argc; i++)
    {
      //Open material file
      FILE* fmat = 0;
      fmat = fopen(argv[i],"r");

      //Check if material has been opened
      if(fmat == 0)
	{
	  printf("Error: Can't open material file '%s'\n",argv[i]);
	  return -1;
	}

      printf("\n\nReading material %d: %s\n\n",i-3,argv[i]);
      
      //Load material
      initStructs initStore;
      pen_elementDataBase elements;
      pen_logGrid grid;
      
      //Initialize grid energy limits
      grid.init(Elow, Emax);

      pen_material* mat = new pen_material;
      mat->load(fmat, stdout, initStore,elements,grid, 4);
      if(penGetError() != PEN_SUCCESS)
	{
	  printf("%s",penErrorString(penGetError()));
	  return -3;
	}
      delete mat;
    }
  
  return 0;
}
