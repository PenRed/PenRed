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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

#include "materialCreator.hh"
#include <stdio.h>

int main (){

  //******************************************************
  //Fix the minimum number of exponent digits in MVS to 2 
#ifdef _MSC_VER
  unsigned int prev_exponent_format =
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
  //******************************************************  

  penred::penMaterialCreator::materialCreator creator;
  
  creator.PEMATW();
  if(creator.IRETRN != 0){
    printf ("%s\n", creator.REASON);
    printf ("IRETRN =%d\n", creator.IRETRN);
    return -1;
  }
  return 0;
}
