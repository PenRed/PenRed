
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

 
#ifndef __PENRED_AUXILIAR_FUNCTIONS__
#define __PENRED_AUXILIAR_FUNCTIONS__
 
#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <cmath>
 
 template <class type>
unsigned seeki(const type* x, const type xf, const unsigned n){

  //Finds the interval (x[i],x[i+1]] containing the value "xf".
  // x is supposed to be continuous ascendent
  
  unsigned low, top;

  if(xf > x[n-1]){
    throw std::out_of_range("Error in seeki, xf is greater than x[n-1]");
  }

  if(xf < x[0]){
    throw std::out_of_range("Error in seeki, xf is lower than x[0]");    
  }

  low = 0;
  top = n-1;

  do{

    unsigned mid = (low+top)/2;
    if(xf > x[mid]){
      low = mid;
    } else {
      top = mid;
    }
    
  }while(top-low > 1);

  return low;
}

//Function to check if a integer is power of 2
template<class T>
//Ensure that T is an unsigned integral type
inline typename std::enable_if<std::is_unsigned<T>::value,bool>::type  
isPowerOf2(T n){
  return n && !(n & (n - 1));
}

//Function to get the log base 2 of an integer
template<class T>
//Ensure that T is an unsigned integral type
inline typename std::enable_if<std::is_unsigned<T>::value,unsigned int>::type  
logb2(T n){
  unsigned int e = 0;
  while(n >>= 1) ++e;
  return e;
}

#endif
