
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

 
#ifndef __GAMMA_POL_STATE__
#define __GAMMA_POL_STATE__

#include "pen_baseState.hh"

struct pen_state_gPol : public pen_particleState{
  //  ****  Photon polarisation.
  //  ----  Polarised photons if IPOL=1, otherwise unpolarised photons.
  int IPOL;
  //  ----  Stokes parameters.
  double SP1, SP2, SP3;
  
  pen_state_gPol() : pen_particleState(), IPOL(0), SP1(0.0), SP2(0.0), SP3(0.0){}
  
  void reset(){

    baseReset();
    IPOL = 0;
    SP1 = 0.0;
    SP2 = 0.0;
    SP3 = 0.0;
  }
 
 std::string stringify() const{
     char text[300];
     sprintf(text,"%s  %12.4E  %12.4E  %12.4E  %d",stringifyBase().c_str(),
             SP1,SP2,SP3,IPOL);
     
     return std::string(text);
 }
 
};

#endif
