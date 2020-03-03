
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

 
 
#ifndef __PEN_PARTICLES_ID__ 
#define __PEN_PARTICLES_ID__ 

#include <cstring>

enum pen_KPAR{
  PEN_ELECTRON,
  PEN_PHOTON,
  PEN_POSITRON,

  ALWAYS_AT_END
};

inline bool isKpar(const unsigned kpar){return (kpar < ALWAYS_AT_END);}

inline const char* particleName(const unsigned kpar){
  switch(kpar){
  case PEN_ELECTRON: return "electron";
  case PEN_PHOTON: return "gamma";
  case PEN_POSITRON: return "positron";
  default: return nullptr;
  }
}

inline unsigned particleID(const char* name){
    
    for(unsigned i = 0; i < ALWAYS_AT_END; i++){
        if(strcmp(name,particleName(i)) == 0)
            return i;
    }
    return ALWAYS_AT_END;
}

#endif
