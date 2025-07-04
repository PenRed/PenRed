
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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

//Interactions enumerations

enum pen_betaE_interact{
  BETAe_HARD_ELASTIC = 0,
  BETAe_HARD_INELASTIC,
  BETAe_HARD_BREMSSTRAHLUNG,
  BETAe_HARD_INNER_SHELL,
  BETAe_DELTA,
  BETAe_SOFT_INTERACTION,
  BETAe_HARD_TOTAL
};

enum pen_gamma_interact{
  GAMMA_RAYLEIGH = 0,
  GAMMA_COMPTON,
  GAMMA_PHOTOELECTRIC,
  GAMMA_PAIR_PRODUCTION,
  GAMMA_DELTA  
};

enum pen_betaP_interact{
  BETAp_HARD_ELASTIC = 0,
  BETAp_HARD_INELASTIC,
  BETAp_HARD_BREMSSTRAHLUNG,
  BETAp_HARD_INNER_SHELL,
  BETAp_ANNIHILATION,
  BETAp_DELTA,
  BETAp_SOFT_INTERACTION,
  BETAp_HARD_TOTAL
};

inline bool isKpar(const unsigned kpar){return (kpar < ALWAYS_AT_END);}

inline bool isRealInteraction(unsigned kpar, unsigned inter) {
  switch(kpar){
  case PEN_ELECTRON: return inter != BETAe_DELTA && inter <= BETAe_HARD_TOTAL;
  case PEN_PHOTON: return inter < GAMMA_DELTA;
  case PEN_POSITRON: return inter != BETAp_DELTA && inter <= BETAp_HARD_TOTAL;
  default: return false;
  }
}

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
