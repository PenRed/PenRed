
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

 
#ifndef __PEN_CONTEXT_
#define __PEN_CONTEXT_

#include <cmath>
#include <cstdlib>


//-------------------
// Context
//-------------------

class pen_context : public abc_context<pen_material>{
  
public:
  
  // ****  Energy grid and interpolation constants.
  pen_logGrid grid;

  // Element data
  pen_elementDataBase& elements;

  // RITA random sampling variables
  CRNDG3 rndg3;

  //  ****  Variance-reduction parameters.

  //  ----  Parameter values defined for each body. NBV is the maximum
  //        number of bodies. When using PENGEOM, NBV must have the same
  //        value as the parameter NB in module PENGEOM_mod.
  const static unsigned int NBV = pen_geoconst::NB;
  //  ----  Forcing factors, FORCE(IBODY,KPAR,ICOL).
  double FORCE[NBV][constants::nParTypes][constants::MAXINTERACTIONS];
  //  ----  Interaction forcing bodies state 
  bool  LFORCE[NBV][constants::nParTypes];
  // WRANGES store interaction forcing weight ranges (wlow,wup)
  // for each particle in each body
  double WRANGES[NBV][2*constants::nParTypes];
  //  ----  Bremsstrahlung splitting numbers, IBRSPL(IBODY).
  unsigned int IBRSPL[NBV];
  //  ----  x-ray splitting numbers, IXRSPL(IBODY).
  unsigned int IXRSPL[NBV];
  bool         LXRSPL[NBV];
  
  pen_context(pen_elementDataBase& inElements) : elements(inElements){
    
    //  ****  Variance-reduction parameters.
    
    // Forcing
    for(unsigned int i = 0; i < NBV; i++)
      for(unsigned int j = 0; j < constants::nParTypes; j++){
	WRANGES[i][2*j]   = 0.0;
	WRANGES[i][2*j+1] = 1.0e6;
	LFORCE[i][j] = false;
	for(unsigned int k = 0; k < constants::MAXINTERACTIONS; k++)
	  FORCE[i][j][k] = 1.0;
      }
    
    // Bremsstrahlung splitting
    for(unsigned int i = 0; i < NBV; i++)
      IBRSPL[i] = 1;

    // x-ray splitting
    for(unsigned int i = 0; i < NBV; i++){
      IXRSPL[i] = 1;
      LXRSPL[i] = false;
    }
  }
  int init(double EMAX, FILE *IWR, int INFO, std::string PMFILE[constants::MAXMAT]);

  inline bool isForcing(const unsigned kpar,
			const unsigned ibody,
			const double wght) const {
    unsigned kpar2 = 2*kpar;
    if(LFORCE[ibody][kpar] &&
       wght >= WRANGES[ibody][kpar2] &&
       wght <= WRANGES[ibody][kpar2+1]){
      return true;
    }
    return false;
  }
  
};

//-----------------------------------------------
// PENELOPE common  functions
//-----------------------------------------------

void DIRECT(const double CDT,
	    const double DF,
	    double &U,
	    double &V,
	    double &W);


void DIRPOL(const double CDT,
	    double &DF,
	    const double CONS,
	    double &SP1,
	    double &SP2,
	    double &SP3,
	    double &U,
	    double &V,
	    double &W,
	    pen_rand& random);

void RELAX(const pen_elementDataBase& elements,
	   const pen_material& mat,
	   pen_particleState& state,
	   const int ICOL,
	   const int MODER,
	   const int IZ,
	   const int IS,
	   int &KS,
	   pen_particleStack<pen_particleState>& stackE,
	   pen_particleStack<pen_state_gPol>& stackG,
	   pen_rand& random);

void PANaR(const pen_particleState& betaState,
	   const double ECUT,
	   pen_particleStack<pen_state_gPol>& stackG,
	   pen_rand& penRand);


#endif
