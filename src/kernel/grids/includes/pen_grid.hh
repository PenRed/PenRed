
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

 
#ifndef __PEN_GRID_
#define __PEN_GRID_

#include <cstdio>
#include <cmath>
#include "pen_errors.hh"
#include "pen_constants.hh"
#include "pen_classes.hh"

//Typedef for SUMGA function argument
typedef double (* FCT_SUMGA)(double, void*);

class pen_logGrid : public abc_grid{

public:
  //  NOTE: To determine the interval KE where the energy E is located, we
  //  do the following,
  //     XEL=LOG(E)
  //     XE=(XEL-DLEMP1)*DLFC
  //     KE=XE
  //     XEK=XE-KE  ! 'fractional' part of XE (used for interpolation).
  //

  pen_logGrid() : abc_grid() {}

  int init(double EMINu, double EMAXu);
  void getInterval(const double E, int& KE, double& XEL, double& XE, double& XEK) const;
  
};

template<size_t gridSize>
class pen_genericLogGrid : public abc_genericGrid<gridSize>{
  using baseType = abc_genericGrid<gridSize>;
public:
  //  NOTE: To determine the interval KE where the energy E is located, we
  //  do the following,
  //     XEL=LOG(E)
  //     XE=(XEL-DLEMP1)*DLFC
  //     KE=XE
  //     XEK=XE-KE  ! 'fractional' part of XE (used for interpolation).
  //

  pen_genericLogGrid() : abc_genericGrid<gridSize>() {}

  int init(double EMINu, double EMAXu){

    //  ****  Consistency of the interval end-points.
  
    if(EMINu >= EMAXu)
      {
	return 1;
      }

    //  ****  Energy grid points.

    baseType::EMIN = EMINu;
    baseType::EL = 0.99999*EMINu;
    baseType::EU = 1.00001*EMAXu;
    baseType::DLFC = log(baseType::EU/baseType::EL)/double(baseType::size-1);
    baseType::DLEMP1 = log(baseType::EL);
    baseType::DLEMP[0] = baseType::DLEMP1;
    baseType::ET[0] = baseType::EL;
    for(unsigned int I = 1; I < baseType::size; I++)
      {
	baseType::DLEMP[I] = baseType::DLEMP[I-1]+baseType::DLFC;
	baseType::ET[I] = exp(baseType::DLEMP[I]);
      }
    baseType::DLFC = (double)1.0/baseType::DLFC;

    baseType::initialized = true;
    return 0;
  }
  void getInterval(const double E, long int& KE,
		   double& XEL,double& XE,
		   double& XEK) const{

    XEL = log(E);
    XE = (XEL-baseType::DLEMP1)*baseType::DLFC;
    KE = (int)XE;
    XEK = XE-KE;
  }
  
};

void SPLINE(double* X,
	    double* Y,
	    double* A,
	    double* B,
	    double* C,
	    double* D,
	    double S1,
	    double SN,
	    int N);

void FINDI(const double* X,
	   const double XC,
	   const int N,
	   unsigned int &I);

double RLMOM(const double* X,
	     double* FCT,
	     double XC,
	     int NPpar,
	     int MOM);


void SLAG6(double &H,
	   double* Y,
	   double* S,
	   int N);

double SUMGA(FCT_SUMGA FCT,
	     void* arg,
	     double XL,
	     double XU,
	     double TOL,
	     const int ISGAW = 0);

void SINTEG(double* X,
	    double* A,
	    double* B,
	    double* C,
	    double* D,
	    double &XL,
	    double &XU,
	    double &SUM,
	    int N);

double RMOMX(const double* X,
	     const double* PDF,
	     const double XD,
	     const double XU,
	     const int NPpar,
	     const int MOM);

void RLPAC(const double* X,
	   double* PDF,
	   double* PAC,
	   int NPpar);

void MERGE2(double* X1,
	    double* Y1,
	    double* X2,
	    double* Y2,
	    double* XM,
	    double* YM,
	    int &N1,
	    int &N2,
	    int &N);

void SORT2(double* X,
	   double* Y,
	   int &N);

#endif
