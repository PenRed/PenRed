
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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

 
#ifndef __PEN_RANDOM_
#define __PEN_RANDOM_

#include <cstdio>
#include <cmath>
#include <cstring>
#include "pen_errors.hh"

//Typedef for RITA function argument
typedef double (* PDF_RITA)(double, void*);

struct CRITAA
{
  static const unsigned int NM=512;
  double XA[NM], AA[NM], BA[NM], FA[NM];
  int IA[NM], NPM1A; 
};

struct CRNDG3
{
  static const unsigned int NR = 128;
  double X[NR], A[NR], B[NR], F[NR];
  int KA[NR], NPM1;
};

struct CRITA{
  const static unsigned int NM=512;
  double XT[NM], PAC[NM], DPAC[NM], A[NM], B[NM];
  int IL[NM], IU[NM], NPM1;  
  //double* XTI=XT;
  //double* QTI=XT;
  //double* PACI=PAC;
  //double* DPACI=DPAC;
  //double* AI=A;
  //double* BI=B;
  //int* ITLI=IL;
  //int* ITUI=IU;
  //int* NPM1I=&NPM1;
};


class pen_rand{

private:

  int ISEED1, ISEED2;
  
public:

  pen_rand() : ISEED1(1), ISEED2(1) {}
  
  inline double rand(){
  //  This is an adapted version of subroutine RANECU written by F. James
  //  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
  //  give a single random number at each call.
  //
  //  The 'seeds' ISEED1 and ISEED2 must be initialised in the main program
  //  and transferred through the named common block /RSEED/.
  //
  //  Some compilers incorporate an intrinsic random number generator with
  //  the same name (but with different argument lists). To avoid conflict,
  //  it is advisable to declare RAND as an external function in all sub-
  //  programs that call it.

  const double USCALE = 1.0E0/2.147483563E9;

  int I1 = ISEED1/53668;
  ISEED1 = 40014*(ISEED1-I1*53668)-I1*12211;
  if(ISEED1 < 0){ ISEED1 = ISEED1+2147483563;}

  int I2 = ISEED2/52774;
  ISEED2 = 40692*(ISEED2-I2*52774)-I2*3791;
  if(ISEED2 < 0){ ISEED2 = ISEED2+2147483399;}

  int IZ = ISEED1-ISEED2;
  if(IZ < 1){ IZ = IZ+2147483562;}
  double RAND_RETURN = IZ*USCALE;
  
  return RAND_RETURN;
}
    
  void rand0(const int N);
  inline void setSeeds(const int newSeed1, const int newSeed2){
    ISEED1 = newSeed1;
    ISEED2 = newSeed2;
  }
  inline void getSeeds(int& gseed1, int& gseed2) const {
    gseed1 = ISEED1;
    gseed2 = ISEED2;
  }
};

void rand0(const unsigned N, int& seed1, int& seed2);


double RNDG3F(double X,
	      void* arg);

double RNDG3(const CRNDG3& rndg3, pen_rand& rand);

void RNDG30(CRNDG3& rndg3);

void RITA0(PDF_RITA PDF,
	   void* arg,
	   const double XLOW,
	   const double XHIGH,
	   const int N,
	   const int NU,
	   double &ERRM,
	   int fileidx,
	   CRITAA& ritaa);

void RITAI0(PDF_RITA PDF,
	    void* arg,
	    const double XLOW,
	    const double XHIGH,
	    const int N,
	    const int NU,
	    double& ERRM,
	    int fileidx,
	    CRITA& rita);

int IRND(double* FA,
	 int* IA,
	 int &N,
	 pen_rand& random);

long int IRND(const double* FA,
            const long int* IA,
            const long int  N,
            pen_rand& random);

void IRND0(double* W,
	   double* F,
	   int* K,
	   int &N);

void IRND0(const double* W,
	   double* F,
	   long int* K,
	   const long int  N);

void RITAM(double XD,
	   double XU,
	   double &XM0,
	   double &XM1,
	   double &XM2,
	   CRITA& rita);

#endif

