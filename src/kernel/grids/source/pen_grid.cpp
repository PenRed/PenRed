
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

 
#include "pen_grid.hh"



//  *********************************************************************
//                       SUBROUTINE EGRID
//  *********************************************************************
int pen_logGrid::init(double EMINu, double EMAXu)
{
  
  //This subroutine sets the energy grid where transport functions are
  //tabulated. The grid is logarithmically spaced and we assume that it
  //is dense enough to permit accurate linear log-log interpolation of
  //he tabulated functions.


  //  ****  Consistency of the interval end-points.
  
  if(EMINu < (double)constants::MINEGRID){ EMINu = (double)constants::MINEGRID;}
  if(EMINu > EMAXu - 1.0) //Mirem si la distancia entre EMINu i EMAXu, es d'1 eV com a minim
    {
      printf("   EMIN =%11.4E eV, EMAX =%11.4E eV\n", EMINu, EMAXu);
      penError(ERR_EGRID_E);
      return 1;
    }

  //  ****  Energy grid points.

  EMIN = EMINu;
  EL = 0.99999*EMINu;
  EU = 1.00001*EMAXu;
  DLFC = log(EU/EL)/double(constants::NEGP-1);
  DLEMP1 = log(EL);
  DLEMP[0] = DLEMP1;
  ET[0] = EL;
  for(unsigned int I = 1; I < constants::NEGP; I++)
  {
    DLEMP[I] = DLEMP[I-1]+DLFC;
    ET[I] = exp(DLEMP[I]);
  }
  DLFC = (double)1.0/DLFC;

  //  NOTE: To determine the interval KE where the energy E is located, we
  //  do the following,
  //     XEL=LOG(E)
  //     XE=1.0D0+(XEL-DLEMP1)*DLFC (+1.0 to avoid first position)
  //     KE=XE
  //     XEK=XE-KE  ! 'fractional' part of XE (used for interpolation).
  //

  initialized = true;
  return 0;
}

void pen_logGrid::getInterval(const double E, int& KE,
				     double& XEL, double& XE,
				     double& XEK) const
{
  XEL = log(E);
  XE = (XEL-DLEMP1)*DLFC;
  KE = (int)XE;
  XEK = XE-KE;
}

//  *********************************************************************
//                       SUBROUTINE SPLINE
//  *********************************************************************
void SPLINE(double* X, double* Y, double* A, double* B, double* C, double* D, double S1, double SN, int N)
{
  //     Cubic spline interpolation of tabulated data.
  //
  //  Input:
  //     X(I) (I=1:N) ... grid points (the X values must be in increasing
  //                      order).
  //     Y(I) (I=1:N) ... corresponding function values.
  //     S1,SN .......... second derivatives at X(1) and X(N). The natural
  //                      spline corresponds to taking S1=SN=0.
  //     N .............. number of grid points.
  //  Output:
  //     A(I),B(I),C(I),D(I) (I=1:N) ... spline coefficients.
  //
  //  The interpolating cubic polynomial in the I-th interval, from X(I) to
  //  X(I+1), is
  //               P(x) = A(I)+x*(B(I)+x*(C(I)+x*D(I)))

  //  Reference: M.J. Maron, 'Numerical Analysis: a Practical Approach',
  //             MacMillan Publ. Co., New York, 1982.


  if(N < 4)
  {
    printf(" *** Error in SPLINE: interpolation cannot be performed with %4d points.\n", N);
    penError(ERR_SPLINE_N_LESS_4);
  }
  int N1 = N-1;
  int N2 = N-2;
  int K;
      //  ****  Auxiliary arrays H(=A) and DELTA(=D).
  for(int I = 0; I < N1; I++)
  {
    A[I] = X[I+1]-X[I];
    if(A[I] < 1.0E-12*( (fabs(X[I]) > fabs(X[I+1])) ? fabs(X[I]) : fabs(X[I+1]) ))
    {
       printf(" *** Error in SPLINE: X values not in increasing order.\n");
       penError(ERR_SPLINE_X_NOT_ORDERED); return;
    }
    D[I] = (Y[I+1]-Y[I])/A[I];
  }
      //  ****  Symmetric coefficient matrix (augmented).
  for(int I = 0; I < N2; I++)
  {
    B[I] = 2.0*(A[I]+A[I+1]);
    K = N1-(I+1)+1;
    D[K-1] = 6.0*(D[K-1]-D[K-2]);
  }
  D[1] = D[1]-A[0]*S1;
  D[N1-1] = D[N1-1]-A[N1-1]*SN;
      //  ****  Gauss solution of the tridiagonal system.
  for(int I = 1; I < N2; I++)
  {
    double R = A[I]/B[I-1];
    B[I] = B[I]-R*A[I];
    D[I+1] = D[I+1]-R*D[I];
  }
      //  ****  The SIGMA coefficients are stored in array D.
  D[N1-1] = D[N1-1]/B[N2-1];
  for(int I = 1; I < N2; I++)
  {
    K = N1-(I+1)+1;
    D[K-1] = (D[K-1]-A[K-1]*D[K+1-1])/B[K-2];
  }
  D[N-1] = SN;
      //  ****  Spline coefficients.
  double SI1 = S1;
  for(int I = 0; I < N1; I++)
  {
    double SI = SI1;
    SI1 = D[I+1];
    double H = A[I];
    double HI = 1.0/H;
    A[I] = (HI/6.0)*(SI*pow(X[I+1],3)-SI1*pow(X[I],3))
      +HI*(Y[I]*X[I+1]-Y[I+1]*X[I])
      +(H/6.0)*(SI1*X[I]-SI*X[I+1]);
    B[I] = (HI/2.0)*(SI1*pow(X[I],2)-SI*pow(X[I+1],2))
      +HI*(Y[I+1]-Y[I])+(H/6.0)*(SI-SI1);
    C[I] = (HI/2.0)*(SI*X[I+1]-SI1*X[I]);
    D[I] = (HI/6.0)*(SI1-SI);
  }
      //  ****  Natural cubic spline for X.GT.X(N).
  double FNP = B[N1-1]+X[N-1]*(2.0*C[N1-1]+X[N-1]*3.0*D[N1-1]);
  A[N-1] = Y[N-1]-X[N-1]*FNP;
  B[N-1] = FNP;
  C[N-1] = 0.0;
  D[N-1] = 0.0;

}

//  *********************************************************************
//                       SUBROUTINE FINDI
//  *********************************************************************
void FINDI(const double* X, const double XC, const int N, unsigned int& I)
{
  //     This subroutine finds the interval (X(I),X(I+1)) that contains the
  //  value XC using binary search.

  //  Input:
  //     X(I) (I=1:N) ... grid points (the X values must be in increasing
  //                      order).
  //     XC ............. point to be located.
  //     N  ............. number of grid points.
  //  Output:
  //     I .............. interval index.

  if(XC > X[N-1])
  {
    I = N-1;
    return;
  }
  if(XC < X[0])
  {
    I = 0;
    return;
  }
  I = 0;
  unsigned int I1 = N-1;
  while(I1-I > 1)
  {
    int IT = (I+I1)/2;
    if(XC > X[IT])
    {
      I = IT;
    }
    else
    {
      I1 = IT;
    }
  }
}

//  *********************************************************************
//                       FUNCTION RLMOM
//  *********************************************************************
double RLMOM(const double* X, double* FCT, double XC, int NPpar, int MOM)
{
  //  Calculation of the integral of (X**MOM)*FCT(X) over the interval from
  //  X(1) to XC, obtained by linear interpolation on a table of FCT.
  //  The independent variable X is assumed to take only positive values.

  //    X ....... array of values of the variable (in increasing order).
  //    FCT ..... corresponding FCT values.
  //    NPpar ...... number of points in the table.
  //    XC ...... upper limit of the integral, X(1).LE.XC.LE.X(NP).
  //    MOM ..... moment order (GE.-1).
  //    RLMOM ... integral of (X**MOM)*FCT(X) over the interval from X(1)
  //              to XC.


  const double EPS = 1.0E-35;

  double RLMOM_RETURN  = 0.0;
  if(MOM < -1){ penError(ERR_RLMOM_0); return RLMOM_RETURN;}
  if(NPpar < 2){ penError(ERR_RLMOM_1); return RLMOM_RETURN;}
  if(X[0] < 0.0){ penError(ERR_RLMOM_2); return RLMOM_RETURN;}
  for(int I = 1; I < NPpar; I++)
  {
    if(X[I] < 0.0){ penError(ERR_RLMOM_3); return RLMOM_RETURN;}
    if(X[I] < X[I-1]){ penError(ERR_RLMOM_4); return RLMOM_RETURN;}
  }

  RLMOM_RETURN  = 0.0;
  if(XC < X[0]){ return RLMOM_RETURN;}
  int IEND = 0;

  double XT;
  if(XC < X[NPpar-1]){ XT = XC;}
  else{ XT = X[NPpar-1];}
  
  double X1, Y1, X2, Y2, XTC, DX, DY, A, B, DS;
  for(int I = 0; I < NPpar-1; I++)
  {
    if(X[I] < EPS){ X1 = EPS;}
    else{ X1 = X[I];}
      
    Y1 = FCT[I];

    if(X[I+1] < EPS){ X2 = EPS;}
    else{ X2 = X[I+1];}

    Y2 = FCT[I+1];
    if(XT < X2)
    {
      XTC = XT;
      IEND = 1;
    }
    else
    {
      XTC = X2;
    }
    DX = X2-X1;
    DY = Y2-Y1;
    if(fabs(DX) > 1.0E-14*fabs(DY))
    {
      B = DY/DX;
      A = Y1-B*X1;
      if(MOM == -1)
      {
        DS = A*log(XTC/X1)+B*(XTC-X1);
      }
      else
      {
        DS = A*(pow(XTC,(MOM+1))-pow(X1,(MOM+1)))/double(MOM+1)+B*(pow(XTC,(MOM+2))-pow(X1,(MOM+2)))/double(MOM+2);
      }
    }
    else
    {
      DS = 0.5*(Y1+Y2)*(XTC-X1)*pow(XTC,MOM);
    }
    RLMOM_RETURN = RLMOM_RETURN+DS;
    if(IEND != 0){ return RLMOM_RETURN;}
  }
  return RLMOM_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE SLAG6
//  *********************************************************************
void SLAG6(double &H, double* Y, double* S, int N)
{
  //     Piecewise six-point Lagrange integration of a uniformly tabulated
  //  function.

  //  Input arguments:
  //     H ............ grid spacing.
  //     Y(I) (1:N) ... array of function values (ordered abscissas).
  //     N ............ number of data points.

  //  Output argument:
  //     S(I) (1:N) ... array of integral values defined as
  //                S(I)=INTEGRAL(Y) from X(1) to X(I)=X(1)+(I-1)*H.


  if(N < 6){ penError(ERR_SLAG6_DP); return;}
  double HR = H/1440.0;
  double Y1 = 0.0;
  double Y2 = Y[0];
  double Y3 = Y[1];
  double Y4 = Y[2];
  S[0] = 0.0;
  S[1] = HR*(475*Y2+1427*Y3-798*Y4+482*Y[3]-173*Y[4]+27*Y[5]);
  S[2] = S[1]+HR*(-27*Y2+637*Y3+1022*Y4-258*Y[3]+77*Y[4]-11*Y[5]);
  for(int I = 3; I < N-2; I++)
  {
    Y1 = Y2;
    Y2 = Y3;
    Y3 = Y4;
    Y4 = Y[I];
    S[I] = S[I-1]+HR*(11*(Y1+Y[I+2])-93*(Y2+Y[I+1])+802*(Y3+Y4));
  }
  double Y5 = Y[N-1-1];
  double Y6 = Y[N-1];
  S[N-1-1] = S[N-2-1]+HR*(-27*Y6+637*Y5+1022*Y4-258*Y3+77*Y2-11*Y1);
  S[N-1] = S[N-1-1]+HR*(475*Y6+1427*Y5-798*Y4+482*Y3-173*Y2+27*Y1);
}

//  *********************************************************************
//                       FUNCTION SUMGA
//  *********************************************************************
double SUMGA(FCT_SUMGA FCT, void* arg, double XL, double XU, double TOL, const int ISGAW)
{
  //  This function calculates the value SUMGA of the integral of the
  //  (external) function FCT over the interval (XL,XU) using the 20-point
  //  Gauss quadrature method with an adaptive-bisection scheme.

  //  TOL is the tolerance, i.e. maximum allowed relative error; it should
  //  not be less than 1.0D-13. A warning message is written in unit 6 when
  //  the required accuracy is not attained. The common block CSUMGA can be
  //  used to transfer the error flag IERGA and the number of calculated
  //  function values to the calling program.

  //                                        Francesc Salvat. 2 April, 2012.

  
  const int NP_S = 10;
  const int NP2_S = 2*NP_S;
  const int NP4 = 4*NP_S;
  const int NOIT = 130;
  const int NOIT5 = NOIT/5;
  const int NCALLT = 100000;
  
  double XM[NP_S], XP[NP_S];
  double S[NOIT], SN[NOIT], XR[NOIT], XRN[NOIT];
  int IERGA = 0;
  
  //  Output error codes:
  //     IERGA = 0, no problem, the calculation has converged.
  //           = 1, too many open subintervals.
  //           = 2, too many function calls.
  //           = 3, subintervals are too narrow.

  //  ****  Gauss 20-point integration formula.
  //  Abscissas.
  double X[NP_S] = {7.6526521133497334E-02,2.2778585114164508E-01,3.7370608871541956E-01,5.1086700195082710E-01,6.3605368072651503E-01,7.4633190646015079E-01,8.3911697182221882E-01,9.1223442825132591E-01,9.6397192727791379E-01,9.9312859918509492E-01};
  //  Weights.
  double W[NP_S] = {1.5275338713072585E-01,1.4917298647260375E-01,1.4209610931838205E-01,1.3168863844917663E-01,1.1819453196151842E-01,1.0193011981724044E-01,8.3276741576704749E-02,6.2672048334109064E-02,4.0601429800386941E-02,1.7614007139152118E-02};

  for(int I = 0; I < NP_S; I++)
  {
    XM[I] = 1.0-X[I];
    XP[I] = 1.0+X[I];
  }
  //  ****  Global and partial tolerances.

  double TOL1; // Global tolerance.

  if(TOL < 1.0E-13){ TOL1 = 1.0E-13;}
  else{ TOL1 = TOL;}

  if( TOL1 > 1.0E-5){ TOL1 = 1.0E-5;}
  
  double TOL2 = TOL1;  // Effective tolerance.
  double TOL3 = 1.0E-13;  // Round-off protection.
  double SUMGA_RETURN = 0.0; //Valor que torna la funcio
  //  ****  Straight integration from XL to XU.
  double H = XU-XL;
  double HH = 0.5*H;
  double X1 = XL;
  double SP = W[0]*(FCT(X1+XM[0]*HH,arg)+FCT(X1+XP[0]*HH,arg));
  for(int J = 1; J < NP_S; J++)
  {
    SP = SP+W[J]*(FCT(X1+XM[J]*HH,arg)+FCT(X1+XP[J]*HH,arg));
  }
  S[0] = SP*HH;
  XR[0] = X1;
  int NCALL = NP2_S;
  int NOI = 1;
  int IDONE = 1;  // To prevent a compilation warning.

    //  ****  Adaptive-bisection scheme.

  bool Eixir = false;
  int NOIP;
  double SUMR;
  while(!Eixir)
  {
    Eixir = true;
    H = HH;  // Subinterval length.
    HH = 0.5*H;
    double AHH = fabs(HH);
    if(TOL2 > 0.01*TOL1){ TOL2 = TOL2*0.5;}
    SUMR = 0.0;
    NOIP = NOI;
    NOI = 0;
    bool Eixir2 = false;
    for(int I = 0; I < NOIP; I++)
    {
      double SI = S[I];  // Bisect the I-th open interval.
      
      X1 = XR[I];
      if(AHH < fabs(X1)*TOL3){ IERGA = 3;}  // The interval is too narrow.
      SP = W[0]*(FCT(X1+XM[0]*HH,arg)+FCT(X1+XP[0]*HH,arg));
      for(int J = 1; J < NP_S; J++)
      {
        SP = SP+W[J]*(FCT(X1+XM[J]*HH,arg)+FCT(X1+XP[J]*HH,arg));
      }
      double S1 = SP*HH;
      
      double X2 = X1+H;
      if(AHH < fabs(X2)*TOL3){ IERGA = 3;}  // The interval is too narrow.
      SP = W[0]*(FCT(X2+XM[0]*HH,arg)+FCT(X2+XP[0]*HH,arg));
      for(int J = 1; J < NP_S; J++)
      {
        SP = SP+W[J]*(FCT(X2+XM[J]*HH,arg)+FCT(X2+XP[J]*HH,arg));
      }
      double S2 = SP*HH;
      
      IDONE = I+1;
      NCALL = NCALL+NP4;
      double S12 = S1+S2;  // Sum of integrals on the two subintervals.
      if(fabs(S12-SI) < ((TOL2*fabs(S12) > 1.0E-35) ? TOL2*fabs(S12) : 1.0E-35))
      {
      //  ****  The integral over the parent interval has converged.
        SUMGA_RETURN = SUMGA_RETURN+S12;
      }
      else
      {
        SUMR = SUMR+S12;
        NOI = NOI+2;
        if(NOI < NOIT)
        {
         //  ****  Store open intervals.
          SN[NOI-2] = S1;
          XRN[NOI-2] = X1;
          SN[NOI-1] = S2;
          XRN[NOI-1] = X2;
        }
        else
        {
        //  ****  Too many open intervals.
          IERGA = 1;
          Eixir2 = true;
          break;
        }
      }
      if(NCALL > NCALLT)
      {
       //  ****  Too many calls to FCT.
        IERGA = 2;
        Eixir2 = true;
        break;
      }
    }
    if(Eixir2){ break;}

  //  ****  Analysis of partial results and error control.

    if(IERGA == 3)  // Intervals are too narrow.
    {
      if(NOI < NOIT5)
      {
        IERGA = 0;  // The result is probably correct.
        SUMGA_RETURN = SUMGA_RETURN+SUMR;
        return SUMGA_RETURN;
      }
      break;
    }

    if(IERGA == 0)
    {
      double Aux_Double = TOL1*fabs(SUMGA_RETURN+SUMR);
      if(Aux_Double < 1.0E-35){ Aux_Double = 1.0E-35;}
      if(fabs(SUMR) < Aux_Double || NOI == 0)
      {
        return SUMGA_RETURN;
      }
      else
      {
        for(int I = 0; I < NOI; I++)
        {
          S[I] = SN[I];
          XR[I] = XRN[I];
        }
        Eixir = false;
        continue;
      }
    }
  }

    //  ****  Warning (low accuracy) message.

  if(IDONE < NOIP)
  {
    for(int I = IDONE; I < NOIP; I++)
    {
      SUMR = SUMR+S[I];
    }
    NOI = NOI+(NOIP-IDONE);
  }
  SUMGA_RETURN = SUMGA_RETURN+SUMR;
  if(ISGAW == 0){ return SUMGA_RETURN;}
  printf("  >>> SUMGA. Gauss adaptive-bisection quadrature.\n");
  printf("  XL =%15.8E, XU =%15.8E, TOL =%8.1E\n", XL, XU, TOL);
  if(fabs(SUMGA_RETURN) > 1.0E-35)
  {
    double RERR = fabs(SUMR)/fabs(SUMGA_RETURN);
    printf("  SUMGA =%22.15E, relative error =%8.1E\n", SUMGA_RETURN, RERR);
  }
  else
  {
    double AERR = fabs(SUMR);
    printf("  SUMGA =%22.15E, absolute error =%8.1E\n", SUMGA_RETURN, AERR);
  }
  printf("  NCALL =%6d, open subintervals =%4d, H =%10.3E\n", NCALL, NOI, HH);
  if(IERGA == 1)
  {
    printf("  IERGA = 1, too many open subintervals.\n");
  }
  else if(IERGA == 2)
  {
    printf("  IERGA = 2, too many function calls.\n");
  }
  else if(IERGA == 3)
  {
    printf("  IERGA = 3, subintervals are too narrow.\n");
  }
  printf("  WARNING: the required accuracy has not been attained.\n");

  return SUMGA_RETURN; //MIRAR
  
}

//  *********************************************************************
//                       SUBROUTINE SINTEG
//  *********************************************************************
void SINTEG(double* X, double* A, double* B, double* C, double* D, double &XL, double &XU, double &SUM, int N)
{

  //  ****  Set integration limits in increasing order.
  double XLL, XUU, SIGN;
  if(XU > XL)
  {
    XLL = XL;
    XUU = XU;
    SIGN = 1.0;
  }
  else
  {
    XLL = XU;
    XUU = XL;
    SIGN = -1.0;
  }
  //  ****  Check integral limits.
  if(XLL < X[0] || XUU > X[N-1])
  {
    printf("\n     Integral limits out of range. Stop.\n");
    penError(ERR_SINTEG_LIMITS); return;
  }
  //  ****  Find involved intervals.
  SUM = 0.0;
  unsigned int IL, IU;
  FINDI(X,XLL,N,IL);
  FINDI(X,XUU,N,IU);

  double X1, X2;
  if(IL == IU)
  {
    //  ****  Only a single interval involved.
    X1 = XLL;
    X2 = XUU;
    SUM = X2*(A[IL]+X2*((B[IL]/2)+X2*((C[IL]/3)+X2*D[IL]/4)))
         -X1*(A[IL]+X1*((B[IL]/2)+X1*((C[IL]/3)+X1*D[IL]/4)));
  }
  else
  {
    //  ****  Contributions from several intervals.
    X1 = XLL;
    X2 = X[IL+1];
    SUM = X2*(A[IL]+X2*((B[IL]/2)+X2*((C[IL]/3)+X2*D[IL]/4)))
         -X1*(A[IL]+X1*((B[IL]/2)+X1*((C[IL]/3)+X1*D[IL]/4)));
    IL = IL+1;
    for(unsigned int I = IL; I <= IU;  I++)
    {
      X1 = X[I];
      X2 = X[I+1];
      if(I == IU){ X2 = XUU;}
      double SUMP = X2*(A[I]+X2*((B[I]/2)+X2*((C[I]/3)+X2*D[I]/4)))
	           -X1*(A[I]+X1*((B[I]/2)+X1*((C[I]/3)+X1*D[I]/4)));
      SUM = SUM+SUMP;
    }
  }
  SUM = SIGN*SUM;     
}

//  *********************************************************************
//                       FUNCTION RMOMX
//  *********************************************************************
double RMOMX(const double* X,
	     const double* PDF,
	     const double XD,
	     const double XU,
	     const int NP,
	     const int MOM)
{
  //  Calculation of momenta of a pdf, PDF(X), obtained from linear log-log
  //  interpolation of the input table. The independent variable X is
  //  assumed to take only positive values.

  //     X ........ array of variable values (in increasing order).
  //     PDF ...... corresponding PDF values (must be non-negative).
  //     NP ....... number of points in the table.
  //     XD, XU ... limits of the integration interval.
  //     MOM ...... moment order.
  //     RMOM = INTEGRAL (X**N)*PDF(X) dX over the interval (XD,XU).

  
  const double EPS = 1.0E-12;
  const double ZERO = 1.0E-35;

  double RMOMX_RETURN=0.0; //Valor que torna la funcio
  
  if(NP < 2){ penError(ERR_RMOMX_1); return RMOMX_RETURN;}
  if(X[0] < 0.0 || PDF[0] < 0.0)
  {
    printf("X(1),PDF(1) = %E %E\n",X[0],PDF[0]);
    penError(ERR_RMOMX_2); return RMOMX_RETURN;
  }
  for(int I = 1; I < NP; I++)
  {
    if(X[I] < 0.0 || PDF[I] < 0.0)
    {
      printf("I,X(I),PDF(I) = %d %E %E\n",I+1,X[I],PDF[I]);
      penError(ERR_RMOMX_3); return RMOMX_RETURN;
    }
    if(X[I] < X[I-1]){ penError(ERR_RMOMX_4); return RMOMX_RETURN;}
  }

  double XLOW = (X[0] > XD ? X[0] : XD);
  if(XLOW < ZERO){ XLOW = ZERO;}
  double XUP = (X[NP-1] < XU ? X[NP-1] : XU);

  if(XLOW > XUP)
  {
    printf("\n WARNING: XLOW is greater than XUP in RMOMX.");
    printf("\n XLOW =%E,   XUP =%E", XLOW, XUP);
    RMOMX_RETURN = 0.0;
    return RMOMX_RETURN;
  }

  int IL = 1;
  int IU = NP-1;
  for(int I = 0; I < NP-1; I++)
  {
    if(X[I] < XLOW){ IL = I+1;}
    if(X[I] < XUP){ IU = I+1;}
  }

  //  ****  A single interval.

  double XIL, XFL, YIL, YFL, X1, X2, DENOM, Y1, Y2, DXL, DYL, DSUM, AP1;
  if(IU == IL)
  {
    XIL = log((X[IL-1] > ZERO ? X[IL-1] : ZERO));
    XFL = log(X[IL+1-1]);
    YIL = log((PDF[IL-1] > ZERO ? PDF[IL-1] : ZERO));
    YFL = log((PDF[IL+1-1] > ZERO ? PDF[IL+1-1] : ZERO));
    X1 = XLOW;
    X2 = XUP;
    DENOM = XFL-XIL;
    if(fabs(DENOM) > ZERO)
    {
      Y1 = exp(YIL+(YFL-YIL)*(log(X1)-XIL)/DENOM)*pow(X1,MOM);
      Y2 = exp(YIL+(YFL-YIL)*(log(X2)-XIL)/DENOM)*pow(X2,MOM);
    }
    else
    {
      Y1 = exp(YIL)*pow(X1,MOM);
      Y2 = exp(YIL)*pow(X2,MOM);
    }
    DXL = log(X2)-log(X1);
    DYL = log((Y2 > ZERO ? Y2 : ZERO))-log((Y1 > ZERO ? Y1 : ZERO));
    if(fabs(DXL) > EPS*fabs(DYL))
    {
      AP1 = 1.0+(DYL/DXL);
      if(fabs(AP1) > EPS)
      {
        DSUM = (Y2*X2-Y1*X1)/AP1;
      }
      else
      {
        DSUM = Y1*X1*DXL;
      }
    }
    else
    {
      DSUM = 0.5*(Y1+Y2)*(X2-X1);
    }
    RMOMX_RETURN = DSUM;
    return RMOMX_RETURN;
  }

      //  ****  Multiple intervals.

  XIL = log((X[IL-1] > ZERO ? X[IL-1] : ZERO));
  XFL = log(X[IL+1-1]);
  YIL = log((PDF[IL-1] > ZERO ? PDF[IL-1] : ZERO));
  YFL = log((PDF[IL+1-1] > ZERO ? PDF[IL+1-1] : ZERO));
  X1 = XLOW;
  DENOM = XFL-XIL;
  if(fabs(DENOM) > ZERO)
  {
    Y1 = exp(YIL+(YFL-YIL)*(log(X1)-XIL)/DENOM)*pow(X1,MOM);
  }
  else
  {
    Y1 = exp(YIL)*pow(X1,MOM);
  }
  X2 = X[IL+1-1];
  Y2 = (PDF[IL+1-1] > ZERO ? PDF[IL+1-1] : ZERO)*pow(X2,MOM);
  DXL = log(X2)-log(X1);
  DYL = log((Y2 > ZERO ? Y2 : ZERO))-log((Y1 > ZERO ? Y1 : ZERO));
  if(fabs(DXL) > EPS*fabs(DYL))
  {
    AP1 = 1.0+(DYL/DXL);
    if(fabs(AP1) > EPS)
    {
      DSUM = (Y2*X2-Y1*X1)/AP1;
    }
    else
    {
      DSUM = Y1*X1*DXL;
    }
  }
  else
  {
    DSUM = 0.5*(Y1+Y2)*(X2-X1);
  }
  RMOMX_RETURN = DSUM;
      
  if(IU > IL+1)
  {
    for(int I = IL; I < IU-1; I++)
    {
      X1 = X[I];
      Y1 = (PDF[I] > ZERO ? PDF[I] : ZERO)*pow(X1,MOM);
      X2 = X[I+1];
      Y2 = (PDF[I+1] > ZERO ? PDF[I+1] : ZERO)*pow(X2,MOM);
      DXL = log(X2)-log(X1);
      DYL = log((Y2 > ZERO ? Y2 : ZERO))-log((Y1 > ZERO ? Y1 : ZERO));
      if(fabs(DXL) > EPS*fabs(DYL))
      {
        AP1 = 1.0+(DYL/DXL);
        if(fabs(AP1) > EPS)
        {
          DSUM = (Y2*X2-Y1*X1)/AP1;
        }
        else
        {
          DSUM = Y1*X1*DXL;
        }
      }
      else
      {
        DSUM = 0.5*(Y1+Y2)*(X2-X1);
      }
      RMOMX_RETURN = RMOMX_RETURN+DSUM;
    }
  }

  X1 = X[IU-1];
  Y1 = (PDF[IU-1] > ZERO ? PDF[IU-1] : ZERO)*pow(X1,MOM);
  XIL = log(X[IU-1]);
  XFL = log(X[IU+1-1]);
  YIL = log((PDF[IU-1] > ZERO ? PDF[IU-1] : ZERO));
  YFL = log((PDF[IU+1-1] > ZERO ? PDF[IU+1-1] : ZERO));
  X2 = XUP;
  DENOM = XFL-XIL;
  if(fabs(DENOM) > ZERO)
  {
    Y2 = exp(YIL+(YFL-YIL)*(log(X2)-XIL)/DENOM)*pow(X2,MOM);
  }
  else
  {
    Y2 = exp(YIL)*pow(X2,MOM);
  }
  DXL = log(X2)-log(X1);
  DYL = log((Y2 > ZERO ? Y2 : ZERO))-log((Y1 > ZERO ? Y1 : ZERO));
  if(fabs(DXL) > EPS*fabs(DYL))
  {
    AP1 = 1.0+(DYL/DXL);
    if(fabs(AP1) > EPS)
    {
      DSUM = (Y2*X2-Y1*X1)/AP1;
    }
    else
    {
      DSUM = Y1*X1*DXL;
    }
  }
  else
  {
    DSUM = 0.5*(Y1+Y2)*(X2-X1);
  }
  RMOMX_RETURN = RMOMX_RETURN+DSUM;

  return RMOMX_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE RLPAC
//  *********************************************************************
void RLPAC(const double* X, double* PDF, double* PAC, int NPpar)
{
  //  Cumulative distribution function of PDF(X)/X, obtained from linear
  //  interpolation on a table of PDF values.
  //  The independent variable X is assumed to take only positive values.

  //    X ..... array of values of the variable (in increasing order).
  //    PDF ... corresponding PDF values.
  //    PAC ... cumulative probability function.
  //    NP .... number of points in the table.


  const double EPS = 1.0E-35;
  PAC[0] = 0.0;
  double X1, Y1, X2, Y2, DX, DY, B, A, DS;
  for(int I = 1; I < NPpar; I++)
  {
    if(X[I-1] < EPS){ X1 = EPS;}
    else{ X1 = X[I-1];}
    Y1 = PDF[I-1];

    if(X[I] < EPS){ X2 = EPS;}
    else{ X2 = X[I];}
      
    Y2 = PDF[I];
    DX = X2-X1;
    DY = Y2-Y1;
    B = DY/DX;
    A = Y1-B*X1;
    DS = A*log(X2/X1)+B*(X2-X1);
    PAC[I] = PAC[I-1]+DS;
  }
     
}

//  *********************************************************************
//                       SUBROUTINE MERGE2
//  *********************************************************************
void MERGE2(double* X1, double* Y1, double* X2, double* Y2, double* XM, double* YM, int &N1, int &N2, int &N)
{
  //  This subroutine merges two tables (X1,Y1), (X2,Y2) of two functions
  //  to produce a table (XM,YM) of the sum of these functions, with abs-
  //  cissas in increasing order. The abscissas and function values are
  //  assumed to be positive. N1, N2 and N are the numbers of grid points
  //  in the input and merged tables. A discontinuity in the function is
  //  described by giving twice the abscissa. Log-log linear interpolation
  //  is used to interpolate the input tables.
  

  const double EPS = 1.0E-10;
  const int NP_S = 12000;
  const int NP2_S = NP_S + NP_S;
  double X[NP2_S], Y12[NP2_S];

  if(N1 > NP_S || N2 > NP_S)
  {
    printf("NP =%7d\n", N1 > N2 ? N1 : N2);
    penError(ERR_MERGE2_NP); return;
  }

  SORT2(X1,Y1,N1);
  if(penGetError() != PEN_SUCCESS){ return;}
  SORT2(X2,Y2,N2);
  if(penGetError() != PEN_SUCCESS){ return;}

  for(int I1 = 0; I1 < N1; I1++)
  {
    X[I1] = X1[I1];
  }
  double XC, TST1, TST2, TST3, TST4, TST12, TST34, TST;
  N = N1;
  for(int I2 = 0; I2 < N2; I2++)
  {
    unsigned int I1; 
    XC = X2[I2];
    FINDI(X1,XC,N1,I1);
    if(I1 == unsigned(N1-1)){I1--;}
    TST1 = fabs(XC-X1[I1]);
    TST2 = fabs(XC-X1[I1+1]);
    TST12 = (TST1 < TST2 ? TST1 : TST2);
    if(I2 > 0)
    {
      TST3 = fabs(XC-X2[I2-1]);
    }
    else
    {
      TST3 = 1.0;
    }
    if(I2+1 < N2)
    {
      TST4 = fabs(XC-X2[I2+1]);
    }
    else
    {
      TST4 = 1.0;
    }
    TST34 = (TST3 < TST4 ? TST3 : TST4);
    TST = EPS*XC;
    if(TST34 > TST)
    {
      if(TST12 > TST)
      {
        N = N+1;
        X[N-1] = XC;
      }
    }
    else
    {
      N = N+1;
      X[N-1] = XC;
    }
  }

  //  ****  Sort and clean the merged grid.

  bool Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    for(int I = 0; I < N-1; I++)
    {
      int IMIN = I+1;
      double XMIN = X[I];
      for(int J = I+1; J < N; J++)
      {
        if(X[J] < XMIN)
        {
          IMIN = J+1;
          XMIN = X[J];
        }
      }
      double SAVE = X[I];
      X[I] = X[IMIN-1];
      X[IMIN-1] = SAVE;
    }

    for(int I = 0; I < N-2; I++)
    {
      if(X[I] > X[I+2]*(1.0E0-EPS))
      {
        X[I+1] = X[N-1];
        N = N-1;
        Eixir = false;
        break;
      }
    }
  }

  for(int I = 0; I < N; I++)
  {
    XC = X[I];
    if(I+1 < N)
    {
      if(X[I] > X[I+1]*(1.0-EPS)){ XC = X[I]*(1.0-EPS);}
    }
    if(I+1 > 1)
    {
      if(X[I] < X[I-1]*(1.0+EPS)){ XC = X[I]*(1.0+EPS);}
    }
    unsigned int J;
    double YI1, YI2;
    FINDI(X1,XC,N1,J);
    if(J == unsigned(N1-1)){ J--;}
    if(X1[J+1] > X1[J]+EPS)
    {
      YI1 = exp(log(Y1[J])+log(XC/X1[J])*log(Y1[J+1]/Y1[J])/log(X1[J+1]/X1[J]));
    }
    else
    {
      YI1 = Y1[J];
    }
    FINDI(X2,XC,N2,J);
    if(J == unsigned(N2-1)){J--;}
    if(X2[J+1] > X2[J]+EPS)
    {
      YI2 = exp(log(Y2[J])+log(XC/X2[J])*log(Y2[J+1]/Y2[J])/log(X2[J+1]/X2[J]));
    }
    else
    {
      YI2 = Y2[J];
    }
    Y12[I] = YI1+YI2;
    if(Y12[I] < 1.0E-75){ Y12[I] = 1.0E-75;}
  }

  if(N > NP_S)
  {
    printf("NP = %7d\n", N);
    penError(ERR_MERGE2_NP); return;
  }
  for(int I = 0; I < N; I++)
  {
    XM[I] = X[I];
    YM[I] = Y12[I];
  }     
}

//  *********************************************************************
//                       SUBROUTINE SORT2
//  *********************************************************************
void SORT2(double* X, double* Y, int &N)
{
  //  This subroutine sorts a table (X,Y) of a function with n data points.
  //  A discontinuity of the function is described by giving twice the abs-
  //  cissa. It is assumed that the function is strictly positive (negative
  //  values of Y are set to zero).


  const int NP_S=12000;
  int IORDER[NP_S];
  
  if(N > NP_S)
  {
    printf("NP =%7d\n", N);
    penError(ERR_SORT2_NP); return;
  }

  if(N == 1){ return;}
  for(int I = 0; I < N; I++)
  {
    IORDER[I] = I+1;
    if(Y[I] < 1.0E-75){ Y[I] = 1.0E-75;}
  }

  int IMIN;
  double XMIN;
  for(int I = 0; I < N-1; I++)
  {
    IMIN = I+1;
    XMIN = X[I];
    for(int J = I+1; J < N; J++)
    {
      if(X[J] < XMIN)
      {
        IMIN = J+1;
        XMIN = X[J];
      }
    }
    double SAVE = X[I];
    X[I] = X[IMIN-1];
    X[IMIN-1] = SAVE;
    SAVE = Y[I];
    Y[I] = Y[IMIN-1];
    Y[IMIN-1] = SAVE;
    int ISAVE = IORDER[I];
    IORDER[I] = IORDER[IMIN-1];
    IORDER[IMIN-1] = ISAVE;
    if(I+1 == 1){ continue;}
    if(IORDER[I] < IORDER[I-1] && fabs(X[I]-X[I-1]) < 1.0E-15)
    {
      SAVE = X[I-1];
      X[I-1] = X[I];
      X[I] = SAVE;
      SAVE = Y[I-1];
      Y[I-1] = Y[I];
      Y[I] = SAVE;
      ISAVE = IORDER[I-1];
      IORDER[I-1] = IORDER[I];
      IORDER[I] = ISAVE;
    }
  }
  int I = N;
  if(IORDER[I-1] < IORDER[I-1-1] && fabs(X[I-1]-X[I-1-1]) < 1.0E-15)
  {
    double SAVE = X[I-1-1];
    X[I-1-1] = X[I-1];
    X[I-1] = SAVE;
    SAVE = Y[I-1-1];
    Y[I-1-1] = Y[I-1];
    Y[I-1] = SAVE;
  }
}
