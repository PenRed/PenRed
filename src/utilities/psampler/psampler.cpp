#include <stdio.h>
#include <math.h>

#include "PenRed.hh"

struct C_PARIN_DCSIN;
struct DOICSC;
struct CADATA_PSAMPLER;
struct CEBRA0;
struct CPAN01;
struct CVERIF;
struct CGCO01;
struct CGPHVS;
struct CGPP02;
struct CGRA00;

void EELaV(const pen_context& context, const pen_material& mat, const double E, const int KPAR, const CEEL00& eel00, const int NSAMP, pen_rand& randoms);
void EELdV(const pen_context& context, const pen_material& mat, const double E, const int KPAR, const CEEL00& eel00, const int NSAMP, CDCSEP& dcsep, pen_rand& randoms);
void EINaV(const pen_context& context, const pen_material& mat, const double E, int NSAMP, int IDIST, CEIN00& ein00, const CESIN& esin, pen_rand& randoms);
void PINaV(const pen_context& context, const pen_material& mat, const double E, CPSIN& psin, int NSAMP, int IDIST, pen_rand& randoms);
void ESIb(const pen_material& mat, const int KE, const double XEL, const double* DLEMP, const double DLFC, int &IZZ, int &ISH, pen_rand& penRand);
void ESIaV(const pen_context& context, const pen_material& mat, const double E, const CESI0& esio, int NSAMP, pen_rand& randoms);
void PSIb(const pen_material& mat, const int KE, const double XEL, const double* DLEMP, const double DLFC, int &IZZ, int &ISH, pen_rand& penRand);
void PSIaV(const pen_context& context, const pen_material& mat, const double E, const CPSI0& psio, int NSAMP, pen_rand& randoms);
void EBRaV(const pen_context& context, const pen_material& mat, const double E, int NSAMP, pen_rand& randoms);
void EBRaAV(const pen_material& mat, const double E, double DE,int NSAMP, pen_rand& randoms);
void PANaV(const pen_material &mat, const double E, int NSAMP, pen_rand& randoms);
void GRAaV(const pen_context& context, pen_material& mat, const double E, int NSAMP, pen_rand& randoms);
void GCOaV(const pen_material& mat, const double E, int NSAMP, pen_rand& randoms);
void GPHaAV(const double E, int NSAMP, pen_rand& randoms);
void GPHaSV(const pen_context& context, const pen_material& mat, const double E, int NSAMP, pen_rand& randoms);
void GPPaV(const pen_context& context, const pen_material& mat, const double E, int NSAMP, pen_rand& randoms);
void RELAXV(const pen_context& context, const pen_material& mat, pen_elementDataBase& elements, const double E, int IZ, int IS, int NSAMP, pen_rand& randoms);
void DOICSS(double WCCM,double UI,double WR,double QI,double WM,int IELEC,double E,double RMU,double& DCSD,double& DCSC, DOICSC& doicsc);
double RIDCS(const pen_material& mat, const int KPAR, double XX, int IE);
double EDCSIW(double W, C_PARIN_DCSIN& parin_dcsin);
double EDCSIA(double RMU, void* arg);
double DODCSD(double RMU, DOICSC& doicsc);
double DODCSC(double WCCM,double RMU, DOICSC& doicsc);
double EBRADX(double T, void* arg);
double EBRDA(const pen_material& mat, double EE, double DEE, double CDT);
double PANDX(double CHI, void* arg);
void GCOTS3(const pen_material& mat, CVERIF& verif, const double E,double& CSIN0,double& CSIN1,double& CSIN2);
double GCOSDT(const pen_material& mat, const double E,double EP);
double GCODDS(double CDT, void* arg);
double SAUTDX(double T, void* arg);
void GPPTX(const pen_elementDataBase& elements, CGPP02& gpp02, const double Z, const double E, double& TX0, double& EPS1AV, double& EPS2AV);
double GPPaD(double EPS, void* arg);

double GRAaD(double CDT, void* arg);
double GRAaF2( double Q2, CGRA00& gra00);
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C    PPPPP    SPSSS    AA    M    M  PPPPP   L       EEEEEE  RRRRR     C
//C    P    P  S        A  A   MM  MM  P    P  L       E       R    R    C
//C    P    P  S       A    A  M MM M  P    P  L       E       R    R    C
//C    PPPPP    SSSS   AAAAAA  M    M  PPPPP   L       EEEE    RRRRR     C
//C    P            S  A    A  M    M  P       L       E       R   R     C
//C    P       SSSSS   A    A  M    M  P       LLLLLL  EEEEEE  R    R    C
//C                                                                      C
//C                                                   (version 2018).    C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C                                                                      C
//C  PENELOPE/PENGEOM (version 2018)                                     C
//C  Copyright (c) 2001-2018                                             C
//C  Universitat de Barcelona                                            C
//C                                                                      C
//C  Permission to use, copy, modify, distribute and sell this software  C
//C  and its documentation for any purpose is hereby granted without     C
//C  fee, provided that the above copyright notice appears in all        C
//C  copies and that both that copyright notice and this permission      C
//C  notice appear in all supporting documentation. The Universitat de   C
//C  Barcelona makes no representations about the suitability of this    C
//C  software for any purpose. It is provided "as is" without express    C
//C  or implied warranty.                                                C
//C                                                                      C
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C
//C     This program verifies the correctness of the sampling algorithms
//C  used in PENELOPE for the different interaction mechanisms and models.
//C  The material data file of the material where radiation propagates
//C  must have been generated previously by running one of the programs
//C  MATERIAL.F or TABLES.F.
//C
//C     The user selects an interaction mechanism and the energy of the
//C  corresponding particle. The calculation may take up to several
//C  minutes, depending on the interaction and the required number of
//C  sampled events. The program generates a file named 'calc.dat' with
//C  the theoretical distribution, and a file 'simul.dat' that contains
//C  the results of the random sampling. During the simulation, the ratios
//C  of the simulated and calculated first and second moments of the
//C  distribution are displayed on the screen, these values should be
//C  close to unity.
//C
//C  Different kinds of interactions are identified by the integer label
//C  ICOL, as indicated in the following table:
//C
//C     +----+-----------------+-----------------+-----------------+
//C     |ICOL|electron (KPAR=1)|photon (KPAR=2)  |positron (KPAR=3)|
//C     +----+-----------------+-----------------+-----------------+
//C     | 1  |hinge            |Rayleigh         |hinge            |
//C     +----+-----------------+-----------------+-----------------+
//C     | 2  |elastic          |Compton          |elastic          |
//C     +----+-----------------+-----------------+-----------------+
//C     | 3  |inelastic        |photoabsorption  |inelastic        |
//C     +----+-----------------+-----------------+-----------------+
//C     | 4  |bremsstrahlung   |pair production  |bremsstrahlung   |
//C     +----+-----------------+-----------------+-----------------+
//C     | 5  |inner-shell ion. |not defined      |inner-shell ion. |
//C     +----+-----------------+-----------------+-----------------+
//C     | 6  |not defined      |not defined      |annihilation     |
//C     +----+-----------------+-----------------+-----------------+
//C     | 7  |delta scattering |delta scattering |delta scattering |
//C     +----+-----------------+-----------------+-----------------+
//C     | 8  |not defined      |not defined      |not defined      |
//C     +----+-----------------+-----------------+-----------------+
//C
//C
//C  *********************************************************************
//C                       MAIN PROGRAM
//C  *********************************************************************

//*******************************************************************
//*  Structs for psampler                                           *
//*******************************************************************

struct C_PARIN_DCSIN
{
  double ZTOTM,EXPOTM;
  double FM[constants::NO],UIM[constants::NO],WRIM[constants::NO],QRIM[constants::NO],WTH[constants::NO],FCIN[constants::NO];
  int NOSCM;

  double EM,DELTA,WCCM;
  int IELEC,MOM;  
};

struct DOICSC
{
  double EI,EE,BETA2,RMU1,RMU2,RMU3,AMOL,BHA1,BHA2,BHA3,BHA4,CDIST,C1,C2,CCLO;
  int    JELEC;
};

struct CADATA_PSAMPLER
{
  double ATW[99],EPX[99],RSCR[99],ETA[99],EB[99][30],ALW[99][30],CP0[99][30];
  int IFI[99][30],IKS[99][30],NSHT[99],LASYMB[99];
};

struct CEBRA0
{
  double EE, DEE;
  int MOM;
  const pen_material* pmat = 0;
};

struct CPAN01
{
  double GAM;
  int MOM;
};

struct CVERIF
{
  static const int NVF=12000;
  double XC1[NVF],XC2[NVF],XC3[NVF],XC4[NVF];
  int KC1[NVF],KC2[NVF],KC3[NVF],KC4[NVF],NX1,NX2,NX3,NX4,NCHAN,NK2,NK3,NK4;
};

struct CGCO01
{
  double EE,EPP;
  const pen_material* pmat;
};

struct CGPHVS
{
  double GAM,BETA;
  int MOM;
};

struct CGPP02
{
  double GAM,FACT,G0,BCBT;
  int MOM;
};

struct CGRA00
{
  double FACTE, Q2MAX;
  int MOM;
  const pen_material* pmat;
};

int main()
{
  //using namespace TRACK_mod;
  //using namespace PENELOPE_mod;
  //using namespace PENERROR_mod;

  //using namespace CEGRID;
  //using namespace COMPOS;
  //using namespace CEIMFP;
  //using namespace CPIMFP;

  //using namespace RSEED;

  //Declare random generator
  pen_rand randoms;
  
  //Create elements data base
  pen_elementDataBase elementsDB;
  
  //Create a context
  pen_context context(elementsDB);

  //Add material to context
  int errmat = context.setMats<pen_material>(1);
  if(errmat != 0)
    {
      printf("Error at context material creation: %d.\n",errmat);
      return -1;
    }

  //Get material
  pen_material& mat = context.getBaseMaterial(0);
  
  char PMFILE[21];
  //char PMFILE[constants::MAXMAT][21];
  std::string PMFILEstr[constants::MAXMAT];

  randoms.rand0(23);

  printf("  \n");
  printf("Enter the name of the material data file...\n");
  scanf("%s",PMFILE);
  printf(" Material data file: %s\n",PMFILE);

  double EMAX=1.0E9;

  printf("  \n");
  printf("Enter C1, C2, WCC, WCR (in eV)...\n");
  double C11, C21, WCC1, WCR1;
  scanf("%lf %lf %lf %lf",&C11, &C21, &WCC1, &WCR1);

  mat.C1=C11;
  mat.C2=C21;
  mat.WCC=WCC1;
  mat.WCR=WCR1;

  if(mat.WCR < 10.0E0){ mat.WCR=-10.0E0;}
  printf("      C1 =%11.4E       C2 =%11.4E\n",mat.C1,mat.C2);
  printf("     WCC =%11.4E eV,   WCR =%11.4E eV\n",mat.WCC,(mat.WCR > 10.0E0 ? mat.WCR : 10.0E0));

  mat.EABS[PEN_ELECTRON] = (mat.WCC > 50.0E0 ? mat.WCC : 50.0E0);
  mat.EABS[PEN_PHOTON]    = (mat.WCR > 50.0E0 ? mat.WCR : 50.0E0);
  mat.EABS[PEN_POSITRON] = (mat.WCC > 50.0E0 ? mat.WCC : 50.0E0);

  //Get minimum energy
  double EMIN = mat.EABS[PEN_ELECTRON] < mat.EABS[PEN_PHOTON] ? mat.EABS[PEN_ELECTRON] : mat.EABS[PEN_PHOTON];
  if(EMIN > mat.EABS[PEN_POSITRON]){ EMIN = mat.EABS[PEN_POSITRON];}

  context.grid.init(EMIN,EMAX);  // Defines energy grid interval.
  elementsDB.GPHa0();  // Initialises photoelectric routines.
  elementsDB.RELAX0();  // Initialises atomic relaxation routines.
  RNDG30(context.rndg3);  // Initialises the Gaussian sampling routine.
    
  printf("  \n");
  printf("Processing material data. Please, wait...\n");
  
  FILE* IWR = fopen("tables.dat","w");
  FILE* IWRR = fopen("tables2.dat","w");
  int INFO=4;

  //fprintf(IWR, "\n Material data file: %-20s\n", PMFILE);

  //Open material data file
  FILE* IRD = nullptr;
  IRD = fopen(PMFILE,"r");
  if(IRD == nullptr)
    {
      fprintf(IWR, "Error: Material file %s could not be opened\n", PMFILE);
      return 1;
    }
  initStructs initStore;
  mat.load(IRD,IWRR,initStore,context.elements,context.grid,INFO);
  PMFILEstr[0].append(PMFILE);
  context.init(EMAX,IWR,INFO,PMFILEstr);
  
  if(penGetError() != PEN_SUCCESS){
    printf("%s\n",penErrorString(penGetError()));
    return -1;
  }
  fclose(IRD);
  fclose(IWR);
  fclose(IWRR);

//
//  ************  Test of the simulation routines.
//

  
  //Declare particle Energy
  double E;

  
  bool Eixir = false;
  while(!Eixir)
  {
//    Eixir = true;

    printf("  \n");
    printf("\n What mechanism do you wish to test?                     ");
    printf("\n Options:  1=EEL,  2=PEL,  3=EIN,  4=PIN,  5=ESI,  6=PSI,");
    printf("\n           7=BRE,  8=BRA,  9=PAN, 10=GRA, 11=GCO, 12=PHA,");
    printf("\n          13=PHS, 14=GPP, 15=RLX, ^C=exit                \n");
    int LMECH;
    scanf("%d",&LMECH);

    if(LMECH==1)
    {
//
//  ****  Elastic electron scattering.
//
      printf("  \n");
      printf("Testing elastic electron scattering.\n");
      printf("Enter electron kinetic energy in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      if(E < 1.0E8)
      {
        printf("  \n");
        printf("Select the DCS model...\n");
        printf("    (1 = Modified Wentzel, or 2 = ELSEPA)\n");
        int IMOD;
        scanf("%d",&IMOD);
        if(IMOD == 1)
        {
	  EELaV(context, mat, E, PEN_ELECTRON, *(initStore.pceel00), NSAMP, randoms);
        }
        else
        {
          EELdV(context, mat, E, PEN_ELECTRON, *(initStore.pceel00), NSAMP, *(initStore.pcdcsep), randoms);
        }
      }
      else
      {
        EELaV(context, mat, E, PEN_ELECTRON, *(initStore.pceel00), NSAMP, randoms);
      }      
    }
    else if(LMECH==2)
    {
//
//  ****  Elastic positron scattering.
//
      printf("  \n");
      printf("Testing elastic positron scattering.\n");
      printf("Enter positron kinetic energy in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      if(E < 1.0E8)
      {
        printf("  \n");
        printf("Select the DCS model...\n");
        printf("    (1 = Modified Wentzel, or 2 = ELSEPA)\n");
        int IMOD;
        scanf("%d",&IMOD);
        if(IMOD == 1)
        {
          EELaV(context, mat, E, PEN_POSITRON, *(initStore.pceel00), NSAMP, randoms);
        }
        else
        {
          EELdV(context, mat, E, PEN_POSITRON, *(initStore.pceel00), NSAMP, *(initStore.pcdcsep), randoms);
        }
      }
      else
      {
        EELaV(context, mat, E, PEN_POSITRON, *(initStore.pceel00), NSAMP, randoms);
      }      
    }
    else if(LMECH==3)
    {
//
//  ****  Inelastic scattering of electrons.
//
      printf("  \n");
      printf("Testing inelastic electron scattering.\n");
      printf("Enter electron kinetic energy in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
//
      printf("  \n");
      printf("Select the distribution to be simulated...\n");
      printf("    1 = angular distribution, or \n");
      printf("    2 = energy distribution \n");
      int IDIST;
      scanf("%d",&IDIST);
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      EINaV(context,mat,E,NSAMP,IDIST,*(initStore.pcein00),*(initStore.pcesin), randoms);
    }
    else if(LMECH==4)
    {
//
//  ****  Inelastic scattering of positrons.
//
      printf("  \n");
      printf("Testing inelastic positron scattering.\n");
      printf("Enter positron kinetic energy in eV...\n");
      scanf("%lf",&E);
//
      printf("  \n");
      printf("Select the distribution to be simulated...\n");
      printf("    1 = angular distribution, or \n");
      printf("    2 = energy distribution \n");
      int IDIST;
      scanf("%d",&IDIST);
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      PINaV(context,mat,E,*(initStore.pcpsin),NSAMP,IDIST,randoms);
    }
    else if(LMECH==5)
    {
//
//  ****  Inner-shell ionization by electron impact.
//
      printf("  \n");
      printf("Testing the distribution of ionized shells.\n");
      printf("Enter kinetic energy of the projectile in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      ESIaV(context,mat,E,*(initStore.pcesi0),NSAMP,randoms);
    }
    else if(LMECH==6)
    {
//
//  ****  Inner-shell ionization by positron impact.
//
      printf("  \n");
      printf("Testing the distribution of ionized shells.\n");
      printf("Enter kinetic energy of the projectile in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      PSIaV(context,mat,E,*(initStore.pcpsi0),NSAMP,randoms);
    }
    else if(LMECH==7)
    {
//
//  ****  Bresstrahlung energy-loss spectrum.
//
      printf("  \n");
      printf("Testing bremss energy-loss spectrum.\n");
      printf("Enter electron kinetic energy in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      EBRaV(context, mat, E, NSAMP,randoms);
    }
    else if(LMECH==8)
    {
//
//  ****  Bresstrahlung angular distribution.
//
      printf("  \n");
      printf("Testing bremss angular distribution.\n");
      printf("Enter electron kinetic energy in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Enter photon energy in eV...\n");
      double DE;
      scanf("%lf",&DE);
      if(DE > E-1.0E0)
      {
        printf(" DE must be less than E.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      EBRaAV(mat,E,DE,NSAMP,randoms);
    }
    else if(LMECH==9)
    {
//
//  ****  Positron annihilation.
//
      printf("  \n");
      printf("Testing positron annihilation.\n");
      printf("Enter positron kinetic energy in eV...\n");
      scanf("%lf",&E);
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      PANaV(mat, E, NSAMP,randoms);
    }
    else if(LMECH==10)
    {
//
//  ****  Rayleigh scattering of photons.
//
      printf("  \n");
      printf("Testing Rayleigh scattering of photons.\n");
      printf("Enter photon energy in eV...\n");
      scanf("%lf",&E);
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      GRAaV(context,mat,E,NSAMP,randoms);
    }
    else if(LMECH==11)
    {
//
//  ****  Compton scattering of photons.
//
      printf("  \n");
      printf("Testing Compton scattering of photons.\n");
      printf("Enter photon energy in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      printf("This calculation may take a few minutes...\n");
      GCOaV(mat, E, NSAMP,randoms);
    }
    else if(LMECH==12)
    {
//
//  ****  Photoelectron angular distribution.
//
      printf("  \n");
      printf("Testing photoelectron angular distribution.\n");
      printf("Enter photon energy in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      GPHaAV(E,NSAMP,randoms);
    }
    else if(LMECH==13)
    {
//
//  ****  Selection of the active shell in photoelectric absorption.
//
      printf("  \n");
      printf("Testing the active shell in photoabsorption.\n");
      printf("Enter photon energy in eV...\n");
      scanf("%lf",&E);
      if(E < context.grid.EL || E > context.grid.EU)
      {
        printf(" Energy outside the allowed range.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      GPHaSV(context, mat, E, NSAMP,randoms);
    }
    else if(LMECH==14)
    {
//
//  ****  Electron-positron pair production.
//
      printf("  \n");
      printf("Testing electron-positron pair production.\n");
      printf("Enter photon energy in eV...\n");
      scanf("%lf",&E);
      if(E < 1.022E6)
      {
        printf(" Energy must be larger than 1.022 MeV.\n");
        Eixir=false; continue;
      }
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      GPPaV(context, mat, E, NSAMP,randoms);
    }
    else if(LMECH==15)
    {
//
//  ****  Atomic relaxation.
//
      int IZZ,IS0;
      printf("  \n");
      printf("Testing atomic relaxation.\n");
      if(mat.NELEM>1)
      {
        printf("Loaded elements:");
        for(int IEL=1; IEL<=mat.NELEM; IEL++)
        {
          printf("%3d",mat.IZ[IEL-1]);
        }
        printf("\nEnter the atomic number of the element to be simulated...\n");
        scanf("%d",&IZZ);
      }
      else
      {
        printf("Loaded element:%3d\n",mat.IZ[1-1]);
        IZZ=mat.IZ[1-1];
      }
      printf("Shell where the initial vacancy is?\n");
      scanf("%d",&IS0);
      printf("Number of values to be sampled?\n");
      double SAMP;
      scanf("%lf",&SAMP);
      if(SAMP > 1.0E9 || SAMP < 1.0)
      {
        printf(" Try again...\n");
        Eixir=false; continue;
      }
      int NSAMP=(int)SAMP;
      RELAXV(context, mat, context.elements, E, IZZ, IS0, NSAMP, randoms);
    }
    else
    {   
    }
  }
}

//  *********************************************************************
//                       SUBROUTINE EELaV
//  *********************************************************************
void EELaV(const pen_context& context, const pen_material& mat, const double E, const int KPAR, const CEEL00& eel00, const int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm for elastic scatter-
//  ing of electrons and positrons. To run this subroutine, load a single
//  material file (since it uses tables that are overwritten when a new
//  material is loaded).  Modified Wentzel (MW) model.
//
  //using namespace TRACK_mod;
  //using namespace PENELOPE_mod;

//  ****  Composition data.
  //using namespace COMPOS;
//  ****  Energy grid and interpolation constants for the current energy.
  //using namespace CEGRID;
//  ****  Electron simulation tables.
  //using namespace CEIMFP;
//  ****  Positron simulation tables.
  //using namespace CPIMFP;
//  ****  Total and transport cross sections.
  //using namespace CEEL00;

//  ****  Auxiliary arrays for random sampling verification.
  double XC1[1000];
  int    KC1[1000];
  double TH1, TH2, SH0, AA, BB, RNDC, XS0; 

  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);

//  ****  Moments of the 'exact' distribution.
//
  if(KPAR == PEN_ELECTRON)
    {
      TH1=exp(log(eel00.XE1[KE])+(log(eel00.XE1[KE+1])-log(eel00.XE1[KE]))*XEK);
      if(eel00.T1E0[KE] > 1.0E-35)
	{ 
	  TH1=TH1-exp(log(eel00.T1E0[KE])+(log(eel00.T1E0[KE+1])-log(eel00.T1E0[KE]))*XEK);
	}  
      TH2=exp(log(eel00.XE2[KE])+(log(eel00.XE2[KE+1])-log(eel00.XE2[KE]))*XEK);
      if(eel00.T2E0[KE] > 1.0E-35)
	{ 
	  TH2=TH2-exp(log(eel00.T2E0[KE])+(log(eel00.T2E0[KE+1])-log(eel00.T2E0[KE]))*XEK);
	}
      SH0=exp(mat.SEHEL[KE]+(mat.SEHEL[KE+1]-mat.SEHEL[KE])*XEK)/mat.VMOL;
      AA=exp(mat.AE[KE]+(mat.AE[KE+1]-mat.AE[KE])*XEK);
      BB=mat.BE[KE]+(mat.BE[KE+1]-mat.BE[KE])*XEK;
      RNDC=mat.RNDCE[KE]+(mat.RNDCE[KE+1]-mat.RNDCE[KE])*XEK;
      printf(" # Electron energy = %13.5E eV\n",E);
      XS0=exp(log(eel00.XE0[KE])+(log(eel00.XE0[KE+1])-log(eel00.XE0[KE]))*XEK);
    }
  else
    {
      TH1=exp(log(eel00.XP1[KE])+(log(eel00.XP1[KE+1])-log(eel00.XP1[KE]))*XEK);
      if(eel00.T1P0[KE] > 1.0E-35)
	{
	  TH1=TH1-exp(log(eel00.T1P0[KE])+(log(eel00.T1P0[KE+1])-log(eel00.T1P0[KE]))*XEK);
	}  
      TH2=exp(log(eel00.XP2[KE])+(log(eel00.XP2[KE+1])-log(eel00.XP2[KE]))*XEK);
      if(eel00.T2P0[KE] > 1.0E-35)
	{
	  TH2=TH2-exp(log(eel00.T2P0[KE])+(log(eel00.T2P0[KE+1])-log(eel00.T2P0[KE]))*XEK);
	}  
      SH0=exp(mat.SPHEL[KE]+(mat.SPHEL[KE+1]-mat.SPHEL[KE])*XEK)/mat.VMOL;
      AA=exp(mat.AP[KE]+(mat.AP[KE+1]-mat.AP[KE])*XEK);
      BB=mat.BP[KE]+(mat.BP[KE+1]-mat.BP[KE])*XEK;
      RNDC=mat.RNDCP[KE]+(mat.RNDCP[KE+1]-mat.RNDCP[KE])*XEK;
      printf(" # Positron energy = %13.5E eV\n",E);
      XS0=exp(log(eel00.XP0[KE])+(log(eel00.XP0[KE+1])-log(eel00.XP0[KE]))*XEK);      
    }

  double X1AVE=TH1/(2.0E0*SH0);
  double X2AVE=(3.0E0*TH1-TH2)/(6.0E0*SH0);
  printf(" # Elastic scattering of electrons/positrons\n");
  printf(" # Hard x-section (cm**2)               = %14.7E\n",SH0);
  printf(" # Hard 1st transport x-section (cm**2) = %14.7E\n",TH1);
  printf(" # Hard 2nd transport x-section (cm**2) = %14.7E\n",TH2);
  printf(" #    A =%14.7E\n",AA);
  printf(" #    B =%14.7E\n",BB);
  printf(" # RNDC =%14.7E\n",RNDC);
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double RENORM=SH0/XS0;
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Elastic scattering of electrons/positrons\n");
  fprintf(IOUT," # Electron energy (eV)  =%14.7E\n",E);
  fprintf(IOUT," # Hard x-section (cm**2)               = %14.7E\n",SH0);
  fprintf(IOUT," # Hard 1st transport x-section (cm**2) = %14.7E\n",TH1);
  fprintf(IOUT," # Hard 2nd transport x-section (cm**2) = %14.7E\n",TH2);
  double DX=1.0E0/200.0E0;
  double PDF;
  for(unsigned int K=1; K<=201; K++)
  {
    double RMU=(K-1)*DX;
    if(BB > 0.0E0)
    {
      PDF=(1.0E0-BB)*AA*(1.0E0+AA)/pow((AA+RMU),2);
    }
    else
    {
      double RB=-BB;
      PDF=(1.0E0-RB)*AA*(1.0E0+AA)/pow((AA+RMU),2);
      if(RMU > 0.5E0){ PDF=PDF+RB*8.0E0*(RMU-0.5E0);}
    }
    fprintf(IOUT,"%14.6E%14.6E\n",RMU,PDF);
  }
  fclose(IOUT);
//
//  ************  Random sampling test.
//
  double X1MAX=1.0E0;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I <= NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=1000000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      double XX;
      EELa(AA,BB,RNDC,XX,randoms);
      X1AV=X1AV+XX;
      double XX2=XX*XX;
      X2AV=X2AV+XX2;
      X4AV=X4AV+XX2*XX2;
      int INDEX=(int)((XX-X1MIN)/DX1+1.0E0);
      if(INDEX > 0 && INDEX <= NCHAN)
      {
        KC1[INDEX-1]=KC1[INDEX-1]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    double X1AVMC=X1AV/FNT;
    double X2AVMC=X2AV/FNT;
    printf(" N,<MU>,<MU**2> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  double X1AVMC=X1AV*FNT;
  double X2AVMC=X2AV*FNT;
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  if(KPAR == PEN_ELECTRON)
    printf(" # Electron energy = %13.5E eV\n",E);
  else
    printf(" # Positron energy = %13.5E eV\n",E);
  
  printf(" # Hard x-section (cm**2)               = %14.7E\n",SH0);
  printf(" # Hard 1st transport x-section (cm**2) = %14.7E\n",TH1);
  printf(" # Hard 2nd transport x-section (cm**2) = %14.7E\n",TH2);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Elastic scattering of electrons/positrons\n");

  if(KPAR == PEN_ELECTRON)
    fprintf(IOUT," # Electron energy (eV)  =%13.5E\n",E);
  else
    fprintf(IOUT," # Positron energy (eV)  =%13.5E\n",E);
    
  
  fprintf(IOUT," # Hard x-section (cm**2)               = %14.7E\n",SH0);
  fprintf(IOUT," # Hard 1st transport x-section (cm**2) = %14.7E\n",TH1);
  fprintf(IOUT," # Hard 2nd transport x-section (cm**2) = %14.7E\n",TH2);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  for(int I=1; I <= NCHAN; I++)
  {
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC*RENORM > 1.0E-35 ? PDFMC*RENORM : 1.0E-35),
     (ERR*RENORM > 1.0E-35 ? ERR*RENORM : 1.0E-35));
  }
  fclose(IOUT);
  return;
}

//  *********************************************************************
//                       SUBROUTINE EELdV
//  *********************************************************************
void EELdV(const pen_context& context, const pen_material& mat, const double E, const int KPAR, const CEEL00& eel00, const int NSAMP, CDCSEP& dcsep, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm for elastic scatter-
//  ing of electrons and positrons. To run this subroutine, load a single
//  material file (since it uses tables that are overwritten when a new
//  material is loaded).  DCSs from the ELSEPA database.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
//  ****  Composition data.
//  using namespace COMPOS;
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Electron simulation tables.
//  using namespace CEIMFP;
//  ****  Positron simulation tables.
//  using namespace CPIMFP;
//  ****  Total and transport cross sections.
//  using namespace CEEL01;
//  ****  Auxiliary arrays for random sampling verification.
  double XC1[1000];
  int    KC1[1000];
//
//  ****  Differential cros sections (ELSEPA).
//
//  using namespace CDCSEP;
//
//  using namespace CRITA;
//
//  using namespace CEELDB;
//  using namespace CPELDB;
//
//  using namespace CELSEP;

  double XS0,SH0,TH1,TH2,RNDC;
  CRITA rita;
  //Get rita pointers alias
  const double* XTI  = rita.XT; 
  const double* PACI = rita.PAC; 
  const double* AI   = rita.A; 
  const double* BI   = rita.B;      
  
//
//  ****  Table of cutoff randon number values for class II simulation.
//        (used to verify that linear interpolation in XEL is accurate).
//
  FILE* IOUT=fopen("rndcs.dat","w");
  for(unsigned int I=0; I < constants::NEGP; I++)
  {
    fprintf(IOUT,"%4d%17.9E%17.9E%17.9E\n",I+1,context.grid.ET[I],mat.RNDCEd[I],mat.RNDCPd[I]);
  }
  fclose(IOUT);

  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);

//
//  ****  Moments of the 'exact' distribution.
//
	  
  if (KPAR == PEN_ELECTRON)
    {
      if(E > mat.EELMAX)
	{
	  printf(" # Modified Wentzel DCS model\n");
	  EELaV(context, mat, E, KPAR, eel00, NSAMP, randoms);
	  return;
	}
      else
	{
	  printf(" # Partial-wave differential cross section\n");
	  RNDC=mat.RNDCEd[KE]+(mat.RNDCEd[KE+1]-mat.RNDCEd[KE])*XEK;
	  DCSEL0(E,-1,dcsep);
	  if(penGetError() != PEN_SUCCESS){
	    printf("%s\n",penErrorString(penGetError()));
	    exit(penGetError());
	  }
	  double ERRM;
	  RITAI0(DCSEL,&dcsep,0.0E0,1.0E0,pen_material::NP_P,pen_material::NP_P/4,ERRM,1,rita);
	  
	  //
	  double RU=RNDC;
	  int I=1;
	  int J=pen_material::NP_P;
	  bool Eixir = false;
	  while(!Eixir)
	    {
	      Eixir = true;
	      int K = (I+J)/2;
	      if(RU > PACI[K-1])
		{
		  I = K;
		}
	      else
		{
		  J = K;
		}
	      if(J-I > 1){ Eixir = false; continue;}
	    }
	  //
	  double RR=RU-PACI[I-1];
	  double DPRO=PACI[I+1-1]-PACI[I-1];
	  double RMUC;
	  if(DPRO < 1.0E-10)
	    {
	      RMUC=XTI[I-1];
	    }
	  else
	    {
	      double CI=(1.0E0+AI[I-1]+BI[I-1])*DPRO;
	      RMUC=XTI[I-1]+(CI*RR/(pow(DPRO,2)+(DPRO*AI[I-1]+BI[I-1]*RR)*RR))*(XTI[I+1-1]-XTI[I-1]);
	    }
	  //
	  double XM0,XM1,XM2;
	  RITAM(RMUC,1.0E0,XM0,XM1,XM2,rita);
	  XS0=dcsep.CSI;
	  double ESH0=dcsep.CSI*XM0;
	  double ETH1=dcsep.CSI*XM1;
	  double ETH2=dcsep.CSI*XM2;	  
	  SH0=ESH0;
	  TH1=2.0E0*ETH1;
	  TH2=6.0E0*(ETH1-ETH2);
	}
      printf("# Electron energy = %13.5E eV\n",E);
    }
  else
    {
      if(E > mat.PELMAX)
	{
	  printf("# Modified Wentzel DCS model\n");
	  EELaV(context, mat, E, KPAR, eel00, NSAMP, randoms);
	  return;
	}
      else
	{
	  printf("# Partial-wave differential cross section\n");
	  RNDC=mat.RNDCPd[KE]+(mat.RNDCPd[KE+1]-mat.RNDCPd[KE])*XEK;
	  DCSEL0(E,+1,dcsep);
	  double ERRM;
	  RITAI0(DCSEL,&dcsep,0.0E0,1.0E0,pen_material::NP_P,pen_material::NP_P/4,ERRM,1,rita);
	  
	  //
	  double RU=RNDC;
	  int I=1;
	  int J=pen_material::NP_P;
	  bool Eixir = false;
	  while(!Eixir)
	    {
	      Eixir = true;
	      int K = (I+J)/2;
	      if(RU > PACI[K-1])
		{
		  I = K;
		}
	      else
		{
		  J = K;
		}
	      if(J-I > 1){ Eixir = false; continue;}
	    }
	  //
	  double RR=RU-PACI[I-1];
	  double DPRO=PACI[I+1-1]-PACI[I-1];
	  double RMUC;
	  if(DPRO < 1.0E-10)
	    {
	      RMUC=XTI[I-1];
	    }
	  else
	    {
	      double CI=(1.0E0+AI[I-1]+BI[I-1])*DPRO;
	      RMUC=XTI[I-1]+(CI*RR/(pow(DPRO,2)+(DPRO*AI[I-1]+BI[I-1]*RR)*RR))*(XTI[I+1-1]-XTI[I-1]);
	    }
	  //
	  double XM0,XM1,XM2;
	  RITAM(RMUC,1.0E0,XM0,XM1,XM2,rita);
	  XS0=dcsep.CSI;
	  double ESH0=dcsep.CSI*XM0;
	  double ETH1=dcsep.CSI*XM1;
	  double ETH2=dcsep.CSI*XM2;
	  SH0=ESH0;
	  TH1=2.0E0*ETH1;
	  TH2=6.0E0*(ETH1-ETH2);
	}
      printf("# Positron energy = %13.5E eV\n",E);
    }
  
//
  double X1AVE=TH1/(2.0E0*SH0);
  double X2AVE=(3.0E0*TH1-TH2)/(6.0E0*SH0);
  printf("# Elastic scattering of electrons/positrons\n");
  printf("# Hard x-section (cm**2)               = %14.7E\n",SH0);
  printf("# Hard 1st transport x-section (cm**2) = %14.7E\n",TH1);
  printf("# Hard 2nd transport x-section (cm**2) = %14.7E\n",TH2);
  printf("# RNDC =%14.7E\n",RNDC);
  printf("#    First moment (numerical) =%14.7E\n",X1AVE);
  printf("#   Second moment (numerical) =%14.7E\n",X2AVE);
  double RENORM=SH0/XS0;
  RENORM=RENORM/(4.0E0*constants::PI);
//
  IOUT = fopen("calc.dat","w");

  fprintf(IOUT," # Elastic scattering of electrons/positrons\n");
  fprintf(IOUT," # Electron energy (eV)  =%14.7E\n",E);
  fprintf(IOUT," # Hard x-section (cm**2)               = %14.7E\n",SH0);
  fprintf(IOUT," # Hard 1st transport x-section (cm**2) = %14.7E\n",TH1);
  fprintf(IOUT," # Hard 2nd transport x-section (cm**2) = %14.7E\n",TH2);

  int NPDEC=40;
  double FACT=pow(10.0E0,(1.0E0/double(NPDEC)));
  double RMU=1.0E-13/FACT;
  for(unsigned int K=1; K <= 10000; K++)
  {
    if(RMU < 0.0999999E0)
    {
      RMU=RMU*FACT;
    }
    else
    {
      RMU=RMU+0.005E0;
    }
    double PDF=DCSEL(RMU,&dcsep)/XS0;
//
    double RMUM=(RMU < 0.99999999999E0 ? RMU : 0.99999999999E0);
    double DCSKE=RIDCS(mat, KPAR, RMUM, KE)/(4.0E0*constants::PI);
    double DCSK1=RIDCS(mat, KPAR, RMUM, KE+1)/(4.0E0*constants::PI);
    double DENOM=log(context.grid.ET[KE+1])-log(context.grid.ET[KE]);
    double EDCSI=((log(context.grid.ET[KE+1])-log(E))*DCSKE+(log(E)-log(context.grid.ET[KE]))*DCSK1)/DENOM;
    double ERR=(EDCSI-PDF)/PDF;
    fprintf(IOUT,"%14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",RMU,PDF,DCSKE,DCSK1,EDCSI,ERR);
    if(RMU > 0.9999999999E0){break;}
  }
  fclose(IOUT);
//
//  ************  Random sampling test.
//
  double X1MAX=1.0E0;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=1000000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double TIMEI=CPUtime();
  double XX;
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      if(KPAR == PEN_ELECTRON)
      {
	
        EELd(mat,KE,XEL,RNDC,XX,context.grid.DLEMP,context.grid.DLFC,randoms);  // Uses the ELSEPA database.
      }
      else
      {
        PELd(mat,KE,XEL,RNDC,XX,context.grid.DLEMP,context.grid.DLFC,randoms);  // Uses the ELSEPA database.
      }

      X1AV=X1AV+XX;
      double XX2=XX*XX;
      X2AV=X2AV+XX2;
      X4AV=X4AV+XX2*XX2;
      int INDEX=(int)((XX-X1MIN)/DX1+1.0E0);
      if(INDEX > 0 && INDEX <= NCHAN)
      {
        KC1[INDEX-1]=KC1[INDEX-1]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    double X1AVMC=X1AV/FNT;
    double X2AVMC=X2AV/FNT;
    printf(" N,<MU>,<MU**2> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  double X1AVMC=X1AV*FNT;
  double X2AVMC=X2AV*FNT;
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  if(KPAR == PEN_ELECTRON)
  {
    printf(" # Electron energy = %13.5E eV\n",E);
  }
  else
  {
    printf(" # Positron energy = %13.5E eV\n",E);
  }
  printf("# Hard x-section (cm**2)               = %14.7E\n",SH0);
  printf("# Hard 1st transport x-section (cm**2) = %14.7E\n",TH1);
  printf("# Hard 2nd transport x-section (cm**2) = %14.7E\n",TH2);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT,"# Elastic scattering of electrons/positrons\n");
  if(KPAR == PEN_ELECTRON)
  {
    fprintf(IOUT,"# Electron energy (eV)  =%13.5E\n",E);
  }
  else
  {
    fprintf(IOUT,"# Positron energy (eV)  =%13.5E\n",E);
  }
  fprintf(IOUT,"# Hard x-section (cm**2)               = %14.7E\n",SH0);
  fprintf(IOUT,"# Hard 1st transport x-section (cm**2) = %14.7E\n",TH1);
  fprintf(IOUT,"# Hard 2nd transport x-section (cm**2) = %14.7E\n",TH2);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  for(int I=1; I <= NCHAN; I++)
  {
    double PDFMC=(double)KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;

    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC*RENORM > 1.0E-35 ? PDFMC*RENORM : 1.0E-35),
     (ERR*RENORM > 1.0E-35 ? ERR*RENORM : 1.0E-35));
  }
  fclose(IOUT);
  return;
}

//  *********************************************************************
//                       FUNCTION RIDCS
//  *********************************************************************
double RIDCS(const pen_material& mat, const int KPAR, double XX, int IE)
{
//
//     Delivers the value of the RITA rational approximation to the pdf
//  at the point X.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;

//  ****  Composition data.
//  using namespace COMPOS;
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//
//  const unsigned int NP_P=128;
//  using namespace CEELDB;
//  using namespace CPELDB;
//
  double RIDCS_RETURN;

  if(KPAR == PEN_ELECTRON)
  {
    if(XX > mat.XSE[pen_material::NP_P-1][IE] || XX < mat.XSE[0][IE])
    { 
      return RIDCS_RETURN=0.0E0;
    }
    int I=0;
    int I1=pen_material::NP_P-1;
    bool Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      int IT = (I+I1)/2;
      if(XX > mat.XSE[IT][IE])
      {
        I = IT;
      }
      else
      {
        I1 = IT;
      }
      if(I1-I > 1){ Eixir = false; continue;}
    }
//
    double TAU=(XX-mat.XSE[I][IE])/(mat.XSE[I+1][IE]-mat.XSE[I][IE]);
    double CON1=2.0E0*mat.BSE[I][IE]*TAU;
    double CI=1.0E0+mat.ASE[I][IE]+mat.BSE[I][IE];
    double CON2=CI-mat.ASE[I][IE]*TAU;
    double ETA;
    if(fabs(CON1) > 1.0E-10*fabs(CON2))
    {
      ETA=CON2*(1.0E0-sqrt(1.0E0-2.0E0*TAU*CON1/pow(CON2,2)))/CON1;
    }
    else
    {
      ETA=TAU/CON2;
    }
    RIDCS_RETURN=(mat.PSE[I+1][IE]-mat.PSE[I][IE])*pow(1.0E0+(mat.ASE[I][IE]+mat.BSE[I][IE]*ETA)*ETA,2)/((1.0E0-mat.BSE[I][IE]*ETA*ETA)*CI*(mat.XSE[I+1][IE]-mat.XSE[I][IE]));
    return RIDCS_RETURN;
  }
  else
  {
//
   if(XX > mat.XSP[pen_material::NP_P-1][IE] || XX < mat.XSP[0][IE])
    { 
      return RIDCS_RETURN=0.0E0;
    }
    int I=0;
    int I1=pen_material::NP_P-1;
    bool Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      int IT = (I+I1)/2;
      if(XX > mat.XSP[IT][IE])
      {
        I = IT;
      }
      else
      {
        I1 = IT;
      }
      if(I1-I > 1){ Eixir = false; continue;}
    }
//
    double TAU=(XX-mat.XSP[I][IE])/(mat.XSP[I+1][IE]-mat.XSP[I][IE]);
    double CON1=2.0E0*mat.BSP[I][IE]*TAU;
    double CI=1.0E0+mat.ASP[I][IE]+mat.BSP[I][IE];
    double CON2=CI-mat.ASP[I][IE]*TAU;
    double ETA;
    if(fabs(CON1) > 1.0E-10*fabs(CON2))
    {
      ETA=CON2*(1.0E0-sqrt(1.0E0-2.0E0*TAU*CON1/pow(CON2,2)))/CON1;
    }
    else
    {
      ETA=TAU/CON2;
    }
    RIDCS_RETURN=(mat.PSP[I+1][IE]-mat.PSP[I][IE])*pow(1.0E0+(mat.ASP[I][IE]+mat.BSP[I][IE]*ETA)*ETA,2)/((1.0E0-mat.BSP[I][IE]*ETA*ETA)*CI*(mat.XSP[I+1][IE]-mat.XSP[I][IE]));
    return RIDCS_RETURN;
  }
}

//  *********************************************************************
//                       SUBROUTINE EINaV
//  *********************************************************************
void EINaV(const pen_context& context, const pen_material& mat, const double E, int NSAMP,int IDIST, CEIN00& ein00, const CESIN& esin, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm of electron inelastic
//  interactions.
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
//
  const double TREV=constants::REV+constants::REV;
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Composition data.
//  using namespace COMPOS;
//  ****  E/P inelastic collisions.
//  using namespace CEIN;
//  ****  Electron simulation tables.
//  using namespace CEIMFP;
//  ****  Partial cross sections of individual shells/oscillators.
//  using namespace CEIN00;
//  ****  Electron inelastic coll. and inner-shell ionization tables.
//  using namespace CEINAC;
//  using namespace CESIAC;
//  using namespace CESIN;
//  ****  Inner-shell ionization by electron and positron impact.
//  using namespace CESI0;
//  ****  Auxiliary arrays for random sampling verification.
  double WT[1000],XC1[1000];
  int    KC1[1000];

//  using namespace CPARIN;
//  using namespace CDCSIN;
//
//  using namespace CSUMGA;

  C_PARIN_DCSIN parin_dcsin;
  //Get struct variables for code clarity
  double& ZTOTM = parin_dcsin.ZTOTM;
  double& EXPOTM = parin_dcsin.EXPOTM;
  double* FM = parin_dcsin.FM;
  double* UIM = parin_dcsin.UIM;
  double* WRIM = parin_dcsin.WRIM;
  double* QRIM = parin_dcsin.QRIM;
  double* WTH = parin_dcsin.WTH;
  double* FCIN = parin_dcsin.FCIN;
  int& NOSCM = parin_dcsin.NOSCM;

  double& EM = parin_dcsin.EM;
  double& DELTA = parin_dcsin.DELTA;
  double& WCCM = parin_dcsin.WCCM;
  int& IELEC = parin_dcsin.IELEC;
  int& MOM = parin_dcsin.MOM;  
  
  EM=E;
  WCCM=mat.WCC;
  IELEC=-1;
  double WU=E*0.5E0;
//
  ZTOTM=mat.ZT;
  EXPOTM=mat.EXPOT;
  NOSCM=mat.NOSC;
  for(int I=0; I<NOSCM; I++)
  {
    FM[I]=mat.F[I];
    double UK=mat.UI[I];
    double WK=mat.WRI[I];
    double WM,WKP,QKP;
    if(UK > 1.0E-3)
    {
      WTH[I]=UK;
      WM=3.0E0*WK-2.0E0*UK;
      if(E > WM)
      {
        WKP=WK;
        QKP=UK;
      }
      else
      {
        WKP=(E+2.0E0*UK)/3.0E0;
        QKP=UK*(E/WM);
      }
    }
    else
    {
      UK=0.0E0;
      WTH[I]=WK;
      WKP=WK;
      QKP=WK;
    }
    UIM[I]=UK;
    WRIM[I]=WKP;
    if(UIM[I] < E){ WU=(WU > 0.5E0*(E+UIM[I]) ? WU : 0.5E0*(E+UIM[I]));} 
    QRIM[I]=QKP;
  }
  
  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);

  //
  //  ****  Calculation of 'exact' moments.
  //
  double XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2;
  EINaT(mat,ein00,E,WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA);

  double XX0=0.0E0;
  double XX1=0.0E0;
  double XX2=0.0E0;
  XT1=0.0E0;
  XT2=0.0E0;
  for(int IO=0; IO<mat.NEIN; IO++)
  {
    int KO=mat.IEIN[IO];
    double TCS0=ein00.SEH0[KO-1];
    if(TCS0 > 1.0E-35)
    {
      FCIN[KO-1]=exp(log((esin.XSEIN[KE][IO] > 1.0E-35 ? esin.XSEIN[KE][IO] : 1.0E-35))
        +(log((esin.XSEIN[KE+1][IO] > 1.0E-35 ? esin.XSEIN[KE+1][IO] : 1.0E-35))
        -log((esin.XSEIN[KE][IO] > 1.0E-35 ? esin.XSEIN[KE][IO] : 1.0E-35)))*XEK)/TCS0;
    }
    else
    {
      FCIN[KO-1]=1.0E0;
    }
    XX0=XX0+ein00.SEH0[KO-1]*FCIN[KO-1];
    XX1=XX1+ein00.SEH1[KO-1]*FCIN[KO-1];
    XX2=XX2+ein00.SEH2[KO-1]*FCIN[KO-1];
    XT1=XT1+ein00.SET1[KO-1]*FCIN[KO-1];
    XT2=XT2+ein00.SET2[KO-1]*FCIN[KO-1];
  }
  double XX0C=XX0;
//
  for(int IO=0; IO<mat.NESI; IO++)
  {
    int KO=mat.IESI[IO];
    double TCS0=ein00.SES0[KO-1]+ein00.SEH0[KO-1];
    double SNUM=exp(log((esin.XSESI[KE][IO] > 1.0E-35 ? esin.XSESI[KE][IO] : 1.0E-35))
             +(log((esin.XSESI[KE+1][IO] > 1.0E-35 ? esin.XSESI[KE+1][IO] : 1.0E-35))
             -log((esin.XSESI[KE][IO] > 1.0E-35 ? esin.XSESI[KE][IO] : 1.0E-35)))*XEK);
    if(TCS0 > 1.0E-35)
    {
      FCIN[KO-1]=SNUM/TCS0;
      XX0=XX0+(ein00.SES0[KO-1]+ein00.SEH0[KO-1])*FCIN[KO-1];
      XX1=XX1+(ein00.SES1[KO-1]+ein00.SEH1[KO-1])*FCIN[KO-1];
      XX2=XX2+(ein00.SES2[KO-1]+ein00.SEH2[KO-1])*FCIN[KO-1];
      XT1=XT1+ein00.SET1[KO-1]*FCIN[KO-1];
      XT2=XT2+ein00.SET2[KO-1]*FCIN[KO-1];
    }
    else if(SNUM > 1.0E-35)
    {
      FCIN[KO-1]=1.0E0;
      XX0=XX0+SNUM;
      XX1=XX1+SNUM*mat.UI[KO-1];
      XX2=XX2+SNUM*pow(mat.UI[KO-1],2);
    }
    else
    {
      FCIN[KO-1]=0.0E0;
    }
  }
//
  printf(" # Inelastic scattering of electrons\n");
  printf(" # Electron energy (eV) =%14.7E\n",E);
  printf(" #   Cutoff energy (eV) =%14.7E\n",WCCM);
  printf(" # Hard x-section (cm**2) =             %14.7E\n",XX0);
  printf(" # Hard stop. x-section (eV*cm**2) =    %14.7E\n",XX1);
  printf(" # Hard strg. x-section (eV**2*cm**2) = %14.7E\n",XX2);
  if(XX0 < 1.0E-40){ return;}
  double W1AVE=XX1/XX0;
  double W2AVE=XX2/XX0;
  double RCUT=XX0C/XX0;
//
  double WCCINF=E+E;
  double XH0T,XH1T,XH2T,XS0T,XS1T,XS2T,XT1T,XT2T;
  EINaT(mat,ein00,E,WCCINF,XH0T,XH1T,XH2T,XS0T,XS1T,XS2T,XT1T,XT2T,DELTA);
  XT1T=0.0E0;
  XT2T=0.0E0;
  for(int IO=0; IO<mat.NEIN; IO++)
  {
    int KO=mat.IEIN[IO];
    XT1T=XT1T+ein00.SET1[KO-1]*FCIN[KO-1];
    XT2T=XT2T+ein00.SET2[KO-1]*FCIN[KO-1];
  }
//
  for(int IO=0; IO<mat.NESI; IO++)
  {
    int KO=mat.IESI[IO];
    XT1T=XT1T+ein00.SET1[KO-1]*FCIN[KO-1];
    XT2T=XT2T+ein00.SET2[KO-1]*FCIN[KO-1];
  }
  double U1AVE=0.5E0*(XT1T-XT1)/XX0;
  double U2AVE=U1AVE-(XT2T-XT2)/(6.0E0*XX0);
//
  double CP2=E*(E+TREV);
  double CP=sqrt(CP2);
  double CPP2=(E-WU)*(E-WU+TREV);
  double CPP=sqrt(CPP2);
  double RMUM;
  if(CPP > 1.0E-9)
  {
    RMUM=(WU*(WU+TREV)-pow(CP-CPP,2))/(4.0E0*CP*CPP);
  }
  else
  {
    RMUM=0.5E0;
  }
//
  printf("  \n");
  MOM=0;
  double SUM0=SUMGA(EDCSIA,&parin_dcsin,0.0E0,RMUM,1.0E-8);
  if(penGetError() != PEN_SUCCESS){ printf(" SUMGA error at 1e\n");}
  MOM=1;
  double U1N=SUMGA(EDCSIA,&parin_dcsin,0.0E0,RMUM,1.0E-8)/SUM0;
  if(penGetError() != PEN_SUCCESS){ printf(" SUMGA error at 2e\n");}
  MOM=2;
  double U2N=SUMGA(EDCSIA,&parin_dcsin,0.0E0,RMUM,1.0E-8)/SUM0;
  if(penGetError() != PEN_SUCCESS){ printf(" SUMGA error at 3e\n");}
  printf(" MU1AV A =%15.7E\n",U1AVE);
  printf(" MU1AV B =%15.7E\n",U1N);
  printf(" MU2AV A =%15.7E\n",U2AVE);
  printf(" MU2AV B =%15.7E\n",U2N);
//
  if(IDIST == 1)
  {

  printf("  \n");
//
//  ************  Angular distribution.
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Inelastic scattering of electrons\n");
  fprintf(IOUT," # Angular distribution\n");
  fprintf(IOUT," # Electron energy (eV) =%14.7E\n",E);
  fprintf(IOUT," #   Cutoff energy (eV) = %14.7E\n",WCCM);
  fprintf(IOUT," # Hard x-section (cm**2) =             %14.7E\n",XX0);
  fprintf(IOUT," # Hard stop. x-section (eV*cm**2) =    %14.7E\n",XX1);
  fprintf(IOUT," # Hard strg. x-section (eV**2*cm**2) = %14.7E\n",XX2);
  MOM=0;
  int NTAB=1000;
  double DX=RMUM*1.05E0/double(NTAB);
//
  for(int K=1; K<=NTAB+1; K++)
  {
    double RMU=(K-1)*DX;
    double PDF=4.0E0*constants::PI*EDCSIA(RMU,&parin_dcsin)/XX0;
    fprintf(IOUT,"%14.6E%14.6E\n",RMU,PDF);
  }
  fclose(IOUT);
//
//  ****  Random sampling test.
//
  double X1MAX=RMUM;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=500000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double W1AV=0.0E0;
  double W2AV=0.0E0;
  double W4AV=0.0E0;
  double U1AV=0.0E0;
  double U2AV=0.0E0;
  double U4AV=0.0E0;
  double W1AVMC=0.0E0,W2AVMC=0.0E0,U1AVMC=0.0E0,U2AVMC=0.0E0;
  double TIMEI=CPUtime();
  int NFAIL=0;
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
//
      double DE,EP,CDT,ES,CDTS;
      int IOSC,IZZ,ISH;
      double RR=randoms.rand();
      if(RR < RCUT)
      {
        EINa(mat,E,KE,XEL,DE,EP,CDT,ES,CDTS,IOSC,DELTA,context.grid.DLEMP,context.grid.DLFC,randoms);
      }
      else
      {
        ESIa(mat,E,KE,XEL,DE,EP,CDT,ES,CDTS,IZZ,ISH,context.grid.DLEMP,context.grid.DLFC,DELTA,randoms);	
      }
      if(DE < E)
      {
        W1AV=W1AV+DE;
        double DE2=DE*DE;
        W2AV=W2AV+DE2;
        W4AV=W4AV+DE2*DE2;
        double RMU=(1.0E0-CDT)*0.5E0;
        U1AV=U1AV+RMU;
        double RMU2=RMU*RMU;
        U2AV=U2AV+RMU2;
        U4AV=U4AV+RMU2*RMU2;
        int INDEX=(int)((RMU-X1MIN)/DX1+1.0E0);
        if(INDEX > 0 && INDEX <= NCHAN)
        {  
          KC1[INDEX-1]=KC1[INDEX-1]+1;
        }
        else
        {
          printf(" Index out of range.\n");
        }
      }
      else
      {
        NFAIL=NFAIL+1;
      }
    }
    double FNT=double(I*LOOP)-NFAIL;
    W1AVMC=W1AV/FNT;
    W2AVMC=W2AV/FNT;
    U1AVMC=U1AV/FNT;
    U2AVMC=U2AV/FNT;
    printf(" N,<W>,<W*W>,<MU>,<MU*MU> =%9d%10.3E%10.3E%10.3E%10.3E\n",I*LOOP,W1AVMC/W1AVE,W2AVMC/W2AVE,U1AVMC/U1AVE,U2AVMC/U2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP-NFAIL);
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Electron energy = %13.5E eV\n",E);
  printf(" #   Cutoff energy = %13.5E eV\n",WCCM);
  printf(" # Total x-section (cm**2) =           %13.5E\n",XX0);
  printf(" # Stopping x-section (eV*cm**2) =     %13.5E\n",XX1);
  printf(" # Straggling x-section (eV**2*cm**2) =%13.5E\n",XX2);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #      <W>   (numerical) =%14.7E\n",W1AVE);
  double ERR1MC=3.0E0*sqrt(((W2AVMC-pow(W1AVMC,2))*FNT > 1.0E-35 ? (W2AVMC-pow(W1AVMC,2))*FNT : 1.0E-35));
  printf(" #      <W> (Monte Carlo) =%14.7E +-%9.2E\n",W1AVMC,ERR1MC);
  printf(" #   <W**2>   (numerical) =%14.7E\n",W2AVE);
  double W4AVMC=W4AV*FNT;
  double ERR2MC=3.0E0*sqrt(((W4AVMC-pow(W2AVMC,2))*FNT > 1.0E-35 ? (W4AVMC-pow(W2AVMC,2))*FNT : 1.0E-35));
  printf(" #   <W**2> (Monte Carlo) =%14.7E +-%9.2E\n",W2AVMC,ERR2MC);
  printf(" #     <MU>   (numerical) =%14.7E\n",U1AVE);
  double ERU1MC=3.0E0*sqrt(((U2AVMC-pow(U1AVMC,2))*FNT > 1.0E-35 ? (U2AVMC-pow(U1AVMC,2))*FNT : 1.0E-35));
  printf(" #     <MU> (Monte Carlo) =%14.7E +-%9.2E\n",U1AVMC,ERU1MC);
  printf(" #  <MU**2>   (numerical) =%14.7E\n",U2AVE);
  double U4AVMC=U4AV*FNT;
  double ERU2MC=3.0E0*sqrt(((U4AVMC-pow(U2AVMC,2))*FNT > 1.0E-35 ? (U4AVMC-pow(U2AVMC,2))*FNT : 1.0E-35));
  printf(" #  <MU**2> (Monte Carlo) =%14.7E +-%9.2E\n",U2AVMC,ERU2MC);

  printf(" # Angular distribution\n");
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Inelastic scattering of electrons\n");
  fprintf(IOUT," # Angular distribution\n");
  fprintf(IOUT," # Electron energy (eV) =%14.7E\n",E);
  fprintf(IOUT," #   Cutoff energy (eV) = %14.7E\n",WCCM);
  fprintf(IOUT," # Total x-section (cm**2) =             %14.7E\n",XX0);
  fprintf(IOUT," # Stopping x-section (eV*cm**2) =       %14.7E\n",XX1);
  fprintf(IOUT," # Straggling x-section (eV**2*cm**2) =%14.7E\n",XX2);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #      <W>   (numerical) =%14.7E\n",W1AVE);
  fprintf(IOUT," #      <W> (Monte Carlo) =%14.7E +-%9.2E\n",W1AVMC,ERR1MC);
  fprintf(IOUT," #   <W**2>   (numerical) =%14.7E\n",W2AVE);
  fprintf(IOUT," #   <W**2> (Monte Carlo) =%14.7E +-%9.2E\n",W2AVMC,ERR2MC);
  fprintf(IOUT," #     <MU>   (numerical) =%14.7E\n",U1AVE);
  fprintf(IOUT," #     <MU> (Monte Carlo) =%14.7E +-%9.2E\n",U1AVMC,ERU1MC);
  fprintf(IOUT," #  <MU**2>   (numerical) =%14.7E\n",U2AVE);
  fprintf(IOUT," #  <MU**2> (Monte Carlo) =%14.7E +-%9.2E\n",U2AVMC,ERU2MC);

//
  for(int I=1; I<=NCHAN; I++)
  { 
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC > 1.0E-35 ? PDFMC : 1.0E-35),
     (ERR > 1.0E-35 ? ERR : 1.0E-35));
  }
  fclose(IOUT);
  return;
  }
//
//  ************  Energy distribution.
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Inelastic scattering of electrons\n");
  fprintf(IOUT," # Energy-loss distribution\n");
  fprintf(IOUT," # Electron energy (eV) =%14.7E\n",E);
  fprintf(IOUT," #   Cutoff energy (eV) = %14.7E\n",WCCM);
  fprintf(IOUT," # Hard x-section (cm**2) =             %14.7E\n",XX0);
  fprintf(IOUT," # Hard stop. x-section (eV*cm**2) =    %14.7E\n",XX1);
  fprintf(IOUT," # Hard strg. x-section (eV**2*cm**2) = %14.7E\n",XX2);
  MOM=0;
  int NTT=500;
  double DW=(E-WCCM)/double(NTT);
  int NTAB;
  for(int I=1; I<=NTT+2; I++)
  {
    NTAB=I;
    WT[I-1]=WCCM+DW*(I-1);
  }
  for(int K=1; K<=mat.NOSC; K++)
  {
    if(mat.UI[K-1] > WCCM && mat.UI[K-1] < E)
    {
      WT[NTAB+1-1]=mat.UI[K-1]-1.0E-2;
      WT[NTAB+2-1]=mat.UI[K-1]+1.0E-2;
      NTAB=NTAB+2;
    }
  }
  for(int I=1; I<=NTAB-1; I++)
  {
    int IMIN=I;
    double WMIN=WT[I-1];
    for(int J=I+1; J<=NTAB; J++)
    {
      if(WT[J-1] < WMIN)
      {
        IMIN=J;
        WMIN=WT[J-1];
      }
    }
    double SAVE=WT[I-1];
    WT[I-1]=WT[IMIN-1];
    WT[IMIN-1]=SAVE;
  }
//
  for(int K=1; K<=NTAB; K++)
  {
    double PDF=EDCSIW(WT[K-1],parin_dcsin)/XX0;
    fprintf(IOUT,"%14.6E%14.6E\n",WT[K-1],PDF);
  }
  fclose(IOUT);
//
//  ****  Random sampling test.
//
  double X1MAX=E;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=500000;
  int NLOOP=(int)(1+(NSAMP/LOOP));
  printf("NLOOP =%d\n",NLOOP);
  double W1AV=0.0E0;
  double W2AV=0.0E0;
  double W4AV=0.0E0;
  double U1AV=0.0E0;
  double U2AV=0.0E0;
  double U4AV=0.0E0;
  double TIMEI=CPUtime();
  int NFAIL=0;
  double W1AVMC=0.0E0,W2AVMC=0.0E0,U1AVMC=0.0E0,U2AVMC=0.0E0;

  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
//
      double RR=randoms.rand();
      double DE,EP,CDT,ES,CDTS;
      if(RR < RCUT)
      {
        int IOSC; 
        EINa(mat,E,KE,XEL,DE,EP,CDT,ES,CDTS,IOSC,DELTA,context.grid.DLEMP,context.grid.DLFC,randoms);
      }
      else
      {
        int IZZ,ISH; 
        ESIa(mat,E,KE,XEL,DE,EP,CDT,ES,CDTS,IZZ,ISH,context.grid.DLEMP,context.grid.DLFC,DELTA,randoms);
      }
      if(DE < E)
      {
        W1AV=W1AV+DE;
        double DE2=DE*DE;
        W2AV=W2AV+DE2;
        W4AV=W4AV+DE2*DE2;
        double RMU=(1.0E0-CDT)*0.5E0;
        U1AV=U1AV+RMU;
        double RMU2=RMU*RMU;
        U2AV=U2AV+RMU2;
        U4AV=U4AV+RMU2*RMU2;
        int INDEX=(int)((DE-X1MIN)/DX1+1.0E0);
        if(INDEX > 0 && INDEX <= NCHAN)
        {
          KC1[INDEX-1]=KC1[INDEX-1]+1;
        }
        else
        {
          printf(" Index out of range.\n");
        }
      }
      else
      {
        NFAIL=NFAIL+1;
      }
    }
    double FNT=double(I*LOOP)-NFAIL;
    W1AVMC=W1AV/FNT;
    W2AVMC=W2AV/FNT;
    U1AVMC=U1AV/FNT;
    U2AVMC=U2AV/FNT;
    printf(" N,<W>,<W*W>,<MU>,<MU*MU> =%9d%10.3E%10.3E%10.3E%10.3E\n",I*LOOP,W1AVMC/W1AVE,W2AVMC/W2AVE,U1AVMC/U1AVE,U2AVMC/U2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP-NFAIL);
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Electron energy = %13.5E eV\n",E);
  printf(" #   Cutoff energy = %13.5E eV\n",WCCM);
  printf(" # Total x-section (cm**2) =           %13.5E\n",XX0);
  printf(" # Stopping x-section (eV*cm**2) =     %13.5E\n",XX1);
  printf(" # Straggling x-section (eV**2*cm**2) =%13.5E\n",XX2);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #      <W>   (numerical) =%14.7E\n",W1AVE);
  double ERR1MC=3.0E0*sqrt(((W2AVMC-pow(W1AVMC,2))*FNT > 1.0E-35 ? (W2AVMC-pow(W1AVMC,2))*FNT : 1.0E-35));
  printf(" #      <W> (Monte Carlo) =%14.7E +-%9.2E\n",W1AVMC,ERR1MC);
  printf(" #   <W**2>   (numerical) =%14.7E\n",W2AVE);
  double W4AVMC=W4AV*FNT;
  double ERR2MC=3.0E0*sqrt(((W4AVMC-pow(W2AVMC,2))*FNT > 1.0E-35 ? (W4AVMC-pow(W2AVMC,2))*FNT : 1.0E-35));
  printf(" #   <W**2> (Monte Carlo) =%14.7E +-%9.2E\n",W2AVMC,ERR2MC);
  printf(" #     <MU>   (numerical) =%14.7E\n",U1AVE);
  double ERU1MC=3.0E0*sqrt(((U2AVMC-pow(U1AVMC,2))*FNT > 1.0E-35 ? (U2AVMC-pow(U1AVMC,2))*FNT : 1.0E-35));
  printf(" #     <MU> (Monte Carlo) =%14.7E +-%9.2E\n",U1AVMC,ERU1MC);
  printf(" #  <MU**2>   (numerical) =%14.7E\n",U2AVE);
  double U4AVMC=U4AV*FNT;
  double ERU2MC=3.0E0*sqrt(((U4AVMC-pow(U2AVMC,2))*FNT > 1.0E-35 ? (U4AVMC-pow(U2AVMC,2))*FNT : 1.0E-35));
  printf(" #  <MU**2> (Monte Carlo) =%14.7E +-%9.2E\n",U2AVMC,ERU2MC);

  printf(" # Energy-loss distribution\n");
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Inelastic scattering of electrons\n");
  fprintf(IOUT," # Energy-loss distribution\n");
  fprintf(IOUT," # Electron energy (eV) =%14.7E\n",E);
  fprintf(IOUT," #   Cutoff energy (eV) = %14.7E\n",WCCM);
  fprintf(IOUT," # Total x-section (cm**2) =           %14.7E\n",XX0);
  fprintf(IOUT," # Stopping x-section (eV*cm**2) =     %14.7E\n",XX1);
  fprintf(IOUT," # Straggling x-section (eV**2*cm**2) =%14.7E\n",XX2);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #      <W>   (numerical) =%14.7E\n",W1AVE);
  ERR1MC=3.0E0*sqrt(((W2AVMC-pow(W1AVMC,2))*FNT > 1.0E-35 ? (W2AVMC-pow(W1AVMC,2))*FNT : 1.0E-35));
  fprintf(IOUT," #      <W> (Monte Carlo) =%14.7E +-%9.2E\n",W1AVMC,ERR1MC);
  fprintf(IOUT," #   <W**2>   (numerical) =%14.7E\n",W2AVE);
  W4AVMC=W4AV*FNT;
  ERR2MC=3.0E0*sqrt(((W4AVMC-pow(W2AVMC,2))*FNT > 1.0E-35 ? (W4AVMC-pow(W2AVMC,2))*FNT : 1.0E-35));
  fprintf(IOUT," #   <W**2> (Monte Carlo) =%14.7E +-%9.2E\n",W2AVMC,ERR2MC);
  fprintf(IOUT," #     <MU>   (numerical) =%14.7E\n",U1AVE);
  ERU1MC=3.0E0*sqrt(((U2AVMC-pow(U1AVMC,2))*FNT > 1.0E-35 ? (U2AVMC-pow(U1AVMC,2))*FNT : 1.0E-35));
  fprintf(IOUT," #     <MU> (Monte Carlo) =%14.7E +-%9.2E\n",U1AVMC,ERU1MC);
  fprintf(IOUT," #  <MU**2>   (numerical) =%14.7E\n",U2AVE);
  U4AVMC=U4AV*FNT;
  ERU2MC=3.0E0*sqrt(((U4AVMC-pow(U2AVMC,2))*FNT > 1.0E-35 ? (U4AVMC-pow(U2AVMC,2))*FNT : 1.0E-35));
  fprintf(IOUT," #  <MU**2> (Monte Carlo) =%14.7E +-%9.2E\n",U2AVMC,ERU2MC);
//
  for(int I=1; I<=NCHAN; I++)
  {
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],(PDFMC > 1.0E-35 ? PDFMC : 1.0E-35),(ERR > 1.0E-35 ? ERR : 1.0E-35));
  }
  fclose(IOUT);
}

//  *********************************************************************
//                       FUNCTION EDCSIW
//  *********************************************************************
double EDCSIW(double W, C_PARIN_DCSIN& parin_dcsin)
{
//
//     Energy-loss DCS for inelastic collisions of electrons or positrons
//  with kinetic energy E. The output value EDCSIW is the DCS, per unit
//  energy loss W. Does not include distant interactions with pure delta
//  oscillators (resonances).
//

//  ****  Fundamental constants and related quantities.
  const double TREV=2.0E0*constants::REV;
  const double RTREV=1.0E0/TREV;
  const double PIELR2=constants::PI*constants::ELRAD*constants::ELRAD;
//
  //using namespace CPARIN;
  //using namespace CDCSIN;

  double  (&F)[constants::NO]=parin_dcsin.FM;
  double  (&UI)[constants::NO]=parin_dcsin.UIM;
  double  (&WRI)[constants::NO]=parin_dcsin.WRIM;
  double  (&QRI)[constants::NO]=parin_dcsin.QRIM;
  int&    NOSC=parin_dcsin.NOSCM;

  double& E=parin_dcsin.EM;
  double& WCCM = parin_dcsin.WCCM;
  double& DELTA = parin_dcsin.DELTA;
  double* WTH = parin_dcsin.WTH;
  double* FCIN = parin_dcsin.FCIN;

  int& IELEC = parin_dcsin.IELEC;
  
  double EDCSIW_RETURN;
//
  if(W < WCCM || W > E)
  {
    EDCSIW_RETURN=0.0E0;
    return EDCSIW_RETURN;
  }
//
//  ****  Constants.
//
  double GAM=1.0E0+E/constants::REV;
  double GAM2=GAM*GAM;
  double BETA2=(GAM2-1.0E0)/GAM2;
  double CONST=PIELR2*TREV/BETA2;
//
  double CPS=E*(E+TREV);
  double CP=sqrt(CPS);
  double AMOL=pow(E/(E+constants::REV),2);
  double G12=pow(GAM+1.0E0,2);
  double BHA1=AMOL*(2.0E0*G12-1.0E0)/(GAM2-1.0E0);
  double BHA2=AMOL*(3.0E0+1.0E0/G12);
  double BHA3=AMOL*2.0E0*GAM*(GAM-1.0E0)/G12;
  double BHA4=AMOL*pow(GAM-1.0E0,2)/G12;
//
  double SDT1=(log(GAM2)-BETA2-DELTA > 0.0E0 ? log(GAM2)-BETA2-DELTA : 0.0E0);
//
  double DCS=0.0E0;
  for(int K=1; K<=NOSC; K++)
  {
    if(W > WTH[K-1])
    {
      if(UI[K-1] > 1.0E-3)
      {
        double WM=3.0E0*WRI[K-1]-2.0E0*UI[K-1];
        double WU;
        if(IELEC == -1)
        {
          WU=(WM < 0.5E0*(E+UI[K-1]) ? WM : 0.5E0*(E+UI[K-1]));
        }
        else
        {
          WU=(WM < E ? WM : E);
        }
        if(W < WU)
        {
          double CPPS=(E-WRI[K-1])*(E-WRI[K-1]+TREV);
          double CPP=sqrt(CPPS);
          double QM;
          if(WRI[K-1] > 1.0E-6*E)
          {
            QM=sqrt(pow(CP-CPP,2)+pow(constants::REV,2))-constants::REV;
          }
          else
          {
            QM=pow(WRI[K-1],2)/(BETA2*TREV);
            QM=QM*(1.0E0-QM*RTREV);
          }
          if(QM < QRI[K-1])
          {
            double SDL1=log(QRI[K-1]*(QM+TREV)/(QM*(QRI[K-1]+TREV)));
            double SD1=(SDL1+SDT1)/WRI[K-1];
            if(W < WU)
            {
              double FW=2.0E0*(WM-W)/pow(WM-UI[K-1],2);
              DCS=DCS+FCIN[K-1]*F[K-1]*FW*SD1;
            }
          }
        }
      }
      if(IELEC == -1)
      {
        double EE=E+UI[K-1];
        if(W < 0.5E0*EE)
        {
          double WE=W/EE;
          double R=1.0E0+pow(WE/(1.0E0-WE),2)-(1.0E0-AMOL)*(WE/(1.0E0-WE))+AMOL*pow(WE,2);
          DCS=DCS+FCIN[K-1]*F[K-1]*R/pow(W,2);
        }
      }
      else
      {
        double WE=W/E;
        double R=1.0E0-WE*(BHA1-WE*(BHA2-WE*(BHA3-WE*BHA4)));
        DCS=DCS+FCIN[K-1]*F[K-1]*R/pow(W,2);
      }
    }
  }
//
  EDCSIW_RETURN=CONST*DCS;
  return EDCSIW_RETURN;
}


//  *********************************************************************
//                       FUNCTION EDCSIA
//  *********************************************************************
double EDCSIA(double RMU, void* arg)
{
//
//     Scattering DCS for inelastic collisions of electrons or positrons
//  with kinetic energy E. The output value EDCSIA is the DCS, per unit
//  solid angle, for the polar angular deflection theta corresponding to
//  RMU=(1-cos(theta))/2.
//
  const double R4PI=1.0E0/(4.0E0*constants::PI);
//
  //using namespace CPARIN;
  //using namespace CDCSIN;

  C_PARIN_DCSIN& parin_dcsin = *((C_PARIN_DCSIN*)arg);  
  
  double  (&F)[constants::NO]=parin_dcsin.FM;
  double  (&UI)[constants::NO]=parin_dcsin.UIM;
  double  (&WRI)[constants::NO]=parin_dcsin.WRIM;
  double  (&QRI)[constants::NO]=parin_dcsin.QRIM;
  double* FCIN = parin_dcsin.FCIN;

  int&    NOSC=parin_dcsin.NOSCM;
  int&    MOM=parin_dcsin.MOM;
  int&    IELEC=parin_dcsin.IELEC;

  double& WCCM = parin_dcsin.WCCM;
  double& E=parin_dcsin.EM;
  double EDCSIA_RETURN;
//
  EDCSIA_RETURN=0.0E0;
  for(int K=1; K<=NOSC; K++)
  {
    if(E > UI[K-1])
    {
      double UII=UI[K-1];
      double WR=WRI[K-1];
      double QI=QRI[K-1];
      double F0,WM;
      if(UII > 1.0E-3)
      {
        double W1=(UII>WCCM ? UII : WCCM);
        WM=3.0E0*WR-2.0E0*UII;
        double W2;
        if(IELEC == -1)
        {
          W2=(WM < 0.5E0*(E+UII) ? WM : 0.5E0*(E+UII));
        }
        else
        {
          W2=(WM < E ? WM : E);
        }
        if(W2 > W1)
        {
          F0=(W2-W1)*(WM+WM-W2-W1)/pow(WM-UII,2);
        }
        else
        {
          F0=0.0E0;
        }
      }
      else
      {
        WM=WR;
        F0=1.0E0;
      }
      double CSD,CSC;
      DOICSC doicsc;
      DOICSS(WCCM,UII,WR,QI,WM,IELEC,E,RMU,CSD,CSC,doicsc);
      EDCSIA_RETURN=EDCSIA_RETURN+FCIN[K-1]*F[K-1]*(CSC+F0*CSD);
    }
  }
  EDCSIA_RETURN=EDCSIA_RETURN*R4PI;
  if(MOM == 0){return EDCSIA_RETURN;}
  EDCSIA_RETURN=EDCSIA_RETURN*pow(RMU,MOM);
  return EDCSIA_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE DOICSS
//  *********************************************************************
void DOICSS(double WCCM,double UI,double WR,double QI,double WM,int IELEC,double E,double RMU,double& DCSD,double& DCSC, DOICSC& doicsc)
{
//
//     Scattering DCS for collisions of electrons and positrons with a
//  unit-strength delta-oscillator with resonance energy WR and cutoff
//  recoil energy QI.
//
//  Input arguments:
//     WCCM .... cutoff energy loss (eV),
//     UI ...... ionization energy (eV),
//     WR ...... resonance energy of the oscillator (eV),
//               it should be larger than 0.01 eV.
//     QI ...... cutoff recoil energy.
//     WM ...... maximum energy loss in distant ionizing interactions,
//     E ....... kinetic energy of the projectile (eV).
//     IELEC ... electron-positron flag,
//               =-1, electrons; =+1, positrons.
//     RMU ..... =(1-COS(THETA))/2, angular deflection.
//  Output arguments:
//     DCSD .... scattering DCS (cm**2) for distant collisions, differen-
//               tial in RMU.
//     DCSC .... scattering DCS (cm**2) for close collisions.
//

//  ****  Fundamental constants and related quantities.
  const double TREV=2.0E0*constants::REV;
  const double PIELR2=constants::PI*constants::ELRAD*constants::ELRAD;

  //using namespace DOICSC;
//
  doicsc.EI=E;
  doicsc.JELEC=IELEC;
  double WU;
  if(IELEC == -1)
  {
    doicsc.EE=E+UI;
    WU=0.5E0*doicsc.EE;
  }
  else
  {
    doicsc.EE=E;
    WU=E;
  }
  DCSD=0.0E0;
  DCSC=0.0E0;
  if(WU < UI){ return;}
//
  double CPS=E*(E+TREV);
  double CP=sqrt(CPS);
  doicsc.BETA2=CPS/pow(E+constants::REV,2);
  double GAM1=E/constants::REV;
  double GAM=1.0E0+GAM1;
  double CONST=PIELR2*TREV/doicsc.BETA2;
//
//  ****  Distant interactions.
//
  if(E > WR)
  {
    double CPP2=(E-WR)*(E-WR+TREV);
    double CPP=sqrt(CPP2);
    double XPP=(CP-CPP)/constants::REV;
    if(XPP > 1.0E-4)
    {
      doicsc.C1=pow(XPP,2);
    }
    else
    {
      double WREG=WR/(E*(GAM+1.0E0));
      double XC=1.0E0+0.5E0*WREG*((1.0E0/GAM)+WREG);
      doicsc.C1=pow(WR*XC,2)/(doicsc.BETA2*pow(constants::REV,2));
    }
//
    doicsc.C2=4.0E0*CP*CPP/pow(constants::REV,2);
    doicsc.CDIST=CONST*doicsc.C2/WR;
    doicsc.RMU1=((QI/constants::REV)*((QI/constants::REV)+2.0E0)-doicsc.C1)/doicsc.C2;
//
    if(WM > WCCM)
    {
      DCSD=(DODCSD(RMU,doicsc) > 0.0E0 ? DODCSD(RMU,doicsc) : 0.0E0);
    }
    else
    {
      DCSD=0.0E0;
    }
  }
  else
  {
    DCSD=0.0E0;
  }
//
//  ****  Close collisions.
//
  doicsc.AMOL=pow((GAM1/GAM),2);
  if(IELEC != -1)
  {
    double G12=pow(GAM+1.0E0,2);
    doicsc.BHA1=doicsc.AMOL*(2.0E0*G12-1.0E0)/(GAM*GAM-1.0E0);
    doicsc.BHA2=doicsc.AMOL*(3.0E0+1.0E0/G12);
    doicsc.BHA3=doicsc.AMOL*2.0E0*GAM*GAM1/G12;
    doicsc.BHA4=doicsc.AMOL*pow(GAM1,2)/G12;
  }
//
  double CPP2=(E-UI)*(E-UI+TREV);
  double CPP=sqrt(CPP2);
  if(CPP > 1.0E-6)
  {
    doicsc.RMU2=(UI*(UI+TREV)-pow(CP-CPP,2))/(4.0E0*CP*CPP);
  }
  else
  {
    doicsc.RMU2=0.5E0;
  }
  CPP2=(E-WU)*(E-WU+TREV);
  CPP=sqrt(CPP2);
  if(CPP > 1.0E-6)
  {
    doicsc.RMU3=(WU*(WU+TREV)-pow(CP-CPP,2))/(4.0E0*CP*CPP);
  }
  else
  {
    doicsc.RMU3=0.5E0;
  }
  doicsc.CCLO=CONST*E*(E+TREV)*TREV;
//
  double WRC=(QI > WCCM ? QI : WCCM);
  if(WRC < WU)
  {
    DCSC=(DODCSC(WRC,RMU,doicsc) > 0.0E0 ? DODCSC(WRC,RMU,doicsc) : 0.0E0);
  }
  else
  {
    DCSC=0.0E0;
  }
}

//  *********************************************************************
double DODCSD(double RMU, DOICSC& doicsc)
{
//  ****  Angular DCS for distant collisions.
  //using namespace DOICSC;
  double DODCSD_RETURN;
//
  if(RMU < 0.0E0 || RMU > doicsc.RMU1)
  {
    DODCSD_RETURN=0.0E0;
    return DODCSD_RETURN;
  }
  double QQM=doicsc.C1+doicsc.C2*RMU;
  double QR;
  if(QQM > 1.0E-4)
  {
    QR=sqrt(QQM+1.0E0)-1.0E0;
  }
  else
  {
    double HQM=0.5E0*QQM;
    QR=HQM*(1.0E0-0.5E0*HQM*(1.0E0-HQM));
  }
  DODCSD_RETURN=doicsc.CDIST/(QQM*(1.0E0+QR));
  return DODCSD_RETURN;
}

//  *********************************************************************
double DODCSC(double WCCM,double RMU, DOICSC& doicsc)
{
//  ****  Angular DCS for close collisions.
  const double TREV=constants::REV+constants::REV;
  //using namespace DOICSC;

  double DODCSC_RETURN;
//
  if(RMU < doicsc.RMU2 || RMU > doicsc.RMU3)
  {
    DODCSC_RETURN=0.0E0;
    return DODCSC_RETURN;
  }
  double EMU=doicsc.EI*2.0E0*RMU*(1.0E0-RMU);
  double WN=(doicsc.EI+TREV)*EMU;
  double WE=WN/(EMU+constants::REV);
  if(WE < WCCM)
  {
    DODCSC_RETURN=0.0E0;
    return DODCSC_RETURN;
  }
  WE=WE/doicsc.EE;
  double R;
  if(doicsc.JELEC == -1)
  {
    R=1.0E0+pow(WE/(1.0E0-WE),2)-(1.0E0-doicsc.AMOL)*(WE/(1.0E0-WE))+doicsc.AMOL*pow(WE,2);
  }
  else
  {
    R=1.0E0-WE*(doicsc.BHA1-WE*(doicsc.BHA2-WE*(doicsc.BHA3-WE*doicsc.BHA4)));
  }
  DODCSC_RETURN=doicsc.CCLO*R*(1.0E0-2.0E0*RMU)/pow(WN,2);
  return DODCSC_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE PINaV
//  *********************************************************************
void PINaV(const pen_context& context, const pen_material& mat, const double E, CPSIN& psin, int NSAMP, int IDIST, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm of positron inelastic
//  interactions.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
//
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Composition data.
//  using namespace COMPOS;
//  ****  E/P inelastic collisions.
//  using namespace CEIN;
//  ****  Positron simulation tables.
//  using namespace CPIMFP;
//  ****  Partial cross sections of individual shells/oscillators.
//  using namespace CPIN00;
//  ****  Electron inelastic coll. and inner-shell ionization tables.
//  using namespace CPINAC;
//  using namespace CPSIAC;
//  using namespace CPSIN;
//  ****  Inner-shell ionization by electron and positron impact.
//  using namespace CESI0;
//  ****  Auxiliary arrays for random sampling verification.
  double WT[1000],XC1[1000];
  int    KC1[1000];

//  using namespace CPARIN;
//  using namespace CDCSIN;
//
//  using namespace CSUMGA;

  CPIN00 pin00;
  
  C_PARIN_DCSIN parin_dcsin;
  //Get struct variables for code clarity
  double& ZTOTM = parin_dcsin.ZTOTM;
  double& EXPOTM = parin_dcsin.EXPOTM;
  double* FM = parin_dcsin.FM;
  double* UIM = parin_dcsin.UIM;
  double* WRIM = parin_dcsin.WRIM;
  double* QRIM = parin_dcsin.QRIM;
  double* WTH = parin_dcsin.WTH;
  double* FCIN = parin_dcsin.FCIN;
  int& NOSCM = parin_dcsin.NOSCM;

  double& EM = parin_dcsin.EM;
  double& DELTA = parin_dcsin.DELTA;
  double& WCCM = parin_dcsin.WCCM;
  int& IELEC = parin_dcsin.IELEC;
  int& MOM = parin_dcsin.MOM;  
  
  EM=E;
  WCCM=mat.WCC;
  IELEC=+1;
//
  ZTOTM=mat.ZT;
  EXPOTM=mat.EXPOT;
  NOSCM=mat.NOSC;
  for(int I=1; I<=NOSCM; I++)
  {
    FM[I-1]=mat.F[I-1];
    double UK=mat.UI[I-1];
    double WK=mat.WRI[I-1];
    double WM,WKP,QKP;
    if(UK > 1.0E-3)
    {
      WTH[I-1]=UK;
      WM=3.0E0*WK-2.0E0*UK;
      if(E > WM)
      {
        WKP=WK;
        QKP=UK;
      }
      else
      {
        WKP=(E+2.0E0*UK)/3.0E0;
        QKP=UK*(E/WM);
      }
    }
    else
    {
      UK=0.0E0;
      WTH[I-1]=WK;
      WKP=WK;
      QKP=WK;
    }
    UIM[I-1]=UK;
    WRIM[I-1]=WKP;
    QRIM[I-1]=QKP;
  }
//

  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);

//
//  ****  Calculation of 'exact' moments.
//
  double XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2;
  PINaT(mat,pin00,E,WCCM,XH0,XH1,XH2,XS0,XS1,XS2,XT1,XT2,DELTA);
//
  double XX0=0.0E0;
  double XX1=0.0E0;
  double XX2=0.0E0;
  XT1=0.0E0;
  XT2=0.0E0;
  for(int IO=1; IO<=mat.NPIN; IO++)
  {
    int KO=mat.IPIN[IO-1];
    double TCS0=pin00.SPH0[KO-1];
    if(TCS0 > 1.0E-35)
    {
      FCIN[KO-1]=exp(log((psin.XSPIN[KE][IO-1] > 1.0E-35 ? psin.XSPIN[KE][IO-1] : 1.0E-35))
        +(log((psin.XSPIN[KE+1][IO-1] > 1.0E-35 ? psin.XSPIN[KE+1][IO-1] : 1.0E-35))
        -log((psin.XSPIN[KE][IO-1] > 1.0E-35 ? psin.XSPIN[KE][IO-1] : 1.0E-35)))*XEK)/TCS0;
    }
    else
    {
      FCIN[KO-1]=1.0E0;
    }
    XX0=XX0+pin00.SPH0[KO-1]*FCIN[KO-1];
    XX1=XX1+pin00.SPH1[KO-1]*FCIN[KO-1];
    XX2=XX2+pin00.SPH2[KO-1]*FCIN[KO-1];
    XT1=XT1+pin00.SPT1[KO-1]*FCIN[KO-1];
    XT2=XT2+pin00.SPT2[KO-1]*FCIN[KO-1];
  }
  double XX0C=XX0;
//
  for(int IO=1; IO<=mat.NPSI; IO++)
  {
    int KO=mat.IPSI[IO-1];
    double TCS0=pin00.SPS0[KO-1]+pin00.SPH0[KO-1];
    double SNUM=exp(log((psin.XSPSI[KE][IO-1] > 1.0E-35 ? psin.XSPSI[KE][IO-1] : 1.0E-35))
             +(log((psin.XSPSI[KE+1][IO-1] > 1.0E-35 ? psin.XSPSI[KE+1][IO-1] : 1.0E-35))
             -log((psin.XSPSI[KE][IO-1] > 1.0E-35 ? psin.XSPSI[KE][IO-1] : 1.0E-35)))*XEK);
    if(TCS0 > 1.0E-35)
    {
      FCIN[KO-1]=SNUM/TCS0;
      XX0=XX0+(pin00.SPS0[KO-1]+pin00.SPH0[KO-1])*FCIN[KO-1];
      XX1=XX1+(pin00.SPS1[KO-1]+pin00.SPH1[KO-1])*FCIN[KO-1];
      XX2=XX2+(pin00.SPS2[KO-1]+pin00.SPH2[KO-1])*FCIN[KO-1];
      XT1=XT1+pin00.SPT1[KO-1]*FCIN[KO-1];
      XT2=XT2+pin00.SPT2[KO-1]*FCIN[KO-1];
    }
    else if(SNUM > 1.0E-35)
    {
      FCIN[KO-1]=1.0E0;
      XX0=XX0+SNUM;
      XX1=XX1+SNUM*mat.UI[KO-1];
      XX2=XX2+SNUM*pow(mat.UI[KO-1],2);
    }
    else
    {
      FCIN[KO-1]=0.0E0;
    }
    printf("si %15.7E%15.7E%15.7E\n",XX0,XX1,XX2);
  }
//
  printf(" # Inelastic scattering of positrons\n");
  printf(" # Positron energy (eV) =%14.7E\n",E);
  printf(" #   Cutoff energy (eV) =%14.7E\n",WCCM);
  printf(" # Hard x-section (cm**2) =             %14.7E\n",XX0);
  printf(" # Hard stop. x-section (eV*cm**2) =    %14.7E\n",XX1);
  printf(" # Hard strg. x-section (eV**2*cm**2) = %14.7E\n",XX2);
  if(XX0 < 1.0E-40){ return;}
  double W1AVE=XX1/XX0;
  double W2AVE=XX2/XX0;
  double RCUT=XX0C/XX0;
//
  double WCCINF=E+E;
  double XH0T,XH1T,XH2T,XS0T,XS1T,XS2T,XT1T,XT2T;
  PINaT(mat,pin00,E,WCCINF,XH0T,XH1T,XH2T,XS0T,XS1T,XS2T,XT1T,XT2T,DELTA);
  XT1T=0.0E0;
  XT2T=0.0E0;
  for(int IO=1; IO<=mat.NPIN; IO++)
  {
    int KO=mat.IPIN[IO-1];
    XT1T=XT1T+pin00.SPT1[KO-1]*FCIN[KO-1];
    XT2T=XT2T+pin00.SPT2[KO-1]*FCIN[KO-1];
  }
//
  for(int IO=1; IO<=mat.NPSI; IO++)
  {
    int KO=mat.IPSI[IO-1];
    XT1T=XT1T+pin00.SPT1[KO-1]*FCIN[KO-1];
    XT2T=XT2T+pin00.SPT2[KO-1]*FCIN[KO-1];
  }
  double U1AVE=0.5E0*(XT1T-XT1)/XX0;
  double U2AVE=U1AVE-(XT2T-XT2)/(6.0E0*XX0);
//
  double RMUM=0.5E0;
//
  printf("  \n");
  MOM=0;
  double SUM0=SUMGA(EDCSIA,&parin_dcsin,0.0E0,RMUM,1.0E-8);
  if(penGetError() != PEN_SUCCESS){ printf(" SUMGA error at 1e\n");}
  MOM=1;
  double U1N=SUMGA(EDCSIA,&parin_dcsin,0.0E0,RMUM,1.0E-8)/SUM0;
  if(penGetError() != PEN_SUCCESS){ printf(" SUMGA error at 2e\n");}
  MOM=2;
  double U2N=SUMGA(EDCSIA,&parin_dcsin,0.0E0,RMUM,1.0E-8)/SUM0;
  if(penGetError() != PEN_SUCCESS){ printf(" SUMGA error at 3e\n");}
  printf(" MU1AV A =%15.7E\n",U1AVE);
  printf(" MU1AV B =%15.7E\n",U1N);
  printf(" MU2AV A =%15.7E\n",U2AVE);
  printf(" MU2AV B =%15.7E\n",U2N);
//
  if(IDIST == 1)
  {

  printf("  \n");
//
//  ************  Angular distribution.
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Inelastic scattering of positrons\n");
  fprintf(IOUT," # Angular distribution\n");
  fprintf(IOUT," # Positron energy (eV) =%14.7E\n",E);
  fprintf(IOUT," #   Cutoff energy (eV) = %14.7E\n",WCCM);
  fprintf(IOUT," # Hard x-section (cm**2) =             %14.7E\n",XX0);
  fprintf(IOUT," # Hard stop. x-section (eV*cm**2) =    %14.7E\n",XX1);
  fprintf(IOUT," # Hard strg. x-section (eV**2*cm**2) = %14.7E\n",XX2);
  MOM=0;
  int NTAB=1000;
  double DX=RMUM*1.05E0/double(NTAB);
//
  for(int K=1; K<=NTAB+1; K++)
  {
    double RMU=(K-1)*DX;
    double PDF=4.0E0*constants::PI*EDCSIA(RMU,&parin_dcsin)/XX0;
    fprintf(IOUT,"%14.6E%14.6E\n",RMU,PDF);
  }
  fclose(IOUT);
//
//  ****  Random sampling test.
//
  double X1MAX=RMUM;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=500000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double W1AV=0.0E0;
  double W2AV=0.0E0;
  double W4AV=0.0E0;
  double U1AV=0.0E0;
  double U2AV=0.0E0;
  double U4AV=0.0E0;
  double W1AVMC=0.0E0,W2AVMC=0.0E0,U1AVMC=0.0E0,U2AVMC=0.0E0;
  double TIMEI=CPUtime();
  int NFAIL=0;
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
//
      double DE,EP,CDT,ES,CDTS;
      int IOSC,IZZ,ISH;
      double RR=randoms.rand();
      if(RR < RCUT)
      {
        PINa(context,mat,E,KE,XEL,DELTA,DE,EP,CDT,ES,CDTS,IOSC,randoms);
      }
      else
      {
        PSIa(context,mat,E,KE,XEL,DELTA,DE,EP,CDT,ES,CDTS,IZZ,ISH,randoms);
      }
      if(DE < E)
      {
        W1AV=W1AV+DE;
        double DE2=DE*DE;
        W2AV=W2AV+DE2;
        W4AV=W4AV+DE2*DE2;
        double RMU=(1.0E0-CDT)*0.5E0;
        U1AV=U1AV+RMU;
        double RMU2=RMU*RMU;
        U2AV=U2AV+RMU2;
        U4AV=U4AV+RMU2*RMU2;
        int INDEX=(int)((RMU-X1MIN)/DX1+1.0E0);
        if(INDEX > 0 && INDEX <= NCHAN)
        {  
          KC1[INDEX-1]=KC1[INDEX-1]+1;
        }
        else
        {
          printf(" Index out of range.\n");
        }
      }
      else
      {
        NFAIL=NFAIL+1;
      }
    }
    double FNT=double(I*LOOP)-NFAIL;
    W1AVMC=W1AV/FNT;
    W2AVMC=W2AV/FNT;
    U1AVMC=U1AV/FNT;
    U2AVMC=U2AV/FNT;
    printf(" N,<W>,<W*W>,<MU>,<MU*MU> =%9d%10.3E%10.3E%10.3E%10.3E\n",I*LOOP,W1AVMC/W1AVE,W2AVMC/W2AVE,U1AVMC/U1AVE,U2AVMC/U2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP-NFAIL);
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Positron energy = %13.5E eV\n",E);
  printf(" #   Cutoff energy = %13.5E eV\n",WCCM);
  printf(" # Total x-section (cm**2) =           %13.5E\n",XX0);
  printf(" # Stopping x-section (eV*cm**2) =     %13.5E\n",XX1);
  printf(" # Straggling x-section (eV**2*cm**2) =%13.5E\n",XX2);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #      <W>   (numerical) =%14.7E\n",W1AVE);
  double ERR1MC=3.0E0*sqrt(((W2AVMC-pow(W1AVMC,2))*FNT > 1.0E-35 ? (W2AVMC-pow(W1AVMC,2))*FNT : 1.0E-35));
  printf(" #      <W> (Monte Carlo) =%14.7E +-%9.2E\n",W1AVMC,ERR1MC);
  printf(" #   <W**2>   (numerical) =%14.7E\n",W2AVE);
  double W4AVMC=W4AV*FNT;
  double ERR2MC=3.0E0*sqrt(((W4AVMC-pow(W2AVMC,2))*FNT > 1.0E-35 ? (W4AVMC-pow(W2AVMC,2))*FNT : 1.0E-35));
  printf(" #   <W**2> (Monte Carlo) =%14.7E +-%9.2E\n",W2AVMC,ERR2MC);
  printf(" #     <MU>   (numerical) =%14.7E\n",U1AVE);
  double ERU1MC=3.0E0*sqrt(((U2AVMC-pow(U1AVMC,2))*FNT > 1.0E-35 ? (U2AVMC-pow(U1AVMC,2))*FNT : 1.0E-35));
  printf(" #     <MU> (Monte Carlo) =%14.7E +-%9.2E\n",U1AVMC,ERU1MC);
  printf(" #  <MU**2>   (numerical) =%14.7E\n",U2AVE);
  double U4AVMC=U4AV*FNT;
  double ERU2MC=3.0E0*sqrt(((U4AVMC-pow(U2AVMC,2))*FNT > 1.0E-35 ? (U4AVMC-pow(U2AVMC,2))*FNT : 1.0E-35));
  printf(" #  <MU**2> (Monte Carlo) =%14.7E +-%9.2E\n",U2AVMC,ERU2MC);

  printf(" # Angular distribution\n");
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Inelastic scattering of positrons\n");
  fprintf(IOUT," # Angular distribution\n");
  fprintf(IOUT," # Positron energy (eV) =%14.7E\n",E);
  fprintf(IOUT," #   Cutoff energy (eV) = %14.7E\n",WCCM);
  fprintf(IOUT," # Total x-section (cm**2) =             %14.7E\n",XX0);
  fprintf(IOUT," # Stopping x-section (eV*cm**2) =       %14.7E\n",XX1);
  fprintf(IOUT," # Straggling x-section (eV**2*cm**2) =%14.7E\n",XX2);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #      <W>   (numerical) =%14.7E\n",W1AVE);
  fprintf(IOUT," #      <W> (Monte Carlo) =%14.7E +-%9.2E\n",W1AVMC,ERR1MC);
  fprintf(IOUT," #   <W**2>   (numerical) =%14.7E\n",W2AVE);
  fprintf(IOUT," #   <W**2> (Monte Carlo) =%14.7E +-%9.2E\n",W2AVMC,ERR2MC);
  fprintf(IOUT," #     <MU>   (numerical) =%14.7E\n",U1AVE);
  fprintf(IOUT," #     <MU> (Monte Carlo) =%14.7E +-%9.2E\n",U1AVMC,ERU1MC);
  fprintf(IOUT," #  <MU**2>   (numerical) =%14.7E\n",U2AVE);
  fprintf(IOUT," #  <MU**2> (Monte Carlo) =%14.7E +-%9.2E\n",U2AVMC,ERU2MC);

//
  for(int I=1; I<=NCHAN; I++)
  { 
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC > 1.0E-35 ? PDFMC : 1.0E-35),
     (ERR > 1.0E-35 ? ERR : 1.0E-35));
  }
  fclose(IOUT);
  return;
  }
//
//  ************  Energy distribution.
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Inelastic scattering of positrons\n");
  fprintf(IOUT," # Energy-loss distribution\n");
  fprintf(IOUT," # Positron energy (eV) =%14.7E\n",E);
  fprintf(IOUT," #   Cutoff energy (eV) = %14.7E\n",WCCM);
  fprintf(IOUT," # Hard x-section (cm**2) =             %14.7E\n",XX0);
  fprintf(IOUT," # Hard stop. x-section (eV*cm**2) =    %14.7E\n",XX1);
  fprintf(IOUT," # Hard strg. x-section (eV**2*cm**2) = %14.7E\n",XX2);
  MOM=0;
  int NTT=500;
  double DW=(E-WCCM)/double(NTT);
  int NTAB;
  for(int I=1; I<=NTT+2; I++)
  {
    NTAB=I;
    WT[I-1]=WCCM+DW*(I-1);
  }
  for(int K=1; K<=mat.NOSC; K++)
  {
    if(mat.UI[K-1] > WCCM && mat.UI[K-1] < E)
    {
      WT[NTAB+1-1]=mat.UI[K-1]-1.0E-2;
      WT[NTAB+2-1]=mat.UI[K-1]+1.0E-2;
      NTAB=NTAB+2;
    }
  }
  for(int I=1; I<=NTAB-1; I++)
  {
    int IMIN=I;
    double WMIN=WT[I-1];
    for(int J=I+1; J<=NTAB; J++)
    {
      if(WT[J-1] < WMIN)
      {
        IMIN=J;
        WMIN=WT[J-1];
      }
    }
    double SAVE=WT[I-1];
    WT[I-1]=WT[IMIN-1];
    WT[IMIN-1]=SAVE;
  }
//
  for(int K=0; K<NTAB; K++)
  {
    double PDF=EDCSIW(WT[K],parin_dcsin)/XX0;
    fprintf(IOUT,"%14.6E%14.6E\n",WT[K],PDF);
  }
  fclose(IOUT);
//
//  ****  Random sampling test.
//
  double X1MAX=E;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=500000;
  int NLOOP=(int)(1+(NSAMP/LOOP));
  printf("NLOOP =%d\n",NLOOP);
  double W1AV=0.0E0;
  double W2AV=0.0E0;
  double W4AV=0.0E0;
  double U1AV=0.0E0;
  double U2AV=0.0E0;
  double U4AV=0.0E0;
  double TIMEI=CPUtime();
  int NFAIL=0;
  double W1AVMC=0.0E0,W2AVMC=0.0E0,U1AVMC=0.0E0,U2AVMC=0.0E0;

  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
//
      double RR=randoms.rand();
      double DE,EP,CDT,ES,CDTS;
      if(RR < RCUT)
      {
        int IOSC; 
        PINa(context,mat,E,KE,XEL,DELTA,DE,EP,CDT,ES,CDTS,IOSC,randoms);
      }
      else
      {
        int IZZ,ISH; 
        PSIa(context,mat,E,KE,XEL,DELTA,DE,EP,CDT,ES,CDTS,IZZ,ISH,randoms);
      }
      if(DE < E)
      {
        W1AV=W1AV+DE;
        double DE2=DE*DE;
        W2AV=W2AV+DE2;
        W4AV=W4AV+DE2*DE2;
        double RMU=(1.0E0-CDT)*0.5E0;
        U1AV=U1AV+RMU;
        double RMU2=RMU*RMU;
        U2AV=U2AV+RMU2;
        U4AV=U4AV+RMU2*RMU2;
        int INDEX=(int)((DE-X1MIN)/DX1+1.0E0);
        if(INDEX > 0 && INDEX <= NCHAN)
        {
          KC1[INDEX-1]=KC1[INDEX-1]+1;
        }
        else
        {
          printf(" Index out of range.\n");
        }
      }
      else
      {
        NFAIL=NFAIL+1;
      }
    }
    double FNT=double(I*LOOP)-NFAIL;
    W1AVMC=W1AV/FNT;
    W2AVMC=W2AV/FNT;
    U1AVMC=U1AV/FNT;
    U2AVMC=U2AV/FNT;
    printf(" N,<W>,<W*W>,<MU>,<MU*MU> =%9d%10.3E%10.3E%10.3E%10.3E\n",I*LOOP,W1AVMC/W1AVE,W2AVMC/W2AVE,U1AVMC/U1AVE,U2AVMC/U2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP-NFAIL);
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Positron energy = %13.5E eV\n",E);
  printf(" #   Cutoff energy = %13.5E eV\n",WCCM);
  printf(" # Total x-section (cm**2) =           %13.5E\n",XX0);
  printf(" # Stopping x-section (eV*cm**2) =     %13.5E\n",XX1);
  printf(" # Straggling x-section (eV**2*cm**2) =%13.5E\n",XX2);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #      <W>   (numerical) =%14.7E\n",W1AVE);
  double ERR1MC=3.0E0*sqrt(((W2AVMC-pow(W1AVMC,2))*FNT > 1.0E-35 ? (W2AVMC-pow(W1AVMC,2))*FNT : 1.0E-35));
  printf(" #      <W> (Monte Carlo) =%14.7E +-%9.2E\n",W1AVMC,ERR1MC);
  printf(" #   <W**2>   (numerical) =%14.7E\n",W2AVE);
  double W4AVMC=W4AV*FNT;
  double ERR2MC=3.0E0*sqrt(((W4AVMC-pow(W2AVMC,2))*FNT > 1.0E-35 ? (W4AVMC-pow(W2AVMC,2))*FNT : 1.0E-35));
  printf(" #   <W**2> (Monte Carlo) =%14.7E +-%9.2E\n",W2AVMC,ERR2MC);
  printf(" #     <MU>   (numerical) =%14.7E\n",U1AVE);
  double ERU1MC=3.0E0*sqrt(((U2AVMC-pow(U1AVMC,2))*FNT > 1.0E-35 ? (U2AVMC-pow(U1AVMC,2))*FNT : 1.0E-35));
  printf(" #     <MU> (Monte Carlo) =%14.7E +-%9.2E\n",U1AVMC,ERU1MC);
  printf(" #  <MU**2>   (numerical) =%14.7E\n",U2AVE);
  double U4AVMC=U4AV*FNT;
  double ERU2MC=3.0E0*sqrt(((U4AVMC-pow(U2AVMC,2))*FNT > 1.0E-35 ? (U4AVMC-pow(U2AVMC,2))*FNT : 1.0E-35));
  printf(" #  <MU**2> (Monte Carlo) =%14.7E +-%9.2E\n",U2AVMC,ERU2MC);

  printf(" # Energy-loss distribution\n");
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Inelastic scattering of positrons\n");
  fprintf(IOUT," # Energy-loss distribution\n");
  fprintf(IOUT," # Positron energy (eV) =%14.7E\n",E);
  fprintf(IOUT," #   Cutoff energy (eV) = %14.7E\n",WCCM);
  fprintf(IOUT," # Total x-section (cm**2) =           %14.7E\n",XX0);
  fprintf(IOUT," # Stopping x-section (eV*cm**2) =     %14.7E\n",XX1);
  fprintf(IOUT," # Straggling x-section (eV**2*cm**2) =%14.7E\n",XX2);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #      <W>   (numerical) =%14.7E\n",W1AVE);
  ERR1MC=3.0E0*sqrt(((W2AVMC-pow(W1AVMC,2))*FNT > 1.0E-35 ? (W2AVMC-pow(W1AVMC,2))*FNT : 1.0E-35));
  fprintf(IOUT," #      <W> (Monte Carlo) =%14.7E +-%9.2E\n",W1AVMC,ERR1MC);
  fprintf(IOUT," #   <W**2>   (numerical) =%14.7E\n",W2AVE);
  W4AVMC=W4AV*FNT;
  ERR2MC=3.0E0*sqrt(((W4AVMC-pow(W2AVMC,2))*FNT > 1.0E-35 ? (W4AVMC-pow(W2AVMC,2))*FNT : 1.0E-35));
  fprintf(IOUT," #   <W**2> (Monte Carlo) =%14.7E +-%9.2E\n",W2AVMC,ERR2MC);
  fprintf(IOUT," #     <MU>   (numerical) =%14.7E\n",U1AVE);
  ERU1MC=3.0E0*sqrt(((U2AVMC-pow(U1AVMC,2))*FNT > 1.0E-35 ? (U2AVMC-pow(U1AVMC,2))*FNT : 1.0E-35));
  fprintf(IOUT," #     <MU> (Monte Carlo) =%14.7E +-%9.2E\n",U1AVMC,ERU1MC);
  fprintf(IOUT," #  <MU**2>   (numerical) =%14.7E\n",U2AVE);
  U4AVMC=U4AV*FNT;
  ERU2MC=3.0E0*sqrt(((U4AVMC-pow(U2AVMC,2))*FNT > 1.0E-35 ? (U4AVMC-pow(U2AVMC,2))*FNT : 1.0E-35));
  fprintf(IOUT," #  <MU**2> (Monte Carlo) =%14.7E +-%9.2E\n",U2AVMC,ERU2MC);
//
  for(int I=1; I<=NCHAN; I++)
  {
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],(PDFMC > 1.0E-35 ? PDFMC : 1.0E-35),(ERR > 1.0E-35 ? ERR : 1.0E-35));
  }
  fclose(IOUT);
}

//  *********************************************************************
//                       SUBROUTINE ESIaV
//  *********************************************************************
void ESIaV(const pen_context& context, const pen_material& mat, const double E, const CESI0& esio, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling of the active atomic electron
//  shell in inner shell ionization by electron impact.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Element data.
//  using namespace CADATA;
//  ****  Composition data.
//  using namespace COMPOS;
//  ****  E/P inelastic collisions.
//  using namespace CEIN;
//  ****  Electron simulation tables.
//  using namespace CEIMFP;
//  ****  Electron inelastic coll. and inner-shell ionization tables.
//  using namespace CEINAC;
//  using namespace CESIAC;
//  using namespace CESIN;
//  using namespace CESI0;
//
  double PROB[30][16];
  int NCOUNT[30][16],ISTORE[99];
//
//  ****  Initialisation.
//
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL];
    for(int ISH=1; ISH<=esio.NSESI[IZZ-1]; ISH++)
    {
      PROB[IEL][ISH-1]=0.0E0;
    }
  }
//
  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);

//
//     DELTA=DEL(M,KE)+(DEL(M,KE+1)-DEL(M,KE))*XEK
//
  double PTOT=0.0E0;
  double WCSI=(mat.EABS[0] < mat.EABS[1] ? mat.EABS[0] : mat.EABS[1]);
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL];
    ISTORE[IZZ-1]=IEL;
    int INDC=esio.IESIF[IZZ-1]-1;
    for(int ISH=1; ISH<=esio.NSESI[IZZ-1]; ISH++)
    {
      double WCUT=context.elements.EB[IZZ-1][ISH-1];
      if(WCUT > WCSI && WCUT < E)
      {
        double PCSI=exp(esio.XESI[INDC+KE][ISH-1]
                +(esio.XESI[INDC+KE+1][ISH-1]-esio.XESI[INDC+KE][ISH-1])*XEK);
        if(PCSI > 1.1E-35)
        {
          PROB[IEL][ISH-1]=PCSI*mat.STF[IEL];
          PTOT=PTOT+PROB[IEL][ISH-1];
          printf(" Z =%3d, ISH =%3d, CS =%13.6E CM**2\n",IZZ,ISH,PCSI); 
        }
      }
    }
  }
  if(PTOT < 1.0E-35)
  {
    printf(" Sorry, this material has no inner shells.\n");
    return;
  }
//
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    for(unsigned int ISH=1; ISH<=16; ISH++)
    {
      NCOUNT[IEL][ISH-1]=0;
    }
  }
//
//  ****  Simulation starts here.
//
  int LOOP=500000;
  int NLOOP=NSAMP/LOOP;
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=0; J<LOOP; J++)
    {
//    ESIa(E,DELTA,DE,EP,CDT,ES,CDTS,M,IZZ,ISH)
      int IZZ, ISH;
      ESIb(mat,KE,XEL,context.grid.DLEMP,context.grid.DLFC,IZZ,ISH,randoms);
      if(IZZ > 0)
      {
        int IEL=ISTORE[IZZ-1];
        if(ISH < 1 || ISH > esio.NSESI[IZZ-1])
        {
          printf("The sampling has given a wrong shell.\n");
        }
        NCOUNT[IEL][ISH-1]=NCOUNT[IEL][ISH-1]+1;
      }
    }
    int NSAMPLE=I*LOOP;
    double RN=double(NSAMPLE);
    printf("N =%9d\n",I*LOOP);
//  ****  The shell ionization probabilities obtained from the simulation
//        are compared with the probabilities obtained from the tabulated
//        ionization x-sections.
    for(int IEL=0; IEL<mat.NELEM; IEL++)
    {
      int IZZ=mat.IZ[IEL];
      for(int ISH=1; ISH<=esio.NSESI[IZZ-1]; ISH++)
      {
        double PMC=NCOUNT[IEL][ISH-1]/RN;
        double ERR=3.0E0*sqrt(PMC*(1.0E0-PMC)/RN);
        double PINP=PROB[IEL][ISH-1]/PTOT;
        if(PINP > 1.0E-15)
        { 
          printf(" %4d%4d%13.5E%13.5E%9.1E%13.5E\n",mat.IZ[IEL],ISH,PINP,PMC,ERR,100.0E0*fabs(1.0E0-PMC/PINP));
        }
      }
    }
  }
  int NSAMPLE=NLOOP*LOOP;
  double RN=double(NSAMPLE);
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Shell ionization probabilities\n");
  fprintf(IOUT," # Projectile energy (eV) = %14.7E\n",E);
  int ITR=0;
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL];
    for(int ISH=1; ISH<=esio.NSESI[IZZ-1]; ISH++)
    {
      ITR=ITR+1;
      double PINP=PROB[IEL][ISH-1]/PTOT;
      if(PINP > 1.0E-15){ fprintf(IOUT,"%4d%13.5E  %4d%4d\n",ITR,PINP,mat.IZ[IEL],ISH);}
    }
  }
  fclose(IOUT);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Shell ionization probabilities\n");
  fprintf(IOUT," # Projectile energy (eV) = %14.7E\n",E);
  ITR=0;
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL];
    for(int ISH=1; ISH<=esio.NSESI[IZZ-1]; ISH++)
    {
      ITR=ITR+1;
      double PMC=NCOUNT[IEL][ISH-1]/RN;
      double ERR=3.0E0*sqrt(PMC*(1.0E0-PMC)/RN);
      if(PMC > 0.0E0){ fprintf(IOUT,"%4d%13.5E%13.5E\n",ITR,PMC,ERR);}
    }
  }
  fclose(IOUT);
}

//  *********************************************************************
//                       SUBROUTINE ESIb
//  *********************************************************************
void ESIb(const pen_material& mat,
    const int KE,
    const double XEL,
    const double* DLEMP,
    const double DLFC,
    int &IZZ,
    int &ISH,
    pen_rand& penRand)
{
  //  Random sampling of the active sell in inner-shell ionisation by
  //  electron impact.
  //
  //  Input arguments:
  //    M ....... material where electrons propagate.
  //  Output arguments:
  //    IZZ ..... atomic number of the element that has been ionized.
  //    ISH ..... atomic electron shell that has been ionized.
  
  //  using namespace PENELOPE_mod;

  //  using namespace CEGRID;
  //  using namespace CEIN;
  //  using namespace CESIAC;

  //  ****  Energy grid point.
  double PK = (XEL-DLEMP[KE])*DLFC;
  int JE;
  if(penRand.rand() < PK)
  {
    JE = KE+1;
  }
  else
  {
    JE = KE;
  }

  //  ************  Selection of the active oscillator.

  double TST = penRand.rand();
  //  ****  Binary search.
  int IO = 0;
  int JO = mat.NESI;
  
  while(JO-IO > 1)
  {
    int IT = (IO+JO)/2;
    if(TST > mat.ESIAC[JE][IT])
    {
      IO = IT;
    }
    else
    {
      JO = IT;
    }
  }
  int IOSC = mat.IESI[IO];

  IZZ = mat.KZ[IOSC-1];
  ISH = mat.KS[IOSC-1];
}

//  *********************************************************************
//                       SUBROUTINE PSIb
//  *********************************************************************
void PSIb(const pen_material& mat,
    const int KE,
    const double XEL,
    const double* DLEMP,
    const double DLFC,
    int &IZZ,
    int &ISH,
    pen_rand& penRand)
{
  //  Random sampling of the active sell in inner-shell ionisation by
  //  positron impact.
  //
  //  Input arguments:
  //    M ....... material where electrons propagate.
  //  Output arguments:
  //    IZZ ..... atomic number of the element that has been ionized.
  //    ISH ..... atomic electron shell that has been ionized.
  
  //  using namespace PENELOPE_mod;

  //  using namespace CEGRID;
  //  using namespace CEIN;
  //  using namespace CESIAC;

  //  ****  Energy grid point.
  double PK = (XEL-DLEMP[KE])*DLFC;
  int JE;
  if(penRand.rand() < PK)
  {
    JE = KE+1;
  }
  else
  {
    JE = KE;
  }

  //  ************  Selection of the active oscillator.

  double TST = penRand.rand();
  //  ****  Binary search.
  int IO = 0;
  int JO = mat.NPSI;
  
  while(JO-IO > 1)
  {
    int IT = (IO+JO)/2;
    if(TST > mat.PSIAC[JE][IT])
    {
      IO = IT;
    }
    else
    {
      JO = IT;
    }
  }
  int IOSC = mat.IPSI[IO];

  IZZ = mat.KZ[IOSC-1];
  ISH = mat.KS[IOSC-1];
}

//  *********************************************************************
//                       SUBROUTINE PSIaV
//  *********************************************************************
void PSIaV(const pen_context& context, const pen_material& mat, const double E, const CPSI0& psio, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling of the active atomic electron
//  shell in inner shell ionization by positron impact.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Element data.
//  using namespace CADATA;
//  ****  Composition data.
//  using namespace COMPOS;
//  ****  E/P inelastic collisions.
//  using namespace CEIN;
//  ****  Electron simulation tables.
//  using namespace CEIMFP;
//  ****  Positron simulation tables.
//  using namespace CPIMFP;
//  ****  Positron inelastic coll. and inner-shell ionization tables.
//  using namespace CPINAC;
//  using namespace CPSIAC;
//  using namespace CPSIN;
//  using namespace CPSI0;
//
  double PROB[30][16];
  int NCOUNT[30][16],ISTORE[99];
//
//  ****  Initialisation.
//
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL];
    for(int ISH=1; ISH<=psio.NSPSI[IZZ-1]; ISH++)
    {
      PROB[IEL][ISH-1]=0.0E0;
    }
  }

  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);


  double DELTA=mat.DEL[KE]+(mat.DEL[KE+1]-mat.DEL[KE])*XEK;

  double PTOT=0.0E0;
  double WCSI=(mat.EABS[0] < mat.EABS[1] ? mat.EABS[0] : mat.EABS[1]);
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL];
    ISTORE[IZZ-1]=IEL;
    int INDC=psio.IPSIF[IZZ-1]-1;
    for(int ISH=1; ISH<=psio.NSPSI[IZZ-1]; ISH++)
    {
      double WCUT=context.elements.EB[IZZ-1][ISH-1];
      if(WCUT > WCSI && WCUT < E)
      {
        double PCSI=exp(psio.XPSI[INDC+KE][ISH-1]
                +(psio.XPSI[INDC+KE+1][ISH-1]-psio.XPSI[INDC+KE][ISH-1])*XEK);
        if(PCSI > 1.1E-35)
        {
          PROB[IEL][ISH-1]=PCSI*mat.STF[IEL];
          PTOT=PTOT+PROB[IEL][ISH-1];
          printf(" Z =%3d, ISH =%3d, CS =%13.6E CM**2\n",IZZ,ISH,PCSI); 
        }
      }
    }
  }
  if(PTOT < 1.0E-35)
  {
    printf(" Sorry, this material has no inner shells.\n");
    return;
  }
//
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    for(int ISH=1; ISH<=16; ISH++)
    {
      NCOUNT[IEL][ISH-1]=0;
    }
  }
//
//  ****  Simulation starts here.
//
  int LOOP=500000;
  int NLOOP=NSAMP/LOOP;
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      double DE,EP,CDT,ES,CDTS;
      int IZZ, ISH;
      PSIa(context,mat,E,KE,XEL,DELTA,DE,EP,CDT,ES,CDTS,IZZ,ISH,randoms);
      if(IZZ > 0)
      {
        int IEL=ISTORE[IZZ-1];
        if(ISH < 1 || ISH > psio.NSPSI[IZZ-1])
        {
          printf("The sampling has given a wrong shell.\n");
        }
        NCOUNT[IEL][ISH-1]=NCOUNT[IEL][ISH-1]+1;
      }
    }
    int NSAMPLE=I*LOOP;
    double RN=double(NSAMPLE);
    printf("N =%9d\n",I*LOOP);
//  ****  The shell ionization probabilities obtained from the simulation
//        are compared with the probabilities obtained from the tabulated
//        ionization x-sections.
    for(int IEL=0; IEL<mat.NELEM; IEL++)
    {
      int IZZ=mat.IZ[IEL];
      for(int ISH=1; ISH<=psio.NSPSI[IZZ-1]; ISH++)
      {
        double PMC=NCOUNT[IEL][ISH-1]/RN;
        double ERR=3.0E0*sqrt(PMC*(1.0E0-PMC)/RN);
        double PINP=PROB[IEL][ISH-1]/PTOT;
        if(PINP > 1.0E-15)
        { 
          printf(" %4d%4d%13.5E%13.5E%9.1E%13.5E\n",mat.IZ[IEL],ISH,PINP,PMC,ERR,100.0E0*fabs(1.0E0-PMC/PINP));
        }
      }
    }
  }
  int NSAMPLE=NLOOP*LOOP;
  double RN=double(NSAMPLE);
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Shell ionization probabilities\n");
  fprintf(IOUT," # Projectile energy (eV) = %14.7E\n",E);
  int ITR=0;
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL];
    for(int ISH=1; ISH<=psio.NSPSI[IZZ-1]; ISH++)
    {
      ITR=ITR+1;
      double PINP=PROB[IEL][ISH-1]/PTOT;
      if(PINP > 1.0E-15){ fprintf(IOUT,"%4d%13.5E  %4d%4d\n",ITR,PINP,mat.IZ[IEL],ISH);}
    }
  }
  fclose(IOUT);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Shell ionization probabilities\n");
  fprintf(IOUT," # Projectile energy (eV) = %14.7E\n",E);
  ITR=0;
  for(int IEL=0; IEL<mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL];
    for(int ISH=1; ISH<=psio.NSPSI[IZZ-1]; ISH++)
    {
      ITR=ITR+1;
      double PMC=NCOUNT[IEL][ISH-1]/RN;
      double ERR=3.0E0*sqrt(PMC*(1.0E0-PMC)/RN);
      if(PMC > 0.0E0){ fprintf(IOUT,"%4d%13.5E%13.5E\n",ITR,PMC,ERR);}
    }
  }
  fclose(IOUT);
}

//  *********************************************************************
//                       SUBROUTINE EBRaV
//  *********************************************************************
void EBRaV(const pen_context& context, const pen_material& mat, const double E, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm for electron bremss
//  emission.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
  const double TREV=2.0E0*constants::REV;
  const unsigned& NBW = constants::NBW;
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Bremsstrahlung emission.
//  using namespace CEBR;
//  using namespace CEBR01;
//  using namespace CEBR02;
//
  double XC1[1000];
  int KC1[1000];
  CEBR01 ebr01;

  double  (&XX)[constants::NBE]=ebr01.X;
  double  (&YY)[constants::NBE]=ebr01.Y;

  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);

//  ****  'Exact' distribution.
//
  double BI1=pow(context.grid.ET[KE]+constants::REV,2)/(context.grid.ET[KE]*(context.grid.ET[KE]+TREV));
  double BI2=pow(context.grid.ET[KE+1]+constants::REV,2)/(context.grid.ET[KE+1]*(context.grid.ET[KE+1]+TREV));
  double FACT=mat.ZBR2*1.0E-27;
  double WCRE=mat.WCR/E;
  for(unsigned IW=0; IW < NBW; IW++)
  {
    XX[IW]=constants::WB[IW];
    YY[IW]=mat.P0[KE][IW];
  }
  double XBR0A=RLMOM(XX,YY,XX[NBW-1],NBW,-1)-RLMOM(XX,YY,WCRE,NBW,-1);
  double XBR1A=RLMOM(XX,YY,XX[NBW-1],NBW,0)-RLMOM(XX,YY,WCRE,NBW,0);
  double XBR2A=RLMOM(XX,YY,XX[NBW-1],NBW,1)-RLMOM(XX,YY,WCRE,NBW,1);
  for(unsigned IW=0; IW<NBW; IW++)
  {
    YY[IW]=mat.P0[KE+1][IW];
  }
  double XBR0B=RLMOM(XX,YY,XX[NBW-1],NBW,-1)-RLMOM(XX,YY,WCRE,NBW,-1);
  double XBR1B=RLMOM(XX,YY,XX[NBW-1],NBW,0)-RLMOM(XX,YY,WCRE,NBW,0);
  double XBR2B=RLMOM(XX,YY,XX[NBW-1],NBW,1)-RLMOM(XX,YY,WCRE,NBW,1);
  double XBR0=((1.0E0-XEK)*BI1*XBR0A+XEK*BI2*XBR0B)*FACT;
  double XBR1=((1.0E0-XEK)*BI1*XBR1A+XEK*BI2*XBR1B)*FACT*E;
  double XBR2=((1.0E0-XEK)*BI1*XBR2A+XEK*BI2*XBR2B)*FACT*E*E;
  printf(" # Bremsstrahlung energy-loss spectrum\n");
  printf(" # Electron energy  (eV) = %14.7E\n",E);
  printf(" # Cutoff energy  (eV) = %14.7E\n",mat.WCR);
  printf(" # Total x-section (cm**2)            =%14.7E\n",XBR0);
  printf(" # Stopping x-section (eV*cm**2)      =%14.7E\n",XBR1);
  printf(" # Straggling x-section (ev**2*cm**2) =%14.7E\n",XBR2);
//
  for(unsigned IW=0; IW<NBW; IW++)
  {
    XX[IW]=constants::WB[IW];
    YY[IW]=FACT*((1.0E0-XEK)*BI1*mat.P0[KE][IW]+XEK*BI2*mat.P0[KE+1][IW]);
  }
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Bremsstrahlung energy-loss spectrum\n");
  fprintf(IOUT," # Electron energy  (eV) = %14.7E\n",E);
  fprintf(IOUT," # Cutoff energy  (eV) = %14.7E\n",mat.WCR);
  fprintf(IOUT," # Total x-section (cm**2)           =%14.7E\n",XBR0);
  fprintf(IOUT," # Stopping x-section (eV*cm**2)     =%14.7E\n",XBR1);
  fprintf(IOUT," # Stragling x-section (ev**2*cm**2) =%14.7E\n",XBR2);

  double XL=(10.0E0>mat.WCR-1.0E0 ? 10.0E0 : mat.WCR-1.0E0)/E;
  int NNP=500;
  double XU=1.0E0;
  double DX=(XU-XL)/double(NNP+1);
  XL=XL-DX;
  for(int I=1; I<=NNP; I++)
  {
    XL=XL+DX;
    unsigned J;
    FINDI(XX,XL,NBW,J);
    J=(J < NBW-2 ? J : NBW-2);
    double X1=XX[J];
    double X2=XX[J+1];
    double Y1=YY[J];
    double Y2=YY[J+1];
    double FC=(Y1+(Y2-Y1)*(XL-X1)/(X2-X1))/XBR0;
    fprintf(IOUT,"%14.6E%14.6E\n",XL*E,FC/(XL*E));
  }
  fclose(IOUT);
  if(NSAMP < 1){ return;}
//  ****  'Exact' moments.
  if(XBR0 < 1.0E-40){ return;}
  double X1AVE=XBR1/XBR0;
  double X2AVE=XBR2/XBR0;
//
//  ************  Random sampling test.
//
  double X1MAX=E;
  double X1MIN=mat.WCR;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=1000000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double TIMEI=CPUtime();
  double X1AVMC,X2AVMC;
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      double DE;
      EBRa(mat,E,KE,XEK,DE,constants::WB,randoms);
      X1AV=X1AV+DE;
      double DE2=DE*DE;
      X2AV=X2AV+DE2;
      X4AV=X4AV+DE2*DE2;
      int INDEX=(int)((DE-X1MIN)/DX1+1.0E0);
      if(INDEX > 0 && INDEX <= NCHAN)
      {
        KC1[INDEX-1]=KC1[INDEX-1]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    X1AVMC=X1AV/FNT;
    X2AVMC=X2AV/FNT;
    printf(" N,<MU>,<MU**2> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Electron energy = %13.5E eV\n",E);
  printf(" #   Cutoff energy = %13.5E eV\n",mat.WCR);
  printf(" # Total x-section (cm**2) =           %13.5E\n",XBR0);
  printf(" # Stopping x-section (eV*cm**2) =     %13.5E\n",XBR1);
  printf(" # Straggling x-section (eV**2*cm**2) =%13.5E\n",XBR2);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Bremsstrahlung energy-loss spectrum\n");
  fprintf(IOUT," # Electron energy (eV)  =%13.5E\n",E);
  fprintf(IOUT," # Cutoff energy = %13.5E eV\n",mat.WCR);
  fprintf(IOUT," # Total x-section (cm**2) =           %13.5E\n",XBR0);
  fprintf(IOUT," # Stopping x-section (eV*cm**2) =     %13.5E\n",XBR1);
  fprintf(IOUT," # Straggling x-section (eV**2*cm**2) =%13.5E\n",XBR2);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  for(int I=1; I <= NCHAN; I++)
  {
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC > 1.0E-35 ? PDFMC : 1.0E-35),
     (ERR > 1.0E-35 ? ERR : 1.0E-35));
  }
  fclose(IOUT);
  return;
}

//  *********************************************************************
//                       SUBROUTINE EBRaAV
//  *********************************************************************
void EBRaAV(const pen_material& mat, const double E, double DE,int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm for the initial dir-
//  ection of bremms photons.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//  ****  Bremsstrahlung emission.
//  using namespace CEBR;
//
//  using namespace CEBRA0;

  const double TREV = 2.0*constants::REV;
  
  double Ebet[7], BET[7];

  Ebet[0] = 1.0E3;
  Ebet[1] = 5.0E3;
  Ebet[2] = 1.0E4;
  Ebet[3] = 5.0E4;
  Ebet[4] = 1.0E5;
  Ebet[5] = 5.0E5;
  Ebet[6] = 1.0E6;

  for(int IE = 0; IE < 7; IE++)
    {
      BET[IE] = sqrt(Ebet[IE]*(Ebet[IE]+TREV))/(Ebet[IE]+constants::REV);
    }  
  
//  ****  Auxiliary arrays for random sampling verification.
  double XC1[1000];
  int KC1[1000];
  CEBRA0 ebra0;
//
  ebra0.EE=E;
  ebra0.DEE=DE;
  ebra0.pmat = &mat;
//
  ebra0.MOM=0;
  double TX0=SUMGA(EBRADX,&ebra0,0.0E0,1.0E0,1.0E-8);
  ebra0.MOM=1;
  double TX1=SUMGA(EBRADX,&ebra0,0.0E0,1.0E0,1.0E-8);
  ebra0.MOM=2;
  double TX2=SUMGA(EBRADX,&ebra0,0.0E0,1.0E0,1.0E-8);
  double X1AVE=TX1/TX0;
  double X2AVE=TX2/TX0;
  printf(" # Electron energy = %13.5E eV\n",E);
  printf(" #   Photon energy = %13.5E eV\n",DE);
  printf(" # Normalisation integral = %13.5E\n",TX0);
  printf(" # Exp(mu)    = %13.5E\n",X1AVE);
  printf(" # Exp(mu**2) = %13.5E\n",X2AVE);
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Bremsstrahlung angular distribution\n");
  fprintf(IOUT," # Electron energy  (eV) =%14.7E\n",E);
  fprintf(IOUT," # Photon energy  (eV) =%14.7E\n",DE);
  fprintf(IOUT," # Normalisation integral = %14.7E\n",TX0);
  fprintf(IOUT," # Exp(mu)    = %14.7E\n",X1AVE);
  fprintf(IOUT," # Exp(mu**2) = %14.7E\n",X2AVE);
  double DX=1.0E0/200.0E0;
  ebra0.MOM=0;
  for(int K=1; K<=201; K++)
  {
    double RMU=(K-1)*DX;
    double CDT=1.0E0-2.0E0*RMU;
    double PDF=2.0E0*EBRDA(*(ebra0.pmat),ebra0.EE,ebra0.DEE,CDT)/TX0;
    double PDF2=PDF/(4.0E0*constants::PI);
    fprintf(IOUT,"%14.6E%14.6E%14.6E%14.6E\n",RMU,PDF,PDF2,acos(CDT)*180/constants::PI);
  }
  fclose(IOUT);
//
//  ************  Random sampling test.
//
  double X1MAX=1.0E0;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=250000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      double CDT;
      EBRaA(*(ebra0.pmat),ebra0.EE,ebra0.DEE,CDT,BET,randoms);
      double RMU=(1.0E0-CDT)*0.5E0;
      X1AV=X1AV+RMU;
      double XX2=RMU*RMU;
      X2AV=X2AV+XX2;
      X4AV=X4AV+XX2*XX2;
      int INDEX=(int)((RMU-X1MIN)/DX1+1.0E0);
      if(INDEX > 0 && INDEX <= NCHAN)
      {
        KC1[INDEX-1]=KC1[INDEX-1]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    double X1AVMC=X1AV/FNT;
    double X2AVMC=X2AV/FNT;
    printf(" N,<MU>,<MU**2> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  double X1AVMC=X1AV*FNT;
  double X2AVMC=X2AV*FNT;
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Electron energy = %13.5E eV\n",E);
  printf(" #   Photon energy = %13.5E eV\n",DE);
  printf(" # Normalisation integral = %13.5E\n",TX0);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Bremsstrahlung angular distribution\n");
  fprintf(IOUT," # Electron energy = %13.5E eV\n",E);
  fprintf(IOUT," #   Photon energy = %13.5E eV\n",DE);
  fprintf(IOUT," # Normalisation integral = %13.5E\n",TX0);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  for(int I=1; I<=NCHAN; I++)
  { 
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC > 1.0E-35 ? PDFMC : 1.0E-35),
     (ERR > 1.0E-35 ? ERR : 1.0E-35));
  }
  fclose(IOUT);
}

//  *********************************************************************
//                       SUBROUTINE EBRADX
//  *********************************************************************
double EBRADX(double T, void* arg)
{
//
//  Angular distribution of bremss photons, relative to the direction of
//  the projectile.   T=(1-CDT)/2
//
//  using namespace CEBRA0;
  CEBRA0* pcebra0 = (CEBRA0*)arg;
  double EBRADX_RETURN;
  //
  double CDT=1.0E0-2.0E0*T;
  EBRADX_RETURN=EBRDA(*(pcebra0->pmat),pcebra0->EE,pcebra0->DEE,CDT)*2.0E0;
  if(pcebra0->MOM > 0){ EBRADX_RETURN=EBRADX_RETURN*pow(T,pcebra0->MOM);}
  return EBRADX_RETURN;
}

//  *********************************************************************
//                       FUNCTION EBRDA
//  *********************************************************************
double EBRDA(const pen_material& mat, double E, double DE, double CDT)
{
//
//  Probability density function of cos(theta) in bremsstrahlung events.
//  Numerical fit/interpolation of partial-wave shape functions generated
//  by the program BREMS of A. Poskus, Comp. Phys. Commun. (2018).
//
//  Input parameters:
//    M ..... material where the projectile moves.
//    E .... kinetic energy of the projectile.
//    DE ... energy of the emitted photon.
//    CDT ... cosine of the polar emission angle.
//  Output value:
//    EBRDA ... normalised PDF of CDT=cos(theta).
//
//  ****  Bremsstrahlung angular distributions.
//  using namespace CBRANG;
  const double TREV = 2.0*constants::REV;
  double EBRDA_RETURN;

  double Ebet[7], BET[7];

  Ebet[0] = 1.0E3;
  Ebet[1] = 5.0E3;
  Ebet[2] = 1.0E4;
  Ebet[3] = 5.0E4;
  Ebet[4] = 1.0E5;
  Ebet[5] = 5.0E5;
  Ebet[6] = 1.0E6;

  for(int IE = 0; IE < 7; IE++)
    {
      BET[IE] = sqrt(Ebet[IE]*(Ebet[IE]+TREV))/(Ebet[IE]+constants::REV);
    }  

//
//  ****  Distribution parameters.
//
  double BETA=sqrt(E*(E+TREV))/(E+constants::REV);
  if(E>1.0E6)
  {
    double RN1=3.0E0*(1.0E0-pow(BETA,2))/8.0E0;
    double XP=(CDT-BETA)/(1.0E0-BETA*CDT);
    EBRDA_RETURN=RN1*(1.0E0+pow(XP,2))/pow(1.0E0-BETA*CDT,2);
  }
  else
  {
    if(BETA < BET[0]){BETA = BET[0];}
    if(BETA > BET[mat.NET-1]){BETA = BET[mat.NET-1];}
//
    unsigned IE;
    FINDI(BET,BETA,mat.NET,IE);
    double RK=1.0E0+40.0E0*DE/E;
    int IK=(int(RK)<40 ? int(RK) : 40);
//
    double P10=mat.BP1[IE][IK-1][1-1]+BETA*(mat.BP1[IE][IK-1][2-1]
           +BETA*(mat.BP1[IE][IK-1][3-1]+BETA*mat.BP1[IE][IK-1][4-1]));
    double P11=mat.BP1[IE][IK+1-1][1-1]+BETA*(mat.BP1[IE][IK+1-1][2-1]
           +BETA*(mat.BP1[IE][IK+1-1][3-1]+BETA*mat.BP1[IE][IK+1-1][4-1]));
    double P1=exp(P10+(RK-IK)*(P11-P10));
    if(P1 > 1.0E0){P1 = 1.0E0;}
//
    double P20=mat.BP2[IE][IK-1][1-1]+BETA*(mat.BP2[IE][IK-1][2-1]
           +BETA*(mat.BP2[IE][IK-1][3-1]+BETA*mat.BP2[IE][IK-1][4-1]));
    double P21=mat.BP2[IE][IK+1-1][1-1]+BETA*(mat.BP2[IE][IK+1-1][2-1]
           +BETA*(mat.BP2[IE][IK+1-1][3-1]+BETA*mat.BP2[IE][IK+1-1][4-1]));
    double P2=P20+(RK-IK)*(P21-P20);
//
//  ****  Normalised probability density function.
//
    double RN1=3.0E0/8.0E0;
    double RN2=3.0E0/4.0E0;
    double XP2=pow((CDT-P2)/(1.0E0-P2*CDT),2);
    EBRDA_RETURN=(P1*RN1*(1.0E0+XP2)+(1.0E0-P1)*RN2*(1.0E0-XP2))
     *(1.0E0-pow(P2,2))/pow(1.0E0-P2*CDT,2);
  }
  return EBRDA_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE PANaV
//  *********************************************************************
void PANaV(const pen_material &mat, const double E, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm for positron annihil-
//  ation.
//
//  using namespace TRACK_mod;
//
  const double TREV=2.0E0*constants::REV;
//
//  using namespace CPAN01;
  CPAN01 pan01;
//  ****  Auxiliary arrays for random sampling verification.
  double XC1[1000];
  int KC1[1000];
//
//
//  ****  Moments of the 'exact' distribution.
//
  pan01.GAM=1.0E0+(E>1.0E0 ? E : 1.0E0)/constants::REV;
  double CHIMIN=1.0E0/(pan01.GAM+1.0E0+sqrt(pow(pan01.GAM,2)-1.0E0));
  double CHIMAX=1.0E0-CHIMIN;
  pan01.MOM=0;
  double TX0N=SUMGA(PANDX,&pan01,CHIMIN,CHIMAX,1.0E-9);
  double TX0;
  PANaT(E,TX0);
  pan01.MOM=1;
  double TX1N=SUMGA(PANDX,&pan01,CHIMIN,CHIMAX,1.0E-9);
  pan01.MOM=2;
  double TX2N=SUMGA(PANDX,&pan01,CHIMIN,CHIMAX,1.0E-9);
  double X1AVE=TX1N/TX0N;
  double X2AVE=TX2N/TX0N;
  printf(" # Positron energy = %13.5E eV\n",E);
  printf(" # Annihilation x-section (numr.,cm**2) = %14.7E\n",TX0N);
  printf(" # Annihilation x-section (anal.,cm**2) = %14.7E\n",TX0);
  printf(" # Exp(chi)    = %14.7E\n",X1AVE);
  printf(" # Exp(chi**2) = %14.7E\n",X2AVE);
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Positron annihilation\n");
  fprintf(IOUT," # Positron energy = %13.5E eV\n",E);
  fprintf(IOUT," # Annihilation x-section (numr.,cm**2) = %14.7E\n",TX0N);
  fprintf(IOUT," # Annihilation x-section (anal.,cm**2) = %14.7E\n",TX0);
  fprintf(IOUT," # Exp(chi)    = %14.7E\n",X1AVE);
  fprintf(IOUT," # Exp(chi**2) = %14.7E\n",X2AVE);
  double DX=1.0E0/200.0E0;
  pan01.MOM=0;
  for(int K=1; K<=201; K++)
  {
    double CHI=(K-1)*DX;
    double PDF;
    if(CHI>CHIMIN && CHI<CHIMAX)
    {
      PDF=PANDX(CHI,&pan01);
    }
    else
    {
      PDF=0.0E0;
    }
    fprintf(IOUT,"%14.6E%14.6E\n",CHI,PDF);
  }
  fclose(IOUT);
//
//  ************  Random sampling test.
//
  double X1MAX=1.0E0;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=1000000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double FACT=1.0E0/(E+TREV);
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      double E1,CDT1,E2,CDT2;
      PANa(mat.EABS[2],E,E1,CDT1,E2,CDT2,randoms);
      double CHI=E1*FACT;
      X1AV=X1AV+CHI;
      double XX2=CHI*CHI;
      X2AV=X2AV+XX2;
      X4AV=X4AV+XX2*XX2;
      int INDEX=(int)((CHI-X1MIN)/DX1+1.0E0);
      if(INDEX>0 && INDEX<=NCHAN)
      {
        KC1[INDEX-1]=KC1[INDEX-1]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    double X1AVMC=X1AV/FNT;
    double X2AVMC=X2AV/FNT;
    printf(" N,<X>,<X*X> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  double X1AVMC=X1AV*FNT;
  double X2AVMC=X2AV*FNT;
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Positron energy = %13.5E eV\n",E);
  printf(" # Annihilation x-section (numr.,cm**2) = %13.5E\n",TX0N);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Positron annihilation\n");
  fprintf(IOUT," # Positron energy (eV)  =%13.5E\n",E);
  fprintf(IOUT," # Annihilation x-section (numr.,cm**2) = %13.5E\n",TX0N);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  double RENORM=TX0;
  for(int I=1; I <= NCHAN; I++)
  {
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC*RENORM > 1.0E-35 ? PDFMC*RENORM : 1.0E-35),
     (ERR*RENORM > 1.0E-35 ? ERR*RENORM : 1.0E-35));
  }
  fclose(IOUT);
  return;
}

//  *********************************************************************
//                       SUBROUTINE PANDX
//  *********************************************************************
double PANDX(double CHI, void* arg)
{
//
//  Differential cross section (per electron and in cm**2) for annihila-
//  tion of positrons with kinetic energy E. CHI=E_-/(E+TREV)
//
//  using namespace CPAN01;

  CPAN01* pan01 = (CPAN01*)arg; 
  
  const double PIELR2=constants::PI*constants::ELRAD*constants::ELRAD;
  double PANDX_RETURN;
//
  double C1=-pow(pan01->GAM+1.0E0,2);
  double C2=pan01->GAM*pan01->GAM+4.0E0*pan01->GAM+1.0E0;
  double S1=C1+C2/CHI-1.0E0/pow(CHI,2);
  PANDX_RETURN=(PIELR2*S1/((pan01->GAM+1.0E0)*(pow(pan01->GAM,2)-1.0E0)))*pow(CHI,pan01->MOM);
  return PANDX_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE GRAaV
//  *********************************************************************
void GRAaV(const pen_context& context, pen_material& mat, const double E, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm for Rayleigh scatter-
//  ing of photons.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
//  ****  Composition data.
//  using namespace COMPOS;
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Photon simulation tables.
//  using namespace CGIMFP;
//  ****  Rayleigh scattering.
//  using namespace CGRA00;
//  using namespace CGRA02;
//  ****  Auxiliary arrays for random sampling verification.
  double XC1[1000];
  int KC1[1000];

  CGRA00 gra00;
  
  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);
//
//  ****  Moments of the 'exact' distribution.
//
  gra00.pmat = &mat;
  double SAVE=mat.SGRA[KE];
  mat.SGRA[KE]=0.0E0;
  gra00.FACTE=2.0E0*pow(E/constants::REV,2);
  double CDTMIN=-1.0E0;
  gra00.MOM=0;
  double TX0=SUMGA(GRAaD,&gra00,CDTMIN,0.90E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.90E0,0.99E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.99E0,0.995E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.995E0,1.00E0,1.0E-5);
  gra00.MOM=1;
  double TX1=SUMGA(GRAaD,&gra00,CDTMIN,0.90E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.90E0,0.99E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.99E0,0.995E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.995E0,1.00E0,1.0E-5);
  gra00.MOM=2;
  double TX2=SUMGA(GRAaD,&gra00,CDTMIN,0.90E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.90E0,0.99E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.99E0,0.995E0,1.0E-5)
         +SUMGA(GRAaD,&gra00,0.995E0,1.00E0,1.0E-5);
  double X1AVE=TX1/TX0;
  double X2AVE=TX2/TX0;
  printf(" # Photon energy = %13.5E eV\n",E);
  printf(" # Normalisation integral = %14.7E\n",TX0);
  printf(" # Exp(mu)    = %14.7E\n",X1AVE);
  printf(" # Exp(mu**2) = %14.7E\n",X2AVE);
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Rayleigh scattering\n");
  fprintf(IOUT," # Photon energy = %13.5E eV\n",E);
  fprintf(IOUT," # Normalisation integral = %14.7E\n",TX0);
  fprintf(IOUT," # Exp(mu)    = %14.7E\n",X1AVE);
  fprintf(IOUT," # Exp(mu**2) = %14.7E\n",X2AVE);
  fprintf(IOUT," # cdtmin     = %14.7E\n",CDTMIN);
  double X1MAX=1.0E0;
  double DX=X1MAX/800.0E0;
  gra00.MOM=0;
  for(int K=1; K<=801; K++)
  {
    double RMU=(K-1)*DX;
    double CDT=1.0E0-2.0E0*RMU;
    double PDF=GRAaD(CDT,&gra00)*2.0E0;
    fprintf(IOUT,"%14.6E%14.6E\n",RMU,PDF);
  }
  fclose(IOUT);
//
//  ************  Random sampling test.
//
  X1MAX=1.0E0;
  int NCHAN;
  double X1MIN=0.0E0;
  if(E<5.0E4)
  {
    NCHAN=100;
  }
  else
  {
    NCHAN=400;
  }
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  int LOOP=1000000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      double CDT;
      int IEFF;
      GRAa(mat,E,KE,XEL,CDT,IEFF,randoms);
      double RMU=(1.0E0-CDT)*0.5E0;
      X1AV=X1AV+RMU;
      double XX2=RMU*RMU;
      X2AV=X2AV+XX2;
      X4AV=X4AV+XX2*XX2;
      int INDEX=(int)((RMU-X1MIN)/DX1+1.0E0);
      if(INDEX>0 && INDEX<=NCHAN)
      {
        KC1[INDEX-1]=KC1[INDEX-1]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    double X1AVMC=X1AV/FNT;
    double X2AVMC=X2AV/FNT;
    printf(" N,<MU>,<MU**2> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  double X1AVMC=X1AV*FNT;
  double X2AVMC=X2AV*FNT;
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Photon energy = %13.5E eV\n",E);
  printf(" # Rayleigh x-section (cm**2) = %13.5E\n",TX0);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Rayleigh scattering\n");
  fprintf(IOUT," # Photon energy (eV)  =%13.5E\n",E);
  fprintf(IOUT," # Rayleigh x-section (cm**2) = %13.5E\n",TX0);
  fprintf(IOUT," # Normalisation integral = %13.5E\n",TX0);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  double RENORM=TX0;
  for(int I=1; I <= NCHAN; I++)
  {
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC*RENORM > 1.0E-35 ? PDFMC*RENORM : 1.0E-35),
     (ERR*RENORM > 1.0E-35 ? ERR*RENORM : 1.0E-35));
  }
  fclose(IOUT);
//
  mat.SGRA[KE]=SAVE; 
  return;
}

//  *********************************************************************
//                       FUNCTION GRAaD
//  *********************************************************************
double GRAaD(double CDT, void* arg)
{
  //  Differential x-section for Rayleigh scattering. Form-factor approx.

  //using namespace PENELOPE_mod;

  //using namespace CGRA00;
  //using namespace CGRA02;

  CGRA00& gra00 = *((CGRA00*)arg); 
  
  double GRAaD_RETURN;
  double Q2 = gra00.FACTE*(1.0-CDT);
  if(Q2 > gra00.pmat->QQM)
  {
    GRAaD_RETURN = 0.0;
    return GRAaD_RETURN;
  }
  GRAaD_RETURN = (1.0+CDT*CDT)*GRAaF2(Q2,gra00);
  if(gra00.MOM > 0){ GRAaD_RETURN = GRAaD_RETURN*pow((1.0-CDT)*0.5,gra00.MOM);}
  return GRAaD_RETURN;
}

//  *********************************************************************
//                       FUNCTION GRAaF2
//  *********************************************************************
double GRAaF2( double Q2, CGRA00& gra00)
{
  //  Squared molecular form factor, as a function of (Q*SL/REV)**2.

  //using namespace PENELOPE_mod;

  //using namespace CGRA00;
  //using namespace CGRA02;

  double GRAaF2_RETURN; //Substitueix a la variable GRAaF2 en fortran. Es el valor que torna la funcio
  if(Q2 < 1.0E-9)
  {
    GRAaF2_RETURN = gra00.pmat->FF0;
  }
  else if(Q2 > gra00.pmat->QQM)
  {
    GRAaF2_RETURN = 0.0;
  }
  else
  {
    double QL = log(Q2);
    unsigned I;
    FINDI(gra00.pmat->QQ,QL,constants::NQ,I);
    double F2 = gra00.pmat->AR[I]+QL*(gra00.pmat->BR[I]+QL*(gra00.pmat->CR[I]+QL*gra00.pmat->DR[I]));
    GRAaF2_RETURN = exp(F2);
  }
  return GRAaF2_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE GCOaV
//  *********************************************************************
void GCOaV(const pen_material& mat, const double E, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm of photon Compton
//  interactions.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//  ****  Compton scattering.
//  using namespace CGCO;
//  ****  Auxiliary arrays for random sampling verification.
//  using namespace CVERIF;
//

  CVERIF verif;
  double CS,CSIN0,CSIN1,CSIN2;
  GCOaT(mat,E,CS);
  GCOTS3(mat,verif,E,CSIN0,CSIN1,CSIN2);
  printf("    \n");
  printf(" # Photon energy = %13.5E eV\n",E);
  printf(" #   Total x-section (approx.) =%13.5E\n",CS);
  printf(" #      Total x-section (num.) =%13.5E\n",CSIN0);
  printf(" # 1st moment x-section (num.) =%13.5E\n",CSIN1);
  printf(" # 2nd moment x-section (num.) =%13.5E\n",CSIN2);
  if(CSIN0<1.0E-35)
  {
    printf(" Sorry, the total cross section is too small.\n");
    return;
  }
//
//  ****  Sort DCS table in increasing energies.
//
  for(int I=1; I<=verif.NX1-1; I++)
  {
    for(int J=I+1; J<=verif.NX1; J++)
    {
      if(verif.XC1[I-1] >= verif.XC1[J-1])
      {
        double SAVE=verif.XC1[I-1];
        verif.XC1[I-1]=verif.XC1[J-1];
        verif.XC1[J-1]=SAVE;
        SAVE=verif.XC2[I-1];
        verif.XC2[I-1]=verif.XC2[J-1];
        verif.XC2[J-1]=SAVE;
      }
    }
  }
  double X1AVE=CSIN1/CSIN0;
  double X2AVE=CSIN2/CSIN0;
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Compton scattering\n");
  fprintf(IOUT," # Photon energy = %13.5E eV\n",E);
  fprintf(IOUT," #   Total x-section (approx.) =%13.5E\n",CS);
  fprintf(IOUT," #      Total x-section (num.) =%13.5E\n",CSIN0);
  fprintf(IOUT," # 1st moment x-section (num.) =%13.5E\n",CSIN1);
  fprintf(IOUT," # 2nd moment x-section (num.) =%13.5E\n",CSIN2);
  for(int I=1; I<=verif.NX1-1; I++)
  {
    fprintf(IOUT," %14.6E%14.6E\n",verif.XC1[I-1],(verif.XC2[I-1]/CSIN0>1.0E-35 ? verif.XC2[I-1]/CSIN0 : 1.0E-35));
  }
  fclose(IOUT);
  if(NSAMP<1){ return;}
//
//  ************  Random sampling test.
//
  double X1MAX=E;
  double X1MIN=0.0E0;
  verif.NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(verif.NCHAN);
  for(int I=1; I<=verif.NCHAN; I++)
  {
    verif.XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    verif.KC1[I-1]=0;
  }
//  ****  Sampling.
  int LOOP=100000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double TIMEI=CPUtime();
      double DE,EP,CDT,ES,CDTS;
      int IZZ,ISH;
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      GCOa(mat,E,DE,EP,CDT,ES,CDTS,IZZ,ISH,randoms);
      double XX=DE;
      X1AV=X1AV+XX;
      double XX2=XX*XX;
      X2AV=X2AV+XX2;
      X4AV=X4AV+XX2*XX2;
      int INDEX=(EP-X1MIN)/DX1+1.0E0;
      if(INDEX>0 && INDEX<=verif.NCHAN)
      {
        verif.KC1[INDEX-1]=verif.KC1[INDEX-1]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    double X1AVMC=X1AV/FNT;
    double X2AVMC=X2AV/FNT;
    printf(" N,<W>,<W**2> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  double X1AVMC=X1AV*FNT;
  double X2AVMC=X2AV*FNT;
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Photon energy = %13.5E eV\n",E);
  printf(" #   Total x-section (approx.) =%13.5E\n",CS);
  printf(" #      Total x-section (num.) =%13.5E\n",CSIN0);
  printf(" # 1st moment x-section (num.) =%13.5E\n",CSIN1);
  printf(" # 2nd moment x-section (num.) =%13.5E\n",CSIN2);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Compton scattering\n");
  fprintf(IOUT," # Photon energy = %13.5E eV\n",E);
  fprintf(IOUT," #   Total x-section (approx.) =%13.5E\n",CS);
  fprintf(IOUT," #      Total x-section (num.) =%13.5E\n",CSIN0);
  fprintf(IOUT," # 1st moment x-section (num.) =%13.5E\n",CSIN1);
  fprintf(IOUT," # 2nd moment x-section (num.) =%13.5E\n",CSIN2);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  for(int I=1; I <= verif.NCHAN; I++)
  {
    double PDFMC=verif.KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",verif.XC1[I-1],
     (PDFMC > 1.0E-35 ? PDFMC : 1.0E-35),
     (ERR > 1.0E-35 ? ERR : 1.0E-35));
  }
  fclose(IOUT);
  return;
}

//  *********************************************************************
//                        SUBROUTINE GCOTS3
//  *********************************************************************
void GCOTS3(const pen_material& mat, CVERIF& verif, const double E,double& CSIN0,double& CSIN1,double& CSIN2)
{
//
//  Total cross section for Compton scattering of photons with energy E
//  in material M.
//
//  using namespace PENELOPE_mod;
//  ****  Compton scattering.
//  using namespace CGCO;
//  ****  Auxiliary arrays for random sampling verification.
//  using namespace CVERIF;
//
  double X[800],Y[800],SUM[800];
//
  double DEM=E/400.0E0;
  CSIN0=0.0E0;
  CSIN1=0.0E0;
  CSIN2=0.0E0;
  int IEND=0;
  verif.NX1=0;
  for(int I=1; I<=mat.NOSCCO+1; I++)
  {	
  	double EPU,EPL;
    if(I==1)
    {
      EPU=E-1.0E-9;
      EPL=E-mat.UICO[I-1]+1.0E-9;
    }
    else if(I==mat.NOSCCO+1)
    {
      EPU=E-mat.UICO[I-1-1]-1.0E-9;
      EPL=1.0E-9;
    }
    else
    {
      EPU=E-mat.UICO[I-1-1]-1.0E-9;
      EPL=E-mat.UICO[I-1]+1.0E-9;
    }
    if(EPL<1.0E-9)
    {
      EPL=1.0E-9;
      IEND=1;
    }
    if(EPL<EPU)
    {
      int NP_P=10+(EPU-EPL)/DEM;
      double DEP=(EPU-EPL)/double(NP_P-1);
      for(int K=1; K<=NP_P; K++)
      {
        double EPA=EPL+double(K-1)*DEP;
        X[K-1]=EPA;
        Y[K-1]=GCOSDT(mat,E,EPA);
        verif.NX1=verif.NX1+1;
        if(verif.NX1>verif.NVF)
        {
          printf("NX1 =%8d\n",verif.NX1);
	        int NERROR = 500; exit(NERROR);
        }
        verif.XC1[verif.NX1-1]=X[K-1];
        verif.XC2[verif.NX1-1]=Y[K-1];
      }
      SLAG6(DEP,Y,SUM,NP_P);
      CSIN0=CSIN0+SUM[NP_P-1];
      for(int K=1; K<=NP_P; K++)
      {
        Y[K-1]=(E-X[K-1])*Y[K-1];
      }
      SLAG6(DEP,Y,SUM,NP_P);
      CSIN1=CSIN1+SUM[NP_P-1];
      for(int K=1; K<=NP_P; K++)
      {
        Y[K-1]=(E-X[K-1])*Y[K-1];
      }
      SLAG6(DEP,Y,SUM,NP_P);
      CSIN2=CSIN2+SUM[NP_P-1];
    }
    if(IEND==1){ break;}
  }
}

//  *********************************************************************
//                        FUNCTION GCOSDT
//  *********************************************************************
double GCOSDT(const pen_material& mat, const double E,double EP)
{
//
//  Single differential cross section for photon Compton scattering, dif-
//  ferential in the energy of the scattered photon. Calculated by num-
//  erical integration of the double differential cross section.
//
//  using namespace PENELOPE_mod;
//  ****  Compton scattering.
//  using namespace CGCO;
//  ****  Auxiliary arrays for random sampling verification.
//  using namespace CVERIF;
//
//  using namespace CGCO01;
//
  double GCOSDT_RETURN;

  CGCO01 gco01;
  
  if(EP<0.0E0 || EP>E)
  {
    GCOSDT_RETURN=0.0E0;
    return GCOSDT_RETURN;
  }
//
  gco01.EE=E;
  gco01.EPP=EP;
  gco01.pmat=&mat;
  double CDTC=1.0E0+(constants::REV/E)-(constants::REV/EP);   //Direction of the Compton line.
  if(CDTC<-1.0E0){ CDTC=-1.0E0+1.0E-10;}
//
  double CS0=0.0E0;
  double DX=1.0E-8;
  double XL=CDTC;
  bool Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    double XU=XL;
    XL=XU-DX;
    if(XL<-1.0E0)
    {
      XL=-1.0E0;
      double CS0P=SUMGA(GCODDS,&gco01,XL,XU,1.0E-7);
      CS0=CS0+CS0P;
    }
    else
    {
      double CS0P=SUMGA(GCODDS,&gco01,XL,XU,1.0E-7);
      CS0=CS0+CS0P;
      DX=DX+DX;
      Eixir=false;
    }
  }
//
  DX=1.0E-8;
  double XU=CDTC;
  Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    XL=XU;
    XU=XL+DX;
    if(XU>1.0E0)
    {
      XU=1.0E0;
      double CS0P=SUMGA(GCODDS,&gco01,XL,XU,1.0E-7);
      CS0=CS0+CS0P;
    }
    else
    {
      double CS0P=SUMGA(GCODDS,&gco01,XL,XU,1.0E-7);
      CS0=CS0+CS0P;
      DX=DX+DX;
      Eixir=false;
    }
  }
  GCOSDT_RETURN=CS0;
  return GCOSDT_RETURN;
}

//  *********************************************************************
//                        FUNCTION GCODDS
//  *********************************************************************
double GCODDS(double CDT, void* arg)
{
//
//  Double differential cross section for photon Compton scattering, dif-
//  differential in the direction and energy of the scattered photon.
//
//  The energies of the primary, E, and the secondary, EP, are entered
//  through common CGCO01.
//
//  using namespace PENELOPE_mod;
//
  const double PIELR2=constants::PI*constants::ELRAD*constants::ELRAD;
  const double D2=1.4142135623731E0;
  const double D1=1.0E0/D2;
  const double D12=0.5E0;
//  ****  Compton scattering.
//  using namespace CGCO;
//
//  using namespace CGCO01;
//
  double GCODDS_RETURN;

  CGCO01& gco01 = *((CGCO01*)arg);
  
  double& E=gco01.EE;
  double &EP=gco01.EPP;
  const pen_material& mat=*(gco01.pmat);

  double DE=E-EP;
//  ****  Momentum transfer.
  double QC2=E*E+EP*EP-2.0E0*E*EP*CDT;
  if(QC2<1.0E-24)
  {
    GCODDS_RETURN=0.0E0;
    return GCODDS_RETURN;
  }
  double QC=sqrt(QC2);
  double CDT1=1.0E0-CDT;
//  ****  Energy of the Compton line.
  double EOEC=1.0E0+(E/constants::REV)*CDT1;
  double ECOE=1.0E0/EOEC;
  double EC=E*ECOE;
//  ****  Electron momentum.
  double PZOMC=E*(EP-EC)/(EC*QC);
  if(fabs(PZOMC)>1.0E0)
  {
    GCODDS_RETURN=0.0E0;
    return GCODDS_RETURN;
  }
//  ****  Klein-Nishina X-factor.
  double XKN=EOEC+ECOE-1.0E0+CDT*CDT;
//  ****  F(Pz)-function.
  double XQC=E*E+EC*EC-2.0E0*E*EC*CDT;
  double aux=(PZOMC<0.2E0 ? PZOMC : 0.2E0);
  double XPZ=(aux > -0.2E0 ? aux : -0.2E0);
  double FPZ=1.0E0+(sqrt(XQC)/E)*(1.0E0+EC*(EC-E*CDT)/XQC)*XPZ;
  if(FPZ<0.0E0)
  {
    GCODDS_RETURN=0.0E0;
    return GCODDS_RETURN;
  }
//  ****  Analytical Compton profile.
  double FJA=0.0E0;
  DE=E-EP;
  for(int I=1; I<=mat.NOSCCO; I++)
  {
    if(DE>mat.UICO[I-1])
    {
      double H=D1+D2*fabs(mat.FJ0[I-1]*PZOMC);
      double FJAP=D2*H*exp(D12-pow(H,2))*mat.FJ0[I-1];
      FJA=FJA+mat.FCO[I-1]*FJAP;
    }
  }
//
  double DPZDEP=(EOEC+PZOMC*(E*CDT-EP)/QC)/QC;
  GCODDS_RETURN=PIELR2*pow(ECOE,2)*XKN*FPZ*FJA*DPZDEP;
  return GCODDS_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE GPHaAV
//  *********************************************************************
void GPHaAV(const double E, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm for the initial dir-
//  ection of photoelectrons. Sauter differential x-section.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
//  ****  Composition data.
//  using namespace COMPOS;
//
//  using namespace CGPHVS;

  CGPHVS gphvs;
  
//  ****  Auxiliary arrays for random sampling verification.
  double XC1[1000];
  int KC1[1000];
//
//  ****  Moments of the 'exact' distribution.
//
  gphvs.GAM=1.0E0+E/constants::REV;
  double GAM2=gphvs.GAM*gphvs.GAM;
  gphvs.BETA=sqrt((GAM2-1.0E0)/GAM2);
  gphvs.MOM=0;
  double TX0=SUMGA(SAUTDX,&gphvs,0.0E0,1.0E0,1.0E-5);
  gphvs.MOM=1;
  double TX1=SUMGA(SAUTDX,&gphvs,0.0E0,1.0E0,1.0E-5);
  gphvs.MOM=2;
  double TX2=SUMGA(SAUTDX,&gphvs,0.0E0,1.0E0,1.0E-5);
  double X1AVE=TX1/TX0;
  double X2AVE=TX2/TX0;
  printf(" # Photon energy = %13.5E eV\n",E);
  printf(" # Sauter x-section = %14.7E\n",TX0);
  printf(" # Exp(mu)    = %14.7E\n",X1AVE);
  printf(" # Exp(mu**2) = %14.7E\n",X2AVE);
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Sauter photoelectron angular distribution\n");
  fprintf(IOUT," # Photon energy = %13.5E eV\n",E);
  fprintf(IOUT," # Sauter x-section = %14.7E\n",TX0);
  fprintf(IOUT," # Exp(mu)    = %14.7E\n",X1AVE);
  fprintf(IOUT," # Exp(mu**2) = %14.7E\n",X2AVE);
  double DX=1.0E0/200.0E0;
  gphvs.MOM=0;
  for(int K=1; K<=201; K++)
  {
    double RMU=(K-1)*DX;
    double PDF=SAUTDX(RMU,&gphvs);
    fprintf(IOUT,"%14.6E%14.6E\n",RMU,PDF);
  }
  fclose(IOUT);
//
//  ************  Random sampling test.
//
  double X1MAX=1.0E0;
  double X1MIN=0.0E0;
  int NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(NCHAN);
  for(int I=1; I<=NCHAN; I++)
  {
    XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    KC1[I-1]=0;
  }
//
  pen_state_gPol gPolstate;

  int LOOP=1000000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      double CDT,DFS;
      SAUTER(E,CDT,DFS,randoms,gPolstate);
      double RMU=(1.0E0-CDT)*0.5E0;
      X1AV=X1AV+RMU;
      double XX2=RMU*RMU;
      X2AV=X2AV+XX2;
      X4AV=X4AV+XX2*XX2;
      int INDEX=(int)((RMU-X1MIN)/DX1+1.0E0);
      if(INDEX>0 && INDEX<=NCHAN)
      {
        KC1[INDEX-1]=KC1[INDEX-1]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    double X1AVMC=X1AV/FNT;
    double X2AVMC=X2AV/FNT;
    printf(" N,<MU>,<MU**2> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  double X1AVMC=X1AV*FNT;
  double X2AVMC=X2AV*FNT;
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Photon energy = %13.5E eV\n",E);
  printf(" # Sauter x-section = %13.5E\n",TX0);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Sauter photoelectron angular distribution\n");
  fprintf(IOUT," # Photon energy (eV)  =%13.5E\n",E);
  fprintf(IOUT," # Sauter x-section = %13.5E\n",TX0);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  double RENORM=TX0;
  for(int I=1; I <= NCHAN; I++)
  {
    double PDFMC=KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",XC1[I-1],
     (PDFMC*RENORM > 1.0E-35 ? PDFMC*RENORM : 1.0E-35),
     (ERR*RENORM > 1.0E-35 ? ERR*RENORM : 1.0E-35));
  }
  fclose(IOUT);
  return;
}

//  *********************************************************************
//                       FUNCTION SAUTDX
//  *********************************************************************
double SAUTDX(double T, void* arg)
{
//
//  Sauter differential x-section.
//
//  using namespace CGPHVS;

  CGPHVS& gphvs = *((CGPHVS*)arg);

  double CTH=1.0E0-2.0E0*T;
  double PROB=((1.0E0-CTH*CTH)/pow(1.0E0-gphvs.BETA*CTH,3))*(0.5E0*gphvs.GAM*(gphvs.GAM-1.0E0)*(gphvs.GAM-2.0E0)+1.0E0/(1.0E0-gphvs.BETA*CTH));
  double SAUTDX_RETURN=PROB*pow(T,gphvs.MOM);
  return SAUTDX_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE GPHaSV
//  *********************************************************************
void GPHaSV(const pen_context& context, const pen_material& mat, const double E, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling of the active atomic electron
//  shell in photoelectric absorption.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//
//  using namespace CECUTR;
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Element data.
//  using namespace CADATA;
//  ****  Composition data.
//  using namespace COMPOS;
//  ****  Photon simulation tables.
//  using namespace CGIMFP;
//  ****  Photoelectric cross sections.
//  using namespace CGPH00;
//
  double PROB[30][17];
  int NCOUNT[30][17],ISTORE[99];
//
//  ****  Initialisation.
//
  for(int IEL=1; IEL<=mat.NELEM; IEL++)
  {
    for(int ISH=1; ISH<=17; ISH++)
    {
      PROB[IEL-1][ISH-1]=0.0E0;
    }
  }

  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  context.grid.getInterval(E,KE,XEL,XE,XEK);
  
  double PTOTA=mat.SGPH[KE]/mat.VMOL;
  printf("PTOTA=%14.7E\n",PTOTA);
  double PTOT=0.0E0;
  for(int IEL=1; IEL<=mat.NELEM; IEL++)
  {
    int IZZ=mat.IZ[IEL-1];
    ISTORE[IZZ-1]=IEL;
//  ****  Binary search.
    int I=context.elements.IPHF[IZZ-1];
    int IU=context.elements.IPHL[IZZ-1];
    bool Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      int IT = (I+IU)/2;
      if(XEL > context.elements.EPH[IT]){
        I = IT;
      }
      else
      {
        IU = IT;
      }
      if(IU-I > 1){ Eixir = false; continue;}
    }
//
    double DEE=context.elements.EPH[I+1]-context.elements.EPH[I];
    double PSHELL=0.0E0;
    for(int ISH=1; ISH<=context.elements.NPHS[IZZ-1]+1; ISH++)
    {
      double PCSL;
      if(DEE>1.0E-15)
      {
        PCSL=context.elements.XPH[I][ISH-1]+(context.elements.XPH[I+1][ISH-1]-context.elements.XPH[I][ISH-1])*(XEL-context.elements.EPH[I])/DEE;
      }
      else
      {
        PCSL=context.elements.XPH[I][ISH-1];
      }
      if(ISH==1)
      {
        PTOT=PTOT+mat.STF[IEL-1]*exp(PCSL);
        PROB[IEL-1][17-1]=mat.STF[IEL-1]*exp(PCSL);    	
      }
      else
      {
        if(context.elements.EB[IZZ-1][ISH-1-1]<E && context.elements.EB[IZZ-1][ISH-1-1]>mat.ECUTR)
        {
          PROB[IEL-1][ISH-1-1]=mat.STF[IEL-1]*exp(PCSL);
          PSHELL=PSHELL+PROB[IEL-1][ISH-1-1];
        }
      }
    }
    PROB[IEL-1][17-1]=PROB[IEL-1][17-1]-PSHELL;
  }
  printf("PTOT=%14.7E\n",PTOT);
//
  for(int IEL=1; IEL<=mat.NELEM; IEL++)
  {
    for(int ISH=1; ISH<=17; ISH++)
    {
      NCOUNT[IEL-1][ISH-1]=0;
    }
  }
//
//  ****  Simulation starts here.
//
  int LOOP=1000000;
  int NLOOP=NSAMP/LOOP;
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      int  IZZ,ISH;
      bool Eixir = false;
      while(!Eixir)
      {
        Eixir = true;
        double ES;
        GPHa(mat,context.elements,E,KE,XEL,ES,IZZ,ISH,randoms);
        if(ISH==0){ Eixir=false; continue;}
      }
      int IEL=ISTORE[IZZ-1];
      NCOUNT[IEL-1][ISH-1]=NCOUNT[IEL-1][ISH-1]+1;
    }
    int NSAMPLE=I*LOOP;
    double RN=double(NSAMPLE);
    printf(" N =%9d\n",I*LOOP);
//  ****  The shell ionization probabilities obtained from the simulation
//        are compared with the probabilities obtained from the photo-
//        electric  x-sections.
    printf("  The total probability may be <1, owing to delta interactions.\n");
    for(int IEL=1; IEL<=mat.NELEM; IEL++)
    {
      for(int ISH=1; ISH<=17; ISH++)
      {
        double PMC=NCOUNT[IEL-1][ISH-1]/RN;
        double ERR=3.0E0*sqrt(PMC*(1.0E0-PMC)/RN);
        double PINP=PROB[IEL-1][ISH-1]/PTOT;
        if(PINP>1.0E-25)
        {
          double TESTE=(PMC-PINP)/(ERR>1.0E-16 ? ERR : 1.0E-16);
          printf(" %4d%4d%13.5E%13.5E%9.1E%9.1E\n",mat.IZ[IEL-1],ISH,PINP,PMC,ERR,TESTE);
        }
      }
    }
  }
  int NSAMPLE=NLOOP*LOOP;
  double RN=double(NSAMPLE);
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Shell photoabsorption probabilities\n");
  fprintf(IOUT," # Photon energy (eV) = %14.7E\n",E);
  int ITR=0;
  for(int IEL=1; IEL<=mat.NELEM; IEL++)
  {
    for(int ISH=1; ISH<=17; ISH++)
    {
      ITR=ITR+1;
      double PINP=PROB[IEL-1][ISH-1]/PTOT;
      if(PINP > 1.0E-15){ fprintf(IOUT,"%4d%13.5E  %4d%4d\n",ITR,PINP,mat.IZ[IEL-1],ISH);}
    }
  }
  fclose(IOUT);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Shell photoabsorption probabilities\n");
  fprintf(IOUT," # Photon energy (eV) = %14.7E\n",E);
  ITR=0;
  for(int IEL=1; IEL<=mat.NELEM; IEL++)
  {
    for(int ISH=1; ISH<=17; ISH++)
    {
      ITR=ITR+1;
      double PMC=NCOUNT[IEL-1][ISH-1]/RN;
      double ERR=3.0E0*sqrt(PMC*(1.0E0-PMC)/RN);
      if(PMC > 0.0E0){ fprintf(IOUT,"%4d%13.5E%13.5E\n",ITR,PMC,ERR);}
    }
  }
  fclose(IOUT);
}

//  *********************************************************************
//                       SUBROUTINE GPPaV
//  *********************************************************************
void GPPaV(const pen_context& context, const pen_material& mat, const double E, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies the sampling algorithm for electron-positron
//  pair production.
//
//  using namespace TRACK_mod;
//  using namespace PENELOPE_mod;
//  ****  Composition data.
//  using namespace COMPOS;
//  ****  Energy grid and interpolation constants for the current energy.
//  using namespace CEGRID;
//  ****  Pair-production cross section parameters.
//  using namespace CGPP00;
//  using namespace CGPP02;
//  ****  Auxiliary arrays for random sampling verification.
//  using namespace CVERIF;
//
//  ****  Moments of the 'exact' distribution.
//

  CGPP02 gpp02;
  CVERIF verif;
  const double REV = constants::REV;
  
  double TX0,X1AVE,X2AVE;
  GPPTX(context.elements,gpp02,mat.ZEQPP,E,TX0,X1AVE,X2AVE);

  printf(" # Electron-positron pair production\n");
  printf(" # Photon energy (eV) = %13.5E\n",E);
  printf(" # Pair production x-section (cm**2) = %14.7E\n",TX0);
  printf(" # Exp(E_electron) (eV)              = %14.7E\n",X1AVE);
  printf(" # Exp(E_electron**2) (eV)           = %14.7E\n",X2AVE);
//
  FILE* IOUT = fopen("calc.dat","w");
  fprintf(IOUT," # Electron-positron pair production\n");
  fprintf(IOUT," # Photon energy (eV) = %13.5E\n",E);
  fprintf(IOUT," # Pair production x-section (cm**2) = %14.7E\n",TX0);
  fprintf(IOUT," # Exp(E_electron) (eV)              = %14.7E\n",X1AVE);
  fprintf(IOUT," # Exp(E_electron**2) (eV)           = %14.7E\n",X2AVE);
  double DX=1.0E0/500.0E0;
  gpp02.MOM=0;
  for(int K=1; K<=500; K++)
  {
    double EPS=(K-0.5E0)*DX;
    double PDF=GPPaD(EPS,&gpp02);
    fprintf(IOUT,"%14.6E%14.6E\n",EPS,PDF);
  }
  fclose(IOUT);
//
//  ************  Random sampling test.
//
  double X1MAX=1.0E0;
  double X1MIN=0.0E0;
  verif.NCHAN=100;
//  ****  Counters initialisation.
  double DX1=(X1MAX-X1MIN)/double(verif.NCHAN);
  for(int I=1; I<=verif.NCHAN; I++)
  {
    verif.XC1[I-1]=X1MIN+(I-0.5E0)*DX1;
    verif.KC1[I-1]=0;
  }
//
  int LOOP=250000;
  int NLOOP=1+(NSAMP/LOOP);
  printf("NLOOP =%d\n",NLOOP);
  double X1AV=0.0E0;
  double X2AV=0.0E0;
  double X4AV=0.0E0;
  double TIMEI=CPUtime();
  double X1AVMC,X2AVMC;

  //double TimeFunc = 0.0;
  
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=0; J<LOOP; J++)
    {
      double EE,CDTE,EP,CDTP;
      int IAZ,ISA;

      //double TimeFunc2=CPUtime();
      GPPa(mat,E,0,0.0,EE,CDTE,EP,CDTP,IAZ,ISA,randoms);

      //TimeFunc += CPUtime()-TimeFunc2;
      
      double EPS=(EE+REV)/E;
      X1AV=X1AV+EPS;
      double XX2=EPS*EPS;
      X2AV=X2AV+XX2;
      X4AV=X4AV+XX2*XX2;
      int INDEX=(int)((EPS-X1MIN)/DX1);
      if(INDEX >= 0 && INDEX < verif.NCHAN)
      {
        verif.KC1[INDEX]=verif.KC1[INDEX]+1;
      }
      else
      {
        printf(" Index out of range.\n");
      }
    }
    double FNT=double(I*LOOP);
    X1AVMC=X1AV/FNT;
    X2AVMC=X2AV/FNT;
    printf(" N,<EPS>,<EPS**2> =%9d%15.7E%15.7E\n",I*LOOP,X1AVMC/X1AVE,X2AVMC/X2AVE);
  }
  double FNT=1.0E0/(double(NLOOP)*LOOP);
  //printf("Function time: %14.5E\n",TimeFunc*FNT);
  X1AVMC=X1AV*FNT;
  X2AVMC=X2AV*FNT;
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  printf("    \n");
  printf(" # Photon energy = %13.5E eV\n",E);
  printf(" # PP x-section  = %13.5E cm**2\n",TX0);
  printf(" # Monte Carlo generated moments.\n");
  printf(" #    First moment (numerical) =%14.7E\n",X1AVE);
  double ERR1MC=3.0E0*sqrt((X2AVMC-pow(X1AVMC,2))*FNT > 1.0E-35 ? (X2AVMC-pow(X1AVMC,2))*FNT : 1.0E-35);
  printf(" #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  printf(" #   Second moment (numerical) =%14.7E\n",X2AVE);
  double X4AVMC=X4AV*FNT;
  double ERR2MC=3.0E0*sqrt((X4AVMC-pow(X2AVMC,2))*FNT > 1.0E-35 ? (X4AVMC-pow(X2AVMC,2))*FNT : 1.0E-35);
  printf(" # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  IOUT = fopen("simul.dat","w");
  fprintf(IOUT," # Electron-positron pair production\n");
  fprintf(IOUT," # Photon energy (eV) = %13.5E\n",E);
  fprintf(IOUT," # Pair production x-section (cm**2) = %13.5E \n",TX0);
  fprintf(IOUT," # Monte Carlo generated moments.\n");
  fprintf(IOUT," #    First moment (numerical) =%14.7E\n",X1AVE);
  fprintf(IOUT," #  First moment (Monte Carlo) =%14.7E +-%9.2E\n", X1AVMC,ERR1MC);
  fprintf(IOUT," #   Second moment (numerical) =%14.7E\n",X2AVE);
  fprintf(IOUT," # Second moment (Monte Carlo) =%14.7E +-%9.2E\n",X2AVMC,ERR2MC);
//
  double RENORM=TX0;
  for(int I=1; I <= verif.NCHAN; I++)
  {
    double PDFMC=verif.KC1[I-1]*FNT;
    double ERR=3.0E0*sqrt(fabs(PDFMC*(1.0E0-PDFMC)*FNT))/DX1;
    PDFMC=PDFMC/DX1;
    fprintf(IOUT," %14.6E%14.6E%14.6E\n",verif.XC1[I-1],
     (PDFMC*RENORM > 1.0E-35 ? PDFMC*RENORM : 1.0E-35),
     (ERR*RENORM > 1.0E-35 ? ERR*RENORM : 1.0E-35));
  }
  fclose(IOUT);
  return;
}

//  *********************************************************************
//                       SUBROUTINE GPPTX
//  *********************************************************************
void GPPTX(const pen_elementDataBase& elements, CGPP02& gpp02, const double Z, const double E, double& TX0, double& EPS1AV, double& EPS2AV)
{
//
//  Total atomic cross section for electron-positron pair production.
//  Bethe-Heitler differential cross section with exponential screening,
//  Coulomb correction and an empirical correction for low energies.
//
//  Input arguments:
//    Z ...... (effective) atomic number,
//    E ...... photon energy (eV).
//
//  Output argument:
//    TX0 .... pair-production cross section (barn/atom).
//    EPS1AV ... expectation value of epsilon.
//    EPS2AV ... expectation value of epsilon**2 squared.
//
  const double ELR2=constants::ELRAD*constants::ELRAD;

//  ****  Element data.
//  using namespace CADATA;
//
//  using namespace CGPP02;

  double &BFACT=gpp02.FACT;
//
  TX0=0.0E0;
  if(E<2.0E0*constants::REV){ return;}
  int IZ=(int)(Z+0.25E0);
  if(IZ <= 0){ IZ=1;}
  if(IZ > 100){ IZ=100;}
  gpp02.GAM=constants::REV/E;
//  ****  DBM COULOMB CORRECTION.
  double ALZ=Z/constants::SL;
  double A=ALZ*ALZ;
  double FC=A*(0.202059E0-A*(0.03693E0-A*(0.00835E0-A*(0.00201E0-A*
         (0.00049E0-A*(0.00012E0-A*0.00003)))))+1.0E0/(A+1.0E0));
//  ****  Low-energy correction.
  double T=sqrt(2.0E0*gpp02.GAM);
  double F00=(-1.774E0-1.210E1*ALZ+1.118E1*A)*T
          +(8.523E0+7.326E1*ALZ-4.441E1*A)*pow(T,2)
          -(1.352E1+1.211E2*ALZ-9.641E1*A)*pow(T,3)
          +(8.946E0+6.205E1*ALZ-6.341E1*A)*pow(T,4);
//  ****  Triplet production.
  double ETAE;
  if(gpp02.GAM<0.25E0)
  {
    T=log(0.25E0/gpp02.GAM);
    double CSF=(2.824E-1-1.909E-1*ALZ)*T
           +(1.095E-1+2.206E-1*ALZ)*pow(T,2)
           -(2.888E-2+4.269E-2*ALZ)*pow(T,3)
           +(2.527E-3+2.623E-3*ALZ)*pow(T,4);
    ETAE=(1.0E0-exp(-CSF))*elements.ETA[IZ-1];
  }
  else
  {
    ETAE=0.0E0;
  }
//
  gpp02.BCBT=2.0E0/elements.RSCR[IZ-1];
  gpp02.G0=4.0E0*log(elements.RSCR[IZ-1])-4.0E0*FC+F00;
  BFACT=1.0093E0*6.66666666666666E-1*ALZ*(Z+ETAE)*ELR2;
//  ****  Total cross section.
  gpp02.MOM=0;
  TX0=SUMGA(GPPaD,&gpp02,gpp02.GAM,1.0E0-gpp02.GAM,1.0E-9);
//  ****  1st and 2nd momenta of the electron kinetic energy
//        distribution.
  gpp02.MOM=1;
  EPS1AV=SUMGA(GPPaD,&gpp02,gpp02.GAM,1.0E0-gpp02.GAM,1.0E-9)/TX0;
  gpp02.MOM=2;
  EPS2AV=SUMGA(GPPaD,&gpp02,gpp02.GAM,1.0E0-gpp02.GAM,1.0E-9)/TX0;
}

//  *********************************************************************
double GPPaD(double EPS, void* arg)
{
  //using namespace CGPP02;

  double GPPaD_RETURN;
  CGPP02& gpp02 = *((CGPP02*)arg);

  double& BFACT=gpp02.FACT;

  if(EPS<gpp02.GAM || EPS > 1.0E0-gpp02.GAM)
  {
    GPPaD_RETURN=0.0E0;
    return GPPaD_RETURN;
  }
  double EPS1=1.0E0-EPS;
  double B=gpp02.GAM/(gpp02.BCBT*EPS*EPS1);
  double G1,G2;
  SCHIFF(B,G1,G2);
  double G10=(G1+gpp02.G0>0.0E0 ? G1+gpp02.G0 : 0.0E0);
  double G20=(G2+gpp02.G0>0.0E0 ? G2+gpp02.G0 : 0.0E0);
  double DCSPP=BFACT*(2.0E0*pow(EPS-0.5E0,2)*G10+G20);
  GPPaD_RETURN=DCSPP*pow(EPS,gpp02.MOM);
  return GPPaD_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE RELAXV
//  *********************************************************************
void RELAXV(const pen_context& context, const pen_material& mat, pen_elementDataBase& elements, const double E, int IZ, int IS, int NSAMP, pen_rand& randoms)
{
//
//  This subroutine verifies that the de-excitation sampling routine does
//  work properly.
//
//  ****  Atomic relaxation data.
  //using namespace CRELAX;

  int NCOUNT[constants::NRX];
//
//  ****  Initialisation.
//
  int KF=elements.IFIRST[IZ-1][IS-1];
  int KL=elements.ILAST[IZ-1][IS-1];
  if(KF==0 || KL==0)
  {          
    printf("psampler:RELAXV:ERROR: KF =%4d KL =%4d\n",KF,KL);
    int NERROR = 501; exit(NERROR);
  }  
  if(KF==KL)
  {
    printf(" Sorry, no transitions for this shell/material.\n");
    return;
  }
//
  double SUM0=0.0E0;
  for(int K=KF; K<=KL; K++)
  {
    SUM0=SUM0+elements.P[K-1];
    printf("%6d%6d%6d%6d%16.8E%16.8E\n",K,IS,elements.IS1[K-1],elements.IS2[K-1],elements.P[K-1],elements.ET[K-1]);
  }
  printf(" SUM0 =%21.17F     \n",SUM0);
//
  for(unsigned I=1; I<=constants::NRX; I++)
  {
    NCOUNT[I-1]=0;
  }
//
//  ****  Here we sample transitions that fill a vacancy in the shell IS.
//
  //Create two dummy particle stacks
  pen_particleStack<pen_particleState> stackE;
  pen_particleStack<pen_state_gPol> stackG;

  //Create auxiliar state
  pen_particleState state;
  state.E = E;
  state.MAT = 1;
  
  double RN=1.0E0;
  int LOOP=2500000;
  int NLOOP=(NSAMP/LOOP>1 ? NSAMP/LOOP : 1);
  double TIMEI=CPUtime();
  for(int I=1; I<=NLOOP; I++)
  {
    for(int J=1; J<=LOOP; J++)
    {
      int KS;
      RELAX(context.elements, mat, state, 1, 0, IZ, IS, KS, stackE, stackG, randoms);
      NCOUNT[KS-1]++;
      //Clear dummy stacks
      stackE.cleans();
      stackG.cleans();
    }
    int NSAMPLE=I*LOOP;
    RN=double(NSAMPLE);
//  ****  The input transition probabilities are compared with the proba-
//        bilities that result from the Monte Carlo simulation. The
//        quantity ERR is the statistical uncertainty (3 sigma).
//        If TESTE.LT.1.0D0 the sampling algorithm works properly.
    for(int J=KF; J<=KL; J++)
    {
      double PMC=double(NCOUNT[J-1])/RN;
      double ERR=3.0E0*sqrt(PMC*(1.0E0-PMC)/RN);
      double PINP=elements.P[J-1]/SUM0;
      double TESTE=(PMC-PINP)/(ERR>1.0E-16 ? ERR : 1.0E-16);
      printf("%4d%4d%4d%13.5E%13.5E%9.1E%9.1E\n",IS,elements.IS1[J-1],elements.IS2[J-1],PINP,PMC,ERR,TESTE);
    }
    printf("     \n");
  }
  double TIMEF=CPUtime();
  double TIMET=TIMEF-TIMEI;
  double SPEED=(double(NLOOP)*LOOP)/TIMET;
  printf(" Simulation speed =%11.3E random_values/sec\n",SPEED);
//
  FILE* IOUT = fopen("calc.dat","w");
  int IX=0;
  for(int J=KF; J<=KL; J++)
  {
    IX=IX+1;
    double PINP=elements.P[J-1]/SUM0;
    fprintf(IOUT,"%4d%14.6E\n",IX,PINP);
  }
  fclose(IOUT);
//
  IOUT = fopen("simul.dat","w");
  IX=0;
  for(int J=KF; J<=KL; J++)
  {
    IX=IX+1;
    double PMC=double(NCOUNT[J-1])/RN;
    double ERR=3.0E0*sqrt(PMC*(1.0E0-PMC)/RN);
    fprintf(IOUT,"%4d%14.6E%14.6E\n",IX,PMC,ERR);
  }
  fclose(IOUT);
}

