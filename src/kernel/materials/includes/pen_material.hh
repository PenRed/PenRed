
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


#ifndef __PEN_MATERIAL_
#define __PEN_MATERIAL_

#include <cmath>
#include <cstdlib>


#include "pen_classes.hh"
#include "pen_grids.hh"
#include "pen_states.hh"


class pen_material;
class pen_elementDataBase;

// Auxiliar initialization structs

struct CEIN00
{
  double SEH0[constants::NO], SEH1[constants::NO], SEH2[constants::NO], SES0[constants::NO], SES1[constants::NO], SES2[constants::NO], SET0[constants::NO], SET1[constants::NO], SET2[constants::NO];  
};

struct CESIN
{
  //  ****  Electron inelastic coll. and inner-shell ionisation tables.
  double XSEIN[constants::NEGP][constants::NO], XSESI[constants::NEGP][constants::NO];
  int ISIE[constants::NO];
};

//  ****  Inner-shell ionisation by electron and positron impact.
struct CESI0
{
  static const unsigned int NRP = 8000; 
  static const unsigned int nElements = 99;
  double XESI[NRP][16];
  int IESIF[nElements], IESIL[nElements], NSESI[nElements], NCURE;

  CESI0()
  {
    //  This subroutine sets all variables in common /CESI0/ to zero.
    //
    //  It has to be invoked before reading the first material definition
    //  file.


    for(unsigned I = 0; I < nElements; I++)
      {
	IESIF[I] = 0;
	IESIL[I] = 0;
	NSESI[I] = 0;
      }

    for(unsigned I = 0; I < NRP; I++)
      {
	for(int J = 0; J < 16; J++)
	  {
	    XESI[I][J] = -80.6;
	  }
      }
    NCURE = 0;
  }  
};

struct CEBR01
{
  double EBT[constants::NBE], XS[constants::NBE][constants::NBW], TXS[constants::NBE], X[constants::NBE], Y[constants::NBE];
};

//  ****  Inner-shell ionisation by electron and positron impact.
struct CPSI0
{
  static const unsigned int NRP = 8000; 
  static const unsigned int nElements = 99;
  double XPSI[NRP][16];
  int IPSIF[nElements], IPSIL[nElements], NSPSI[nElements], NCURP;

  CPSI0()
  {
    //  This subroutine sets all variables in common /CPSI0/ to zero.
    //
    //  It has to be invoked before reading the first material definition
    //  file.


    for(unsigned I = 0; I < nElements; I++)
      {
	IPSIF[I] = 0;
	IPSIL[I] = 0;
	NSPSI[I] = 0;
      }

    for(unsigned I = 0; I < NRP; I++)
      {
	for(int J = 0; J < 16; J++)
	  {
	    XPSI[I][J] = -80.6;
	  }
      }
    NCURP = 0;
  }  
};

//  ****  Elastic scattering of electrons and positrons.
struct CEEL00{
  
  double EJT[constants::NEGP], XE0[constants::NEGP], XE1[constants::NEGP], XE2[constants::NEGP], XP0[constants::NEGP],
    XP1[constants::NEGP], XP2[constants::NEGP], T1E0[constants::NEGP], T2E0[constants::NEGP], T1P0[constants::NEGP],
    T2P0[constants::NEGP], EJTL[constants::NEGP], FJL[constants::NEGP], A[constants::NEGP], B[constants::NEGP], C[constants::NEGP], D[constants::NEGP];  
};

//  ****  E/P inelastic collisions.
struct CEINTF{
  
  double T1EI[constants::NEGP], T2EI[constants::NEGP], T1PI[constants::NEGP], T2PI[constants::NEGP];
};

struct CPIN00
{
  //  ****  Partial cross sections of individual shells/oscillators.
  // Associat a les routines PEMATR i PINaT solament. COMPTE !! Les variables del COMMON tene noms differents en les dues.
  double SPH0[constants::NO], SPH1[constants::NO], SPH2[constants::NO], SPS0[constants::NO], SPS1[constants::NO], SPS2[constants::NO], SPT0[constants::NO], SPT1[constants::NO], SPT2[constants::NO];  
};

struct CPIN01
{
  double EI,CPS,BHA1,BHA2,BHA3,BHA4;
  int MOM;
};

struct CPSIN
{
  //  ****  Positron inelastic coll. and inner-shell ionisation tables.
  double XSPIN[constants::NEGP][constants::NO], XSPSI[constants::NEGP][constants::NO];
  int ISIP[constants::NO];
};

//  ****  Elastic scattering simulation tables. (CDCSEP)
struct CDCSEP{

  const static unsigned int NE = 96;
  const static unsigned int NA = 606;
  
  double ETS[NE], ETL[NE], TH[NA], THR[NA], XMU[NA], XMUL[NA],
    ECS[NE], ETCS1[NE], ETCS2[NE], EDCS[NE][NA], PCS[NE], PTCS1[NE],
    PTCS2[NE], PDCS[NE][NA], DCSI[NA], DCSIL[NA], CSI, TCS1I, TCS2I;
};

struct CGPH01
{
  //  ****  Photon simulation tables.
  static const unsigned NDIM = 12000;
  double ER[NDIM], XSR[NDIM];
  int NPHD;  
};

struct CEIN01
{
  double EI, EE, CPS, AMOL, MOM;
};

//******************************************************************
//Create a struct to encapsulate all auxiliar initalization structs
struct initStructs{
    CEIN00* pcein00;
    CESIN* pcesin;
    CESI0* pcesi0;
    CEBR01* pcebr01;
    CPSI0* pcpsi0;
    CEEL00* pceel00;
    CEINTF* pceintf;
    CPIN00* pcpin00;
    CPSIN* pcpsin;
    CPIN01* pcpin01;
    CDCSEP* pcdcsep;
    CRITA* pcrita;
    CGPH01* pcgph01;
    
  initStructs(){
    allocate();
  }
  
  void clear(){
    free();
    allocate();
  }

  void free(){
    delete pcein00;
    delete pcesin;
    delete pcesi0;
    delete pcebr01;
    delete pcpsi0;
    delete pceel00;
    delete pceintf;
    delete pcpin00;
    delete pcpsin;
    delete pcpin01;
    delete pcdcsep;
    delete pcrita;
    delete pcgph01;            
  }
  //******************************************************************

    ~initStructs(){free();}

private:
  void allocate(){
    pcein00 = new CEIN00;
    pcesin = new CESIN;
    pcesi0 = new CESI0;
    pcebr01 = new CEBR01;
    pcpsi0 = new CPSI0;
    pceel00 = new CEEL00;
    pceintf = new CEINTF;
    pcpin00 = new CPIN00;
    pcpsin = new CPSIN;
    pcpin01 = new CPIN01;
    pcdcsep = new CDCSEP;
    pcrita = new CRITA;
    pcgph01 = new CGPH01;          
  }
};


//-------------------
// Element database
//-------------------

class pen_elementDataBase{

 public:

  //  ****  Element data.

  static const unsigned int nElements = 99;
  
  //  EB  -> Correspond to the Ui value of specified element and shell
  //       (Shell ionization energy).
  //
  //  IKS -> Code number of the shell, oredere as,
  //        "  ","K ","L1","L2","L3","M1","M2","M3","M4","M5","N1","N2",
  //        "N3","N4","N5","N6","N7","O1","O2","O3","O4","O5","O6","O7",
  //        "P1","P2","P3","P4","P5","Q1","X "
  //        in ascendant order.
  //
  //  ALW -> Vacanci life time (s) of specified shell in each element.
  // 
  // NSHT -> Number of shells in each element.
  //
  //  IFI -> Shell occupation number (f_i).
  //
  //  CP0 -> One-electron Compton profile at P_z=0 (atomic units).
  
  double EB[nElements][30], ALW[nElements][30], CP0[nElements][30];
  int IFI[nElements][30], IKS[nElements][30], NSHT[nElements];  

//  Physical data for the elements Z=1-99.
  
//  ************  Chemical symbols of the elements.  
  char LASYMB[nElements][3];
  
  //  ************  Atomic weights (mean relative atomic masses).

  double ATW[nElements];

  //  ************  Mean excitation energies of the elements (eV).

  double EPX[nElements];

  //  ************  Pair-production cross section parameters.

  //  ****  Screening parameter (R mc/hbar).
  double RSCR[nElements];
  
  //  ****  Asymptotic triplet contribution (eta).
  double ETA[nElements];

  //  ****  Photoelectric cross sections.
  // EPH -> Energy grid point for photoelectric cross section for loaded
  //        elements. The position where the data begins and end for each
  //        element will be specified by IPHF and IPHL respectively.
  //
  // XPH -> Cross section for each grid point for photoelectric cross section for loaded
  //        elements. The position where the data begins and end for each
  //        element will be specified by IPHF and IPHL respectively.
  //
  // NCUR-> Total number of elements stored in EPH and XPH.
  //
  double EPH[constants::NTP], XPH[constants::NTP][17];
  unsigned int NCUR;

  // Positions for each loaded element.
  // NPHS -> Number of shells in each element
  // IPHF -> First position where Photoelectric cross section
  //         data begins for each element in arrays EPH and XPH.
  //
  // IPHL -> Last position of Photoelectric cross section
  //         data for each element in arrays EPH and XPH.
  int IPHF[nElements], IPHL[nElements], NPHS[nElements];

  // **** Atomic relaxation
  //
  //  De-excitation data for the loaded elements are stored in the common
  //  block /CRELAX/, in a form designed to minimize the amount of memory
  //  and to facilitate the random sampling. The quantities in the common
  //  block are the following:
  //  IFIRST(99,16) ... de-excitation data for a vacancy in the shell IS of
  //     the element IZ start at the position K=IFIRST(IZ,IS) in the
  //     storage arrays. The allowed values for IS are 1 to 16 (K shell
  //     and L, M and N subshells).
  //  ILAST(99,16) ... the de-excitation data for a vacancy in the shell
  //     IS of the element IZ end at the position K=ILAST(IZ,IS) in the
  //     storage arrays.  
  //  IS1(K), IS2(K) -> shells that are active in the transition (see the
  //     shell label code below). For radiative transitions, IS2(K)=0.
  //  P(K) -> relative probability for the transition IS-IS1(K)-IS2(K).
  //  ET(K) -> energy of the secondary particle emitted in the transition.
  //  F(K), IAL(K) -> cutoff and alias values (Walker's sampling method).  
  //  ---------------------------------------------------------------------
  //  Label code IS for electron shells:
  //      1 = K  (1s1/2),     11 = N2 (4p1/2),     21 = O5 (5d5/2),
  //      2 = L1 (2s1/2),     12 = N3 (4p3/2),     22 = O6 (5f5/2),
  //      3 = L2 (2p1/2),     13 = N4 (4d3/2),     23 = O7 (5f7/2),
  //      4 = L3 (2p3/2),     14 = N5 (4d5/2),     24 = P1 (6s1/2),
  //      5 = M1 (3s1/2),     15 = N6 (4f5/2),     25 = P2 (6p1/2),
  //      6 = M2 (3p1/2),     16 = N7 (4f7/2),     26 = P3 (6p3/2),
  //      7 = M3 (3p3/2),     17 = O1 (5s1/2),     27 = P4 (6d3/2),
  //      8 = M4 (3d3/2),     18 = O2 (5p1/2),     28 = P5 (6d5/2),
  //      9 = M5 (3d5/2),     19 = O3 (5p3/2),     29 = Q1 (7s1/2),
  //     10 = N1 (4s1/2),     20 = O4 (5d3/2),     30 = outer shells.
  //  ---------------------------------------------------------------------
  
  double P[constants::NRX], ET[constants::NRX], F[constants::NRX];
  int IS0[constants::NRX], IS1[constants::NRX], IS2[constants::NRX], NRELAX; //NREAX is a substitute of NCUR to avoid double NCUR declaration
  int* IAL=IS0;
  int IFIRST[nElements][16], ILAST[nElements][16];
 
  pen_elementDataBase();
  unsigned int getPHposition(const unsigned int IZ, const double XEL);
  void GPHa0();
  void RELAX0();
  
};


//-------------------
// Materials
//-------------------

class pen_material : public abc_material{

public:
  
  //  ****  Simulation parameters (must be defined before calling the
  //        initialisation subroutine PEINIT).

  
  //  ----  Electron/positron transport parameters.
  double C1;
  double C2;
  double WCC;
  double WCR;
  
  //  ****  Rayleigh scattering.
  double ERA[constants::NEX];
  int IED[constants::NEGP], IEU[constants::NEGP], NE;  
  
  //  ****  Global information on the material system (defined by
  //        subroutine PEINIT).
  
  // Minimum cut energy for each material
  double ECUTR;
  
  //  ****  Composition data.

  /*
   *   STF -> Store the ratio atmos/molecule of each element for each material.
   *   IZ  -> Store the atomic number of each element for each material.
   *   ZT  -> Store the 'total' atomic number of each material. This value will
   *         be calculated as follows: 
   *               ZT[M-1] = ZT[M-1]+STF[M-1][I]*IZ[M-1][I]
   *   AT  -> Store the 'total' atomic weight calculated as follows,
   *               AT[M-1] = AT[M-1]+ATW[IZZ-1]*STF[M-1][I]
   *         where ATW is the atomic weight (mean relative atomic masses) of
   *         corresponding element.
   *   RHO -> Store the density of each material in g/cm^3
   *  VMOL -> Store the molecular density (Na*RHO/AT) where NA is the avogadro number.
   * NELEM -> Store the number of elements in each material.
   */
  
  double STF[30], ZT, AT, RHO, VMOL;
  int IZ[30], NELEM;


  //  ****  Electron/positron inelastic collisions. Uses an osocillator model
  //
  //
  //  EXPOT ->  Mean excitation energy
  //  OP2   ->  (Plasma energy)^2
  //  NOSC  ->  Number of oscillators
  //  F     -> ?   
  //  UI    ->  oscillator energy?   
  //  WRI   -> ?   
  //  KZ    -> ?   
  //  KS    -> ?   
  
  double EXPOT, OP2, F[constants::NO], UI[constants::NO], WRI[constants::NO];
  int KZ[constants::NO], KS[constants::NO], NOSC;

  //  ****  Electron inelastic coll. and inner-shell ionisation tables.
  //

  // hard inelastic collisions of electrons.
  
  // EINAC -> Oscillators collision probability for each energy bin ?
  // IEIN  -> Oscillators indexes?
  // NEIN  -> Number of oscillators?
  
  double EINAC[constants::NEGP][constants::NO];
  int IEIN[constants::NO], NEIN;  

  // inner-shell ionisation by electron impact.

  // ESIAC -> Oscillators ionization probability for each energy bin ?
  // IESI  -> Oscillators indexes?
  // NESI  -> Number of oscillators?
  
  double ESIAC[constants::NEGP][constants::NO];
  int IESI[constants::NO], NESI;  

  // hard inelastic collisions of positrons.
  
  // PINAC -> Oscillators collision probability for each energy bin ?
  // IPIN  -> Oscillators indexes?
  // NPIN  -> Number of oscillators?
  
  double PINAC[constants::NEGP][constants::NO];
  int IPIN[constants::NO], NPIN;  

  // inner-shell ionisation by positron impact.

  // PSIAC -> Oscillators ionization probability for each energy bin ?
  // IPSI  -> Oscillators indexes?
  // NPSI  -> Number of oscillators?
  
  double PSIAC[constants::NEGP][constants::NO];
  int IPSI[constants::NO], NPSI;  

  //  ****  Compton scattering.

  // NOSCCO -> Number of shells
  // UICO   -> Shell binding energy?
  // FOC    ->
  // FJ0    ->
  // PTRSH  ->
  // KZCO   ->
  // KSCO   ->
  
  double FCO[constants::NOCO], UICO[constants::NOCO], FJ0[constants::NOCO], PTRSH[constants::NOCO];
  int KZCO[constants::NOCO], KSCO[constants::NOCO], NOSCCO;    

  
  //  ****  Electron simulation tables.

  double TSTPE[constants::NEGP], TSTRE[constants::NEGP], TRL1E[constants::NEGP], TRL2E[constants::NEGP];
  
  //
  // Hard interaction variables
  //
  // define "l" as mean free path in cm
  //
  // SEHEL -> log(1/l) for elastic collisions
  // SEHIN -> log(1/l) for inelastic collisions
  // SEISI -> log(1/l) for inner-shell ionisation (isi)
  // SEHBR -> log(1/l) for bremsstrahlung emision
  //
  // SETOT -> log(1/l) total
  //
  double SEHEL[constants::NEGP], SEHIN[constants::NEGP],
    SEISI[constants::NEGP], SEHBR[constants::NEGP],
    SETOT[constants::NEGP];

  // Precalculate distance between grid points
  //
  // Example: DSEHEL[i] = SEHEL[i+1]-SEHEL[i]
  //
  double DSEHEL[constants::NEGP], DSEHIN[constants::NEGP],
    DSEISI[constants::NEGP], DSEHBR[constants::NEGP],
    DSETOT[constants::NEGP];
  
  //
  // Soft interaction variables
  //  
  double CSTPE[constants::NEGP], RSTPE[constants::NEGP],
    DEL[constants::NEGP], W1E[constants::NEGP],
    W2E[constants::NEGP], DW1EL[constants::NEGP],
    DW2EL[constants::NEGP], RNDCE[constants::NEGP],
    AE[constants::NEGP], BE[constants::NEGP],
    T1E[constants::NEGP], T2E[constants::NEGP];

  //  ****  Positron simulation tables.

  double TSTPP[constants::NEGP], TSTRP[constants::NEGP], TRL1P[constants::NEGP], TRL2P[constants::NEGP];  
  //
  // Hard interaction variables
  //
  // define "l" as mean free path in cm
  //
  // SPHEL -> log(1/l) for elastic collisions
  // SPHIN -> log(1/l) for inelastic collisions
  // SPISI -> log(1/l) for inner-shell ionisation (isi)
  // SPHBR -> log(1/l) for bremsstrahlung emision
  // SPAN  -> log(1/l) for positron annihilation
  //
  // SPTOT -> log(1/l) total
  //
  double SPHEL[constants::NEGP], SPHIN[constants::NEGP],
    SPISI[constants::NEGP], SPHBR[constants::NEGP],
    SPAN[constants::NEGP], SPTOT[constants::NEGP];

  // Precalculate distance between grid points
  //
  // Example: DSEHEL[i] = SEHEL[i+1]-SEHEL[i]
  //
  double DSPHEL[constants::NEGP], DSPHIN[constants::NEGP],
    DSPISI[constants::NEGP], DSPHBR[constants::NEGP],
    DSPAN[constants::NEGP], DSPTOT[constants::NEGP];

  
  //
  // Soft interaction variables
  //
  double CSTPP[constants::NEGP], RSTPP[constants::NEGP],
    W1P[constants::NEGP], W2P[constants::NEGP], DW1PL[constants::NEGP],
    DW2PL[constants::NEGP], RNDCP[constants::NEGP], AP[constants::NEGP],
    BP[constants::NEGP], T1P[constants::NEGP], T2P[constants::NEGP];  

  
  //  ****  Photon simulation tables.
  //
  // define "l" as mean free path in cm
  //
  // SGRA  -> 1/l for Rayleigh interaction
  // SGCO  -> log(1/l) for Compton interaction
  // SGPH  -> 1/l for photoelectric interaction
  // SGPP  -> log(1/l) for pair production
  //
  double SGRA[constants::NEGP], SGCO[constants::NEGP],
    SGPH[constants::NEGP], SGPP[constants::NEGP];
  
  //  ****  Photon simulation tables.
  double TRIP[constants::NEGP];

  // Precalculate distance between grid points
  //
  // Example: DSGRA[i] = SGRA[i+1]-SGRA[i]
  //
  double DSGRA[constants::NEGP], DSGCO[constants::NEGP],
    DSGPH[constants::NEGP], DSGPP[constants::NEGP],
    DTRIP[constants::NEGP];
  
  //  ****  Bremsstrahlung emission.
  double PBCUT[constants::NEGP], WBCUT[constants::NEGP],
    PDFB[constants::NEGP][constants::NBW], DPDFB[constants::NEGP][constants::NBW],
    PACB[constants::NEGP][constants::NBW], ZBR2;

  //  ****  Rayleigh scattering.
  double FF[constants::NQ], XSRA[constants::NEX];

  // **** Bremsstrahlung angular distribution parameters
  static const unsigned int NET = 7;
  static const unsigned int NKT = 41;
  double BP1[NET][NKT][4], BP2[NET][NKT][4], ZBEQ;

  static const unsigned int NP_RSC = 150;
  const unsigned int NP2_RSC = 2*NP_RSC; //(NP2 = NP+NP)
  const unsigned int NPM1_RSC = NP_RSC-1; //(NPM1 = NP-1)
  
  // **** coherent (Rayleigh) scattering
  double QRA[NP_RSC], PRA[NP_RSC], DPRA[NP_RSC], ARA[NP_RSC], BRA[NP_RSC], PMAX[constants::NEGP];
  int ITLRA[NP_RSC], ITURA[NP_RSC];    

  // **** electron-positron pair and triplet production
  double ZEQPP, F0[2], BCB;  


  // **** electron hard elastic events

  static const unsigned int NP_P = 128;
  static const unsigned int NPM1_P = NP_P-1;
  
  double XSE[NP_P][constants::NEGP], PSE[NP_P][constants::NEGP], ASE[NP_P][constants::NEGP], BSE[NP_P][constants::NEGP];
  int ITLE[NP_P][constants::NEGP], ITUE[NP_P][constants::NEGP];
  
  // **** positron hard elastic events
  double XSP[NP_P][constants::NEGP], PSP[NP_P][constants::NEGP], ASP[NP_P][constants::NEGP], BSP[NP_P][constants::NEGP];
  int ITLP[NP_P][constants::NEGP], ITUP[NP_P][constants::NEGP];  

  // **** CELSEP
  double EELMAX, PELMAX, RNDCEd[constants::NEGP], RNDCPd[constants::NEGP];  

  //  ****  Particle range in material.
  double RANGE[constants::nParTypes][constants::NEGP], RANGEL[constants::nParTypes][constants::NEGP];      

  //  **** electron bremss cross section
  double P0[constants::NEGP][constants::NBW];


  //  ****  Electron and positron radiative yields.
  double EBRY[constants::NEGP], PBRY[constants::NEGP];  

  
  //  ****  Rayleigh scattering of photons.
  double FACTE, Q2MAX;
  int MOM;

  //  **** CGRA02
  double QQ[constants::NQ], AR[constants::NQ], BR[constants::NQ],
    CR[constants::NQ], DR[constants::NQ], FF0, QQM;    

  
  pen_material();
  ~pen_material(){}
  
  void load(FILE* IRD,
	    FILE* IWR,
	    initStructs& initStore,
	    pen_elementDataBase& elements,
	    pen_logGrid& grid,
	    int INFO = 0);
  
};

//-----------------------------------------------
// PENELOPE initialization functions
//-----------------------------------------------


void RELAXR(pen_elementDataBase& elements,
	    FILE* IRD,
	    FILE* IWR,
	    int INFO);

void ESIaR(pen_material& mat,
	   pen_elementDataBase& elemDB,
	   CESI0& esi0,
	   FILE* IRD,
	   FILE* IWR,
	   int INFO,
	   const double* DLEMP);

void PSIaR(pen_elementDataBase& elements,
	   pen_material& mat,
	   pen_logGrid& grid,
	   CPSI0& cpsi0,
	   FILE* IRD,
	   FILE* IWR,
	   int INFO);

void EBRaR(pen_material& mat,
	   CEBR01& cebr01,
	   double &WCRM,
	   FILE* IRD,
	   FILE* IWR,
	   int INFO,
	   const double ET[constants::NEGP],
	   const double DLEMP[constants::NEGP]);

void BRaAR(pen_material& mat,
	   FILE* IRD,
	   FILE* IWR,
	   int INFO);

void EINaT(const pen_material& mat,
	   CEIN00& cein00,
	   const double E,
	   double &WCCM,
	   double &XH0,
	   double &XH1,
	   double &XH2,
	   double &XS0,
	   double &XS1,
	   double &XS2,
	   double &XT1,
	   double &XT2,
	   double &DELTA);

void EINaT1(const double E,
	    double &UK,
	    double& WK,
	    double DELTA,
	    double &WCCM,
	    double &H0,
	    double &H1,
	    double &H2,
	    double &S0,
	    double &S1,
	    double &S2,
	    double &R0,
	    double &R1,
	    double &R2);

void PINaT(const pen_material& mat,
	   CPIN00& cpin00,
	   const double E,
	   double &WCCM,
	   double &XH0,
	   double &XH1,
	   double &XH2,
	   double &XS0,
	   double &XS1,
	   double &XS2,
	   double &XT1,
	   double &XT2,
	   double &DELTA);

void PINaT1(CPIN01& cpin01,
	    const double E,
	    double &UK,
	    double &WK,
	    double DELTA,
	    double &WCCM,
	    double &H0,
	    double &H1,
	    double &H2,
	    double &S0,
	    double &S1,
	    double &S2,
	    double &R0,
	    double &R1,
	    double &R2);

void EBRaT(const pen_material& mat,
	   const double &E,
	   double &WCRM,
	   double &XH0,
	   double &XH1,
	   double &XH2,
	   double &XS1,
	   double &XS2,
	   const double DLEMP1,
	   const double DLFC,
	   const double* WB);

void PBRaT(const pen_material& mat,
	   const double &E,
	   const double &WCRM,
	   double &XH0,
	   double &XH1,
	   double &XH2,
	   double &XS1,
	   double &XS2,
	   const double DLEMP1,
	   const double DLFC,
	   const double* WB);

void PANaT(const double &E,
	   double &TXS);

void EELaR(pen_material& mat,
	   FILE* IRD,
	   FILE* IWR,
	   int INFO,
	   CEEL00& eel00,
	   CEINTF& eintf,
	   const double* DLEMP,
	   const double* ET);

void EELdR(pen_material& mat,
	   FILE* IRD,
	   FILE* IWR,
	   int &INFO,
	   CDCSEP& dcsep,
	   CRITA& rita,
	   CEINTF& eintf,
	   double* ET);

int GRAaR(pen_logGrid& grid,
	  pen_material& mat,
	  CRITA& rita,
	  FILE* IRD,
	  FILE* IWR,
	  int &INFO);

void GPHaR(pen_material& mat,
	   pen_elementDataBase& elemDB,
	   FILE* IRD,
	   FILE* IWR,
	   int INFO,
	   CGPH01& gph01,
	   const double* DLEMP);

void GRAaTI(const pen_material& mat,
	    const double DLEMP1,
	    const double DLFC,
	    double &E,
	    double &ECS);

void GPPa0(pen_material& mat,
	   const pen_elementDataBase& elemDB);

double EINaDS(double RMU,
	      void* arg);

double PINaDS(double RMU,
	      void* arg);

void EELa0(double &XS0,
	   double &XS1,
	   double &XS2,
	   double &XS0H,
	   double &A,
	   double &B,
	   double &RNDC,
	   double &XS1S,
	   double &XS2S);

void DCSEL0(const double E,
	    int IELEC,
	    CDCSEP& dcsep);

double DCSEL(double RMU,
	     void* arg);

double GRAaF2( double Q2,
	       void* arg);


#endif
