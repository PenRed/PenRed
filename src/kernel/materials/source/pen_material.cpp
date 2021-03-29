
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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


#include "pen_material.hh"

//-------------------
// Materials
//-------------------

pen_material::pen_material() : C1(0.01), C2(0.01), WCC(0.01), WCR(0.01), ECUTR(50.0)
{
  //Init default values
  for(unsigned int i = 0 ; i < constants::nParTypes; i++)
    EABS[i] = 50.0;
  
}

void pen_material::load(FILE* IRD,
			FILE* IWR,
			initStructs& initStore,
			pen_elementDataBase& elements,
			pen_logGrid& grid,
			int INFO)
{

  //  This subroutine reads the definition file of material M (unit IRD)
  //  and initialises the simulation routines for this material. Informa-
  //  tion is written on unit IWR, the amount of written information is
  //  determined by the value of INFO.
      
  //Mantain "mat" variable for possible
  //future changes.
  
  pen_material& mat = *this;
  
  char NAME[63], LNAME[63];
  const double FOURPI=4.0*constants::PI;
  
  //  ****  Partial cross sections of individual shells/oscillators.
  CEIN00& cein00 = *initStore.pcein00;
  
  //  ****  Electron inelastic coll. and inner-shell ionisation tables.
  CESIN& cesin = *initStore.pcesin;
  CESI0& cesi0 = *initStore.pcesi0;
  CEBR01& cebr01 = *initStore.pcebr01;
  CPSI0& cpsi0 = *initStore.pcpsi0;
  CEEL00& ceel00 = *initStore.pceel00;
  CEINTF& ceintf = *initStore.pceintf;
  CPIN00& cpin00 = *initStore.pcpin00;
  CPSIN& cpsin = *initStore.pcpsin;
  CPIN01& cpin01 = *initStore.pcpin01;
  CDCSEP& cdcsep = *initStore.pcdcsep;
  CRITA& crita = *initStore.pcrita;
  CGPH01& cgph01 = *initStore.pcgph01;
  
  ///////////////////////////////////////
  
  //  ****  Positron inelastic coll. and inner-shell ionisation tables.
  double STFI[constants::NO], STFO[constants::NO];
  int INOUT[constants::NO];

  //  ****  Auxiliary arrays.
  double EIT[constants::NEGP], EITL[constants::NEGP], FL[constants::NEGP], F1[constants::NEGP], F2[constants::NEGP], F3[constants::NEGP], F4[constants::NEGP], RADY[constants::NEGP], RADN[constants::NEGP];
  
  //  ****  Rescaling of calculated MFPs.
  //  When ISCALE is set equal to 1, the program rescales the calculated
  //  total cross sections and MFPs to reproduce the cross sections and
  //  stopping powers read from the input material data file. To avoid this
  //  rescaling, set ISCALE=0.

  int ISCALE=1;
    
  //  ************  Material characteristics.

  // Calculate ECUTR
  if(EABS[PEN_ELECTRON] < EABS[PEN_PHOTON]){ ECUTR = EABS[PEN_ELECTRON];}
  else{ ECUTR = EABS[PEN_PHOTON];}  
  
  strcpy(LNAME, " PENELOPE (v. 2018)  Material data file ...............");
  fscanf(IRD,"%55c%*[^\n]",NAME);
  //Append end of string chars
  NAME[55] = '\0';
  getc(IRD);
  
  if(strcmp(NAME, LNAME) != 0)
  {
    fprintf(IWR, "\n I/O error. Corrupt material data file.\n      The first line is: %55s\n      ... and should be: %55s\n",NAME, LNAME);
    penError(ERR_PEMATR_CORRUPT_MAT_FILE);
    return;
  }
  fprintf(IWR, "%55s\n", NAME);

  fscanf(IRD,"%*11c%62[^\n]%*[^\n]",LNAME);
  getc(IRD);
  fprintf(IWR, " Material: %62s\n", LNAME);
  
  fscanf(IRD, "%*15c%lf%*[^\n]", &(mat.RHO));
  getc(IRD);
  
  fprintf(IWR, " Mass density =%15.8E g/cm**3\n", mat.RHO);
  
  fscanf(IRD, "%*37c%3d%*[^\n]", &mat.NELEM);
  getc(IRD);
  fprintf(IWR, " Number of elements in the molecule = %2d\n", mat.NELEM);fflush(IWR);
  if(mat.NELEM > 30){ fprintf(IWR, "To many elements in material %s\n", NAME); penError(ERR_PEMATR_2_MANY_ELEMENTS);
    return;}
  
  mat.ZT = 0.0;
  mat.AT = 0.0;
  
  for(int I = 0; I < mat.NELEM; I++)
  {
    fscanf(IRD, "%*18c%d,%*19c%lf%*[^\n]", &mat.IZ[I], &mat.STF[I]);
    getc(IRD);
    fprintf(IWR, "    Element: %s (Z=%2d), atoms/molecule =%15.8E\n", elements.LASYMB[mat.IZ[I]-1], mat.IZ[I], mat.STF[I]);fflush(IWR);
      
    mat.ZT = mat.ZT+mat.STF[I]*mat.IZ[I];
    int IZZ = mat.IZ[I];
    mat.AT = mat.AT+elements.ATW[IZZ-1]*mat.STF[I];
  }
  mat.VMOL = constants::AVOG*mat.RHO/mat.AT;

  if(INFO >= 2){ fprintf(IWR, "\n Molecular density = %15.8E 1/cm**3\n", mat.VMOL);}
  
  mat.OP2 = FOURPI*mat.ZT*mat.VMOL*pow(constants::A0B,3)*pow(constants::HREV,2);
  double OMEGA = sqrt(mat.OP2);

  mat.DEN = mat.RHO;
  mat.RDEN=1.0/mat.DEN;

  if(INFO >= 2)
  {
    fprintf(IWR, "\n *** Electron/positron inelastic scattering.\n");
    fprintf(IWR, " Plasma energy = %15.8E eV\n", OMEGA);
  }
  
  fscanf(IRD, "%*25c%lf%*[^\n]", &mat.EXPOT);
  getc(IRD);
  fprintf(IWR, " Mean excitation energy =%15.8E eV\n", mat.EXPOT);

  //  ****  E/P inelastic collisions.

  fscanf(IRD, "%*24c%3d%*[^\n]", &mat.NOSC);
  getc(IRD);
  if(INFO >= 2 || (unsigned)mat.NOSC > constants::NO){ fprintf(IWR, " Number of oscillators =%4d\n", mat.NOSC);}
  if((unsigned)mat.NOSC > constants::NO){ penError(ERR_PEMATR_2_MANY_OSCILLATORS);
    return;}
  if(INFO >= 2){ fprintf(IWR, "\n           Fi            Ui (eV)         Wi (eV)      KZ  KS\n ------------------------------------------------------------\n");}
  double EXPT = 0.0;
  for(int I = 0; I < mat.NOSC; I++)
  {
    fscanf(IRD , "%*4c%lf %lf %lf %d %d%*[^\n]", &mat.F[I], &mat.UI[I], &mat.WRI[I], &mat.KZ[I], &mat.KS[I]);
    getc(IRD);
      
    if(mat.UI[I] < 1.0E-3)
    {
      mat.UI[I] = 0.0;
    }
    if(INFO >= 2){ fprintf(IWR, "%4d%16.8E%16.8E%16.8E%4d%4d\n", I+1, mat.F[I], mat.UI[I], mat.WRI[I], mat.KZ[I], mat.KS[I]);}
    EXPT = EXPT+mat.F[I]*log(mat.WRI[I]);
  }
  EXPT = exp(EXPT/mat.ZT);

  if(fabs(EXPT-mat.EXPOT) > 1.0E-6*mat.EXPOT)
  {
    fprintf(IWR, "EXPT      =%.5E\nEXPOT (M) =%.5E\nInconsistent oscillator data.\n", EXPT, mat.EXPOT);
    penError(ERR_PEMATR_INCONSISTENT_OSCILLATOR);
    return;
  }

  //  ****  Compton scattering.

  fscanf(IRD,"%*19c%3d%*[^\n]", &mat.NOSCCO);
  getc(IRD);
  if(INFO >= 2 || (unsigned)mat.NOSCCO > constants::NOCO)
  {
    fprintf(IWR, "\n *** Compton scattering (Impulse Approximation).\n");
    fprintf(IWR, " Number of shells =%4d\n", mat.NOSCCO);
    if((unsigned)mat.NOSCCO > constants::NOCO){ penError(ERR_PEMATR_2_MANY_SHELLS);
      return;}
  }
  if(INFO >= 2){ fprintf(IWR, "\n           Fi            Ui (eV)         Ji(0)        KZ  KS\n ------------------------------------------------------------\n");}

  for(int I = 0; I < mat.NOSCCO; I++)
  {
    fscanf(IRD, "%*4c%lf %lf %lf %d %d%*[^\n]", &mat.FCO[I], &mat.UICO[I], &mat.FJ0[I], &mat.KZCO[I], &mat.KSCO[I]);
    getc(IRD);
    if(INFO >= 2){ fprintf(IWR, "%4d %15.8E %15.8E %15.8E%4d%4d\n", I+1, mat.FCO[I], mat.UICO[I], mat.FJ0[I], mat.KZCO[I], mat.KSCO[I]);}
    mat.FJ0[I] = mat.FJ0[I]*constants::SL;
  }

  //  ************  Atomic relaxation data.

  for(int I = 0; I < mat.NELEM; I++)
  {
    RELAXR(elements, IRD, IWR, INFO);
    if(penGetError() != PEN_SUCCESS){
      return;}
  }

  //  ****  Inner-shell ionisation by electron and positron impact.

  ESIaR(mat,elements,cesi0,IRD,IWR,INFO,grid.DLEMP);
  if(penGetError() != PEN_SUCCESS){
    return;}
  PSIaR(elements,mat,grid,cpsi0,IRD,IWR,INFO);
  if(penGetError() != PEN_SUCCESS){
    return;}

  //  ****  Electron and positron interaction properties.

  //  ****  Scaled bremsstrahlung x-section.
  int IBREMS;
  double WCRM;
  if(mat.WCR >= 0.0)
  {
    if(mat.WCR < 10.0){mat.WCR = 10.0;}     
    WCRM = mat.WCR;
    IBREMS = 1;
  }
  else
  {
    WCRM = 10.0;
    mat.WCR = 10.0;
    IBREMS = 0;
  }
  EBRaR(mat,cebr01,WCRM, IRD, IWR, INFO, grid.ET, grid.DLEMP);
  if(penGetError() != 0){
    return;}
  //  ****  Bremsstrahlung angular distribution.
  BRaAR(mat, IRD, IWR, INFO);
  if(penGetError() != 0){
    return;}

  //  ****  Stopping powers.
  int NDATA;
  fscanf(IRD, "%*58c%4d%*[^\n]", &NDATA);
  getc(IRD);
  
  if(INFO >= 2){ fprintf(IWR, "\n *** Stopping powers for electrons and positrons,  NDATA =%4d\n", NDATA);}
  if((unsigned)NDATA > constants::NEGP){ penError(ERR_PEMATR_2_MANY_DP_1);
    return;}
  if(INFO >= 2){ fprintf(IWR, "\n  Energy     Scol,e-     Srad,e-     Scol,e+     Srad,e+\n   (eV)     (MeV/mtu)   (MeV/mtu)   (MeV/mtu)   (MeV/mtu)\n ----------------------------------------------------------\n");}

  for(int I = 0; I < NDATA; I++)
  {
    fscanf(IRD, "%lf %lf %lf %lf %lf%*[^\n]", &EIT[I], &F1[I], &F2[I], &F3[I], &F4[I]);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR, "%10.3E%12.5E%12.5E%12.5E%12.5E\n", EIT[I], F1[I], F2[I], F3[I], F4[I]);}
    EITL[I] = log(EIT[I]);
  }

  //  ************  Inelastic x-section for electrons.

  double WCCM = mat.WCC;
  double XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2;
  double DELTA;
  double DIFMAX=0.0;
  double STPI;
  double STPC;
  
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] >= grid.EL && EIT[I] <= grid.EU)
    {
      STPI = F1[I]*mat.RHO*1.0E6;  // Collision stopping power.
      EINaT(mat, cein00, EIT[I], WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA);
      STPC = (XS1+XH1)*mat.VMOL;
      double AuxDouble = 100.0*fabs(STPC-STPI)/STPI;
      if(DIFMAX < AuxDouble){DIFMAX = AuxDouble;}
    }
  }

  int ICOR;
  if(DIFMAX > 1.0E-3 && ISCALE == 1)
  {
    ICOR=1;
    for(int I = 0; I < NDATA; I++)
    {
      FL[I] = log(F1[I]*mat.RHO*1.0E6);
    }
    SPLINE(EITL, FL, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, NDATA);
    if(penGetError() != 0){
      return;}
  }
  else
  {
    ICOR=0;
  }

  //  **** Classify inner and outer shells.

  int NI = 0;
  int NS = 0;
  double EBIN;
  
  for(unsigned int IO = 0; IO < constants::NO; IO++)
  {
    STFI[IO]=0.0;
    STFO[IO]=0.0;
  }

  for(int KO = 0; KO < mat.NOSC; KO++)
  {
    int IZZ = mat.KZ[KO];
    int ISH = mat.KS[KO];
    if(IZZ < 3 || ISH > 16 || mat.UI[KO] < 50.0)
    {
      NI = NI+1;
      mat.IEIN[NI-1] = KO+1;
      cesin.ISIE[KO] = -NI;
      INOUT[NI-1] = 0;
      for(int IEL = 0; IEL < mat.NELEM; IEL++)
      {
        if(IZZ == mat.IZ[IEL]){ STFO[NI-1] = mat.STF[IEL];}
      }
    }
    else if(ISH <= cesi0.NSESI[IZZ-1])
    {
      EBIN = elements.EB[IZZ-1][ISH-1];
      if(EBIN > mat.ECUTR)
      {
        NS = NS+1;
        mat.IESI[NS-1] = KO+1;
        cesin.ISIE[KO] = NS;
        for(int IEL = 0; IEL < mat.NELEM; IEL++)
        {
          if(IZZ == mat.IZ[IEL]){ STFI[NS-1] = mat.STF[IEL];}
      
        }
      }   
      else
      {
        NI = NI+1;
        mat.IEIN[NI-1] = KO+1;
        cesin.ISIE[KO] = -NI;
        INOUT[NI-1] = 1;
        for(int IEL = 0; IEL < mat.NELEM; IEL++)
        {
          if(IZZ == mat.IZ[IEL]){ STFO[NI-1] = mat.STF[IEL];}
        }
      }
    }
    else
    {
      NI = NI+1;
      mat.IEIN[NI-1] = KO+1;
      cesin.ISIE[KO] = -NI;
      INOUT[NI-1] = 0;
      for(int IEL = 0; IEL < mat.NELEM; IEL++)
      {
        if(IZZ == mat.IZ[IEL]){ STFO[NI-1] = mat.STF[IEL];}
      }
    }
  }
  mat.NEIN = NI;
  mat.NESI = NS;

  //  ****  Simulation arrays.
  double EE;
  double FACT;
  double EC;
  
  double PCSI;
  double STOT;
  double DFERMI;
  
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    EE = grid.ET[I];       
    EINaT(mat, cein00, EE, WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA);
    STPC = (XS1+XH1)*mat.VMOL;
    if(ICOR == 1)
    {
      unsigned int J;    
      EC = grid.DLEMP[I];
      FINDI(EITL, EC, NDATA, J);
      STPI = exp(ceel00.A[J]+EC*(ceel00.B[J]+EC*(ceel00.C[J]+EC*ceel00.D[J])));
      FACT = STPI/STPC;
    }
    else
    {
      FACT = 1.0;
    }
    mat.DEL[I] = DELTA;
    mat.CSTPE[I] = STPC*FACT;

    double XS1SI = 0.0;
    double XS2SI = 0.0;
    double XT1SI = 0.0;
    double XT2SI = 0.0;
    double STPSI = 0.0;
    double STRSI = 0.0;
    double XS1IN = 0.0;
    double XS2IN = 0.0;
    double XT1IN = 0.0;
    double XT2IN = 0.0;
    double STPIN = 0.0;
    double STRIN = 0.0;

    for(int IO = 0; IO < mat.NEIN; IO++)
    {
      cesin.XSEIN[I][IO] = 0.0;
    }
    for(int IO = 0; IO < mat.NESI; IO++)
    {
      cesin.XSESI[I][IO] = 0.0;
    }

    for(int KO = 0; KO < mat.NOSC; KO++)
    {
      int IO = cesin.ISIE[KO];
    //  ****  Inner shells
      if(IO > 0)
      {
        int IZZ = mat.KZ[KO];
        int ISH = mat.KS[KO];
        EBIN = elements.EB[IZZ-1][ISH-1];
        if(EBIN < EE)
        {
          int INDC = cesi0.IESIF[IZZ-1]-1;
          PCSI = exp(cesi0.XESI[INDC+I][ISH-1])*STFI[IO-1];
          if(PCSI > 1.0E-35)
          {
            double H0, H1, H2, S0, S1, S2, R0, R1, R2;
            EINaT1(EE, mat.UI[KO], mat.WRI[KO], 0.0, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
            STOT = mat.F[KO]*(S0+H0);
            if(STOT > 1.0E-35)
            {
              DFERMI = (cein00.SES0[KO]+cein00.SEH0[KO])/STOT;
              PCSI = PCSI*DFERMI;
              cesin.XSESI[I][IO-1] = PCSI;
              double RNFO = PCSI/(cein00.SES0[KO]+cein00.SEH0[KO]);
              STPSI = STPSI+(cein00.SES1[KO]+cein00.SEH1[KO])*RNFO;
              STRSI = STRSI+(cein00.SES2[KO]+cein00.SEH2[KO])*RNFO;
            }
            else
            {
              cesin.XSESI[I][IO-1]=PCSI;
              STPSI = STPSI+PCSI*EBIN;
              STRSI = STRSI+PCSI*pow(EBIN,2);
            }
          }
        }
      }
      else
      {
        //  ****  Outer shells
        IO = -IO;
        if(INOUT[IO-1] == 1)
        {
          double H0, H1, H2, S0, S1, S2, R0, R1, R2;
          int IZZ = mat.KZ[KO];
          int ISH = mat.KS[KO];
          EBIN = elements.EB[IZZ-1][ISH-1];
          int INDC = cesi0.IESIF[IZZ-1]-1;
          PCSI = exp(cesi0.XESI[INDC+I][ISH-1])*STFO[IO-1];
          EINaT1( EE, mat.UI[KO], mat.WRI[KO], 0.0, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
          STOT = mat.F[KO]*(S0+H0);
          if(STOT > 1.0E-35)
          {
            DFERMI = (cein00.SES0[KO]+cein00.SEH0[KO])/STOT;
            PCSI = PCSI*DFERMI;
            double RNFO = PCSI/(cein00.SES0[KO]+cein00.SEH0[KO]);
            XS1SI = XS1SI+cein00.SES1[KO]*RNFO;
            XS2SI = XS2SI+cein00.SES2[KO]*RNFO;
            XT1SI = XT1SI+cein00.SET1[KO]*RNFO;
            XT2SI = XT2SI+cein00.SET2[KO]*RNFO;
            cesin.XSEIN[I][IO-1] = cein00.SEH0[KO]*RNFO;
            STPSI = STPSI+(cein00.SES1[KO]+cein00.SEH1[KO])*RNFO;
            STRSI = STRSI+(cein00.SES2[KO]+cein00.SEH2[KO])*RNFO;
          }
          else
          {
            XS1SI = XS1SI+cein00.SES1[KO];
            XS2SI = XS2SI+cein00.SES2[KO];
            XT1SI = XT1SI+cein00.SET1[KO];
            XT2SI = XT2SI+cein00.SET2[KO];
            cesin.XSEIN[I][IO-1] = cein00.SEH0[KO];
            STPSI = STPSI+(cein00.SES1[KO]+cein00.SEH1[KO]);
            STRSI = STRSI+(cein00.SES2[KO]+cein00.SEH2[KO]);
          }
        }
        else
        {
          XS1IN = XS1IN+cein00.SES1[KO];
          XS2IN = XS2IN+cein00.SES2[KO];
          XT1IN = XT1IN+cein00.SET1[KO];
          XT2IN = XT2IN+cein00.SET2[KO];
          cesin.XSEIN[I][IO-1] = cein00.SEH0[KO];
          STPIN = STPIN+(cein00.SES1[KO]+cein00.SEH1[KO]);
          STRIN = STRIN+(cein00.SES2[KO]+cein00.SEH2[KO]);
        }
      }
    }

    double SSIT = 0.0;
    if(mat.NESI > 0)
    {
      for(int IO = 0; IO < mat.NESI; IO++)
      {
        mat.ESIAC[I][IO] = SSIT;
        SSIT = SSIT+cesin.XSESI[I][IO];
      }
      double Aux_Double = SSIT*mat.VMOL;
      if(Aux_Double < 1.0E-35){ Aux_Double = 1.0E-35;}
      
      mat.SEISI[I] = log(Aux_Double);
    }
    else
    {
      mat.SEISI[I] = log(1.0E-35);
    }
    if(SSIT > 1.0E-35)
    {
      for(int IO = 0; IO < mat.NESI; IO++)
      {
        mat.ESIAC[I][IO] = mat.ESIAC[I][IO]/SSIT;
      }
    }
    else
    {
      for(int IO = 0; IO < mat.NESI; IO++)
      {
        mat.ESIAC[I][IO] = 1.0;
      }
    }

    STPSI = STPSI*mat.VMOL;
    STPIN = STPIN*mat.VMOL;
    double FNORM = (mat.CSTPE[I]-STPSI)/STPIN;
 
    double SINT = 0.0;
    if(mat.NEIN == 0){ penError(ERR_PEMATR_PHMART);
      return;}
    for(int IO = 0; IO < mat.NEIN; IO++)
    {
      if(INOUT[IO] == 0){ cesin.XSEIN[I][IO] = FNORM*cesin.XSEIN[I][IO];}
      mat.EINAC[I][IO] = SINT;
      SINT = SINT+cesin.XSEIN[I][IO];
    }
    if(SINT > 1.0E-35)
    {
      for(int IO = 0; IO < mat.NEIN; IO++)
      {
        mat.EINAC[I][IO] = mat.EINAC[I][IO]/SINT;
      }
    }
    else
    {
      for(int IO = 0; IO < mat.NEIN; IO++)
      {
        mat.EINAC[I][IO] = 1.0;
      }
    }
        
    mat.SEHIN[I] = SINT*mat.VMOL;
    if(mat.SEHIN[I] < 1.0E-35){ mat.SEHIN[I] = 1.0E-35;}
    mat.SEHIN[I] = log(mat.SEHIN[I]); 
    mat.TSTRE[I] = (STRSI+FNORM*STRIN)*mat.VMOL;

    mat.W1E[I] = (XS1SI+FNORM*XS1IN)*mat.VMOL;
    mat.W2E[I] = (XS2SI+FNORM*XS2IN)*mat.VMOL;
    ceintf.T1EI[I] = (XT1SI+FNORM*XT1IN)*mat.VMOL;
    ceintf.T2EI[I] = (XT2SI+FNORM*XT2IN)*mat.VMOL;
  }

  //  ************  Inelastic x-section for positrons.

  DIFMAX = 0.0;
    
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] >= grid.EL && EIT[I] <= grid.EU)
    {
      STPI = F3[I]*mat.RHO*1.0E6;  // Collision stopping power.
      PINaT(mat,cpin00,EIT[I], WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA);
      STPC = (XS1+XH1)*mat.VMOL;
      if(DIFMAX < 100.0*fabs(STPC-STPI)/STPI){ DIFMAX = 100.0*fabs(STPC-STPI)/STPI;}
    }
  }

  if(DIFMAX > 1.0E-3 && ISCALE == 1)
  {
    ICOR = 1;
    for(int I = 0; I < NDATA; I++)
    {
      FL[I] = log(F3[I]*mat.RHO*1.0E6);
    }
    SPLINE( EITL, FL, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, NDATA);
    if(penGetError() != 0){
      return;}
  }
  else
  {
    ICOR = 0;
  }

  //  **** Classify inner and outer shells.

  NI = 0;
  NS = 0;
  for(unsigned int IO = 0; IO < constants::NO; IO++)
  {
    STFI[IO]=0.0;
    STFO[IO]=0.0;
  }
    
  for(int KO = 0; KO < mat.NOSC; KO++)
  {
    int IZZ = mat.KZ[KO];
    int ISH = mat.KS[KO];
    if(IZZ < 3 || ISH > 16 || mat.UI[KO] < 50.0)
    {
      NI = NI+1;
      mat.IPIN[NI-1] = KO+1;
      cpsin.ISIP[KO] = -NI;
      INOUT[NI-1] = 0;
      for(int IEL = 0; IEL < mat.NELEM; IEL++)
      {
        if(IZZ == mat.IZ[IEL]){ STFO[NI-1] = mat.STF[IEL];}
      }
    }
    else if(ISH <= cpsi0.NSPSI[IZZ-1])
    { 
      EBIN = elements.EB[IZZ-1][ISH-1];
      if(EBIN > mat.ECUTR)
      {
        NS = NS+1;
        mat.IPSI[NS-1] = KO+1;
        cpsin.ISIP[KO] = NS;
        for(int IEL = 0; IEL < mat.NELEM; IEL++)
        {
          if(IZZ == mat.IZ[IEL]){ STFI[NS-1] = mat.STF[IEL];}
        }
      }
      else
      {
        NI = NI+1;
        mat.IPIN[NI-1] = KO+1;
        cpsin.ISIP[KO] = -NI;
        INOUT[NI-1] = 1;
        for(int IEL = 0; IEL < mat.NELEM; IEL++)
        {
          if(IZZ == mat.IZ[IEL]){ STFO[NI-1] = mat.STF[IEL];}
        }
      }
    }
    else
    {
      NI = NI+1;
      mat.IPIN[NI-1] = KO+1;
      cpsin.ISIP[KO] = -NI;
      INOUT[NI-1] = 0;
      for(int IEL = 0; IEL < mat.NELEM; IEL++)
      {
        if(IZZ == mat.IZ[IEL]){ STFO[NI-1] = mat.STF[IEL];}
      }
    }
  }
  mat.NPIN = NI;
  mat.NPSI = NS;

  //  ****  Simulation arrays.

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    EE = grid.ET[I];
    PINaT(mat, cpin00, EE, WCCM, XH0, XH1, XH2, XS0, XS1, XS2, XT1, XT2, DELTA);
    STPC = (XS1+XH1)*mat.VMOL;
    if(ICOR == 1)
    {
      unsigned int J;
      EC = grid.DLEMP[I];
      FINDI( EITL, EC, NDATA, J);
      STPI = exp(ceel00.A[J]+EC*(ceel00.B[J]+EC*(ceel00.C[J]+EC*ceel00.D[J])));
      FACT = STPI/STPC;
    }
    else
    {
      FACT=1.0;
    }
    mat.CSTPP[I] = STPC*FACT;

    double XS1SI=0.0;
    double XS2SI=0.0;
    double XT1SI=0.0;
    double XT2SI=0.0;
    double STPSI=0.0;
    double STRSI=0.0;
    double XS1IN=0.0;
    double XS2IN=0.0;
    double XT1IN=0.0;
    double XT2IN=0.0;
    double STPIN=0.0;
    double STRIN=0.0;

    for(int IO = 0; IO < mat.NPIN; IO++)
    {
      cpsin.XSPIN[I][IO] = 0.0;
    }
    for(int IO = 0; IO < mat.NPSI; IO++)
    {
      cpsin.XSPSI[I][IO] = 0.0;
    }

    for(int KO = 0; KO < mat.NOSC; KO++)
    {
      int IO = cpsin.ISIP[KO];
    //  ****  Inner shells
      if(IO > 0)
      {       
        int IZZ = mat.KZ[KO];
        int ISH = mat.KS[KO];
        EBIN = elements.EB[IZZ-1][ISH-1];
        if(EBIN < EE)
        {
          int INDC = cpsi0.IPSIF[IZZ-1]-1;
          PCSI = exp(cpsi0.XPSI[INDC+I][ISH-1])*STFI[IO-1];
          if(PCSI > 1.0E-35)
          {
            double H0, H1, H2, S0, S1, S2, R0, R1, R2;
            PINaT1(cpin01, EE, mat.UI[KO], mat.WRI[KO], 0.0, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
            STOT = mat.F[KO]*(S0+H0);
            if(STOT > 1.0E-35)
            {
              DFERMI = (cpin00.SPS0[KO]+cpin00.SPH0[KO])/STOT;
              PCSI = PCSI*DFERMI;
              cpsin.XSPSI[I][IO-1] = PCSI;
              double RNFO = PCSI/(cpin00.SPS0[KO]+cpin00.SPH0[KO]);
              STPSI = STPSI+(cpin00.SPS1[KO]+cpin00.SPH1[KO])*RNFO;
              STRSI = STRSI+(cpin00.SPS2[KO]+cpin00.SPH2[KO])*RNFO;
            }
            else
            {
              cpsin.XSPSI[I][IO-1] = PCSI;
              STPSI = STPSI+PCSI*EBIN;
              STRSI = STRSI+PCSI*pow(EBIN,2);
            }
          }
        }
      }
      else
      {
          //  ****  Outer shells
        IO = -IO;
        if(INOUT[IO-1] == 1)
        {
          int IZZ = mat.KZ[KO];
          int ISH = mat.KS[KO];
          double H0, H1, H2, S0, S1, S2, R0, R1, R2;
          EBIN = elements.EB[IZZ-1][ISH-1];
          int INDC = cpsi0.IPSIF[IZZ-1]-1;
          PCSI = exp(cpsi0.XPSI[INDC+I][ISH-1])*STFO[IO-1];
          PINaT1( cpin01, EE, mat.UI[KO], mat.WRI[KO], 0.0, WCCM, H0, H1, H2, S0, S1, S2, R0, R1, R2);
          STOT = mat.F[KO]*(S0+H0);
          if(STOT > 1.0E-35)
          {
            DFERMI = (cpin00.SPS0[KO]+cpin00.SPH0[KO])/STOT;
            PCSI = PCSI*DFERMI;
            double RNFO = PCSI/(cpin00.SPS0[KO]+cpin00.SPH0[KO]);
            XS1SI = XS1SI+cpin00.SPS1[KO]*RNFO;
            XS2SI = XS2SI+cpin00.SPS2[KO]*RNFO;
            XT1SI = XT1SI+cpin00.SPT1[KO]*RNFO;
            XT2SI = XT2SI+cpin00.SPT2[KO]*RNFO;
            cpsin.XSPIN[I][IO-1] = cpin00.SPH0[KO]*RNFO;
            STPSI = STPSI+(cpin00.SPS1[KO]+cpin00.SPH1[KO])*RNFO;
            STRSI = STRSI+(cpin00.SPS2[KO]+cpin00.SPH2[KO])*RNFO;
          }
          else
          {
            XS1SI = XS1SI+cpin00.SPS1[KO];
            XS2SI = XS2SI+cpin00.SPS2[KO];
            XT1SI = XT1SI+cpin00.SPT1[KO];
            XT2SI = XT2SI+cpin00.SPT2[KO];
            cpsin.XSPIN[I][IO-1] = cpin00.SPH0[KO];
            STPSI = STPSI+(cpin00.SPS1[KO]+cpin00.SPH1[KO]);
            STRSI = STRSI+(cpin00.SPS2[KO]+cpin00.SPH2[KO]);
          }
        }
        else
        {
          XS1IN = XS1IN+cpin00.SPS1[KO];
          XS2IN = XS2IN+cpin00.SPS2[KO];
          XT1IN = XT1IN+cpin00.SPT1[KO];
          XT2IN = XT2IN+cpin00.SPT2[KO];
          cpsin.XSPIN[I][IO-1] = cpin00.SPH0[KO];
          STPIN = STPIN+(cpin00.SPS1[KO]+cpin00.SPH1[KO]);
          STRIN = STRIN+(cpin00.SPS2[KO]+cpin00.SPH2[KO]);
        }
      }
    }

    double SSIT = 0.0;
    if(mat.NPSI > 0)
    {
      for(int IO = 0; IO < mat.NPSI; IO++)
      {
        mat.PSIAC[I][IO] = SSIT;
        SSIT = SSIT+cpsin.XSPSI[I][IO];
      }
      mat.SPISI[I] = SSIT*mat.VMOL;
      if(mat.SPISI[I] < 1.0E-35){ mat.SPISI[I] = 1.0E-35;}
      mat.SPISI[I] = log(mat.SPISI[I]);
    }
    else
    {
      mat.SPISI[I] = log(1.0E-35);
    }
    if(SSIT > 1.0E-35)
    {
      for(int IO = 0; IO < mat.NPSI; IO++)
      {
        mat.PSIAC[I][IO] = mat.PSIAC[I][IO]/SSIT;
      }
    }
    else
    {
      for(int IO = 0; IO < mat.NPSI; IO++)
      {
        mat.PSIAC[I][IO] = 1.0;
      }
    }

    STPSI = STPSI*mat.VMOL;
    STPIN = STPIN*mat.VMOL;
    double FNORM = (mat.CSTPP[I]-STPSI)/STPIN;

    double SINT = 0.0;
    if(mat.NPIN == 0){ penError(ERR_PEMATR_PHMART);
      return;}
    
    for(int IO = 0; IO < mat.NPIN; IO++)
    {
      if(INOUT[IO] == 0){ cpsin.XSPIN[I][IO] = FNORM*cpsin.XSPIN[I][IO];}
      mat.PINAC[I][IO] = SINT;
      SINT = SINT+cpsin.XSPIN[I][IO];
    }
    if(SINT > 1.0E-35)
    {
      for(int IO = 0; IO < mat.NPIN; IO++)
      {
        mat.PINAC[I][IO] = mat.PINAC[I][IO]/SINT;
      }
    }
    else
    {
      for(int IO = 0;  IO < mat.NPIN; IO++)
      {
        mat.PINAC[I][IO] = 1.0;
      }
    }
        
    mat.SPHIN[I] = SINT*mat.VMOL;
    if(mat.SPHIN[I] < 1.0E-35){ mat.SPHIN[I] = 1.0E-35;}
    mat.SPHIN[I] = log(mat.SPHIN[I]);
    mat.TSTRP[I] = (STRSI+FNORM*STRIN)*mat.VMOL;
    mat.W1P[I] = (XS1SI+FNORM*XS1IN)*mat.VMOL;
    mat.W2P[I] = (XS2SI+FNORM*XS2IN)*mat.VMOL;
    ceintf.T1PI[I] = (XT1SI+FNORM*XT1IN)*mat.VMOL;
    ceintf.T2PI[I] = (XT2SI+FNORM*XT2IN)*mat.VMOL;
  }

  //  ****  Bremsstrahlung x-section for electrons.

  double STPR;
  DIFMAX=0.0;
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] >= grid.EL && EIT[I] < grid.EU)
    {
      STPI = F2[I]*mat.RHO*1.0E6;  // Radiative stopping power.
      EBRaT(mat,EIT[I], WCRM, XH0, XH1, XH2, XS1, XS2, grid.DLEMP1, grid.DLFC, constants::WB);
      if(penGetError() != 0){
	return;}
      STPR=(XS1+XH1)*mat.VMOL;
      if(DIFMAX < fabs(STPR-STPI)/(0.01*STPI)){ DIFMAX = fabs(STPR-STPI)/(0.01*STPI);}
    }
  }

  if(DIFMAX > 1.0E-3 && ISCALE == 1)
  {
    ICOR = 1;
    for(int I = 0; I < NDATA; I++)
    {
      FL[I] = log(F2[I]*mat.RHO*1.0E6);
    }
    SPLINE( EITL, FL, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, NDATA);
    if(penGetError() != 0){
      return;}
  }
  else
  {
    ICOR = 0;
  }

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    EBRaT(mat, grid.ET[I], WCRM, XH0, XH1, XH2, XS1, XS2, grid.DLEMP1, grid.DLFC, constants::WB);
    if(penGetError() != 0){
      return;}
    STPR = (XS1+XH1)*mat.VMOL;
    if(ICOR == 1)
    {
      unsigned int J;
      EC = grid.DLEMP[I];
      FINDI( EITL, EC, NDATA, J);
      STPI = exp(ceel00.A[J]+EC*(ceel00.B[J]+EC*(ceel00.C[J]+EC*ceel00.D[J])));
      FACT = STPI/STPR;
    }
    else
    {
      FACT=1.0;
    }
    mat.RSTPE[I] = STPR*FACT;
    mat.TSTPE[I] = log(mat.CSTPE[I]+mat.RSTPE[I]);
    mat.TSTRE[I] = log(mat.TSTRE[I]+(XS2+XH2)*mat.VMOL*FACT);

    mat.SEHBR[I] = XH0*mat.VMOL*FACT;
    if(mat.SEHBR[I] < 1.0E-35){mat.SEHBR[I] = 1.0E-35;}
    mat.SEHBR[I] = log(mat.SEHBR[I]);
    if(IBREMS == 1)
    {
      mat.W1E[I] = mat.W1E[I]+XS1*mat.VMOL*FACT;
      mat.W2E[I] = mat.W2E[I]+XS2*mat.VMOL*FACT;
    }
  }

  //  ****  Electron range as a function of energy.

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    F1[I]=1.0/(mat.CSTPE[I]+mat.RSTPE[I]);
  }
  SPLINE( grid.ET, F1, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, constants::NEGP);
  if(penGetError() != 0){
    return;}
  mat.RANGE[PEN_ELECTRON][0] = 1.0E-8;
  mat.RANGEL[PEN_ELECTRON][0] = log(mat.RANGE[PEN_ELECTRON][0]);

  double XL;
  double XU;
  double DRaux;
  for(unsigned int I = 1; I < constants::NEGP; I++)
  {
    XL = grid.ET[I-1];
    XU = grid.ET[I];
    SINTEG( grid.ET, ceel00.A, ceel00.B, ceel00.C, ceel00.D, XL, XU, DRaux, constants::NEGP);
    if(penGetError() != 0){
      return;}
    mat.RANGE[PEN_ELECTRON][I] = mat.RANGE[PEN_ELECTRON][I-1]+DRaux;
    mat.RANGEL[PEN_ELECTRON][I] = log(mat.RANGE[PEN_ELECTRON][I]);
  }

  //  ****  Electron radiative yield as a function of energy.

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    F1[I] = mat.RSTPE[I]/(mat.CSTPE[I]+mat.RSTPE[I]);
  }
  RADY[0] = 1.0E-35;
  mat.EBRY[0] = log(RADY[0]/grid.ET[0]);
  for(unsigned int I = 1; I < constants::NEGP; I++)
  {
    XL = grid.ET[I-1];
    XU = grid.ET[I];
    RADY[I] = RADY[I-1]+RMOMX( grid.ET, F1, XL, XU, constants::NEGP, 0);
    if(penGetError() != 0){
      return;}
    mat.EBRY[I] = log(RADY[I]/grid.ET[I]);
  }

  //  ****  Electron bremss. photon number yield as a function of energy.

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    F1[I] = exp(mat.SEHBR[I])/(mat.CSTPE[I]+mat.RSTPE[I]);
  }
  RADN[0] = 0.0;
  for(unsigned int I = 1; I < constants::NEGP; I++)
  {
    XL = grid.ET[I-1];
    XU = grid.ET[I];
    RADN[I] = RADN[I-1]+RMOMX( grid.ET, F1, XL, XU, constants::NEGP, 0);
    if(penGetError() != 0){
      return;}
  }
    
  //  ****  Print electron stopping power tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Stopping powers for electrons\n");
    fprintf(IWR, "\n   Energy        Scol         Srad         range     Rad. Yield   PhotonYield    delta\n    (eV)       (eV/mtu)     (eV/mtu)       (mtu)                    (W>WCR)\n ------------------------------------------------------------------------------------------\n");
    for(unsigned int I = 0; I < constants::NEGP; I++)
    {
      fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", grid.ET[I], mat.CSTPE[I]/mat.RHO, mat.RSTPE[I]/mat.RHO, mat.RANGE[PEN_ELECTRON][I]*mat.RHO, RADY[I]/grid.ET[I], RADN[I], mat.DEL[I]);
    }
  }
    
  //  ****  Bremss x-section for positrons.

  DIFMAX = 0.0;
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] >= grid.EL && EIT[I]  < grid.EU)
    {
      STPI = F4[I]*mat.RHO*1.0E6;  // Radiative stopping power.
      PBRaT(mat, EIT[I], WCRM, XH0, XH1, XH2, XS1, XS2, grid.DLEMP1, grid.DLFC, constants::WB); 
      if(penGetError() != 0){
	return;}
      STPR = (XS1+XH1)*mat.VMOL;
      if(DIFMAX < fabs(STPR-STPI)/(0.01*STPI)) { DIFMAX = fabs(STPR-STPI)/(0.01*STPI);}
    }
  }

  if(DIFMAX > 1.0E-3 && ISCALE == 1)
  {
    ICOR = 1;
    for(int I = 0; I < NDATA; I++)
    {
      FL[I] = log(F4[I]*mat.RHO*1.0E6);
    }
    SPLINE( EITL, FL, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, NDATA);
    if(penGetError() != 0){
      return;}
  }
  else
  {
    ICOR = 0;
  }

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    PBRaT(mat, grid.ET[I], WCRM, XH0, XH1, XH2, XS1, XS2, grid.DLEMP1, grid.DLFC, constants::WB);
    if(penGetError() != 0){
      return;}
    STPR = (XS1+XH1)*mat.VMOL;
    if(ICOR == 1)
    {
      unsigned int J;
      EC=grid.DLEMP[I];
      FINDI(EITL, EC, NDATA, J);
      STPI = exp(ceel00.A[J]+EC*(ceel00.B[J]+EC*(ceel00.C[J]+EC*ceel00.D[J])));
      FACT = STPI/STPR;
    }
    else
    {
      FACT = 1.0;
    }
    mat.RSTPP[I] = STPR*FACT;
    mat.TSTPP[I] = log(mat.CSTPP[I]+mat.RSTPP[I]);
    mat.TSTRP[I] = log(mat.TSTRP[I]+(XS2+XH2)*mat.VMOL*FACT);
    mat.SPHBR[I] = XH0*mat.VMOL*FACT;
    if(mat.SPHBR[I] < 1.0E-35){ mat.SPHBR[I] = 1.0E-35;}
    mat.SPHBR[I] = log(mat.SPHBR[I]);
    if(IBREMS == 1)
    {
      mat.W1P[I] = mat.W1P[I]+XS1*mat.VMOL*FACT;
      mat.W2P[I] = mat.W2P[I]+XS2*mat.VMOL*FACT;
    }
  }

  //  ****  Positron annihilation.

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    double TXAN;
    PANaT( grid.ET[I], TXAN);
    mat.SPAN[I] = log(TXAN*mat.ZT*mat.VMOL);
  }

  //  ****  Positron range as a function of energy.

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    F3[I] = 1.0/(mat.CSTPP[I]+mat.RSTPP[I]);
  }
  SPLINE( grid.ET, F3, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, constants::NEGP);
  if(penGetError() != 0){
    return;}
  mat.RANGE[PEN_POSITRON][0] = 1.0E-8;
  mat.RANGEL[PEN_POSITRON][0] = log(mat.RANGE[PEN_POSITRON][0]);
  for(unsigned int I = 1; I < constants::NEGP; I++)
  {
    XL = grid.ET[I-1];
    XU = grid.ET[I];
    SINTEG( grid.ET, ceel00.A, ceel00.B, ceel00.C, ceel00.D, XL, XU, DRaux, constants::NEGP);
    if(penGetError() != 0){
      return;}
    mat.RANGE[PEN_POSITRON][I] = mat.RANGE[PEN_POSITRON][I-1]+DRaux;
    mat.RANGEL[PEN_POSITRON][I] = log(mat.RANGE[PEN_POSITRON][I]);
  }

  //  ****  Positron radiative yield as a function of energy.

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    F3[I] = mat.RSTPP[I]/(mat.CSTPP[I]+mat.RSTPP[I]);
  }
  RADY[0] = 1.0E-35;
  mat.PBRY[0] = log(RADY[0]/grid.ET[0]);
  for(unsigned int I = 1; I < constants::NEGP; I++)
  {
    XL = grid.ET[I-1];
    XU = grid.ET[I];
    RADY[I] = RADY[I-1]+RMOMX( grid.ET, F3, XL, XU, constants::NEGP, 0);
    if(penGetError() != 0){
      return;}
    mat.PBRY[I] = log(RADY[I]/grid.ET[I]);
  }

  //  ****  Positron bremss. photon number yield as a function of energy.

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    F3[I] = exp(mat.SPHBR[I])/(mat.CSTPP[I]+mat.RSTPP[I]);
  }
  RADN[0]=0.0;
  for(unsigned int I = 1; I < constants::NEGP; I++)
  {
    XL = grid.ET[I-1];
    XU = grid.ET[I];
    RADN[I] = RADN[I-1]+RMOMX( grid.ET, F3, XL, XU, constants::NEGP, 0);
    if(penGetError() != 0){
      return;}
  }

  //  ****  Print positron stopping power tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Stopping powers for positrons\n");
    fprintf(IWR, "\n   Energy        Scol         Srad         range     Rad. Yield   PhotonYield  annih. mfp\n    (eV)       (eV/mtu)     (eV/mtu)       (mtu)                    (W>WCR)      (mtu)\n ------------------------------------------------------------------------------------------\n");
    for(unsigned int I = 0; I < constants::NEGP; I++)
    {
      fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", grid.ET[I], mat.CSTPP[I]/mat.RHO, mat.RSTPP[I]/mat.RHO, mat.RANGE[PEN_POSITRON][I]*mat.RHO, RADY[I]/grid.ET[I], RADN[I], mat.RHO/exp(mat.SPAN[I]));
    }
  }

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Soft stopping power and energy straggling\n");
    fprintf(IWR, "\n   Energy       SStp,e-      SStr,e-     STP,e-     SStp,e+      SStr,e+     Stp,e+\n    (eV)       (eV/mtu)    (eV**2/mtu)  soft/tot   (eV/mtu)    (eV**2/mtu)  soft/tot\n ------------------------------------------------------------------------------------\n");
  }
  double FMTU = 1.0/mat.RHO;
  double FSOFTE;
  double FSOFTP;
  double W1EW;
  double W2EW;
  double W1PW;
  double W2PW;
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
      //  ****  Soft energy-loss events are switched off when W1 is too small.
    FSOFTE = mat.W1E[I]/(mat.CSTPE[I]+mat.RSTPE[I]);
    FSOFTP = mat.W1P[I]/(mat.CSTPP[I]+mat.RSTPP[I]);
    if(mat.W1E[I] > 1.0E-5*(mat.CSTPE[I]+mat.RSTPE[I]))
    {
      W1EW = mat.W1E[I];
      W2EW = mat.W2E[I];
      if(mat.W1E[I] < 1.0E-35){ mat.W1E[I] = 1.0E-35;}
      mat.W1E[I] = log(mat.W1E[I]);
      if(mat.W2E[I] < 1.0E-35){ mat.W2E[I] = 1.0E-35;}
      mat.W2E[I] = log(mat.W2E[I]);
    }
    else
    {
      W1EW = 0.0;
      W2EW = 0.0;
      mat.W1E[I] = -100.0;
      mat.W2E[I] = -100.0;
    }
    if(mat.W1P[I] > 1.0E-5*(mat.CSTPP[I]+mat.RSTPP[I]))
    {
      W1PW = mat.W1P[I];
      W2PW = mat.W2P[I];
      if(mat.W1P[I] < 1.0E-35){ mat.W1P[I] = 1.0E-35;}
      mat.W1P[I] = log(mat.W1P[I]);
      if(mat.W2P[I] < 1.0E-35){ mat.W2P[I] = 1.0E-35;}
      mat.W2P[I] = log(mat.W2P[I]);
    }
    else
    {
      W1PW = 0.0;
      W2PW = 0.0;
      mat.W1P[I] = -100.0;
      mat.W2P[I] = -100.0;
    }
    if(INFO >= 3){ fprintf(IWR, "%12.5E%13.5E%13.5E%10.2E%13.5E%13.5E%10.2E\n", grid.ET[I], W1EW*FMTU, W2EW*FMTU, FSOFTE, W1PW*FMTU, W2PW*FMTU, FSOFTP);}
  }

  //  ****  Elastic scattering of electrons and positrons.
  EELaR(mat, IRD, IWR, INFO, ceel00, ceintf, grid.DLEMP, grid.ET);
  if(penGetError() != PEN_SUCCESS){
    return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    mat.TRL1E[I] = -log(ceel00.XE1[I]*mat.VMOL);
    mat.TRL2E[I] = -log(ceel00.XE2[I]*mat.VMOL);
    mat.TRL1P[I] = -log(ceel00.XP1[I]*mat.VMOL);
    mat.TRL2P[I] = -log(ceel00.XP2[I]*mat.VMOL);
  }
  EELdR(mat, IRD, IWR, INFO, cdcsep, crita, ceintf, grid.ET);  // Uses the ELSEPA database.
  if(penGetError() != 0){
    return;}
      
  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Soft angular transport coefficients\n");
    fprintf(IWR, "\n   Energy      SITMFP1,e-   SITMFP2,e-   SITMFP1,e+   SITMFP2,e+\n    (eV)        (1/mtu)      (1/mtu)      (1/mtu)      (1/mtu)\n ----------------------------------------------------------------\n");
    for(unsigned int I = 0; I < constants::NEGP; I++)
    {
      fprintf(IWR, "%12.5E%13.5E%13.5E%13.5E%13.5E\n", grid.ET[I], exp(mat.T1E[I])*FMTU, exp(mat.T2E[I])*FMTU, exp(mat.T1P[I])*FMTU, exp(mat.T2P[I])*FMTU);
    }
  }
    
  //  ****  Print electron mean free path tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Electron mean free paths (hard events)\n *** NOTE: The MFP for inner-shell ionisation (isi) is listed only for\n           completeness. The MFP for inelastic collisions (in) has been\n           calculated by considering all inelastic events, including isi.\n");
    fprintf(IWR, "\n   Energy        MFPel        MFPin        MFPbr       MFPtot       MFPisi\n    (eV)         (mtu)        (mtu)        (mtu)        (mtu)        (mtu)\n -----------------------------------------------------------------------------\n");
  }

  double FPEL;
  double FPSI;
  double FPIN;
  double FPBR;
  double FPTOT;
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    FPEL = exp(mat.SEHEL[I]);
    FPSI = exp(mat.SEISI[I]);
    FPIN = exp(mat.SEHIN[I])+FPSI;
    FPBR = exp(mat.SEHBR[I]);
    FPTOT = FPEL+FPIN+FPBR;
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", grid.ET[I], mat.RHO/FPEL, mat.RHO/FPIN, mat.RHO/FPBR, mat.RHO/FPTOT, mat.RHO/FPSI);}
    mat.SETOT[I] = log(FPTOT);
  }

  for(unsigned int I = 1; I < constants::NEGP-1; I++)
  {
    if(exp(mat.SETOT[I]) > 1.005*exp(mat.SETOT[I-1]) && exp(mat.SETOT[I]) > 1.005*exp(mat.SETOT[I+1]) && grid.ET[I] > mat.EABS[PEN_ELECTRON] && grid.ET[I] < 1.0E6)
    {
      fprintf(IWR, "\n WARNING: The electron hard IMFP has a maximum at E = %13.6E eV\n", grid.ET[I]);
    }
  }

  //  ****  Print positron mean free path tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Positron mean free paths (hard events)\n *** NOTE: The MFP for inner-shell ionisation (isi) is listed only for\n           completeness. The MFP for inelastic collisions (in) has been\n           calculated by considering all inelastic events, including isi.\n");
    fprintf(IWR, "\n   Energy        MFPel        MFPin        MFPbr        MFPan       MFPtot       MFPisi\n    (eV)         (mtu)        (mtu)        (mtu)        (mtu)        (mtu)        (mtu)\n ------------------------------------------------------------------------------------------\n");
  }
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    FPEL = exp(mat.SPHEL[I]);
    FPSI = exp(mat.SPISI[I]);
    FPIN = exp(mat.SPHIN[I])+FPSI;
    FPBR = exp(mat.SPHBR[I]);
    double FPAN = exp(mat.SPAN[I]);
    FPTOT = FPEL+FPIN+FPBR+FPAN;
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", grid.ET[I], mat.RHO/FPEL, mat.RHO/FPIN, mat.RHO/FPBR, mat.RHO/FPAN, mat.RHO/FPTOT, mat.RHO/FPSI);}
    mat.SPTOT[I] = log(FPTOT);
  }

  for(unsigned int I = 1; I < constants::NEGP-1; I++)
  {
    if(exp(mat.SPTOT[I]) > 1.005*exp(mat.SPTOT[I-1]) && exp(mat.SPTOT[I]) > 1.005*exp(mat.SPTOT[I+1]) && grid.ET[I] > mat.EABS[PEN_POSITRON] && grid.ET[I] < 1.0E6)
    {
      fprintf(IWR, "\n WARNING: The positron hard IMFP has a maximum at E = %13.6E eV\n", grid.ET[I]);
    }
  }

  //  ****  Correction for the energy dependence of the stopping power and
  //  the energy straggling parameter of soft interactions.

  for(unsigned int I = 0; I < constants::NEGP-1; I++)
  {
    mat.DW1EL[I] = (mat.W1E[I+1]-mat.W1E[I])/(grid.ET[I+1]-grid.ET[I]);
    mat.DW2EL[I] = 0.5*(mat.W2E[I+1]-mat.W2E[I])/(grid.ET[I+1]-grid.ET[I])+mat.DW1EL[I];
    mat.DW1EL[I] = 0.5*mat.DW1EL[I];
    mat.DW1PL[I] = (mat.W1P[I+1]-mat.W1P[I])/(grid.ET[I+1]-grid.ET[I]);
    mat.DW2PL[I] = 0.5*(mat.W2P[I+1]-mat.W2P[I])/(grid.ET[I+1]-grid.ET[I])+mat.DW1PL[I];
    mat.DW1PL[I] = 0.5*mat.DW1PL[I];
  }
  mat.DW1EL[constants::NEGP-1] = 0.0;
  mat.DW2EL[constants::NEGP-1] = 0.0;
  mat.DW1PL[constants::NEGP-1] = 0.0;
  mat.DW2PL[constants::NEGP-1] = 0.0;

  //  ************  Photon interaction properties.

  //  ****  Rayleigh scattering.
  GRAaR(grid, mat, crita, IRD, IWR, INFO);
  if(penGetError() != 0){
    return;}
  //  ****  Compton scattering and pair production.
  fscanf(IRD, "%*57c%4d%*[^\n]", &NDATA);
  getc(IRD);
  
  if(INFO >= 2){ fprintf(IWR, "\n *** Compton and pair-production cross sections,  NDATA =%4d\n", NDATA);}
  if((unsigned)NDATA > constants::NEGP){ penError(ERR_PEMATR_2_MANY_DP_2);
    return;}
  if(INFO >= 2){ fprintf(IWR, "\n  Energy     CS-Comp     CS-pair   CS-triplet\n   (eV)      (cm**2)     (cm**2)     (cm**2)\n ----------------------------------------------\n");}
  for(int I = 0; I < NDATA; I++)
  {
    fscanf(IRD, "%lf %lf %lf %lf%*[^\n]", &EIT[I], &F2[I], &F3[I], &F4[I]);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR, "%10.3E%12.5E%12.5E%12.5E\n", EIT[I], F2[I], F3[I], F4[I]);}
    EITL[I] = log(EIT[I]);
  }

  //  ****  Compton scattering.

  double VMOLL = log(mat.VMOL);
  for(int I = 0; I < NDATA; I++)
  {
    FL[I] = F2[I];
    if(FL[I] < 1.0E-35){ FL[I] = 1.0E-35;}
    FL[I] = log(FL[I]);
  }
  SPLINE( EITL, FL, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, NDATA);
  if(penGetError() != 0){
    return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  { 
    unsigned int J;     
    EC = grid.DLEMP[I];
    FINDI( EITL, EC, NDATA, J);
    mat.SGCO[I] = ceel00.A[J]+EC*(ceel00.B[J]+EC*(ceel00.C[J]+EC*ceel00.D[J]))+VMOLL;
  }

  //  ****  Electron-positron pair and triplet production.

  //  ****  Pairs.
  int NP_S = 0;
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] < 1.023E6) { continue;}
    NP_S = NP_S+1;
    F1[NP_S-1] = EITL[I];
    FACT = pow(EIT[I]/(EIT[I]-1.022E6),3);
        
    FL[NP_S-1] = F3[I];
    if(FL[NP_S-1] < 1.0E-35){ FL[NP_S-1] = 1.0E-35;}
    FL[NP_S-1] = log(FACT*FL[NP_S-1]);
  }
  SPLINE( F1, FL, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, NP_S);
  if(penGetError() != 0){
    return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    if(grid.ET[I] < 1.023E6)
    {
      mat.SGPP[I] = -80.6;
    }
    else
    {
      unsigned int J;  
      FACT = pow(grid.ET[I]/(grid.ET[I]-1.022E6),3);
      EC = grid.DLEMP[I];
      FINDI( F1, EC, NP_S, J);
      mat.SGPP[I] = ceel00.A[J]+EC*(ceel00.B[J]+EC*(ceel00.C[J]+EC*ceel00.D[J]))+VMOLL-log(FACT);
    }
  }
  //  ****  Triplets.
  NP_S = 0;
  for(int I = 0; I < NDATA; I++)
  {
    if(EIT[I] < 2.045E6){ continue;}
    NP_S = NP_S+1;
    FACT = pow(EIT[I]/(EIT[I]-2.044E6),3);
    F1[NP_S-1] = EITL[I];

    FL[NP_S-1] = F4[I];
    if(FL[NP_S-1] < 1.0E-35){ FL[NP_S-1] = 1.0E-35;}
    FL[NP_S-1] = log(FACT*FL[NP_S-1]);
  }
  SPLINE( F1, FL, ceel00.A, ceel00.B, ceel00.C, ceel00.D, 0.0, 0.0, NP_S);
  if(penGetError() != 0){
    return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    if(grid.ET[I] < 2.045E6)
    {
      mat.TRIP[I] = -80.6;
    }
    else
    {
      unsigned int J;
      FACT = pow(grid.ET[I]/(grid.ET[I]-2.044E6),3);
      EC = grid.DLEMP[I];
      FINDI( F1, EC, NP_S, J);
      double TRIPL = exp(ceel00.A[J]+EC*(ceel00.B[J]+EC*(ceel00.C[J]+EC*ceel00.D[J]))+VMOLL-log(FACT));
      double PAIR = exp(mat.SGPP[I]);
      mat.SGPP[I] = log(PAIR+TRIPL);
      mat.TRIP[I] = TRIPL/(PAIR+TRIPL);
    }
  }

  if(mat.NOSCCO > 1)
  {
    mat.PTRSH[0] = mat.FCO[0];
    for(int I = 1; I < mat.NOSCCO; I++)
    {
      mat.PTRSH[I] = mat.PTRSH[I-1]+mat.FCO[I];
    }
    for(int I = 0; I < mat.NOSCCO; I++)
    {
      mat.PTRSH[I] = mat.PTRSH[I]/mat.PTRSH[mat.NOSCCO-1];
    }
  }
  else
  {
    mat.PTRSH[0] = 1.0;
  }

  //  ****  Photoelectric absorption.

  GPHaR( mat, elements, IRD, IWR, INFO, cgph01, grid.DLEMP);
  if(penGetError() != 0){
    return;}

  //  ****  Photon 'range' (= mean free path, slightly underestimated).

  double PRAY;
  double PCO;
  double PPP;
  double PPH;
  double PT;
  for(unsigned KE = 0; KE < constants::NEGP; KE++)
  { 
    double ECS;
    GRAaTI(mat, grid.DLEMP1, grid.DLFC, grid.ET[KE], ECS);
    PRAY = ECS*mat.VMOL;
    PCO = exp(mat.SGCO[KE]);
    if(grid.ET[KE] < 1.023E6)
    {
      PPP = 1.0E-35;
    }
    else
    {
      PPP = exp(mat.SGPP[KE]);
    }
    PPH = mat.SGPH[KE];
    PT = (PRAY+PCO+PPP+PPH);
    mat.RANGE[PEN_PHOTON][KE] = 1.0/PT;
    mat.RANGEL[PEN_PHOTON][KE] = log(mat.RANGE[PEN_PHOTON][KE]);
  }

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Photon mass attenuation coefficients\n");
    fprintf(IWR, "\n   Energy      Rayleigh      Compton    Photoelect.     Pair         Total\n    (eV)        (1/mtu)      (1/mtu)      (1/mtu)      (1/mtu)      (1/mtu)\n -----------------------------------------------------------------------------\n");
  }
  for(int I = 0; I < cgph01.NPHD; I++)
  {
    if(cgph01.ER[I] >= grid.EL && cgph01.ER[I] <= grid.EU)
    {
      double ECS;
      double XEL = log(cgph01.ER[I]);
      double XE = 1.0+(XEL-grid.DLEMP1)*grid.DLFC;
      int KE = (int)XE;
      double XEK = XE-(double)KE;
      GRAaTI(mat, grid.DLEMP1, grid.DLFC, cgph01.ER[I], ECS);
      PRAY = ECS*mat.VMOL/mat.RHO;
      PCO = exp(mat.SGCO[KE-1]+(mat.SGCO[KE]-mat.SGCO[KE-1])*XEK)/mat.RHO;
      if(cgph01.ER[I] < 1.022E6)
      {
        PPP = 1.0E-35;
      }
      else
      {
        PPP = exp(mat.SGPP[KE-1]+(mat.SGPP[KE]-mat.SGPP[KE-1])*XEK)/mat.RHO;
      }
      PPH = cgph01.XSR[I]*mat.VMOL/mat.RHO;
      PT = PRAY+PCO+PPP+PPH;
      if(INFO >= 3)
      {
        fprintf(IWR, "%12.5E%13.5E%13.5E%13.5E%13.5E%13.5E\n", cgph01.ER[I], PRAY, PCO, PPH, PPP, PT);
      }
    }
  }

  //  ****  Pair production. Initialisation routine.

  GPPa0(mat, elements);

  //Precalculate distance between the logarithm of inverse mean
  //free path points.
  const unsigned NEGPm1 = constants::NEGP-1;

  //Electron interactions
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSEHEL[I] = mat.SEHEL[I+1] - mat.SEHEL[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSEHIN[I] = mat.SEHIN[I+1] - mat.SEHIN[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSEISI[I] = mat.SEISI[I+1] - mat.SEISI[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSEHBR[I] = mat.SEHBR[I+1] - mat.SEHBR[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSETOT[I] = mat.SETOT[I+1] - mat.SETOT[I];
  }

  //Gamma interactions
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSGRA[I] = mat.SGRA[I+1] - mat.SGRA[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSGCO[I] = mat.SGCO[I+1] - mat.SGCO[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSGPH[I] = mat.SGPH[I+1] - mat.SGPH[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSGPP[I] = mat.SGPP[I+1] - mat.SGPP[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DTRIP[I] = mat.TRIP[I+1] - mat.TRIP[I];
  }

  //Positron interactions
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSPHEL[I] = mat.SPHEL[I+1] - mat.SPHEL[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSPHIN[I] = mat.SPHIN[I+1] - mat.SPHIN[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSPISI[I] = mat.SPISI[I+1] - mat.SPISI[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSPHBR[I] = mat.SPHBR[I+1] - mat.SPHBR[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSPAN[I]  = mat.SPAN[I+1]  - mat.SPAN[I];
  }
  for(unsigned int I = 0; I < NEGPm1; ++I){
    mat.DSPTOT[I] = mat.SPTOT[I+1] - mat.SPTOT[I];
  }
  
  
  strcpy(LNAME," PENELOPE (v. 2018)  End of material data file ........");
  fscanf(IRD,"%55[^\n]%*[^\n]",NAME);
  getc(IRD);
  if(strcmp(NAME,LNAME) != 0)
  {
    fprintf(IWR, "\n I/O error. Corrupt material data file.\n");
    fprintf(IWR, "       The last line is: [%s]\n", NAME);
    fprintf(IWR, "      ... and should be: [%s]\n", LNAME);
    penError(ERR_PEMATR_CORRUPT_MAT_FILE);
    return;
  }
  fprintf(IWR, "%55s\n", NAME);

}

//-------------------
// Element database
//-------------------

pen_elementDataBase::pen_elementDataBase()
{
  
  //  ************  Chemical symbols of the elements.
  const char LASYMB_1[nElements][3] = {"H ","He","Li","Be","B ","C ",
				       "N ","O ","F ","Ne","Na","Mg",
				       "Al","Si","P ","S ","Cl","Ar",
				       "K ","Ca","Sc","Ti","V ","Cr",
				       "Mn","Fe","Co","Ni","Cu","Zn",
				       "Ga","Ge","As","Se","Br","Kr",
				       "Rb","Sr","Y ","Zr","Nb","Mo",
				       "Tc","Ru","Rh","Pd","Ag","Cd",
				       "In","Sn","Sb","Te","I ","Xe",
				       "Cs","Ba","La","Ce","Pr","Nd",
				       "Pm","Sm","Eu","Gd","Tb","Dy",
				       "Ho","Er","Tm","Yb","Lu","Hf",
				       "Ta","W ","Re","Os","Ir","Pt",
				       "Au","Hg","Tl","Pb","Bi","Po",
				       "At","Rn","Fr","Ra","Ac","Th",
				       "Pa","U ","Np","Pu","Am","Cm",
				       "Bk","Cf","Es"};

  for(unsigned int i = 0; i < nElements; i++)
    strcpy(LASYMB[i],LASYMB_1[i]);

  //  ************  Atomic weights (mean relative atomic masses).
  
  const double ATW_1[nElements] = {1.0079, 4.0026, 6.9410, 9.0122, 1.0811E1,
				   1.2011E1, 1.4007E1, 1.5999E1, 1.8998E1,
				   2.0179E1, 2.2990E1, 2.4305E1, 2.6982E1,
				   2.8086E1, 3.0974E1, 3.2066E1, 3.5453E1,
				   3.9948E1, 3.9098E1, 4.0078E1, 4.4956E1,
				   4.7880E1, 5.0942E1, 5.1996E1, 5.4938E1,
				   5.5847E1, 5.8933E1, 5.8690E1, 6.3546E1,
				   6.5390E1, 6.9723E1, 7.2610E1, 7.4922E1,
				   7.8960E1, 7.9904E1, 8.3800E1, 8.5468E1,
				   8.7620E1, 8.8906E1, 9.1224E1, 9.2906E1,
				   9.5940E1, 9.7907E1, 1.0107E2, 1.0291E2,
				   1.0642E2, 1.0787E2, 1.1241E2, 1.1482E2,
				   1.1871E2, 1.2175E2, 1.2760E2, 1.2690E2,
				   1.3129E2, 1.3291E2, 1.3733E2, 1.3891E2,
				   1.4012E2, 1.4091E2, 1.4424E2, 1.4491E2,
				   1.5036E2, 1.5196E2, 1.5725E2, 1.5893E2,
				   1.6250E2, 1.6493E2, 1.6726E2, 1.6893E2,
				   1.7304E2, 1.7497E2, 1.7849E2, 1.8095E2,
				   1.8385E2, 1.8621E2, 1.9020E2, 1.9222E2,
				   1.9508E2, 1.9697E2, 2.0059E2, 2.0438E2,
				   2.0720E2, 2.0898E2, 2.0898E2, 2.0999E2,
				   2.2202E2, 2.2302E2, 2.2603E2, 2.2703E2,
				   2.3204E2, 2.3104E2, 2.3803E2, 2.3705E2,
				   2.3905E2, 2.4306E2, 2.4707E2, 2.4707E2,
				   2.5108E2, 2.5208E2};
  
  memcpy(ATW,ATW_1,sizeof(double)*nElements);

  //  ************  Mean excitation energies of the elements (eV).
  
  const double EPX_1[nElements] ={19.2, 41.8, 40.0, 63.7, 76.0, 81.0, 82.0, 95.0,
				  115.0, 137.0, 149.0, 156.0, 166.0, 173.0, 173.0,
				  180.0, 174.0, 188.0, 190.0, 191.0, 216.0, 233.0,
				  245.0, 257.0, 272.0, 286.0, 297.0, 311.0, 322.0,
				  330.0, 334.0, 350.0, 347.0, 348.0, 343.0, 352.0,
				  363.0, 366.0, 379.0, 393.0, 417.0, 424.0, 428.0,
				  441.0, 449.0, 470.0, 470.0, 469.0, 488.0, 488.0,
				  487.0, 485.0, 491.0, 482.0, 488.0, 491.0, 501.0,
				  523.0, 535.0, 546.0, 560.0, 574.0, 580.0, 591.0,
				  614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 694.0,
				  705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0,
				  790.0, 800.0, 810.0, 823.0, 823.0, 830.0, 825.0,
				  794.0, 827.0, 826.0, 841.0, 847.0, 878.0, 890.0,
				  902.0, 921.0, 934.0, 939.0, 952.0, 966.0, 980.0};

  memcpy(EPX,EPX_1,sizeof(double)*nElements);


  //  ************  Pair-production cross section parameters.
  
  //  ****  Screening parameter (R mc/hbar).
  const double RSCR_1[nElements] = {1.2281E2, 7.3167E1, 6.9228E1, 6.7301E1, 6.4696E1,
				    6.1228E1, 5.7524E1, 5.4033E1, 5.0787E1, 4.7851E1,
				    4.6373E1, 4.5401E1, 4.4503E1, 4.3815E1, 4.3074E1,
				    4.2321E1, 4.1586E1, 4.0953E1, 4.0524E1, 4.0256E1,
				    3.9756E1, 3.9144E1, 3.8462E1, 3.7778E1, 3.7174E1,
				    3.6663E1, 3.5986E1, 3.5317E1, 3.4688E1, 3.4197E1,
				    3.3786E1, 3.3422E1, 3.3068E1, 3.2740E1, 3.2438E1,
				    3.2143E1, 3.1884E1, 3.1622E1, 3.1438E1, 3.1142E1,
				    3.0950E1, 3.0758E1, 3.0561E1, 3.0285E1, 3.0097E1,
				    2.9832E1, 2.9581E1, 2.9411E1, 2.9247E1, 2.9085E1,
				    2.8930E1, 2.8721E1, 2.8580E1, 2.8442E1, 2.8312E1,
				    2.8139E1, 2.7973E1, 2.7819E1, 2.7675E1, 2.7496E1,
				    2.7285E1, 2.7093E1, 2.6911E1, 2.6705E1, 2.6516E1,
				    2.6304E1, 2.6108E1, 2.5929E1, 2.5730E1, 2.5577E1,
				    2.5403E1, 2.5245E1, 2.5100E1, 2.4941E1, 2.4790E1,
				    2.4655E1, 2.4506E1, 2.4391E1, 2.4262E1, 2.4145E1,
				    2.4039E1, 2.3922E1, 2.3813E1, 2.3712E1, 2.3621E1,
				    2.3523E1, 2.3430E1, 2.3331E1, 2.3238E1, 2.3139E1,
				    2.3048E1, 2.2967E1, 2.2833E1, 2.2694E1, 2.2624E1,
				    2.2545E1, 2.2446E1, 2.2358E1, 2.2264E1};
  
  memcpy(RSCR,RSCR_1,sizeof(double)*nElements);

  //  ****  Asymptotic triplet contribution (eta).
  const double ETA_1[nElements] = {1.1570, 1.1690, 1.2190, 1.2010, 1.1890, 1.1740, 1.1760,
				 1.1690, 1.1630, 1.1570, 1.1740, 1.1830, 1.1860, 1.1840,
				 1.1800, 1.1780, 1.1750, 1.1700, 1.1800, 1.1870, 1.1840,
				 1.1800, 1.1770, 1.1660, 1.1690, 1.1660, 1.1640, 1.1620,
				 1.1540, 1.1560, 1.1570, 1.1580, 1.1570, 1.1580, 1.1580,
				 1.1580, 1.1660, 1.1730, 1.1740, 1.1750, 1.1700, 1.1690,
				 1.1720, 1.1690, 1.1680, 1.1640, 1.1670, 1.1700, 1.1720,
				 1.1740, 1.1750, 1.1780, 1.1790, 1.1800, 1.1870, 1.1940,
				 1.1970, 1.1960, 1.1940, 1.1940, 1.1940, 1.1940, 1.1940,
				 1.1960, 1.1970, 1.1960, 1.1970, 1.1970, 1.1980, 1.1980,
				 1.2000, 1.2010, 1.2020, 1.2040, 1.2050, 1.2060, 1.2080,
				 1.2070, 1.2080, 1.2120, 1.2150, 1.2180, 1.2210, 1.2240,
				 1.2270, 1.2300, 1.2370, 1.2430, 1.2470, 1.2500, 1.2510,
				 1.2520, 1.2550, 1.2560, 1.2570, 1.2590, 1.2620, 1.2620,
				 1.2650};

  memcpy(ETA,ETA_1,sizeof(double)*nElements);

  GPHa0();
  RELAX0();
}

//  *********************************************************************
//                       SUBROUTINE GPHa0
//  *********************************************************************
void pen_elementDataBase::GPHa0()
{
  //  This subroutine sets all variables in common /CGPH00/ to zero.
  //  It has to be invoked before reading the first material definition
  //  file.

  
  for(unsigned I = 0; I < nElements; I++)
  {
    NPHS[I] = 0;
    IPHF[I] = 0;
    IPHL[I] = 0;
  }

  for(unsigned I = 0; I < constants::NTP; I++)
  {
    EPH[I] = 0.0;
    for(int J = 0; J < 17; J++)
    {
      XPH[I][J] = 1.0E-35;
    }
  }
  NCUR = 0;
}

//  *********************************************************************
//                       SUBROUTINE RELAX0
//  *********************************************************************
void pen_elementDataBase::RELAX0()
{
  //  This subroutine sets all variables in common /CRELAX/ to zero.

  
  for(unsigned I = 0; I < nElements; I++)
  {
    NSHT[I] = 0;
    for(int J = 0; J < 16; J++)
    {
      IFIRST[I][J] = 0;
      ILAST[I][J] = 0;
    }
    for(int J = 0; J < 30; J++)
    {
      IFI[I][J] = 0;
      EB[I][J] = 0.0;
      IKS[I][J] = 0;
    }
  }

  for(unsigned I = 0; I < constants::NRX; I++)
  {
    IS0[I] = 0;
    IS1[I] = 0;
    IS2[I] = 0;
    ET[I] = 0.0;
    P[I] = 0.0;
    F[I] = 0.0;
  }
  NRELAX = 0;
}


unsigned int pen_elementDataBase::getPHposition(const unsigned int IZ, const double XEL)
{
  //Search first position in photoelectric cross sections arrays
  //for element with atomic number Z = IZ.

  //Perform a binary search
  unsigned int I = IPHF[IZ-1];
  unsigned int IU = IPHL[IZ-1];

  unsigned int IT;
  while(IU - I > 1)
    {
      IT = (I+IU)/2;
      if(XEL > EPH[IT])
	I = IT;
      else
	IU = IT;
    };

  return I;
  
}

//-----------------------------------------------
// PENELOPE initialization functions
//-----------------------------------------------

//  *********************************************************************
//                       SUBROUTINE RELAXR
//  *********************************************************************
void RELAXR(pen_elementDataBase& elements, FILE* IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads atomic relaxation data, for a single element,
  //  from the material definition file (unit IRD) and initialises the
  //  algorithm for random sampling of de-excitation cascades of this
  //  element. (See heading comments in subroutine RELAX).

  
  char CH5[6];

  const int NTRAN = 2500;
  
  double PR[2500], ER[2500], WW[2500], FF[2500];
  int JS0[2500], JS1[2500], JS2[2500], KK[2500], IORD[2500], ISR[2500], IQQ[30];
  double EE[30],ALWR[30],CP0P[30];

  const char LSHELL[31][3] = {"  ","K ","L1","L2","L3","M1","M2","M3","M4","M5","N1","N2","N3","N4","N5","N6","N7","O1","O2","O3","O4","O5","O6","O7","P1","P2","P3","P4","P5","Q1","X "};

  //  ****  Input transition data.

  int IZ, NSHR, NT;
  fscanf(IRD, "%*16c%3d%*18c%3d%*23c%5d%*[^\n]", &IZ, &NSHR, &NT);
  getc(IRD);
  if(INFO >= 2){ fprintf(IWR, "\n *** RELAX:  Z =%3d,  no. of shells =%3d,  no. of transitions =%5d", IZ, NSHR, NT);}

  if(NT > NTRAN){ penError(ERR_RELAXR_NTRAN); return;}
  if(elements.NRELAX+NT > int(constants::NRX))
  {
    fprintf(IWR, "Insufficient memory storage in RELAXR.\n");
    fprintf(IWR, "Increase the value of the parameter NRX to %d\n", elements.NRELAX+NT);
    penError(ERR_RELAXR_MEM); return;
  }
  
  if(INFO >= 2){ fprintf(IWR, "\n\n  i   Shell    f    Ui (eV)    Gamma(1/eV)  lifetime (s)    Ji(0)\n --------------------------------------------------------------------\n");}

  for(int IS = 0; IS < NSHR; IS++)
  {
    fscanf(IRD, "%d %5s %d %lf %lf %lf%*[^\n]", &ISR[IS], CH5, &IQQ[IS], &EE[IS], &ALWR[IS], &CP0P[IS]);
    getc(IRD);
    double ALTIME;
    if(ALWR[IS] > 0.0)
    {
      ALTIME = constants::HBAR/ALWR[IS];
    }
    else
    {
      ALTIME = 0.0;
    }
    if(INFO >= 2){ fprintf(IWR, " %3d %5s %2s  %1d %12.5E %12.5E %12.5E %12.5E\n", ISR[IS], CH5, LSHELL[ISR[IS]], IQQ[IS], EE[IS], ALWR[IS], ALTIME, CP0P[IS]);}
  }

  if(NT > 0)
  {
    if(INFO >= 2){ fprintf(IWR, "\n  S0 S1 S2   Probability     Energy (eV)\n ----------------------------------------\n");}
    for(int I = 0; I < NT; I++)
    {
      fscanf(IRD, "%d %d %d %lf %lf%*[^\n]", &JS0[I], &JS1[I], &JS2[I], &PR[I], &ER[I]);
      getc(IRD);
      if(INFO >= 2){ fprintf(IWR, "  %2s %2s %2s %15.8E %15.8E\n", LSHELL[JS0[I]], LSHELL[JS1[I]], LSHELL[JS2[I]], PR[I], ER[I]);}
      if(PR[I] < 1.0E-35)
      {
        if(INFO < 2){ fprintf(IWR, "  %2s %2s %2s %15.8E %15.8E\n", LSHELL[JS0[I]],LSHELL[JS1[I]], LSHELL[JS2[I]], PR[I], ER[I]);}
        penError(ERR_RELAXR_NEG_TRANS); return;
      }
    }
  }

  //  ****  Check if this element's data have already been loaded.

  if(elements.IFIRST[IZ-1][0] != 0){ return;}

  elements.NSHT[IZ-1] = NSHR;
  for(int IS = 0; IS < NSHR; IS++)
  {
    elements.IKS[IZ-1][IS] = ISR[IS];
    elements.IFI[IZ-1][ISR[IS]-1] = IQQ[IS];
    elements.EB[IZ-1][ISR[IS]-1] = EE[IS];
    if(ALWR[IS] > 0.0)
    {
      elements.ALW[IZ-1][ISR[IS]-1] = constants::HBAR/ALWR[IS];
    }
    else
    {
      elements.ALW[IZ-1][ISR[IS]-1] = 0.0;
    }
    elements.CP0[IZ][ISR[IS]-1] = CP0P[IS];
  }
  if(NT == 0)
  {
    for(int IS = 0; IS < 16; IS++)
    {
      elements.IFIRST[IZ-1][IS] = elements.NRELAX+1;
      elements.ILAST[IZ-1][IS] = elements.NRELAX+1;
    //  ****  The array IS0 contains the alias values.
      elements.IS0[elements.NRELAX] = elements.NRELAX+1;
      elements.P[elements.NRELAX] = 1.0;
      elements.F[elements.NRELAX] = 1.0;
      elements.ET[elements.NRELAX] = 0.0;
      elements.IS1[elements.NRELAX] = 1;
      elements.IS2[elements.NRELAX] = 1;
      elements.NRELAX = elements.NRELAX+1;
    }
    return;
  }

  //  ****  Walker's aliasing.

  for(int IS = 0; IS < 16; IS++)
  {
    int N = 0;
    for(int J = 0; J < NT; J++)
    {
      if(JS0[J] == IS+1)
      {
        N = N+1;
        IORD[N-1] = J+1;
        WW[N-1] = PR[J];
      }
    }
    if(N > 1)
    {
      IRND0(WW,FF,KK,N);
      elements.IFIRST[IZ-1][IS] = elements.NRELAX+1;
      elements.ILAST[IZ-1][IS] = elements.NRELAX+N;
      for(int L = 0; L < N; L++)
      {
        elements.P[elements.NRELAX+L] = WW[L];
        elements.F[elements.NRELAX+L] = FF[L];
        elements.ET[elements.NRELAX+L] = ER[IORD[L]-1];
        //  ****  The array IS0 contains the alias values.
        elements.IS0[elements.NRELAX+L] = elements.IFIRST[IZ-1][IS]+KK[L]-1;
        elements.IS1[elements.NRELAX+L] = JS1[IORD[L]-1];
        elements.IS2[elements.NRELAX+L] = JS2[IORD[L]-1];
      }
      elements.NRELAX = elements.NRELAX+N;
    }
    else
    {
      elements.NRELAX = elements.NRELAX+1;
      elements.IFIRST[IZ-1][IS] = elements.NRELAX;
      elements.ILAST[IZ-1][IS] = elements.NRELAX;
      elements.IS0[elements.NRELAX] = elements.NRELAX;
      elements.P[elements.NRELAX-1] = 1.0;
      elements.F[elements.NRELAX-1] = 1.0;
      elements.ET[elements.NRELAX-1] = ER[0];
      elements.IS1[elements.NRELAX-1] = JS1[0];
      elements.IS2[elements.NRELAX-1] = JS2[0];
    }
  }

  //  ****  Verify that transition probabilities are correctly reproduced.

  double TST = 0.0;
  double PT, PPI;
  int IO, IN;
  for(int IS = 0; IS < 16; IS++)
  {
    IO = elements.IFIRST[IZ-1][IS];
    IN = elements.ILAST[IZ-1][IS];
    PT = 0.0;
    for(int I = IO-1; I < IN; I++)
    {
      PT = PT+elements.P[I];
    }
    for(int I = IO-1; I < IN; I++)
    {
      PPI = 0.0;
      for(int J = IO-1; J < IN; J++)
      {
        if(elements.IS0[J] == I+1){ PPI = PPI+(1.0-elements.F[J]);}
      }
      PPI = (PPI+elements.F[I])/double(IN-IO+1);
      if(TST < fabs(1.0-PPI*PT/elements.P[I])){ TST = fabs(1.0-PPI*PT/elements.P[I]);}
    }
  }
  if(TST > 1.0E-12)
  {
    printf("\nTST =%13.6E",TST);
    penError(ERR_RELAXR_ROUND); return;
  }
}

//  *********************************************************************
//                       SUBROUTINE ESIaR
//  *********************************************************************
void ESIaR(pen_material& mat, pen_elementDataBase& elemDB, CESI0& esi0, FILE* IRD, FILE* IWR, int INFO, const double* DLEMP)
{
  //  This subroutine reads cross sections for inner-shell ionisation by
  //  electron impact of the elements in material M and prepares simulation
  //  tables.


  char CS5[16][6];

  const int NDIN=850;
  double E[NDIN], XESIR[NDIN][16], X[NDIN], Y[NDIN];

  strcpy(CS5[0], "CS-K ");
  strcpy(CS5[1], "CS-L1");
  strcpy(CS5[2], "CS-L2");
  strcpy(CS5[3], "CS-L3");
  strcpy(CS5[4], "CS-M1");
  strcpy(CS5[5], "CS-M2");
  strcpy(CS5[6], "CS-M3");
  strcpy(CS5[7], "CS-M4");
  strcpy(CS5[8], "CS-M5");
  strcpy(CS5[9], "CS-N1");
  strcpy(CS5[10], "CS-N2");
  strcpy(CS5[11], "CS-N3");
  strcpy(CS5[12], "CS-N4");
  strcpy(CS5[13], "CS-N5");
  strcpy(CS5[14], "CS-N6");
  strcpy(CS5[15], "CS-N7");

  //  ************  Read element x-section tables

  int NSHR, NDATA, IZZ;
  for(int IEL = 0; IEL < mat.NELEM; IEL++)
  {
    fscanf(IRD, "%*46c%3d%*11c%3d%*10c%4d%*[^\n]", &IZZ, &NSHR, &NDATA);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR, "\n *** Electron impact ionisation cross sections,  IZ =%3d,  NSHELL =%3d,  NDATA =%4d", IZZ, NSHR, NDATA);}

    if(IZZ != mat.IZ[IEL]){ penError(ERR_ESIaR_MATERIAL_DF); return;}
    if(NDATA > NDIN){ penError(ERR_ESIaR_DP); return;}
    if(NSHR > 16){ penError(ERR_ESIaR_SHELLS); return;}
    for(int IE = 0; IE < NDATA; IE++)
    {
      fscanf(IRD, "%lf", &E[IE]);
      for(int IS = 0; IS < NSHR; IS++)
      {
        if(IS<NSHR-1){fscanf(IRD, "%lf", &XESIR[IE][IS]);}
        else{
          fscanf(IRD, "%lf%*[^\n]",&XESIR[IE][IS]);
          getc(IRD);
        }
      }   
    }

      //  ****  Remove shells with ionisation energies less than 50 eV.

    int NSHA;
    if(NSHR > 1)
    {
      bool Eixir = false;
      NSHA = NSHR;
      for(int IS = NSHA-1; IS >= 0; IS--)
      {
        if(elemDB.EB[IZZ-1][IS] < 50.0)
        {
          NSHR = NSHR-1;
        }
        else
        {
          Eixir = true;
          break;
        }
      }
      if(NSHR < 1 && !Eixir){ NSHR = 1;}
    }

    double TCS;
    if(INFO >= 2)
    {
      fprintf(IWR, "\n\n   Energy");
      for(int IS = 0; IS < NSHR; IS++)
      {
        fprintf(IWR, "       %s", CS5[IS]);
      }
      fprintf(IWR,"\n");
      for(int IE = 0; IE < NDATA; IE++)
      {
        TCS = 0.0;
        for(int IS = 0; IS < NSHR; IS++)
        {
          TCS = TCS+XESIR[IE][IS];
        }
        
        fprintf(IWR, " %11.5E", E[IE]);
        for(int IS = 0; IS < NSHR; IS++)
        {
          fprintf(IWR, " %11.5E", XESIR[IE][IS]);
        }
        fprintf(IWR, " %11.5E\n", TCS);
      }
    }

    double XC, DX;
    int IC;
    esi0.NSESI[IZZ-1] = NSHR;
    if(esi0.IESIF[IZZ-1] == 0)
    {
      esi0.IESIF[IZZ-1] = esi0.NCURE+1;
      if(esi0.NCURE+constants::NEGP > esi0.NRP)
      {
        fprintf(IWR, "\nInsufficient memory storage in ESIaR. \nIncrease the value of the parameter NRP to %d\n", esi0.NCURE+constants::NEGP);
        penError(ERR_ESIaR_MEMORY); return;
      }
      for(int IS = 0; IS < NSHR; IS++)
      {
        int N = 0;
        for(int I = 0; I < NDATA; I++)
        {
          if(XESIR[I][IS] > 1.0E-35)
          {
            N = N+1;
            X[N-1] = log(E[I]);
            if(N > 1)
            {
              if(X[N-1] < X[N-2]+1.0E-6){ X[N-1] = X[N-2]+1.0E-6;}
            }
            Y[N-1] = log(XESIR[I][IS]);
          }
        }
        if(N > 4)
        {
          for(unsigned int I = 0; I < constants::NEGP; I++)
          {
            IC = esi0.NCURE+I;
            XC = DLEMP[I];
            if(XC > X[0])
            {
              unsigned int J;
              FINDI(X,XC,N,J);
              if(J == unsigned(N-1)){ J = N-2;}
              DX = X[J+1]-X[J];
              if(DX > 1.0E-6)
              {
                esi0.XESI[IC][IS] = Y[J]+(XC-X[J])*(Y[J+1]-Y[J])/DX;
              }
              else
              {
                esi0.XESI[IC][IS] = (Y[J+1]+Y[J])/2.0;
              }
            }
            else
            {
              esi0.XESI[IC][IS] = -80.6;
            }
          }
        }
        else
        {
          for(unsigned int I = 0; I < constants::NEGP; I++)
          {
            IC = esi0.NCURE+I;
            esi0.XESI[IC][IS] = -80.6;
          }
        }
      }
      esi0.NCURE += constants::NEGP;
      esi0.IESIL[IZZ-1] = esi0.NCURE;
    }
  }
     
}

//  *********************************************************************
//                       SUBROUTINE PSIaR
//  *********************************************************************
void PSIaR(pen_elementDataBase& elements, pen_material& mat, pen_logGrid& grid, CPSI0& cpsi0, FILE *IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads cross sections for inner-shell ionisation by
  //  positron impact of the elements in material M and prepares simulation
  //  tables.


  char CS5[16][6];
  const int NDIN=800;
  double E[NDIN], XPSIR[NDIN][16], X[NDIN], Y[NDIN];

  strcpy(CS5[0], "CS-K ");
  strcpy(CS5[1], "CS-L1");
  strcpy(CS5[2], "CS-L2");
  strcpy(CS5[3], "CS-L3");
  strcpy(CS5[4], "CS-M1");
  strcpy(CS5[5], "CS-M2");
  strcpy(CS5[6], "CS-M3");
  strcpy(CS5[7], "CS-M4");
  strcpy(CS5[8], "CS-M5");
  strcpy(CS5[9], "CS-N1");
  strcpy(CS5[10], "CS-N2");
  strcpy(CS5[11], "CS-N3");
  strcpy(CS5[12], "CS-N4");
  strcpy(CS5[13], "CS-N5");
  strcpy(CS5[14], "CS-N6");
  strcpy(CS5[15], "CS-N7");

  //  ************  Read element x-section tables

  int NSHR, IZZ, NDATA;
  for(int IEL = 0; IEL < mat.NELEM; IEL++)
  {
    fscanf(IRD, "%*46c%3d%*11c%3d%*10c%4d%*[^\n]", &IZZ, &NSHR, &NDATA);
    getc(IRD);

    if(INFO >= 2){ fprintf(IWR, "\n *** Positron impact ionisation cross sections,  IZ =%3d,  NSHELL =%3d,  NDATA =%4d\n", IZZ, NSHR, NDATA);}

    if(IZZ != mat.IZ[IEL]){ penError(ERR_PSIaR_CORRUPTED_FILE); return;}
    if(NDATA > NDIN){ penError(ERR_PSIaR_DP); return;}
    if(NSHR > 16){ penError(ERR_PSIaR_SHELLS); return;}
    for(int IE = 0; IE < NDATA; IE++)
    {
      fscanf(IRD, "%lf", &E[IE]);
      for(int IS = 0; IS < NSHR; IS++)
      {
        if(IS<NSHR-1){fscanf(IRD, "%lf", &XPSIR[IE][IS]);}else{
        fscanf(IRD, "%lf%*[^\n]",&XPSIR[IE][IS]);
        getc(IRD);
      }
      }
    }

      //  ****  Remove shells with ionisation energies less than 50 eV.

    int NSHA;
    if(NSHR > 1)
    {
      NSHA = NSHR;
      bool Eixir = false;
      for(int IS = NSHA-1; IS >= 0; IS--)
      {
        if(elements.EB[IZZ-1][IS] < 50.0)
        {
          NSHR = NSHR-1;
        }
        else
        {
          Eixir = true;
          break;
        }
      }
      if(NSHR < 1 && !Eixir){ NSHR = 1;}
    }

    double TCS;
    if(INFO >= 2)
    {
      fprintf(IWR, "\n   Energy");
      for(int IS = 0; IS < NSHR; IS++)
      {
        fprintf(IWR, "       %s", CS5[IS]);
      }
      fprintf(IWR,"\n");
      for(int IE = 0; IE < NDATA; IE++)
      { 
        TCS = 0.0;
        for(int IS = 0; IS < NSHR; IS++)
        {
          TCS = TCS+XPSIR[IE][IS];
        }
        fprintf(IWR, " %11.5E", E[IE]);
        for(int IS = 0; IS < NSHR; IS++)
        {
          fprintf(IWR, " %11.5E", XPSIR[IE][IS]);
        }
        fprintf(IWR, " %11.5E\n", TCS);
      }
    }

    double XC, DX;
    int N, IC;
    cpsi0.NSPSI[IZZ-1] = NSHR;
    if(cpsi0.IPSIF[IZZ-1] == 0)
    {
      cpsi0.IPSIF[IZZ-1] = cpsi0.NCURP+1;
      if(cpsi0.NCURP+constants::NEGP > CESI0::NRP)
      {
        fprintf(IWR, "\nInsufficient memory storage in PSIaR.");
        fprintf(IWR,"\nIncrease the value of the parameter NRP to %d\n", cpsi0.NCURP+constants::NEGP);
        penError(ERR_PSIaR_MEM); return;
      }
      for(int IS = 0; IS < NSHR; IS++)
      {
        N = 0;
        for(int I = 0; I < NDATA; I++)
        {
          if(XPSIR[I][IS] > 1.0E-35)
          {
            N = N+1;
            X[N-1] = log(E[I]);
            if(N > 1)
            {
              if(X[N-1] < X[N-1-1]+1.0E-6){ X[N-1] = X[N-1-1]+1.0E-6;}
            }
            Y[N-1] = log(XPSIR[I][IS]);
          }
        }
        if(N > 4)
        {
          for(unsigned int I = 0; I < constants::NEGP; I++)
          {
            IC = cpsi0.NCURP+I+1;
            XC = grid.DLEMP[I];
            if(XC > X[0])
            {
              unsigned int J;
              FINDI(X,XC,N,J);
              if(int(J) == N-1){ J = N-2;}
              DX = X[J+1]-X[J];
              if(DX > 1.0E-6)
              {
                cpsi0.XPSI[IC-1][IS] = Y[J]+(XC-X[J])*(Y[J+1]-Y[J])/DX;
              }
              else
              {
                cpsi0.XPSI[IC-1][IS] = (Y[J+1]+Y[J])/2.0;
              }
            }
            else
            {
              cpsi0.XPSI[IC-1][IS] = -80.6;
            }
          }
        }
        else
        {
          for(unsigned int I = 0; I < constants::NEGP; I++)
          {
            IC = cpsi0.NCURP+I+1;
            cpsi0.XPSI[IC-1][IS] = -80.6;
          }
        }
      }
      cpsi0.NCURP = cpsi0.NCURP+constants::NEGP;
      cpsi0.IPSIL[IZZ-1] = cpsi0.NCURP;
    }
  }
}

//  *********************************************************************
//                       SUBROUTINE EBRaR
//  *********************************************************************
void EBRaR(pen_material& mat, CEBR01& cebr01, double &WCRM, FILE* IRD, FILE* IWR, int INFO, const double ET[constants::NEGP], const double DLEMP[constants::NEGP])
{
  //  This subroutine reads the bremss scaled cross section for electrons
  //  in material M from the material data file. It computes restricted
  //  integrated cross sections and initialises the algorithm for simula-
  //  tion of bremss emission by electrons and positrons.


  double A[constants::NEGP], B[constants::NEGP], C[constants::NEGP], D[constants::NEGP], PAC[constants::NEGP], PDF[constants::NEGP], ZBR;
  int    NBER;

  //  ****  Reading the scaled cross section table.

  fscanf(IRD, "%*45c%lf%*10c%4d%*[^\n]", &ZBR, &NBER);
  getc(IRD);
  
  if(INFO >= 2){ fprintf(IWR, "\n *** Electron scaled bremss x-section,  ZEQ =%12.5E,  NDATA =%4d\n", ZBR, NBER);fflush(IWR);}
  if(NBER != constants::NBE){ penError(ERR_EBRaR_FORMAT); return;}
  mat.ZBR2 = ZBR*ZBR;

  for(unsigned int IE = 0; IE < constants::NBE; IE++)
  {
    fscanf(IRD, "%lf", &cebr01.EBT[IE]);
    for(unsigned int IW = 0; IW < constants::NBW; IW++)
    { 
      if((IW+1)%5==0){
      fscanf(IRD, "%lf%*[^\n]", &cebr01.XS[IE][IW]);
      getc(IRD);
      }
      else{fscanf(IRD, "%lf", &cebr01.XS[IE][IW]);}
    }
    fscanf(IRD, "%lf%*[^\n]", &cebr01.TXS[IE]);
    getc(IRD);
    if(INFO >= 2)
    {
      fprintf(IWR, "%9.2E", cebr01.EBT[IE]);
      int contador_salt = 0;
      for(unsigned int IW = 0; IW < constants::NBW; IW++)
      {
        fprintf(IWR, " %11.5E", cebr01.XS[IE][IW]);
        contador_salt++;
        if(contador_salt == 5)
        {
          fprintf(IWR, "\n         ");
          contador_salt = 0;
        }
      }
      fprintf(IWR, "                                    %10.3E\n", cebr01.TXS[IE]);
    }
    cebr01.X[IE] = log(cebr01.EBT[IE]);
  }

  //  ****  Compute the scaled energy loss distribution and sampling
  //        parameters for the energies in the simulation grid.
  //
  //  ****  Interpolation in E.

  double ELL;
  for(unsigned int IW = 0; IW < constants::NBW; IW++)
  {
    for(unsigned int IE = 0; IE < constants::NBE; IE++)
    {
      cebr01.Y[IE] = log(cebr01.XS[IE][IW]);
    }
    SPLINE(cebr01.X,cebr01.Y,A,B,C,D,0.0,0.0,constants::NBE);
    if(penGetError() != 0){ return;}
    for(unsigned int I = 0; I < constants::NEGP; I++)
    {
      unsigned int J;
      ELL = DLEMP[I];
      FINDI(cebr01.X,ELL,constants::NBE,J);
      if(ELL > cebr01.X[0])
      {
        mat.P0[I][IW] = exp(A[J]+ELL*(B[J]+ELL*(C[J]+ELL*D[J])));
      }
      else
      {
        double F1 = A[0]+cebr01.X[0]*(B[0]+cebr01.X[0]*(C[0]+cebr01.X[0]*D[0]));
        double FP1 = B[0]+cebr01.X[0]*(2.0*C[0]+cebr01.X[0]*3.0*D[0]);
        mat.P0[I][IW] = exp(F1+FP1*(ELL-cebr01.X[0]));
      }
    }
  }
  for(unsigned int IE = 0; IE < constants::NEGP; IE++)
  {
    for(unsigned int IW = 0; IW < constants::NBW; IW++)
    {
      PDF[IW] = mat.P0[IE][IW];
    }

    RLPAC(constants::WB,PDF,PAC,constants::NBW);
    for(unsigned int IW = 0; IW < constants::NBW; IW++)
    {
      mat.PDFB[IE][IW] = PDF[IW];
      mat.PACB[IE][IW] = PAC[IW];
    }
    for(unsigned int IW = 0; IW < constants::NBW-1; IW++)
    {
      mat.DPDFB[IE][IW] = mat.PDFB[IE][IW+1]-mat.PDFB[IE][IW];
    }
    mat.DPDFB[IE][constants::NBW-1] = 0.0;
    //  ****  The cutoff scaled energy loss is slightly modified to ensure
    // that the sampling routine EBR covers the allowed energy loss interval.
    double XC;
    if((unsigned)IE+1 < constants::NEGP)
    {
      XC = WCRM/ET[IE+1];
    }
    else
    {
      XC = WCRM/ET[constants::NEGP-1];
    }
    if(XC < 1.0)
    {
      mat.PBCUT[IE] = RLMOM(constants::WB,PDF,XC,constants::NBW,-1);
      if(penGetError() != 0){ return;}
      mat.WBCUT[IE] = XC;
    }
    else
    {
      mat.PBCUT[IE] = RLMOM(constants::WB,PDF,1.0,constants::NBW,-1);
      if(penGetError() != 0){ return;}
      mat.WBCUT[IE] = 1.0;
    }
  }

}

//  *********************************************************************
//                       SUBROUTINE BRaAR
//  *********************************************************************
void BRaAR(pen_material& mat, FILE* IRD, FILE* IWR, int INFO)
{
  //  This subroutine reads bremsstrahlung angular distribution parameters
  //  of material M from the material data file. It also initialises the
  //  algorithm for generation of the initial direction of bremsstrahlung
  //  photons.


  const double TREV = 2.0*constants::REV;
  
  const int NE = 7;
  const int NK = 10;
  double E[NE], BET[mat.NET], XK[NK], Q1R[NE][NK], Q2R[NE][NK], Q1[NE][mat.NKT], Q2[NE][mat.NKT];
  double X[NK], A[NK], B[NK], C[NK], D[NK];

  double ZEQ;
  int NDATA;

  E[0] = 1.0E3;
  E[1] = 5.0E3;
  E[2] = 1.0E4;
  E[3] = 5.0E4;
  E[4] = 1.0E5;
  E[5] = 5.0E5;
  E[6] = 1.0E6;

  XK[0] = 0.0E0;
  XK[1] = 0.1E0;
  XK[2] = 0.2E0;
  XK[3] = 0.3E0;
  XK[4] = 0.4E0;
  XK[5] = 0.5E0;
  XK[6] = 0.6E0;
  XK[7] = 0.7E0;
  XK[8] = 0.8E0;
  XK[9] = 0.95E0;

  for (unsigned int IE = 0; IE < mat.NET; IE++)
    {
      BET[IE] = sqrt (E[IE] * (E[IE] + TREV)) / (E[IE] + constants::REV);
    }

  //  ****  Grid of reduced photon energies.

  double BK[mat.NKT];
  double DDK = 1.0E0 / double (mat.NKT - 1);
  for(unsigned int IK = 0; IK < mat.NKT; IK++)
  {
    BK[IK] = IK*DDK;
  }

  //  ****  Read angular distribution parameters from file (unit IRD).

  fscanf(IRD, "%*40c%lf%*10c%4d%*[^\n]", &ZEQ, &NDATA);
  getc(IRD);
  
  if(INFO >= 2){ fprintf(IWR,"\n *** Bremss angular distribution,  ZEQ =%12.5E,  NDATA =%4d\n", ZEQ, NDATA);}
  if(NDATA != 70){ penError(ERR_BRAR_INCONSISTENT_DATA); return;}
  if(ZEQ > 1.0){ mat.ZBEQ = ZEQ;}
  else{ mat.ZBEQ = 1.0;}
  
  if(mat.ZBEQ > 99.0){ mat.ZBEQ = 99.0;}

  for(int IE1 = 0; IE1 < NE; IE1++)
  {
    for(int IK1 = 0; IK1 < NK; IK1++)
    {
      int IE, IK;
      double ER,RKR,Q1RR,Q2RR;
      fscanf(IRD, "%d %d %lf %lf %lf %lf%*[^\n]", &IE, &IK, &ER, &RKR, &Q1RR, &Q2RR);
      getc(IRD);
      if((fabs(ER-E[IE-1]) < 1.0E-6) && (fabs(RKR-XK[IK-1]) < 1.0E-6))
      {
        Q1R[IE-1][IK-1] = Q1RR;
        Q2R[IE-1][IK-1] = Q2RR;
      }
      else
      {
        printf("\nCorrupt data file (pdbrang.p08).\n");
        penError(ERR_BRAR_CORRUPT_DATA_FILE); return;
      }
    }
  }

  if(INFO >= 2)
  {
    for(int IE = 0; IE < NE; IE++)
    {
      for(int IK = 0; IK < NK; IK++)
      {
        fprintf(IWR, "%10.3E %10.3E %14.7E %14.7E\n", E[IE], XK[IK], Q1R[IE][IK], Q2R[IE][IK]);
      }
    }
  }

  //  ****  Expanded K-grid of distribution parameters.

  for(int IE = 0; IE < NE; IE++)
  {
    for(int IK = 0; IK < NK; IK++)
    {
      X[IK] = log(Q1R[IE][IK]);
    }
    SPLINE(XK,X,A,B,C,D,0.0,0.0,NK);
    if(penGetError() != PEN_SUCCESS){ return;}
    for(unsigned int IK = 0; IK < mat.NKT; IK++)
    {
      unsigned int J;
      FINDI(XK,BK[IK],NK,J);
      Q1[IE][IK] = A[J]+BK[IK]*(B[J]+BK[IK]*(C[J]+BK[IK]*D[J]));
    }
    for(int IK = 0; IK < NK; IK++)
    {
      X[IK] = Q2R[IE][IK];
    }
    SPLINE(XK,X,A,B,C,D,0.0,0.0,NK);
    if(penGetError() != PEN_SUCCESS){ return;}
    for(unsigned int IK = 0; IK < mat.NKT; IK++)
    {
      unsigned int J;
      FINDI(XK,BK[IK],NK,J);
      Q2[IE][IK] = A[J]+BK[IK]*(B[J]+BK[IK]*(C[J]+BK[IK]*D[J]));
    }
  }

  //  ****  ... and natural cubic spline interpolations.

  for(unsigned int IK = 0; IK < mat.NKT; IK++)
  {
    for(int IE = 0; IE < NE; IE++)
    {
      X[IE] = Q1[IE][IK];
    }
    SPLINE(BET,X,A,B,C,D,0.0,0.0,NE);
    if(penGetError() != PEN_SUCCESS){ return;}
    for(int IE = 0; IE < NE; IE++)
    {
      mat.BP1[IE][IK][0] = A[IE];
      mat.BP1[IE][IK][1] = B[IE];
      mat.BP1[IE][IK][2] = C[IE];
      mat.BP1[IE][IK][3] = D[IE];
    }
    for(int IE = 0; IE < NE; IE++)
    {
      X[IE] = Q2[IE][IK];
    }
    SPLINE(BET,X,A,B,C,D,0.0,0.0,NE);
    if(penGetError() != PEN_SUCCESS){ return;}
    for(int IE = 0; IE < NE; IE++)
    {
      mat.BP2[IE][IK][0] = A[IE];
      mat.BP2[IE][IK][1] = B[IE];
      mat.BP2[IE][IK][2] = C[IE];
      mat.BP2[IE][IK][3] = D[IE];
    }
  }
}

//  *********************************************************************
//                       SUBROUTINE EINaT
//  *********************************************************************
void EINaT(const pen_material& mat, CEIN00& cein00, const double E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA)
{
  //  Integrated cross sections for inelastic collisions of electrons of
  //  energy E in material M, restricted to energy losses larger than and
  //  less than the cutoff energy WCCM.
  //
  //  Sternheimer-Liljequist GOS model.
  //
  //  Output arguments:
  //    XH0 ... total cross section for hard colls. (cm**2).
  //    XH1 ... stopping cross section for hard colls. (eV*cm**2).
  //    XH2 ... straggling cross section for hard colls. (eV**2*cm**2).
  //    XS0 ... total cross section for soft colls. (cm**2).
  //    XS1 ... stopping cross section for soft colls. (eV*cm**2)
  //    XS2 ... straggling cross section for soft colls. (eV**2*cm**2).
  //    XT1 ... 1st transport cross section for soft colls. (cm**2).
  //    XT2 ... 2nd transport cross section for soft colls. (cm**2).
  //    DELTA ... Fermi's density effect correction.

  
  //  ****  Constants.
  
  const double GAM = 1.0+E/constants::REV;
  const double GAM2 = GAM*GAM;

  //  ************  Density effect.

  //  ****  Sternheimer's resonance energy (WL2=L**2).
  double TST = mat.ZT/(GAM2*mat.OP2);
  double WL2 = 0.0;
  double WL2U, WL2L;
  double FDEL = 0.0;
  bool Eixir = false;
  for(int I = 0; I < mat.NOSC; I++)
  {
    FDEL = FDEL+mat.F[I]/(pow(mat.WRI[I],2)+WL2);
  }
  if(FDEL < TST)
  {
    DELTA = 0.0;
  }
  else
  {
    WL2 = pow(mat.WRI[mat.NOSC-1],2);
    Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      WL2 = WL2+WL2;
      FDEL = 0.0;
      for(int I = 0; I < mat.NOSC; I++)
      {
        FDEL = FDEL+mat.F[I]/(pow(mat.WRI[I],2)+WL2);
      }
      if(FDEL > TST){ Eixir = false; continue;}
    }
    WL2L = 0.0;
    WL2U = WL2;
    Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      WL2 = 0.5*(WL2L+WL2U);
      FDEL = 0.0;
      for(int I = 0; I < mat.NOSC; I++)
      {
        FDEL = FDEL+mat.F[I]/(pow(mat.WRI[I],2)+WL2);
      }
      if(FDEL > TST)
      {
        WL2L = WL2;
      }
      else
      {
        WL2U = WL2;
      }
      if(WL2U-WL2L > 1.0E-12*WL2){ Eixir = false; continue;}
    }
      //  ****  Density effect correction (delta).
    DELTA = 0.0;
    for(int I = 0; I < mat.NOSC; I++)
    {
      DELTA = DELTA + mat.F[I]*log(1.0+WL2/pow(mat.WRI[I],2));
    }
    DELTA = (DELTA/mat.ZT)-WL2/(GAM2*mat.OP2);
  }

  //  ****  Shell-oscillator cross sections.

  for(int I = 0; I < mat.NOSC; I++)
  {
    cein00.SEH0[I] = 0.0;
    cein00.SEH1[I] = 0.0;
    cein00.SEH2[I] = 0.0;
    cein00.SES0[I] = 0.0;
    cein00.SES1[I] = 0.0;
    cein00.SES2[I] = 0.0;
    cein00.SET0[I] = 0.0;
    cein00.SET1[I] = 0.0;
    cein00.SET2[I] = 0.0;
  }
  XH0 = 0.0;
  XH1 = 0.0;
  XH2 = 0.0;
  XS0 = 0.0;
  XS1 = 0.0;
  XS2 = 0.0;
  double XT0 = 0.0;
  XT1 = 0.0;
  XT2 = 0.0;

  double UK, WK;
  double H0, H1, H2, S0, S1, S2, R0, R1, R2;
  for(int K = 0; K < mat.NOSC; K++)
  {
    UK = mat.UI[K];
    WK = mat.WRI[K];
    EINaT1(E,UK,WK,DELTA,WCCM,H0,H1,H2,S0,S1,S2,R0,R1,R2);
    cein00.SEH0[K] = mat.F[K]*H0;
    cein00.SEH1[K] = mat.F[K]*H1;
    cein00.SEH2[K] = mat.F[K]*H2;
    cein00.SES0[K] = mat.F[K]*S0;
    cein00.SES1[K] = mat.F[K]*S1;
    cein00.SES2[K] = mat.F[K]*S2;
    cein00.SET0[K] = mat.F[K]*R0;
    cein00.SET1[K] = mat.F[K]*2.0*R1;
    cein00.SET2[K] = mat.F[K]*6.0*(R1-R2);
    XH0 += cein00.SEH0[K];
    XH1 += cein00.SEH1[K];
    XH2 += cein00.SEH2[K];
    XS0 += cein00.SES0[K];
    XS1 += cein00.SES1[K];
    XS2 += cein00.SES2[K];
    XT0 += cein00.SET0[K];
    XT1 += cein00.SET1[K];
    XT2 += cein00.SET2[K];
  }
  
}

//  *********************************************************************
//                       SUBROUTINE EINaT1
//  *********************************************************************
void EINaT1(const double E, double &UK, double& WK, double DELTA, double &WCCM, double &H0, double &H1, double &H2, double &S0, double &S1, double &S2, double &R0, double &R1, double &R2)
{
  //  Integrated cross sections for inelastic collisions of electrons with
  //  a single-shell oscillator, restricted to energy losses larger than,
  //  and smaller than, the cutoff energy loss WCCM.
  //
  //  Sternheimer-Liljequist oscillator model.
  //
  //  Input arguments:
  //    E ..... kinetic energy (eV).
  //    UK .... ionisation energy (eV).
  //    WK .... resonance energy (eV).
  //    DELTA ... Fermi's density effect correction.
  //    WCCM ... cutoff energy loss (eV).
  //
  //  Output arguments:
  //    H0 .... total cross section for hard colls. (cm**2).
  //    H1 .... stopping cross section for hard colls. (eV*cm**2).
  //    H2 .... straggling cross section for hard colls. (eV**2*cm**2).
  //    S0 .... total cross section for soft colls. (cm**2).
  //    S1 .... stopping cross section for soft colls. (eV*cm**2).
  //    S2 .... straggling cross section for soft colls. (eV**2*cm**2).
  //    R0 .... total cross section for soft colls. (cm**2).
  //    R1 .... 1st transport cross section for soft colls. (cm**2).
  //    R2 .... 2nd transport cross section for soft colls. (cm**2).

  
  const double TREV = 2.0*constants::REV;
  const double RTREV = 1.0/TREV;
  const double PIELR2 = constants::PI*constants::ELRAD*constants::ELRAD;

  CEIN01 cein01;
  
  H0 = 0.0;
  H1 = 0.0;
  H2 = 0.0;
  S0 = 0.0;
  S1 = 0.0;
  S2 = 0.0;
  R0 = 0.0;
  R1 = 0.0;
  R2 = 0.0;

  double WTHR;
  if(UK > 1.0E-3)
  {
    WTHR = UK;
  }
  else
  {
    WTHR = WK;
  }
  if(E < WTHR+1.0E-6){return;}

  //  ****  Constants.

  cein01.EI = E;
  const double GAM = 1.0+E/constants::REV;
  const double GAM2 = GAM*GAM;
  const double BETA2 = (GAM2-1.0)/GAM2;
  const double CONST = PIELR2*TREV/BETA2;

  cein01.CPS = E*(E+TREV);
  double CP = sqrt(cein01.CPS);
  cein01.AMOL = pow(E/(E+constants::REV),2);

  //  ****  Trick: The resonance energy and the cutoff recoil energy of
  //        inner shells are varied to yield a smooth threshold.

  double WM, WKP, QKP, WCMAX, WDMAX;
  if(UK > 1.0E-3)
  {
    WM = 3.0*WK-2.0*UK;
    if(E > WM)
    {
      WKP = WK;
      QKP = UK;
    }
    else
    {
      WKP = (E+2.0*UK)/3.0;
      QKP = UK*(E/WM);
      WM = E;
    }
    cein01.EE = E+UK;
    WCMAX = 0.5*cein01.EE;
    if(WCMAX < WM){ WDMAX = WCMAX;}
    else{ WDMAX = WM;}
  }

  else
  {
    WM = E;
    WKP = WK;
    QKP = WK;
    cein01.EE = E;
    WCMAX = 0.5*cein01.EE;
    WDMAX = WKP+1.0;
  }

  //  ****  Distant interactions.

  double SDL1 = 0.0;
  double SDT1 = 0.0;
  double CPPS, CPP;
  double A, B, BA;
  double QM;
  double RMU1;
  if(WDMAX > WTHR+1.0E-6)
  {
    CPPS = (E-WKP)*(E-WKP+TREV);
    CPP = sqrt(CPPS);
    A = 4.0*CP*CPP;
    B = pow(CP-CPP,2);

    if(WKP > 1.0E-6*E)
    {
      QM = sqrt(pow(CP-CPP,2)+pow(constants::REV,2))-constants::REV;
    }
    else
    {
      QM = WKP*WKP/(BETA2*TREV);
      QM = QM*(1.0-QM*RTREV);
    }
    if(QM < QKP)
    {
      SDL1 = log(QKP*(QM+TREV)/(QM*(QKP+TREV)));
      SDT1 = log(GAM2)-BETA2-DELTA;
      if(SDT1 < 0.0){ SDT1 = 0.0;}
    
    //  ****  Soft distant transport moments of orders 0-2.
      if(WCCM > WTHR)
      {
        BA = B/A;
        RMU1 = (QKP*(QKP+TREV)-B)/A;
        R0 = log((RMU1+BA)/BA);
        R1 = RMU1-BA*R0;
        R2 = pow(BA,2)*R0+0.5*RMU1*(RMU1-2.0*BA);
        R0 = R0/WKP;
        R1 = R1/WKP;
        R2 = R2/WKP;
        R0 = R0+SDT1/WKP;
      }
    }
  }

  double F0, F1, WL, WU;
  double SD1 = SDL1+SDT1;
  if(SD1 > 0.0)
  {
    if(UK > 1.0E-3)
    {
    //  ****  Inner-shell excitations (triangle distribution).
      F0 = 1.0/pow(WM-UK,2);
      F1 = 2.0*F0*SD1/WKP;
      if(WCCM < UK)
      {
        WL = UK;
        WU = WDMAX;
        H0 = F1*(WM*(WU-WL)-(pow(WU,2)-pow(WL,2))/2.0);
        H1 = F1*(WM*(pow(WU,2)-pow(WL,2))/2.0-(pow(WU,3)-pow(WL,3))/3.0);
        H2 = F1*(WM*(pow(WU,3)-pow(WL,3))/3.0-(pow(WU,4)-pow(WL,4))/4.0);
      }
      else
      {
        if(WCCM > WDMAX)
        {
          WL = UK;
          WU = WDMAX;
          S0 = F1*(WM*(WU-WL)-(pow(WU,2)-pow(WL,2))/2.0);
          S1 = F1*(WM*(pow(WU,2)-pow(WL,2))/2.0-(pow(WU,3)-pow(WL,3))/3.0);
          S2 = F1*(WM*(pow(WU,3)-pow(WL,3))/3.0-(pow(WU,4)-pow(WL,4))/4.0);
        }
        else
        {
          WL = WCCM;
          WU = WDMAX;
          H0 = F1*(WM*(WU-WL)-(pow(WU,2)-pow(WL,2))/2.0);
          H1 = F1*(WM*(pow(WU,2)-pow(WL,2))/2.0-(pow(WU,3)-pow(WL,3))/3.0);
          H2 = F1*(WM*(pow(WU,3)-pow(WL,3))/3.0-(pow(WU,4)-pow(WL,4))/4.0);
          WL = UK;
          WU = WCCM;
          S0 = F1*(WM*(WU-WL)-(pow(WU,2)-pow(WL,2))/2.0);
          S1 = F1*(WM*(pow(WU,2)-pow(WL,2))/2.0-(pow(WU,3)-pow(WL,3))/3.0);
          S2 = F1*(WM*(pow(WU,3)-pow(WL,3))/3.0-(pow(WU,4)-pow(WL,4))/4.0);
        }
        double F2 = F0*(2.0*WM*(WU-WL)-(pow(WU,2)-pow(WL,2)));
        R0 = F2*R0;
        R1 = F2*R1;
        R2 = F2*R2;
      }
    }
    else
    {
    //  ****  Outer-shell excitations (delta oscillator).
      if(WCCM < WKP)
      {
        H1 = SD1;
        H0 = SD1/WKP;
        H2 = SD1*WKP;
      }
      else
      {
        S1 = SD1;
        S0 = SD1/WKP;
        S2 = SD1*WKP;
      }
    }
  }

  //  ****  Close collisions (Moller's cross section).

  if(WCMAX < WTHR+1.0E-6){}
  else
  {
    if(WCCM < WTHR)   // No soft interactions.
    {
      WL = WTHR;
      WU = WCMAX;
      H0 = H0+(1.0/(cein01.EE-WU))-(1.0/(cein01.EE-WL))-(1.0/WU)+(1.0/WL)+(1.0-cein01.AMOL)*log(((cein01.EE-WU)*WL)/((cein01.EE-WL)*WU))/cein01.EE+cein01.AMOL*(WU-WL)/pow(cein01.EE,2);
    
      H1 = H1+log(WU/WL)+(cein01.EE/(cein01.EE-WU))-(cein01.EE/(cein01.EE-WL))+(2.0-cein01.AMOL)*log((cein01.EE-WU)/(cein01.EE-WL))+cein01.AMOL*(pow(WU,2)-pow(WL,2))/(2.0*pow(cein01.EE,2));
    
      H2 = H2+(2.0-cein01.AMOL)*(WU-WL)+(WU*(2.0*cein01.EE-WU)/(cein01.EE-WU))-(WL*(2.0*cein01.EE-WL)/(cein01.EE-WL))+(3.0-cein01.AMOL)*cein01.EE*log((cein01.EE-WU)/(cein01.EE-WL))+cein01.AMOL*(pow(WU,3)-pow(WL,3))/(3.0*pow(cein01.EE,2));
    }
    else
    {
      if(WCCM > WCMAX)
      {
        WL = WTHR;
        WU = WCMAX;
        S0 = S0+(1.0/(cein01.EE-WU))-(1.0/(cein01.EE-WL))-(1.0/WU)+(1.0/WL)+(1.0-cein01.AMOL)*log(((cein01.EE-WU)*WL)/((cein01.EE-WL)*WU))/cein01.EE+cein01.AMOL*(WU-WL)/pow(cein01.EE,2);
        
        S1 = S1+log(WU/WL)+(cein01.EE/(cein01.EE-WU))-(cein01.EE/(cein01.EE-WL))+(2.0-cein01.AMOL)*log((cein01.EE-WU)/(cein01.EE-WL))+cein01.AMOL*(pow(WU,2)-pow(WL,2))/(2.0*pow(cein01.EE,2));
        
        S2 = S2+(2.0-cein01.AMOL)*(WU-WL)+(WU*(2.0*cein01.EE-WU)/(cein01.EE-WU))-(WL*(2.0*cein01.EE-WL)/(cein01.EE-WL))+(3.0-cein01.AMOL)*cein01.EE*log((cein01.EE-WU)/(cein01.EE-WL))+cein01.AMOL*(pow(WU,3)-pow(WL,3))/(3.0*pow(cein01.EE,2));
      }
      else
      {
        WL = WCCM;
        WU = WCMAX;
        H0 = H0+(1.0/(cein01.EE-WU))-(1.0/(cein01.EE-WL))-(1.0/WU)+(1.0/WL)+(1.0-cein01.AMOL)*log(((cein01.EE-WU)*WL)/((cein01.EE-WL)*WU))/cein01.EE+cein01.AMOL*(WU-WL)/pow(cein01.EE,2);
        
        H1 = H1+log(WU/WL)+(cein01.EE/(cein01.EE-WU))-(cein01.EE/(cein01.EE-WL))+(2.0-cein01.AMOL)*log((cein01.EE-WU)/(cein01.EE-WL))+cein01.AMOL*(pow(WU,2)-pow(WL,2))/(2.0*pow(cein01.EE,2));
        
        H2 = H2+(2.0-cein01.AMOL)*(WU-WL)+(WU*(2.0*cein01.EE-WU)/(cein01.EE-WU))-(WL*(2.0*cein01.EE-WL)/(cein01.EE-WL))+(3.0-cein01.AMOL)*cein01.EE*log((cein01.EE-WU)/(cein01.EE-WL))+cein01.AMOL*(pow(WU,3)-pow(WL,3))/(3.0*pow(cein01.EE,2));
        WL = WTHR;
        WU = WCCM;
        S0 = S0+(1.0/(cein01.EE-WU))-(1.0/(cein01.EE-WL))-(1.0/WU)+(1.0/WL)+(1.0-cein01.AMOL)*log(((cein01.EE-WU)*WL)/((cein01.EE-WL)*WU))/cein01.EE+cein01.AMOL*(WU-WL)/pow(cein01.EE,2);
        
        S1 = S1+log(WU/WL)+(cein01.EE/(cein01.EE-WU))-(cein01.EE/(cein01.EE-WL))+(2.0-cein01.AMOL)*log((cein01.EE-WU)/(cein01.EE-WL))+cein01.AMOL*(pow(WU,2)-pow(WL,2))/(2.0*pow(cein01.EE,2));
        
        S2 = S2+(2.0-cein01.AMOL)*(WU-WL)+(WU*(2.0*cein01.EE-WU)/(cein01.EE-WU))-(WL*(2.0*cein01.EE-WL)/(cein01.EE-WL))+(3.0-cein01.AMOL)*cein01.EE*log((cein01.EE-WU)/(cein01.EE-WL))+cein01.AMOL*(pow(WU,3)-pow(WL,3))/(3.0*pow(cein01.EE,2));
      }
      //  ****  Soft close transport moments of orders 0-2.
      double CP2S = (E-WL)*(E-WL+TREV);
      double CP2 = sqrt(CP2S);
      double RMU2 = (WL*(WL+TREV)-pow(CP-CP2,2))/(4.0*CP*CP2);
      double CP3S = (E-WU)*(E-WU+TREV);
      double CP3 = sqrt(CP3S);
      double RMU3 = (WU*(WU+TREV)-pow(CP-CP3,2))/(4.0*CP*CP3);
      cein01.MOM = 0;
      R0 = R0+SUMGA(EINaDS,&cein01,RMU2,RMU3,1.0E-7);
      cein01.MOM = 1;
      R1 = R1+SUMGA(EINaDS,&cein01,RMU2,RMU3,1.0E-7);
      cein01.MOM = 2;
      R2 = R2+SUMGA(EINaDS,&cein01,RMU2,RMU3,1.0E-7);
    }
  }
      
  H0 = CONST*H0;
  H1 = CONST*H1;
  H2 = CONST*H2;
  S0 = CONST*S0;
  S1 = CONST*S1;
  S2 = CONST*S2;
  R0 = CONST*R0;
  R1 = CONST*R1;
  R2 = CONST*R2;      
}

void PINaT(const pen_material& mat, CPIN00& cpin00, const double E, double &WCCM, double &XH0, double &XH1, double &XH2, double &XS0, double &XS1, double &XS2, double &XT1, double &XT2, double &DELTA)
{
  //  Integrated cross sections for inelastic collisions of positrons of
  //  energy E in material M, restricted to energy losses larger than and
  //  less than the cutoff energy WCCM.
  //
  //  Sternheimer-Liljequist GOS model.
  //
  //  Output arguments:
  //    XH0 ... total cross section for hard colls. (cm**2).
  //    XH1 ... stopping cross section for hard colls. (eV*cm**2).
  //    XH2 ... straggling cross section for hard colls. (eV**2*cm**2).
  //    XS0 ... total cross section for soft colls. (cm**2).
  //    XS1 ... stopping cross section for soft colls. (eV*cm**2)
  //    XS2 ... straggling cross section for soft colls. (eV**2*cm**2).
  //    XT1 ... 1st transport cross section for soft colls. (cm**2).
  //    XT2 ... 2nd transport cross section for soft colls. (cm**2).
  //    DELTA ... Fermi's density effect correction.


  //  ****  Constants.

  const double GAM = 1.0+E/constants::REV;
  const double GAM2 = GAM*GAM;

  CPIN01 cpin01;
  
  //  ************  Density effect.

  //  ****  Sternheimer's resonance energy (WL2=L**2).
  bool Eixir = false;
  double TST, WL2, FDEL, WL2L, WL2U;
  TST = mat.ZT/(GAM2*mat.OP2);
  WL2 = 0.0;
  FDEL = 0.0;
  for(int I = 0; I < mat.NOSC; I++)
  {
    FDEL = FDEL+mat.F[I]/(pow(mat.WRI[I],2)+WL2);
  }
  if(FDEL < TST)
  {
    DELTA = 0.0;
  }
  else
  {
    WL2 = pow(mat.WRI[mat.NOSC-1],2);
    Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      WL2 = WL2+WL2;
      FDEL = 0.0;
      for(int I = 0; I < mat.NOSC; I++)
      {
        FDEL = FDEL+mat.F[I]/(pow(mat.WRI[I],2)+WL2);
      }
      if(FDEL > TST){ Eixir = false; continue;}
    }
    WL2L = 0.0;
    WL2U = WL2;
    Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      WL2 = 0.5*(WL2L+WL2U);
      FDEL = 0.0;
      for(int I = 0; I < mat.NOSC; I++)
      {
        FDEL = FDEL+mat.F[I]/(pow(mat.WRI[I],2)+WL2);
      }
      if(FDEL > TST)
      {
        WL2L = WL2;
      }
      else
      {
        WL2U = WL2;
      }
      if(WL2U-WL2L > 1.0E-12*WL2){ Eixir = false; continue;}
    }
      //  ****  Density effect correction (delta).
    DELTA = 0.0;
    for(int I = 0; I < mat.NOSC; I++)
    {
      DELTA = DELTA+mat.F[I]*log(1.0+WL2/pow(mat.WRI[I],2));
    }
    DELTA = (DELTA/mat.ZT)-WL2/(GAM2*mat.OP2);
  }

  //  ****  Shell-oscillator cross sections.

  for(int I = 0; I < mat.NOSC; I++)
  {
    cpin00.SPH0[I] = 0.0;
    cpin00.SPH1[I] = 0.0;
    cpin00.SPH2[I] = 0.0;
    cpin00.SPS0[I] = 0.0;
    cpin00.SPS1[I] = 0.0;
    cpin00.SPS2[I] = 0.0;
    cpin00.SPT0[I] = 0.0;
    cpin00.SPT1[I] = 0.0;
    cpin00.SPT2[I] = 0.0;
  }
  XH0 = 0.0;
  XH1 = 0.0;
  XH2 = 0.0;
  XS0 = 0.0;
  XS1 = 0.0;
  XS2 = 0.0;
  double XT0 = 0.0;
  XT1 = 0.0;
  XT2 = 0.0;

  double UK, WK, H0, H1, H2, S0, S1, S2, R0, R1, R2;
  for(int K = 0; K < mat.NOSC; K++)
  {
    UK = mat.UI[K];
    WK = mat.WRI[K];
    PINaT1(cpin01,E,UK,WK,DELTA,WCCM,H0,H1,H2,S0,S1,S2,R0,R1,R2);
    cpin00.SPH0[K] = mat.F[K]*H0;
    cpin00.SPH1[K] = mat.F[K]*H1;
    cpin00.SPH2[K] = mat.F[K]*H2;
    cpin00.SPS0[K] = mat.F[K]*S0;
    cpin00.SPS1[K] = mat.F[K]*S1;
    cpin00.SPS2[K] = mat.F[K]*S2;
    cpin00.SPT0[K] = mat.F[K]*R0;
    cpin00.SPT1[K] = mat.F[K]*2.0*R1;
    cpin00.SPT2[K] = mat.F[K]*6.0*(R1-R2);
    XH0 = XH0+cpin00.SPH0[K];
    XH1 = XH1+cpin00.SPH1[K];
    XH2 = XH2+cpin00.SPH2[K];
    XS0 = XS0+cpin00.SPS0[K];
    XS1 = XS1+cpin00.SPS1[K];
    XS2 = XS2+cpin00.SPS2[K];
    XT0 = XT0+cpin00.SPT0[K];
    XT1 = XT1+cpin00.SPT1[K];
    XT2 = XT2+cpin00.SPT2[K];
  }

}

//  *********************************************************************
//                       SUBROUTINE PINaT1
//  *********************************************************************
void PINaT1(CPIN01& cpin01, const double E, double &UK, double &WK, double DELTA, double &WCCM, double &H0, double &H1, double &H2, double &S0, double &S1, double &S2, double &R0, double &R1, double &R2)
{
//  Integrated cross sections for inelastic collisions of positrons with
//  a single-shell oscillator, restricted to energy losses larger than,
//  and smaller than, the cutoff energy loss WCCM.
//
//  Sternheimer-Liljequist oscillator model.
//
//  Input arguments:
//    E ..... kinetic energy (eV).
//    UK .... ionisation energy (eV).
//    WK .... resonance energy (eV).
//    DELTA ... Fermi's density effect correction.
//    WCCM ... cutoff energy loss (eV).
//
//  Output arguments:
//    H0 .... total cross section for hard colls. (cm**2).
//    H1 .... stopping cross section for hard colls. (eV*cm**2).
//    H2 .... straggling cross section for hard colls. (eV**2*cm**2).
//    S0 .... total cross section for soft colls. (cm**2).
//    S1 .... stopping cross section for soft colls. (eV*cm**2).
//    S2 .... straggling cross section for soft colls. (eV**2*cm**2).
//    R0 .... total cross section for soft colls. (cm**2).
//    R1 .... 1st transport cross section for soft colls. (cm**2).
//    R2 .... 2nd transport cross section for soft colls. (cm**2).
//

  
  const double TREV = 2.0*constants::REV;
  const double RTREV = 1.0/TREV;
  const double PIELR2 = constants::PI*constants::ELRAD*constants::ELRAD;
  

  H0 = 0.0;
  H1 = 0.0;
  H2 = 0.0;
  S0 = 0.0;
  S1 = 0.0;
  S2 = 0.0;
  R0 = 0.0;
  R1 = 0.0;
  R2 = 0.0;

  double WTHR;
  if(UK > 1.0E-3)
  {
    WTHR = UK;
  }
  else
  {
    WTHR = WK;
  }
  if(E < WTHR+1.0E-6){ return;}

  //  ****  Constants.

  cpin01.EI = E;
  const double GAM = 1.0+E/constants::REV;
  const double GAM2 = GAM*GAM;
  const double BETA2 = (GAM2-1.0)/GAM2;
  const double CONST = PIELR2*TREV/BETA2;

  cpin01.CPS = E*(E+TREV);
  double CP = sqrt(cpin01.CPS);
  double AMOL = pow(E/(E+constants::REV),2);
  double G12 = pow(GAM+1.0,2);
  cpin01.BHA1 = AMOL*(2.0*G12-1.0)/(GAM2-1.0);
  cpin01.BHA2 = AMOL*(3.0+1.0/G12);
  cpin01.BHA3 = AMOL*2.0*GAM*(GAM-1.0)/G12;
  cpin01.BHA4 = AMOL*pow(GAM-1.0,2)/G12;

  //  ****  Trick: The resonance energy and the cutoff recoil energy of
  //        inner shells are varied to yield a smooth threshold. 

  double WM, WKP, QKP, WCMAX, WDMAX;
  if(UK > 1.0E-3)
  {
    WM = 3.0*WK-2.0*UK;
    if(E > WM)
    {
      WKP = WK;
      QKP = UK;
    }
    else
    {
      WKP = (E+2.0*UK)/3.0;
      QKP = UK*(E/WM);
      WM = E;
    }
    WCMAX = E;
    if(WM < WCMAX){ WDMAX = WM;}
    else{ WDMAX = WCMAX;}
  }
  else
  {
    WM = E;
    WKP = WK;
    QKP = WK;
    WCMAX = E;
    WDMAX = WKP+1.0;
  }     

  //  ****  Distant interactions.

  
  double SDL1 = 0.0;
  double SDT1 = 0.0;
  double CPPS, CPP, A, B, QM, BA, RMU1;
  if(WDMAX > WTHR+1.0E-6)
  {
    CPPS = (E-WKP)*(E-WKP+TREV);
    CPP = sqrt(CPPS);
    A = 4.0*CP*CPP;
    B = pow(CP-CPP,2);
      
    if(WKP > 1.0E-6*E)
    {
      QM = sqrt(pow(CP-CPP,2)+pow(constants::REV,2))-constants::REV;
    }
    else
    {
      QM = WKP*WKP/(BETA2*TREV);
      QM = QM*(1.0-QM*RTREV);
    }
    if(QM < QKP)
    {
      SDL1 = log(QKP*(QM+TREV)/(QM*(QKP+TREV)));

      SDT1 = log(GAM2)-BETA2-DELTA;
      if(SDT1 < 0.0){ SDT1 = 0.0;}
    
    //  ****  Soft distant transport moments of orders 0-2.
      if(WCCM > WTHR)
      {
        BA = B/A;
        RMU1 = (QKP*(QKP+TREV)-B)/A;
        R0 = log((RMU1+BA)/BA);
        R1 = RMU1-BA*R0;
        R2 = pow(BA,2)*R0+0.5*RMU1*(RMU1-2.0*BA);
        R0 = R0/WKP;
        R1 = R1/WKP;
        R2 = R2/WKP;
        R0 = R0+SDT1/WKP;
      }
    }
  }

  double SD1 = SDL1+SDT1;
  double F0, F1, WL, WU;
  if(SD1 > 0.0)
  {
    if(UK > 1.0E-3)
    {
    //  ****  Inner-shell excitations (triangle distribution).
      F0 = 1.0/pow(WM-UK,2);
      F1 = 2.0*F0*SD1/WKP;
      if(WCCM < UK)
      {
        WL = UK;
        WU = WDMAX;
        H0 = F1*(WM*(WU-WL)-(pow(WU,2)-pow(WL,2))/2.0);
        H1 = F1*(WM*(pow(WU,2)-pow(WL,2))/2.0-(pow(WU,3)-pow(WL,3))/3.0);
        H2 = F1*(WM*(pow(WU,3)-pow(WL,3))/3.0-(pow(WU,4)-pow(WL,4))/4.0);
      }
      else
      {
        if(WCCM > WDMAX)
        {
          WL = UK;
          WU = WDMAX;
          S0 = F1*(WM*(WU-WL)-(pow(WU,2)-pow(WL,2))/2.0);
          S1 = F1*(WM*(pow(WU,2)-pow(WL,2))/2.0-(pow(WU,3)-pow(WL,3))/3.0);
          S2 = F1*(WM*(pow(WU,3)-pow(WL,3))/3.0-(pow(WU,4)-pow(WL,4))/4.0);
        }
        else
        {
          WL = WCCM;
          WU = WDMAX;
          H0 = F1*(WM*(WU-WL)-(pow(WU,2)-pow(WL,2))/2.0);
          H1 = F1*(WM*(pow(WU,2)-pow(WL,2))/2.0-(pow(WU,3)-pow(WL,3))/3.0);
          H2 = F1*(WM*(pow(WU,3)-pow(WL,3))/3.0-(pow(WU,4)-pow(WL,4))/4.0);
          WL = UK;
          WU = WCCM;
          S0 = F1*(WM*(WU-WL)-(pow(WU,2)-pow(WL,2))/2.0);
          S1 = F1*(WM*(pow(WU,2)-pow(WL,2))/2.0-(pow(WU,3)-pow(WL,3))/3.0);
          S2 = F1*(WM*(pow(WU,3)-pow(WL,3))/3.0-(pow(WU,4)-pow(WL,4))/4.0);
        }
        double F2 = F0*(2.0*WM*(WU-WL)-(pow(WU,2)-pow(WL,2)));
        R0 = F2*R0;
        R1 = F2*R1;
        R2 = F2*R2;
      }
    }
    else
    {
    //  ****  Outer-shell excitations (delta oscillator).
      if(WCCM < WKP)
      {
        H1 = SD1;
        H0 = SD1/WKP;
        H2 = SD1*WKP;
      }
      else
      {
        S1 = SD1;
        S0 = SD1/WKP;
        S2 = SD1*WKP;
      }
    }
  }

  //  ****  Close collisions (Bhabha's cross section).

  if(WCMAX < WTHR+1.0E-6){}
  else
  {
    if(WCCM < WTHR)  // No soft interactions.
    {
      WL = WTHR;
      WU = WCMAX;
      H0 = H0+(1.0/WL)-(1.0/WU)-cpin01.BHA1*log(WU/WL)/E+cpin01.BHA2*(WU-WL)/pow(E,2)-cpin01.BHA3*(pow(WU,2)-pow(WL,2))/(2.0*pow(E,3))+cpin01.BHA4*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,4));

      H1 = H1+log(WU/WL)-cpin01.BHA1*(WU-WL)/E+cpin01.BHA2*(pow(WU,2)-pow(WL,2))/(2.0*pow(E,2))-cpin01.BHA3*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,3))+cpin01.BHA4*(pow(WU,4)-pow(WL,4))/(4.0*pow(E,4));
    
      H2 = H2+WU-WL-cpin01.BHA1*(pow(WU,2)-pow(WL,2))/(2.0*E)+cpin01.BHA2*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,2))-cpin01.BHA3*(pow(WU,4)-pow(WL,4))/(4.0*pow(E,3))+cpin01.BHA4*(pow(WU,5)-pow(WL,5))/(5.0*pow(E,4));
    }
    else
    {
      if(WCCM > WCMAX)
      {
        WL = WTHR;
        WU = WCMAX;
        S0 = S0+(1.0/WL)-(1.0/WU)-cpin01.BHA1*log(WU/WL)/E+cpin01.BHA2*(WU-WL)/pow(E,2)-cpin01.BHA3*(pow(WU,2)-pow(WL,2))/(2.0*pow(E,3))+cpin01.BHA4*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,4));
        
        S1 = S1+log(WU/WL)-cpin01.BHA1*(WU-WL)/E+cpin01.BHA2*(pow(WU,2)-pow(WL,2))/(2.0*pow(E,2))-cpin01.BHA3*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,3))+cpin01.BHA4*(pow(WU,4)-pow(WL,4))/(4.0*pow(E,4));

        S2 = S2+WU-WL-cpin01.BHA1*(pow(WU,2)-pow(WL,2))/(2.0*E)+cpin01.BHA2*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,2))-cpin01.BHA3*(pow(WU,4)-pow(WL,4))/(4.0*pow(E,3))+cpin01.BHA4*(pow(WU,5)-pow(WL,5))/(5.0*pow(E,4));
      }
      else
      {
        WL = WCCM;
        WU = WCMAX;
        H0 = H0+(1.0/WL)-(1.0/WU)-cpin01.BHA1*log(WU/WL)/E+cpin01.BHA2*(WU-WL)/pow(E,2)-cpin01.BHA3*(pow(WU,2)-pow(WL,2))/(2.0*pow(E,3))+cpin01.BHA4*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,4));
        
        H1 = H1+log(WU/WL)-cpin01.BHA1*(WU-WL)/E+cpin01.BHA2*(pow(WU,2)-pow(WL,2))/(2.0*pow(E,2))-cpin01.BHA3*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,3))+cpin01.BHA4*(pow(WU,4)-pow(WL,4))/(4.0*pow(E,4));
        
        H2 = H2+WU-WL-cpin01.BHA1*(pow(WU,2)-pow(WL,2))/(2.0*E)+cpin01.BHA2*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,2))-cpin01.BHA3*(pow(WU,4)-pow(WL,4))/(4.0*pow(E,3))+cpin01.BHA4*(pow(WU,5)-pow(WL,5))/(5.0*pow(E,4));
        WL = WTHR;
        WU = WCCM;
        S0 = S0+(1.0/WL)-(1.0/WU)-cpin01.BHA1*log(WU/WL)/E+cpin01.BHA2*(WU-WL)/pow(E,2)-cpin01.BHA3*(pow(WU,2)-pow(WL,2))/(2.0*pow(E,3))+cpin01.BHA4*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,4));
        
        S1 = S1+log(WU/WL)-cpin01.BHA1*(WU-WL)/E+cpin01.BHA2*(pow(WU,2)-pow(WL,2))/(2.0*pow(E,2))-cpin01.BHA3*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,3))+cpin01.BHA4*(pow(WU,4)-pow(WL,4))/(4.0*pow(E,4));
        
        S2 = S2+WU-WL-cpin01.BHA1*(pow(WU,2)-pow(WL,2))/(2.0*E)+cpin01.BHA2*(pow(WU,3)-pow(WL,3))/(3.0*pow(E,2))-cpin01.BHA3*(pow(WU,4)-pow(WL,4))/(4.0*pow(E,3))+cpin01.BHA4*(pow(WU,5)-pow(WL,5))/(5.0*pow(E,4));
      }
      //  ****  Soft close transport moments of orders 0-2.
      double CP2S = (E-WL)*(E-WL+TREV);
      double CP2 = sqrt(CP2S);
      double RMU2 = (WL*(WL+TREV)-pow(CP-CP2,2))/(4.0*CP*CP2);
      double CP3S, CP3, RMU3;
      if(WU < E-1.0)
      {
        CP3S = (E-WU)*(E-WU+TREV);
        CP3 = sqrt(CP3S);
        RMU3 = (WU*(WU+TREV)-pow(CP-CP3,2))/(4.0*CP*CP3);
      }
      else
      {
        RMU3 = 0.5;
      }
      cpin01.MOM = 0;
      R0 = R0+SUMGA(PINaDS,&cpin01,RMU2,RMU3,1.0E-7);
      cpin01.MOM = 1;
      R1 = R1+SUMGA(PINaDS,&cpin01,RMU2,RMU3,1.0E-7);
      cpin01.MOM = 2;
      R2 = R2+SUMGA(PINaDS,&cpin01,RMU2,RMU3,1.0E-7);
    }
  }

  H0 = CONST*H0;
  H1 = CONST*H1;
  H2 = CONST*H2;
  S0 = CONST*S0;
  S1 = CONST*S1;
  S2 = CONST*S2;
  R0 = CONST*R0;
  R1 = CONST*R1;
  R2 = CONST*R2;

}

//  *********************************************************************
//                       SUBROUTINE EBRaT
//  *********************************************************************
void EBRaT(const pen_material& mat, const double &E, double &WCRM, double &XH0, double &XH1, double &XH2,double &XS1, double &XS2, const double DLEMP1, const double DLFC, const double* WB)
{
  //  Integrated cross sections for bremss emission by electrons of energy
  //  E in material M, restricted to energy losses larger than and less
  //  than the cutoff energy WCRM.

  //  Output arguments:
  //    XH0 ... total cross section for hard emission (cm**2).
  //    XH1 ... stopping cross section for hard emission (eV*cm**2).
  //    XH2 ... straggling cross section for hard emission (eV**2*cm**2).
  //    XS1 ... stopping cross section for soft emission (eV*cm**2).
  //    XS2 ... straggling cross section for soft emission (eV**2*cm**2).


  const double TREV = 2.0*constants::REV;

  double X[constants::NBE], Y[constants::NBE];
  
  double XEL = log(E);
  if(XEL < DLEMP1){ XEL = DLEMP1;}
  
  double XE = 1.0+(XEL-DLEMP1)*DLFC;
  int KE = (int)XE;
  double XEK = XE-(double)KE;
  //  ****  Global x-section factor.
  double FACT = mat.ZBR2*(pow(E+constants::REV,2)/(E*(E+TREV)))*1.0E-27;

  //  ****  Moments of the scaled bremss x-section.

  double WCRE = WCRM/E;
  for(unsigned int IW = 0; IW < constants::NBW; IW++)
  {
    X[IW] = WB[IW];
    Y[IW] = mat.P0[KE-1][IW];
  }
  double XH0A = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,-1)-RLMOM(X,Y,WCRE,constants::NBW,-1);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XS1A = RLMOM(X,Y,WCRE,constants::NBW,0);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XS2A = RLMOM(X,Y,WCRE,constants::NBW,1);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XH1A = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,0)-XS1A;
  if(penGetError() != PEN_SUCCESS){ return;}
  double XH2A = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,1)-XS2A;
  if(penGetError() != PEN_SUCCESS){ return;}
  for(unsigned int IW = 0; IW < constants::NBW; IW++)
  {
    if(KE+1 < (int)constants::NEGP){ Y[IW] = mat.P0[KE][IW];}
    else{ Y[IW] = mat.P0[constants::NEGP-1][IW];}
  }
  double XH0B = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,-1)-RLMOM(X,Y,WCRE,constants::NBW,-1);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XS1B = RLMOM(X,Y,WCRE,constants::NBW,0);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XS2B = RLMOM(X,Y,WCRE,constants::NBW,1);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XH1B = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,0)-XS1B;
  if(penGetError() != PEN_SUCCESS){ return;}
  double XH2B = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,1)-XS2B;
  if(penGetError() != PEN_SUCCESS){ return;}

  XH0 = ((1.0-XEK)*XH0A+XEK*XH0B)*FACT;
  XS1 = ((1.0-XEK)*XS1A+XEK*XS1B)*FACT*E;
  XH1 = ((1.0-XEK)*XH1A+XEK*XH1B)*FACT*E;
  XS2 = ((1.0-XEK)*XS2A+XEK*XS2B)*FACT*E*E;
  XH2 = ((1.0-XEK)*XH2A+XEK*XH2B)*FACT*E*E;
}

//  *********************************************************************
//                       SUBROUTINE PBRaT
//  *********************************************************************
void PBRaT(const pen_material& mat, const double &E, const double &WCRM, double &XH0, double &XH1, double &XH2, double &XS1, double &XS2, const double DLEMP1, const double DLFC, const double* WB)
{
  //  Integrated cross sections for bremss emission by positrons of energy
  //  E in material M, restricted to energy losses larger than and less
  //  than the cutoff energy WCRM.

  //  Output arguments:
  //    XH0 ... total cross section for hard emission (cm**2).
  //    XH1 ... stopping cross section for hard emission (eV*cm**2).
  //    XH2 ... straggling cross section for hard emission (eV**2*cm**2).
  //    XS1 ... stopping cross section for soft emission (eV*cm**2).
  //    XS2 ... straggling cross section for soft emission (eV**2*cm**2).


  const double TREV = 2.0*constants::REV;

  double X[constants::NBE], Y[constants::NBE];
  
  double XEL = log(E);
  if(XEL < DLEMP1){ XEL = DLEMP1;}
  
  double XE = 1.0+(XEL-DLEMP1)*DLFC;
  int KE = (int)XE;
  double XEK = XE-(double)KE;
  //  ****  Global x-section factor.
  double FACT = mat.ZBR2*(pow(E+constants::REV,2)/(E*(E+TREV)))*1.0E-27;
  //  ****  Positron correction factor.
  double T = log(1.0+1.0E6*E/(constants::REV*mat.ZBR2));
  double FPOS = 1.0-exp(-T*(1.2359E-1-T*(6.1274E-2-T*(3.1516E-2-T*(7.7446E-3-T*(1.0595E-3-T*(7.0568E-5-T*1.8080E-6)))))));
  FACT = FACT*FPOS;

  //  ****  Moments of the scaled bremss x-section.

  double WCRE = WCRM/E;
  for(unsigned int IW = 0; IW < constants::NBW; IW++)
  {
    X[IW] = WB[IW];
    Y[IW] = mat.P0[KE-1][IW];
  }
  double XH0A = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,-1)-RLMOM(X,Y,WCRE,constants::NBW,-1);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XS1A = RLMOM(X,Y,WCRE,constants::NBW,0);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XS2A = RLMOM(X,Y,WCRE,constants::NBW,1);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XH1A = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,0)-XS1A;
  if(penGetError() != PEN_SUCCESS){ return;}
  double XH2A = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,1)-XS2A;
  if(penGetError() != PEN_SUCCESS){ return;}
  for(unsigned int IW = 0; IW < constants::NBW; IW++)
  {
    if(KE+1 < (int)constants::NEGP){ Y[IW] = mat.P0[KE][IW];}
    else{ Y[IW] = mat.P0[constants::NEGP-1][IW];}
  }
  double XH0B = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,-1)-RLMOM(X,Y,WCRE,constants::NBW,-1);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XS1B = RLMOM(X,Y,WCRE,constants::NBW,0);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XS2B = RLMOM(X,Y,WCRE,constants::NBW,1);
  if(penGetError() != PEN_SUCCESS){ return;}
  double XH1B = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,0)-XS1B;
  if(penGetError() != PEN_SUCCESS){ return;}
  double XH2B = RLMOM(X,Y,X[constants::NBW-1],constants::NBW,1)-XS2B;
  if(penGetError() != PEN_SUCCESS){ return;}
  
  XH0 = ((1.0-XEK)*XH0A+XEK*XH0B)*FACT;
  XS1 = ((1.0-XEK)*XS1A+XEK*XS1B)*FACT*E;
  XH1 = ((1.0-XEK)*XH1A+XEK*XH1B)*FACT*E;
  XS2 = ((1.0-XEK)*XS2A+XEK*XS2B)*FACT*E*E;
  XH2 = ((1.0-XEK)*XH2A+XEK*XH2B)*FACT*E*E;
}

//  *********************************************************************
//                       SUBROUTINE PANaT
//  *********************************************************************
void PANaT(const double &E, double &TXS)
{
  //  Total cross section (per electron) for annihilation of positrons with
  //  kinetic energy E. Computed from Heitler's dcs formula for annihila-
  //  tion with free electrons at rest.

  //  Output argument:
  //    XST ... total annihilation cross section (cm**2).


  const double PIELR2 = constants::PI*constants::ELRAD*constants::ELRAD;

  double GAM;
  if(E < 1.0){ GAM = 1.0+1.0/constants::REV;}
  else{ GAM = 1.0+E/constants::REV;}
  
  double GAM2 = GAM*GAM;
  double F2 = GAM2-1.0;
  double F1 = sqrt(F2);
  TXS = PIELR2*((GAM2+4.0*GAM+1.0)*log(GAM+F1)/F2-(GAM+3.0)/F1)/(GAM+1.0);
}

//  *********************************************************************
//                       SUBROUTINE EELaR
//  *********************************************************************
void EELaR(pen_material& mat, FILE* IRD, FILE* IWR, int INFO, CEEL00& eel00, CEINTF& eintf, const double* DLEMP, const double* ET)
{
  //  This subroutine reads elastic cross sections for electrons and posi-
  //  trons in material M from the material data file. It also initialises
  //  the algorithm for simulation of elastic scattering of electrons and
  //  positrons.



  //Create pointers to structure data
  double* EJT = eel00.EJT;
  
  double* XE0 = eel00.XE0;
  double* XE1 = eel00.XE1;
  double* XE2 = eel00.XE2;

  double* XP0 = eel00.XP0;
  double* XP1 = eel00.XP1;
  double* XP2 = eel00.XP2;

  double* EJTL = eel00.EJTL;
  double* FJL  = eel00.FJL;
  
  double* A  = eel00.A;
  double* B  = eel00.B;
  double* C  = eel00.C;
  double* D  = eel00.D;

  double* T1E0 = eel00.T1E0;
  double* T2E0 = eel00.T2E0;

  double* T1P0 = eel00.T1E0;
  double* T2P0 = eel00.T2E0;
  
  double* T1EI = eintf.T1EI;
  double* T2EI = eintf.T2EI;

  double* T1PI = eintf.T1EI;
  double* T2PI = eintf.T2EI;
  
  int NDATA;

  //  ****  Reading the input cross section table.
  
  fscanf(IRD, "%*59c%d%*[^\n]", &NDATA);
  getc(IRD);
  
  if(INFO >= 2)
  {
    fprintf(IWR, "\n *** Electron and positron elastic cross sections,  NDATA =%4d\n", NDATA);
    fprintf(IWR, "\n  Energy       CS0,e-      CS1,e-      CS2,e-      CS0,e+      CS1,e+      CS2,e+\n   (eV)        (cm**2)     (cm**2)     (cm**2)     (cm**2)     (cm**2)     (cm**2)\n ------------------------------------------------------------------------------------\n");
  }
  for(int I = 0; I < NDATA; I++)
  {
    fscanf(IRD, "%lf %lf %lf %lf %lf %lf %lf%*[^\n]", &EJT[I], &XE0[I], &XE1[I], &XE2[I], &XP0[I], &XP1[I], &XP2[I]);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR,"%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n", EJT[I], XE0[I], XE1[I], XE2[I], XP0[I], XP1[I], XP2[I]);}
    EJTL[I] = log(EJT[I]);
  }

  //  ****  Elastic scattering of electrons.

  double EC;
  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XE0[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(penGetError() != PEN_SUCCESS){ return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    unsigned int J;
    EC = DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XE0[I] = exp(A[J]+EC*(B[J]+EC*(C[J]+EC*D[J])));
  }
  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XE1[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(penGetError() != PEN_SUCCESS){ return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    unsigned int J;
    EC = DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XE1[I] = exp(A[J]+EC*(B[J]+EC*(C[J]+EC*D[J])));
  }

  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XE2[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(penGetError() != PEN_SUCCESS){ return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    unsigned int J;
    EC = DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XE2[I] = exp(A[J]+EC*(B[J]+EC*(C[J]+EC*D[J])));
  }

  double XS0, XS1, XS2, XS0H, XS1S, XS2S;
  double AAA, BBB;
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    XS0 = XE0[I];
    XS1 = XE1[I];
    XS2 = XE2[I];
    double FPEL = 1.0/(XS0*mat.VMOL);
    double FPT1 = 1.0/(XS1*mat.VMOL);
    double FPST = ET[I]/(mat.CSTPE[I]+mat.RSTPE[I]);
    
    XS0H = mat.C1*FPT1;
    if(XS0H > mat.C2*FPST){ XS0H = mat.C2*FPST;}
    if(XS0H < FPEL){ XS0H = FPEL;}
    XS0H = 1.0/(mat.VMOL*XS0H);
    
    double RNDC;
    EELa0(XS0,XS1,XS2,XS0H,AAA,BBB,RNDC,XS1S,XS2S);
    mat.SEHEL[I] = XS0H*mat.VMOL;
    
    mat.RNDCE[I] = RNDC;
    mat.AE[I] = AAA;
    mat.BE[I] = BBB;
    T1E0[I] = XS1S;
    mat.T1E[I] = T1EI[I]+XS1S*mat.VMOL;
    T2E0[I] = XS2S;
    mat.T2E[I] = T2EI[I]+XS2S*mat.VMOL;
  }

  //  ****  Print electron elastic scattering tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of electrons\n");
    fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)        A           B           RNDC\n ------------------------------------------------------------------------------------------\n");
  }
  double FP0, FP1;
  double HMFP;
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    FP0 = mat.RHO/(XE0[I]*mat.VMOL);
    FP1 = mat.RHO/(XE1[I]*mat.VMOL);
    HMFP = mat.RHO/mat.SEHEL[I];
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", ET[I], FP0, FP1, HMFP, mat.AE[I], mat.BE[I], mat.RNDCE[I]);}
    mat.SEHEL[I] = log(mat.SEHEL[I]);
    mat.AE[I] = log(mat.AE[I]);
      //  ****  Soft scattering events are switched off when T1E is too small.
    if(mat.T1E[I] > 1.0E-6*XE1[I]*mat.VMOL)
    {
      if(mat.T1E[I] < 1.0E-35){ mat.T1E[I] = 1.0E-35;}
      mat.T1E[I] = log(mat.T1E[I]);
      if(mat.T2E[I] < 1.0E-35){ mat.T2E[I] = 1.0E-35;}
      mat.T2E[I] = log(mat.T2E[I]);
    }
    else
    {
      mat.T1E[I] = -100.0;
      mat.T2E[I] = -100.0;
    }
  }

  //  ****  Elastic scattering of positrons.

  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XP0[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(penGetError() != PEN_SUCCESS){ return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    unsigned int J;
    EC = DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XP0[I] = exp(A[J]+EC*(B[J]+EC*(C[J]+EC*D[J])));
  }

  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XP1[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(penGetError() != PEN_SUCCESS){ return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    unsigned int J;
    EC = DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XP1[I] = exp(A[J]+EC*(B[J]+EC*(C[J]+EC*D[J])));
  }

  for(int I = 0; I < NDATA; I++)
  {
    FJL[I] = log(XP2[I]);
  }
  SPLINE(EJTL,FJL,A,B,C,D,0.0,0.0,NDATA);
  if(penGetError() != PEN_SUCCESS){ return;}
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    unsigned int J;
    EC = DLEMP[I];
    FINDI(EJTL,EC,NDATA,J);
    XP2[I] = exp(A[J]+EC*(B[J]+EC*(C[J]+EC*D[J])));
  }

  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    XS0 = XP0[I];
    XS1 = XP1[I];
    XS2 = XP2[I];
    double FPEL = 1.0/(XS0*mat.VMOL);
    double FPT1 = 1.0/(XS1*mat.VMOL);
    double FPST = ET[I]/(mat.CSTPP[I]+mat.RSTPP[I]);

    XS0H = mat.C1*FPT1;
    if(XS0H > mat.C2*FPST){ XS0H = mat.C2*FPST;}
    if(XS0H < FPEL){ XS0H = 1.0/(mat.VMOL*FPEL);}
    else{ XS0H = 1.0/(mat.VMOL*XS0H);}
      
    double RNDC;
    EELa0(XS0,XS1,XS2,XS0H,AAA,BBB,RNDC,XS1S,XS2S);
    mat.SPHEL[I] = XS0H*mat.VMOL;
    mat.RNDCP[I] = RNDC;
    mat.AP[I] = AAA;
    mat.BP[I] = BBB;
    T1P0[I] = XS1S;
    mat.T1P[I] = T1PI[I]+XS1S*mat.VMOL;
    T2P0[I] = XS2S;
    mat.T2P[I] = T2PI[I]+XS2S*mat.VMOL;
  }

  //  ****  Print positron elastic scattering tables.

  if(INFO >= 3)
  {
    fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of positrons\n");
    fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)        A           B           RNDC\n ------------------------------------------------------------------------------------------\n");
  }
  for(unsigned int I = 0; I < constants::NEGP; I++)
  {
    FP0 = mat.RHO/(XP0[I]*mat.VMOL);
    FP1 = mat.RHO/(XP1[I]*mat.VMOL);
    HMFP = mat.RHO/mat.SPHEL[I];
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n", ET[I], FP0, FP1, HMFP, mat.AP[I], mat.BP[I], mat.RNDCP[I]);}
    mat.SPHEL[I] = log(mat.SPHEL[I]);
    mat.AP[I] = log(mat.AP[I]);
      
      //  ****  Soft scattering events are switched off when T1P is too small.
    if(mat.T1P[I] > 1.0E-6*XP1[I]*mat.VMOL)
    {
      if(mat.T1P[I] < 1.0E-35){ mat.T1P[I] = log(1.0E-35);}
      else{ mat.T1P[I] = log(mat.T1P[I]);}
    
      if(mat.T2P[I] < 1.0E-35){ mat.T2P[I] = log(1.0E-35);}
      else{ mat.T2P[I] = log(mat.T2P[I]);}
    }
    else
    {
      mat.T1P[I] = -100.0;
      mat.T2P[I] = -100.0;
    }
  }
}

//  *********************************************************************
//                       SUBROUTINE EELdR
//  *********************************************************************
void EELdR(pen_material& mat, FILE* IRD, FILE* IWR, int &INFO, CDCSEP& dcsep, CRITA& rita, CEINTF& eintf, double* ET)
{
  //     This subroutine reads elastic cross sections for electrons and
  //  positrons in material M from the elastic scattering database. It also
  //  initialises the algorithm for simulation of elastic collisions of
  //  electrons and positrons.


  const unsigned int NP_P = mat.NP_P;
  const unsigned int NE = dcsep.NE;
  const unsigned int NA = dcsep.NA;
  
  //Create pointers to structure variables
  double* ETS   = dcsep.ETS;
  double* ETL   = dcsep.ETL;
  double* TH    = dcsep.TH;
  double* THR   = dcsep.THR;
  double* XMU   = dcsep.XMU;
  double* XMUL  = dcsep.XMUL;
  double* ECS   = dcsep.ECS;
  double* ETCS1 = dcsep.ETCS1;
  double* ETCS2 = dcsep.ETCS2;
  double* DCSI  = dcsep.DCSI;
  double* PCS   = dcsep.PCS;
  double* PTCS1 = dcsep.PTCS1;
  double* PTCS2 = dcsep.PTCS2;

  double* XTI  = rita.XT; 
  double* PACI = rita.PAC; 
  double* AI   = rita.A; 
  double* BI   = rita.B;
  
  int* ITLI = rita.IL; 
  int* ITUI = rita.IU; 

  double XE0[constants::NEGP];
  double XE1[constants::NEGP];
  //double XE2[constants::NEGP];   //For debug
  //double T1E0[constants::NEGP];  //For debug
  //double T2E0[constants::NEGP];  //For debug

  double XP0[constants::NEGP];
  double XP1[constants::NEGP];
  //double XP2[constants::NEGP];   //For debug
  //double T1P0[constants::NEGP];  //For debug
  //double T2P0[constants::NEGP];  //For debug

  double* T1EI = eintf.T1EI;
  double* T2EI = eintf.T2EI;
  double* T1PI = eintf.T1PI;
  double* T2PI = eintf.T2PI;
  
  char CTEXT[51];

  double EGRD[16] = {1.0, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 1.00E1};

  //  ****  Energy mesh points (in eV).
  
  unsigned int IE = 0;
  int IGRID = 10;
  double FGRID = 10.0;
  bool Eixir = false;
  double EV;
  while(!Eixir)
  {
    Eixir = true;
    IGRID = IGRID+1;
    EV = EGRD[IGRID-1]*FGRID;
    if(IGRID == 16)
    {
      IGRID = 1;
      FGRID = 10.0*FGRID;
    }
    IE = IE+1;
    ETS[IE-1] = EV;
    ETL[IE-1] = log(ETS[IE-1]);
    if(IE < NE){ Eixir = false; continue;}
  }

  //  ****  Angular grid (TH in deg, XMU=(1.0D0-COS(TH))/2).

  unsigned int I = 0;
  TH[I] = 0.0;
  THR[I] = TH[I]*constants::PI/180.0;
  XMU[I] = (1.0-cos(THR[I]))/2.0;
  XMUL[I] = log(1.0E-35);
  I = 1;
  TH[I] = 1.0E-4;
  THR[I] = TH[I]*constants::PI/180.0;
  XMU[I] = (1.0-cos(THR[I]))/2.0;

  XMUL[I] = XMU[I];
  if(XMUL[I] < 1.0E-35){ XMUL[I] = 1.0E-35;}
  XMUL[I]=log(XMUL[I]);

  Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    I = I+1;
    if(TH[I-1] < 0.9999E-3)
    {
      TH[I] = TH[I-1]+2.5E-5;
    }
    else if(TH[I-1] < 0.9999E-2)
    {
      TH[I] = TH[I-1]+2.5E-4;
    }
    else if(TH[I-1] < 0.9999E-1)
    {
      TH[I] = TH[I-1]+2.5E-3;
    }
    else if(TH[I-1] < 0.9999)
    {
      TH[I] = TH[I-1]+2.5E-2;
    }
    else if(TH[I-1] < 0.9999E+1)
    {
      TH[I] = TH[I-1]+1.0E-1;
    }
    else if(TH[I-1] < 2.4999E+1)
    {
      TH[I] = TH[I-1]+2.5E-1;
    }
    else
    {
      TH[I] = TH[I-1]+5.0E-1;
    }
    THR[I] = TH[I]*constants::PI/180.0;
    XMU[I] = (1.0-cos(THR[I]))/2.0;
    if(XMU[I] < 1.0E-35){ XMU[I] = 1.0E-35;}
    if( XMU[I] > 1.0E-35){ XMUL[I] = log(XMU[I]);}
    else{ XMUL[I] = log(1.0E-35);}
      
    if((unsigned)(I+1) < NA){ Eixir = false; continue;}
  }

  //  ****  Read elastic DCS tables.

  if(INFO > 3){ fprintf(IWR, "\n *** Electron elastic differential cross sections\n");}
  int IELEC = -1;
  fscanf(IRD, "%50[^\n]%*[^\n]", CTEXT);
  getc(IRD);
  for(IE = 0; IE < NE; IE++)
  {
    double ETSIE;
    fscanf(IRD, "%d %lf %lf %lf %lf%*[^\n]", &IELEC, &ETSIE, &ECS[IE], &ETCS1[IE], &ETCS2[IE]);
    getc(IRD);
    if(INFO > 3){ fprintf(IWR, "%3d%10.3E%12.5E%12.5E%12.5E\n", IELEC, ETS[IE], ECS[IE], ETCS1[IE], ETCS2[IE]);}
    fflush(IWR);
    if(IELEC != -1 || fabs(ETSIE-ETS[IE]) > 0.1)
    {
      fprintf(IWR, "\n Error reading electron elastic DCS data.\n");
      penError(ERR_EELdR_em_EL_DCS); return;
    }
    for(unsigned int K = 0; K < NA; K++)
    {
      fscanf(IRD, "%lf", &dcsep.EDCS[IE][K]);
    }
    fscanf(IRD,"%*[^\n]");
    getc(IRD);

    if(INFO > 3)
    {
      for(unsigned int K = 0; K < NA; K++)
      {
        if((K+1)%10==0 || K==NA-1){fprintf(IWR, " %11.5E\n", dcsep.EDCS[IE][K]);}
        else{fprintf(IWR, " %11.5E", dcsep.EDCS[IE][K]);}
      }
    }
  //  ****  Consistency test.
    for(unsigned int K = 0; K < NA; K++)
    {
      DCSI[K] = dcsep.EDCS[IE][K];
    }
    double ECS0 = 4.0*constants::PI*RMOMX(XMU,DCSI,0.0,1.0,NA,0);
    if(penGetError() != PEN_SUCCESS){ return;}
    double ECS1 = 4.0*constants::PI*RMOMX(XMU,DCSI,0.0,1.0,NA,1);
    double ECS2 = 4.0*constants::PI*RMOMX(XMU,DCSI,0.0,1.0,NA,2);
    double TCS1 = 2.0*ECS1;
    double TCS2 = 6.0*(ECS1-ECS2);
    double TS0 = (ECS0-ECS[IE])/ECS[IE];
    double TS1 = (TCS1-ETCS1[IE])/ETCS1[IE];
    double TS2 = (TCS2-ETCS2[IE])/ETCS2[IE];
    double TSTE = fabs(TS0);
    if(TSTE < fabs(TS1)){ TSTE = fabs(TS1);}
    if(TSTE < fabs(TS2)){ TSTE = fabs(TS2);}
    if(TSTE > 1.0E-4)
    {
      fprintf(IWR, "\n E=%12.5E\n", ETS[IE]);
      fprintf(IWR, " Electron cross section data are corrupt.\n");
      penError(ERR_EELdR_em_CS_DB); return;
    }
  }

  if(INFO > 3){ fprintf(IWR, "\n *** Positron elastic differential cross sections\n");}
  IELEC = +1;
  fscanf(IRD, "%50[^\n]%*[^\n]", CTEXT);
  getc(IRD);
  for(IE = 0; IE < NE; IE++)
  {
    double ETSIE;
    fscanf(IRD, "%d %lf %lf %lf %lf%*[^\n]", &IELEC, &ETSIE, &PCS[IE], &PTCS1[IE], &PTCS2[IE]);
    getc(IRD);
    if(INFO > 3){ fprintf(IWR, "%3d%10.3E%12.5E%12.5E%12.5E\n",IELEC,ETS[IE],PCS[IE],PTCS1[IE],PTCS2[IE]);}
    if(IELEC != +1 || fabs(ETSIE-ETS[IE]) > 0.1)
    {
      fprintf(IWR, "\n Error reading positron elastic DCS data.\n");
      penError(ERR_EELdR_ep_EL_DCS); return;
    }
    for(unsigned int K = 0; K < NA; K++)
    {
      fscanf(IRD, "%lf", &dcsep.PDCS[IE][K]);
    }
    fscanf(IRD,"%*[^\n]");
    getc(IRD);
      
    if(INFO > 3)
    {
      for(unsigned int K = 0; K < NA; K++)
      {
        if((K+1)%10==0 || K==NA-1){fprintf(IWR, " %11.5E\n", dcsep.PDCS[IE][K]);}
        else{fprintf(IWR, " %11.5E", dcsep.PDCS[IE][K]);}
      }
    }
      //  ****  Consistency test.
    for(unsigned int K = 0; K < NA; K++)
    {
      DCSI[K] = dcsep.PDCS[IE][K];
    }
    double ECS0 = 4.0*constants::PI*RMOMX(XMU,DCSI,0.0,1.0,NA,0);
    if(penGetError() != PEN_SUCCESS){ return;}
    double ECS1 = 4.0*constants::PI*RMOMX(XMU,DCSI,0.0,1.0,NA,1);
    double ECS2 = 4.0*constants::PI*RMOMX(XMU,DCSI,0.0,1.0,NA,2);
    double TCS1 = 2.0*ECS1;
    double TCS2 = 6.0*(ECS1-ECS2);
    double TS0 = (ECS0-PCS[IE])/PCS[IE];
    double TS1 = (TCS1-PTCS1[IE])/PTCS1[IE];
    double TS2 = (TCS2-PTCS2[IE])/PTCS2[IE];

    double TSTE = fabs(TS0);
    if(TSTE < fabs(TS1)){ TSTE = fabs(TS1);}
    if(TSTE < fabs(TS2)){ TSTE = fabs(TS2);}

    if(TSTE > 1.0E-4)
    {
      fprintf(IWR, "\n E=%11.5E\n", ETS[IE]);
      fprintf(IWR, " Positron cross section data are corrupt.\n");
      penError(ERR_EELdR_ep_CS_DB); return;
    }
  }
  int NPP = NP_P;
  int NU = NPP/4;

  //  ************  Electrons.

  unsigned int IEME = 0;
  for(unsigned int KE = 0; KE < constants::NEGP; KE++)
  {
    double ERRM;
    if(ET[KE] > 0.999999E8){ break;}
    DCSEL0(ET[KE],-1,dcsep);
    if(penGetError() != PEN_SUCCESS){ return;}
    RITAI0(DCSEL,&dcsep,0.0,1.0,NPP,NU,ERRM,0,rita);
    for(I = 0; I < NP_P; I++)
    {
      mat.XSE[I][KE] = XTI[I];
      mat.PSE[I][KE] = PACI[I];
      mat.ASE[I][KE] = AI[I];
      mat.BSE[I][KE] = BI[I];
      mat.ITLE[I][KE] = ITLI[I];
      mat.ITUE[I][KE] = ITUI[I];
    }
    double XM0A, XM1, XM2;
    RITAM(0.0,1.0,XM0A,XM1,XM2,rita);
    double ECS0 = dcsep.CSI;
    double ECS1 = dcsep.CSI*XM1/XM0A;
    double ECS2 = dcsep.CSI*XM2/XM0A;
    XE0[KE] = ECS0;
    XE1[KE] = 2.0*ECS1;
    //XE2[KE] = 6.0*(ECS1-ECS2);

    double FPEL = 1.0/(XE0[KE]*mat.VMOL);
    double FPT1 = 1.0/(XE1[KE]*mat.VMOL);
    double FPST = ET[KE]/(mat.CSTPE[KE]+mat.RSTPE[KE]);

    double XS0H = mat.C1*FPT1;
    if(XS0H > mat.C2*FPST){ XS0H = mat.C2*FPST;}
    if(XS0H < FPEL){ XS0H = FPEL;}

    XS0H = 1.0/(mat.VMOL*XS0H);

    double RNDC = 1.0-XS0H/XE0[KE];
    if(RNDC < 1.0E-10){ RNDC = 1.0E-10;}
      
    if(RNDC < 1.0E-6){ RNDC = 0.0;}
    mat.RNDCEd[KE] = RNDC;
      
    double RU = RNDC;
    I = 1;
    int J = NP_P;
    bool Eixir2 = false;
    while(!Eixir2)
    {
      Eixir2 = true;
      int K = (I+J)/2;   
      if(RU > mat.PSE[K-1][KE])
      {
        I = K;
      }
      else
      {
        J = K;
      }
      if(J-I > 1){ Eixir2 = false; continue;}
    }

    double RR = RU-mat.PSE[I-1][KE];
    double DPRO = mat.PSE[I+1-1][KE]-mat.PSE[I-1][KE];
    double RMUC;
    if(DPRO < 1.0E-10)
    {
      RMUC = mat.XSE[I-1][KE];
    }
    else
    {
      double CI = (1.0 + mat.ASE[I-1][KE] + mat.BSE[I-1][KE])*DPRO;
      RMUC = mat.XSE[I-1][KE] + (CI*RR/(pow(DPRO,2)+(DPRO*mat.ASE[I-1][KE] + mat.BSE[I-1][KE]*RR)*RR))*(mat.XSE[I+1-1][KE] - mat.XSE[I-1][KE]);
    }

      //  ****  Moments of the PDF on the restricted interval (0,RMUC).
      //        Total and transport cross sections for soft interactions.
      
    double XM0; 
    RITAM(0.0,RMUC,XM0,XM1,XM2,rita);
    ECS1 = dcsep.CSI*XM1/XM0A;
    ECS2 = dcsep.CSI*XM2/XM0A;
    double TCS1 = 2.0*ECS1;
    double TCS2 = 6.0*(ECS1-ECS2);
    mat.SEHEL[KE] = XS0H*mat.VMOL;
    //T1E0[KE] = TCS1;
    mat.T1E[KE] = T1EI[KE]+TCS1*mat.VMOL;
    //T2E0[KE] = TCS2;
    mat.T2E[KE] = T2EI[KE]+TCS2*mat.VMOL;
    IEME = KE+1;
  }

  mat.EELMAX = ET[IEME-1]-1.0;
  if(mat.EELMAX > 0.999999E8){ mat.EELMAX = 0.999999E8;}
  
  //  ****  Print electron elastic scattering tables.

  if(INFO >= 3){ fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of electrons (ELSEPA database)\n");}
  if(INFO >= 3){ fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)\n --------------------------------------------------\n");}
  for(I = 0; I < IEME; I++)
  {
    double FP0 = mat.RHO/(XE0[I]*mat.VMOL);
    double FP1 = mat.RHO/(XE1[I]*mat.VMOL);
    double HMFP = mat.RHO/mat.SEHEL[I];
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E\n", ET[I], FP0, FP1, HMFP);}
    mat.SEHEL[I] = log(mat.SEHEL[I]);
      //  ****  Soft scattering events are switched off when T1E is too small.
    if(mat.T1E[I] > 1.0E-6*XE1[I]*mat.VMOL)
    {
      if(mat.T1E[I] < 1.0E-35){ mat.T1E[I] = log(1.0E-35);}
      else{ mat.T1E[I] = log(mat.T1E[I]);}

      if(mat.T2E[I] < 1.0E-35){ mat.T2E[I] = log(1.0E-35);}
      else{ mat.T2E[I] = log(mat.T2E[I]);}
    }
    else
    {
      mat.T1E[I] = -100.0;
      mat.T2E[I] = -100.0;
    }
  }

  //  ************  Positrons.

  unsigned int IEMP = 0;
  for(unsigned int KE = 0; KE < constants::NEGP; KE++)
  {
    double ERRM;
    if(ET[KE] > 0.999999E8){break;}
    DCSEL0(ET[KE],+1,dcsep);
    if(penGetError() != PEN_SUCCESS){ return;}
    RITAI0(DCSEL,&dcsep,0.0,1.0,NPP,NU,ERRM,0,rita);
    for(I = 0; I < NP_P; I++)
    {
      mat.XSP[I][KE] = XTI[I];
      mat.PSP[I][KE] = PACI[I];
      mat.ASP[I][KE] = AI[I];
      mat.BSP[I][KE] = BI[I];
      mat.ITLP[I][KE] = ITLI[I];
      mat.ITUP[I][KE] = ITUI[I];
    }
    double XM0A,XM1,XM2;
    RITAM(0.0,1.0,XM0A,XM1,XM2,rita);
    double ECS0 = dcsep.CSI;
    double ECS1 = dcsep.CSI*XM1/XM0A;
    double ECS2 = dcsep.CSI*XM2/XM0A;
    XP0[KE] = ECS0;
    XP1[KE] = 2.0*ECS1;
    //XP2[KE] = 6.0*(ECS1-ECS2);

    double FPEL = 1.0/(XP0[KE]*mat.VMOL);
    double FPT1 = 1.0/(XP1[KE]*mat.VMOL);
    double FPST = ET[KE]/(mat.CSTPP[KE]+mat.RSTPP[KE]);

    double XS0H = mat.C1*FPT1;
    if(XS0H > mat.C2*FPST){ XS0H = mat.C2*FPST;}
    if(XS0H < FPEL){ XS0H = FPEL;}
    XS0H=1.0/(mat.VMOL*XS0H);

    double RNDC = 1.0-XS0H/XP0[KE];
    if(RNDC < 1.0E-10){ RNDC = 1.0E-10;}

    if(RNDC < 1.0E-6){ RNDC = 0.0;}
    mat.RNDCPd[KE] = RNDC;

    double RU = RNDC;
    I = 1;
    int J = NP_P;
    bool Eixir2 = false;
    while(!Eixir2)
    {
      Eixir2 = true;
      int K = (I+J)/2;
      if(RU > mat.PSP[K-1][KE])
      {
        I = K;
      }
      else
      {
        J = K;
      }
      if(J-I > 1){ Eixir2 = false; continue;}
    }

    double RMUC, CI;
    double RR = RU-mat.PSP[I-1][KE];
    double DPRO = mat.PSP[I][KE] - mat.PSP[I-1][KE];
    if(DPRO < 1.0E-10)
    {
      RMUC = mat.XSP[I-1][KE];
    }
    else
    {
      CI = (1.0 + mat.ASP[I-1][KE] + mat.BSP[I-1][KE])*DPRO;
      RMUC = mat.XSP[I-1][KE] + (CI*RR/(pow(DPRO,2)+(DPRO*mat.ASP[I-1][KE] + mat.BSP[I-1][KE]*RR)*RR))*(mat.XSP[I][KE] - mat.XSP[I-1][KE]);
    }

      //  ****  Moments of the PDF on the restricted interval (0,RMUC).
      //        Total and transport cross sections for soft interactions.

    double XM0;
    RITAM(0.0,RMUC,XM0,XM1,XM2,rita);
    ECS1 = dcsep.CSI*XM1/XM0A;
    ECS2 = dcsep.CSI*XM2/XM0A;
    double TCS1 = 2.0*ECS1;
    double TCS2 = 6.0*(ECS1-ECS2);
    mat.SPHEL[KE] = XS0H*mat.VMOL;
    //T1P0[KE] = TCS1;
    mat.T1P[KE] = T1PI[KE]+TCS1*mat.VMOL;
    //T2P0[KE] = TCS2;
    mat.T2P[KE] = T2PI[KE]+TCS2*mat.VMOL;
    IEMP = KE+1;
  }

  mat.PELMAX = ET[IEMP-1]-1.0;
  if(mat.PELMAX > 0.999999E8){ mat.PELMAX = 0.999999E8;}

  //  ****  Print positron elastic scattering tables.

  if(INFO >= 3){ fprintf(IWR, "\n PENELOPE >>>  Elastic scattering of positrons (ELSEPA database)\n");}
  if(INFO >= 3){ fprintf(IWR, "\n   E (eV)      MFP (mtu)   TMFP1 (mtu)  MFPh (mtu)\n --------------------------------------------------\n");}
  for(I = 0; I < IEMP; I++)
  {
    double FP0 = mat.RHO/(XP0[I]*mat.VMOL);
    double FP1 = mat.RHO/(XP1[I]*mat.VMOL);
    double HMFP = mat.RHO/mat.SPHEL[I];
    if(INFO >= 3){ fprintf(IWR, "%12.5E %12.5E %12.5E %12.5E\n", ET[I], FP0, FP1, HMFP);}
    mat.SPHEL[I] = log(mat.SPHEL[I]);
      //  ****  Soft scattering events are switched off when T1P is too small.
    if(mat.T1P[I] > 1.0E-6*XP1[I]*mat.VMOL)
    {
      if(mat.T1P[I] < 1.0E-35){ mat.T1P[I] = log(1.0E-35);}
      else{ mat.T1P[I] = log(mat.T1P[I]);}

      if(mat.T2P[I] < 1.0E-35){ mat.T2P[I] = log(1.0E-35);}
      else{ mat.T2P[I] = log(mat.T2P[I]);}
    }
    else
    {
      mat.T1P[I] = -100.0;
      mat.T2P[I] = -100.0;
    }
  } 
}

//  *********************************************************************
//                       SUBROUTINE GRAaR
//  *********************************************************************
int GRAaR(pen_logGrid& grid, pen_material& mat, CRITA& rita, FILE* IRD, FILE* IWR, int &INFO)
{
  //  This subroutine reads the squared molecular form factor and the DCS
  //  for Rayleigh scattering of photons in material M. These two functions
  //  are tabulated using the same grids for all materials.
  //     The random sampling of the scattering angle is performed using the
  //  RITA algorithm.


  const double RREV = 1.0/constants::REV;

  double Q[constants::NQ], F[constants::NQ], FFI[constants::NQ], ER[constants::NEX], AA[constants::NQ], BB[constants::NQ], CC[constants::NQ], DD[constants::NQ];
  const int NIP =51;
  double QI[51], FUN[51], SUM[51];

  int NQQ;
  
  fscanf(IRD, "%*32c%d%*8c%d%*[^\n]", &NQQ, &mat.NE);
  getc(IRD);
  if(INFO >= 2){ fprintf(IWR, "\n *** Rayleigh scattering.  NQ = %3d,  NE = %4d\n", constants::NQ, mat.NE);}
  for(unsigned int I = 0; I < constants::NQ; I++)
  {
    fscanf(IRD, "%lf %lf%*[^\n]", &Q[I], &FFI[I]);
    getc(IRD);
    if(I == 0){ mat.FF0 = pow(FFI[I],2);}
    F[I] = log(pow(FFI[I],2));
    mat.FF[I] = FFI[I];
  }
  for(int I = 0; I < mat.NE; I++)
  {
    fscanf(IRD, "%lf %lf%*[^\n]", &ER[I], &mat.XSRA[I]);
    getc(IRD);
    mat.XSRA[I] = log(mat.XSRA[I]*mat.VMOL);
  }
  double XS1;
  if(ER[mat.NE-1] < grid.EU)
  {
    XS1 = mat.XSRA[mat.NE-1]+(mat.XSRA[mat.NE-1]-mat.XSRA[mat.NE-1-1])*(log(grid.EU)-log(ER[mat.NE-1-1]))/(log(ER[mat.NE-1])-log(ER[mat.NE-1-1]));
    mat.XSRA[mat.NE-1] = XS1;
    ER[mat.NE-1] = grid.EU;
  }

  mat.QQM = pow(Q[constants::NQ-1],2);
  for(unsigned int I = 0; I < constants::NQ; I++)
    {
      mat.QQ[I] = log(pow(Q[I],2));
    }
  for(int I = 0; I < mat.NE; I++)
    {
      mat.ERA[I] = log(ER[I]);
    }
  for(unsigned int I = 0; I < constants::NEGP; I++)
    {
      unsigned int J;
      FINDI(mat.ERA,grid.DLEMP[I],mat.NE,J);
      mat.IED[I] = J;
    }
  for(unsigned int I = 0; I < constants::NEGP-1; I++)
    {
      mat.IEU[I] = mat.IED[I+1]+1;
      if(mat.IEU[I] >= mat.NE){ mat.IEU[I] = mat.NE-1;}  
    }
  mat.IEU[constants::NEGP-1] = mat.IED[constants::NEGP-1]+1;
  if(mat.IEU[constants::NEGP-1] >= mat.NE){ mat.IEU[constants::NEGP-1] = mat.NE-1;}

  SPLINE(mat.QQ,F,AA,BB,CC,DD,0.0,0.0,constants::NQ);
  if(penGetError() != PEN_SUCCESS){ return 0;}
  for(unsigned int I = 0; I < constants::NQ; I++)
  {
    mat.AR[I] = AA[I];
    mat.BR[I] = BB[I];
    mat.CR[I] = CC[I];
    mat.DR[I] = DD[I];
  }

  //  ****  Total cross section at the simulation grid points (slightly
  //        increased to simplify interpolation).

  double EE = grid.DLEMP[0];
  unsigned int J;
  FINDI(mat.ERA,EE,mat.NE,J);
  XS1 = mat.XSRA[J]+(mat.XSRA[J+1]-mat.XSRA[J])*(EE-mat.ERA[J])/(mat.ERA[J+1]-mat.ERA[J]);
  double XSMAX;
  for(unsigned int IE = 0; IE < constants::NEGP-1; IE++)
  {
    XSMAX = XS1;
    unsigned int J1 = J+1;
    EE = grid.DLEMP[IE+1];
    FINDI(mat.ERA,EE,mat.NE,J);
    XS1 = mat.XSRA[J]+(mat.XSRA[J+1]-mat.XSRA[J])*(EE-mat.ERA[J])/(mat.ERA[J+1]-mat.ERA[J]);
    if(XSMAX < XS1){ XSMAX = XS1;}

    if(J1 < J)
    {
      for(unsigned int I = J1-1; I < J; I++)
      {
        if(XSMAX < mat.XSRA[I]){ XSMAX = mat.XSRA[I];}            
      }
    }
    mat.SGRA[IE] = exp(XSMAX);
  }
  mat.SGRA[constants::NEGP-1] = mat.SGRA[constants::NEGP-2];

  if(INFO >= 2)
  {
    fprintf(IWR, "\n   Q/me*c     Form factor\n -------------------------\n");
    for(unsigned int I = 0; I < constants::NQ; I++)
    {
      fprintf(IWR, "%12.5E%13.5E\n", Q[I], FFI[I]);
    }
    fprintf(IWR, "\n   Energy       CS-Rayl\n    (eV)        (cm**2)\n -------------------------\n");
    for(int I = 0; I < mat.NE; I++)
    {
      fprintf(IWR, "%12.5E%13.5E\n", ER[I], exp(mat.XSRA[I])/mat.VMOL);
    }
  }
  //
  //  ****  Initialisation of the RITA algorithm for random sampling of the
  //  squared momentum transfer from the squared molecular form factor.

  double Q2MIN = 0.0;
  mat.Q2MAX = 0.0;
  int NPT = pen_material::NP_RSC;
  int NU = NPT/4;
  for(unsigned int I = 1; I < constants::NQ; I++)
  {
    if(GRAaF2(pow(Q[I],2),&mat) > 1.0E-35){ mat.Q2MAX = pow(Q[I-1],2);}
  }
  double ERRM;
  RITAI0(GRAaF2,&mat,Q2MIN,mat.Q2MAX,NPT,NU,ERRM,0,rita);
  int NPI = rita.NPM1+1;
  if(NPI != mat.NP_RSC)
  {
    printf("GRAaR. RITA initialisation error.\n");
    printf("The number of fixed grid points is %d\n", NPI);
    printf("The required number of grid points was %d\n", mat.NP_RSC);
    penError(ERR_GRAaR_INIT); return ERR_GRAaR_INIT;
  }
  if(ERRM > 1.0E-5)
  {
    printf("GRAaR. RITA interpolation error is too large.\n");
    printf("The interpolation error is %E\n", ERRM);
    penError(ERR_GRAaR_INTERPOLATION); return ERR_GRAaR_INTERPOLATION;
  }

  //  ****  Upper limit of the X2 interval for the PENELOPE grid energies.

  for(unsigned int IE = 0; IE < constants::NEGP; IE++)
  {
    double QM = 2.0*grid.ET[IE]*RREV;
    double Q2M = QM*QM;
    int I;
    if(Q2M > rita.XT[0])
    {
      if(Q2M < rita.XT[mat.NP_RSC-1])
      {
        I = 1;
        J = NPI;
        bool Eixir = false;
        while(!Eixir)
        {
          Eixir = true;
          int IT = (I+J)/2;
          if(Q2M > rita.XT[IT-1])
          {
            I = IT;
          }
          else
          {
            J = IT;
          }
          if(J-I > 1){ Eixir = false; continue;}
        }

        double Q1 = rita.XT[I-1];
        double Q2 = Q2M;
        double DQ = (Q2-Q1)/double(NIP-1);
        double TAU, CON1, CI, CON2;
        for(int K = 0; K < NIP; K++)
        {
          QI[K] = Q1+double(K+1-1)*DQ;
          TAU = (QI[K]-rita.XT[I-1])/(rita.XT[I]-rita.XT[I-1]);
          CON1 = 2.0*rita.B[I-1]*TAU;
          CI = 1.0+rita.A[I-1]+rita.B[I-1];
          CON2 = CI-rita.A[I-1]*TAU;
          double ETAP;
          if(fabs(CON1) > 1.0E-16*fabs(CON2))
          {
            ETAP = CON2*(1.0-sqrt(1.0-2.0*TAU*CON1/pow(CON2,2)))/CON1;
          }
          else
          {
            ETAP = TAU/CON2;
          }
          FUN[K] = rita.DPAC[I-1]*pow(1.0+(rita.A[I-1]+rita.B[I-1]*ETAP)*ETAP,2)/((1.0-rita.B[I-1]*ETAP*ETAP)*CI*(rita.XT[I]-rita.XT[I-1]));
        }
        SLAG6(DQ,FUN,SUM,NIP);
        if(penGetError() != PEN_SUCCESS){ return -1;}
        mat.PMAX[IE] = rita.PAC[I-1]+SUM[NIP-1];
      }
      else
      {
        mat.PMAX[IE] = 1.0;
      }
    }
    else
    {
      mat.PMAX[IE] = rita.PAC[0];
    }
  }

  for(unsigned int I = 0; I < mat.NP_RSC; I++)
  {
    mat.QRA[I] = rita.XT[I];
    mat.PRA[I] = rita.PAC[I];
    mat.DPRA[I]= rita.DPAC[I];
    mat.ARA[I] = rita.A[I];
    mat.BRA[I] = rita.B[I];
    mat.ITLRA[I] = rita.IL[I];
    mat.ITURA[I] = rita.IU[I];
  }
  return PEN_SUCCESS;
}

//  *********************************************************************
//                       SUBROUTINE GPHaR
//  *********************************************************************
void GPHaR(pen_material& mat, pen_elementDataBase& elemDB, FILE* IRD, FILE* IWR, int INFO, CGPH01& gph01, const double* DLEMP)
{
  //  This subroutine reads photoelectric cross sections of the elements in
  //  material M and prepares simulation tables.
  
  //  NOTE: The array SGPH(M,IE) defines a piecewise constant function that
  //  is larger than the actual photoelectric cross section. SGPH(M,IE) is
  //  defined as the largest value of the photoelectric x-section in the
  //  energy interval from ET(IE) to ET(IE+1). The photon mean free path is
  //  sampled by using this 'augmented' cross section and, to compensate
  //  for this, the photon survives (i.e., it is not absorbed) with a prob-
  //  ability such that the 'exact' photoelectric attenuation coefficient
  //  is reproduced. This trick allows subroutine JUMP to disregard the
  //  existence of absorption edges in the photoelectric x-section and to
  //  perform the tracking of photons faster.
  
  
  char CS5[17][6];

  double XGPHR[gph01.NDIM][17], X1[gph01.NDIM], Y1[gph01.NDIM], X2[gph01.NDIM], Y2[gph01.NDIM];
  int ISH[17], IZZ, NSHR, NDATA;
  
  strcpy(CS5[0], "total");
  strcpy(CS5[1], "CS-K ");
  strcpy(CS5[2], "CS-L1");
  strcpy(CS5[3], "CS-L2");
  strcpy(CS5[4], "CS-L3");
  strcpy(CS5[5], "CS-M1");
  strcpy(CS5[6], "CS-M2");
  strcpy(CS5[7], "CS-M3");
  strcpy(CS5[8], "CS-M4");
  strcpy(CS5[9], "CS-M5");
  strcpy(CS5[10], "CS-N1");
  strcpy(CS5[11], "CS-N2");
  strcpy(CS5[12], "CS-N3");
  strcpy(CS5[13], "CS-N4");
  strcpy(CS5[14], "CS-N5");
  strcpy(CS5[15], "CS-N6");
  strcpy(CS5[16], "CS-N7");

  //  ************  Read element x-section tables

  for(int IEL = 0; IEL < mat.NELEM; IEL++)
  {
    fscanf(IRD, "%*40c%d%*11c%d%*10c%d%*[^\n]", &IZZ, &NSHR, &NDATA);
    getc(IRD);
    
    if(INFO >= 2){ fprintf(IWR, "\n *** Photoelectric cross sections,  IZ =%3d,  NSHELL =%3d,  NDATA =%5d\n", IZZ, NSHR, NDATA);}

    if(IZZ != mat.IZ[IEL]){ penError(ERR_GPHaR_MAT_DF); return;}
    if(NDATA > (int)gph01.NDIM){ penError(ERR_GPHaR_DP); return;}
    if(NSHR > 16){ penError(ERR_GPHaR_SHELLS); return;}

    for(int IS = 0; IS < NSHR+1; IS++)
    {
      fscanf(IRD, "%d", &ISH[IS]);
    }
    fscanf(IRD,"%*[^\n]");
    getc(IRD);
    
    for(int IE = 0; IE < NDATA; IE++)
    {
      fscanf(IRD, "%lf", &gph01.ER[IE]);
      for(int IS = 0; IS < NSHR+1; IS++)
      {
        fscanf(IRD, "%lf", &XGPHR[IE][IS]);
      }
      fscanf(IRD,"%*[^\n]");
      getc(IRD);
    }

  //  ****  Remove shells with ionisation energies less than 50 eV.

    int NSHA;
    if(NSHR > 1)
    {
      bool Eixir = false;
      NSHA = NSHR;
      for(int IS = NSHA; IS >= 1; IS--)
      {
        if(elemDB.EB[IZZ-1][ISH[IS]-1] < 50.0)
        {
          NSHR = NSHR-1;
        }
        else
        {
          Eixir = true;
          break;
        }
      }
      if(NSHR < 0 && !Eixir){ NSHR = 0;}
    }
  
    if(INFO >= 2)
    {
      fprintf(IWR, "\n   Energy       ");
      for(int IS = 0; IS < NSHR+1; IS++)
      {
        if(IS<NSHR+1-1){fprintf(IWR, "%s       ", CS5[ISH[IS]]);}
        else{fprintf(IWR, "%s\n", CS5[ISH[IS]]);}
      }
      for(int IE = 0; IE < NDATA; IE++)
      {
        fprintf(IWR, "%12.5E", gph01.ER[IE]);
        for(int IS = 0; IS < NSHR+1; IS++)
        {
          fprintf(IWR, "%12.5E", XGPHR[IE][IS]);
        }
        fprintf(IWR, "\n");
      }
    }

    //Check if data for this element has already been loaded
    int IC;
    if(elemDB.NPHS[IZZ-1] == 0)
    {
      elemDB.IPHF[IZZ-1] = elemDB.NCUR;
      if(elemDB.NCUR+NDATA > constants::NTP)
      {
        fprintf(IWR, "Insufficient memory storage in GPHaR.\n");
        fprintf(IWR, "Increase the value of the parameter NTP to %d\n", elemDB.NCUR+NDATA);
        penError(ERR_GPHaR_MEM); return;
      }
      for(int IE = 0; IE < NDATA; IE++)
      {
        IC = elemDB.NCUR+IE;
        elemDB.EPH[IC] = log(gph01.ER[IE]);
        for(int IS = 0; IS < NSHR+1; IS++)
        {
          if( XGPHR[IE][ISH[IS]] > 1.0E-35){ elemDB.XPH[IC][ISH[IS]] = log(XGPHR[IE][ISH[IS]]);}
          else{ elemDB.XPH[IC][ISH[IS]] = log(1.0E-35);}
        }
      }
      elemDB.NCUR += NDATA;
      elemDB.IPHL[IZZ-1] = elemDB.NCUR-1;
      elemDB.NPHS[IZZ-1] = NSHR;
    }
  }

  //  ****  Total photoelectric attenuation coefficient.

  IZZ = mat.IZ[0];
  int IST = elemDB.IPHF[IZZ-1];
  int LST = elemDB.IPHL[IZZ-1];
  gph01.NPHD = 0;
  for(int I = IST; I <= LST; I++)
  {
    gph01.NPHD = gph01.NPHD+1;
    gph01.ER[gph01.NPHD-1] = exp(elemDB.EPH[I]);
    gph01.XSR[gph01.NPHD-1] = mat.STF[0]*exp(elemDB.XPH[I][0]);
  }
  int N1, N2;
  if(mat.NELEM > 1)
  {
    for(int IEL = 1; IEL < mat.NELEM; IEL++)
    {
      N1 = gph01.NPHD;
      for(int I = 0; I < N1; I++)
      {
        X1[I] = gph01.ER[I];
        Y1[I] = gph01.XSR[I];
      }
      IZZ = mat.IZ[IEL];
      IST = elemDB.IPHF[IZZ-1];
      LST = elemDB.IPHL[IZZ-1];
      N2 = 0;
      for(int I = IST; I <= LST; I++)
      {
        N2 = N2+1;
        X2[N2-1] = exp(elemDB.EPH[I]);
        Y2[N2-1] = mat.STF[IEL]*exp(elemDB.XPH[I][0]);
      }
      MERGE2(X1,Y1,X2,Y2,gph01.ER,gph01.XSR,N1,N2,gph01.NPHD);
      if(penGetError() != PEN_SUCCESS){ return;}
    }
  }

  //  ****  Total photoelectric cross section at the simulation grid points
  //        (slightly increased to simplify the interpolation).

  for(int I = 0; I < gph01.NPHD; I++)
  {
    X1[I] = log(gph01.ER[I]);
    Y1[I] = log(gph01.XSR[I]*mat.VMOL);
  }

  double EG1, EG2, DX, F1, F2;
  for(unsigned int IE = 0; IE < constants::NEGP-1; IE++)
  {
    unsigned int I1, I2;
    EG1 = DLEMP[IE];
    FINDI(X1,EG1,gph01.NPHD,I1);
    if(I1 == unsigned(gph01.NPHD-1)){I1--;}
    DX = X1[I1+1]-X1[I1];
    if(DX > 1.0E-15)
    {
      F1 = Y1[I1]+(Y1[I1+1]-Y1[I1])*(EG1-X1[I1])/DX;
    }
    else
    {
      F1 = Y1[I1];
    }
    EG2 = DLEMP[IE+1];
    FINDI(X1,EG2,gph01.NPHD,I2);
    if(I2 == unsigned(gph01.NPHD-1)){I2--;}
    DX = X1[I2+1]-X1[I2];
    if(DX > 1.0E-15)
    {
      F2 = Y1[I2]+(Y1[I2+1]-Y1[I2])*(EG2-X1[I2])/DX;
    }
    else
    {
      F2 = Y1[I2];
    }
      //  ****  To avoid interpolating the attenuation coefficient tables, we
      //        replace the photoelectric inverse mean free path (imfp) in each
      //        energy interval by its upper bound. The increase of the imfp
      //        is interpreted as the imfp of delta interactions.
      
    double FM;

    if(F1 > F2){ FM = F1;}
    else{ FM = F2;}
    if(I1+1 <= I2)
    {
      for(unsigned int I = I1+1; I <= I2; I++)
      {
        if(FM < Y1[I]){ FM = Y1[I];}        
      }
    }
    mat.SGPH[IE] = exp(FM);
  }
  mat.SGPH[constants::NEGP-1] = mat.SGPH[constants::NEGP-2];
}

//  *********************************************************************
//                       SUBROUTINE GRAaTI
//  *********************************************************************
void GRAaTI(const pen_material& mat, const double DLEMP1, const double DLFC, double &E, double &ECS)
{
  //  Total cross section for Rayleigh (coherent) photon scattering, in
  //  cm**2. Interpolated from input data.

  
  double XELN = log(E);
  double XEN = 1.0+(XELN-DLEMP1)*DLFC;
  int    KEN = (int)XEN;
  if(KEN == 0){ KEN = 1;}

  //  ****  Binary search.

  int II = mat.IED[KEN-1];
  int IU = mat.IEU[KEN-1];
  bool Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    int IT=(II+IU)/2;
    if(XELN > mat.ERA[IT])
    {
      II = IT;
    }
    else
    {
      IU = IT;
    }
    if(IU-II > 1){ Eixir = false; continue;}
  }
  ECS = exp(mat.XSRA[II]+(mat.XSRA[II+1]-mat.XSRA[II])
	    *(XELN-mat.ERA[II])/(mat.ERA[II+1]-mat.ERA[II]))/mat.VMOL;
}

//  *********************************************************************
//                       SUBROUTINE GPPa0
//  *********************************************************************
void GPPa0(pen_material& mat, const pen_elementDataBase& elemDB)
{
  //  Initialisation of the sampling algorithm for electron-positron pair
  //  production by photons in material M. Bethe-Heitler differential cross
  //  section.

  
  //  ***  Effective atomic number.
  
  double FACT = 0.0;
  int IZZ;
  for(int I = 0; I < mat.NELEM; I++)
  {
    IZZ = mat.IZ[I];
    FACT = FACT+IZZ*elemDB.ATW[IZZ-1]*mat.STF[I];
  }
  mat.ZEQPP = FACT/mat.AT;
  IZZ = (int)(mat.ZEQPP+0.25);
  if(IZZ <= 0){ IZZ = 1;}
  if(IZZ > 99){ IZZ = 99;}
  //  ****  DBM Coulomb correction.
  double ALZ = mat.ZEQPP/constants::SL;
  double A = ALZ*ALZ;
  double FC = A*(0.202059-A*(0.03693-A*(0.00835-A*(0.00201-A*(0.00049-A*(0.00012-A*0.00003)))))+1.0/(A+1.0));
  //  ****  Screening functions and low-energy correction.
  mat.BCB = 2.0/elemDB.RSCR[IZZ-1];
  mat.F0[0] = 4.0*log(elemDB.RSCR[IZZ-1]);
  mat.F0[1] = mat.F0[0]-4.0*FC;
}

//  *********************************************************************
//                       FUNCTION EINaDS
//  *********************************************************************
double EINaDS(double RMU, void* arg)
{
  //  Angular differential cross section for soft close inelastic colli-
  //  sions of electrons.
  //

  CEIN01* cein01 = (CEIN01*)(arg);
  

  double EINaDS_RETURN;
  double AUX = 2.0*RMU*(1.0-RMU);
  double DENOM = cein01->EI*AUX+constants::REV;
  double W = cein01->CPS*AUX/DENOM;
  double DWDMU = cein01->CPS*constants::REV*(2.0-4.0*RMU)/pow(DENOM,2);
  EINaDS_RETURN = (1.0+pow(W/(cein01->EE-W),2)-
		   (1.0-cein01->AMOL)*(W/(cein01->EE-W))+
		   cein01->AMOL*pow(W/cein01->EE,2))*DWDMU*pow(RMU,cein01->MOM)/pow(W,2);
  return EINaDS_RETURN;
}

//  *********************************************************************
//                       FUNCTION PINaDS
//  *********************************************************************
double PINaDS(double RMU, void* arg)
{
  //  Angular differential cross section for soft close inelastic colli-
  //  sions of positrons.


  CPIN01& cpin01 = *(CPIN01*)arg;
  
  double AUX = 2.0*RMU*(1.0-RMU);
  double DENOM = cpin01.EI*AUX+constants::REV;
  double W = cpin01.CPS*AUX/DENOM;
  double DWDMU = cpin01.CPS*constants::REV*(2.0-4.0*RMU)/pow(DENOM,2);
  double WE = W/cpin01.EI;
  double PINaDS_RETURN = (1.0-WE*(cpin01.BHA1-WE*(cpin01.BHA2-WE*(cpin01.BHA3-WE*cpin01.BHA4))))*DWDMU*pow(RMU,cpin01.MOM)/pow(W,2);
  return PINaDS_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE EELa0
//  *********************************************************************
void EELa0(double &XS0, double &XS1, double &XS2, double &XS0H, double &A, double &B, double &RNDC, double &XS1S, double &XS2S)
{

  //  This subroutine determines the parameters of the MW model for elastic
  //  scattering of electrons and positrons and initialises the mixed
  //  simulation algorithm (for particles with a given energy).

  //  Input arguments:
  //    XS0 .... total x-section (cm**2).
  //    XS1 .... 1st transport x-section (cm**2).
  //    XS2 .... 2nd transport x-section (cm**2).
  //    XS0H ... suggested x-section for hard events (cm**2).

  //  Output values:
  //    A, B ... angular distribution parameters.
  //    RNDC ... cutoff probability.
  //    XS0H ... adopted x-section for hard events (cm**2).
  //    XS1S ... 1st transport x-section for soft events (cm**2).
  //    XS2S ... 2nd transport x-section for soft events (cm**2). 


  double RMU1, RMU2;
  if(XS0 < 0.0){ penError(ERR_EELa0_NEG_CS); return;}
  RMU1 = XS1/(2.0*XS0);
  if(RMU1 > 0.48){ RMU1 = 0.48;} // Ensures numerical consistency.
  
  RMU2 = (3.0*XS1-XS2)/(6.0*XS0);
  if(RMU2 > 0.32){ RMU2 = 0.32;}
  
  if(RMU1 < 0.0 || RMU1 < RMU2)
  {
    printf("\n\n   *** The arguments in subroutine EELa0 are inconsistent.\n       XS0 = %14.7E, XS1 = %14.7E", XS0, XS1);
    penError(ERR_EELa0_ARG); return;
  }

  //  ****  Wentzel screening parameter.

  double TST;
  double AU, AL;
  A = 1.0;
  bool Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    A = A+A;
    TST = A*(A+1.0)*log((1.0+A)/A)-A-RMU1;
    if(TST < 0.0){ Eixir = false; continue;}
  }
  AU = A;
  AL = 0.0;
  Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    A = 0.5*(AL+AU);
    TST = A*(A+1.0)*log((1.0+A)/A)-A-RMU1;
    if(TST > 0.0)
    {
      AU = A;
    }
    else
    {
      AL = A;
    }
    if(fabs(TST) > 1.0E-15){ Eixir = false; continue;}
  }
  //  ****  At high energies, when truncation errors in the input tables
  //  are significant, we use delta scattering.
  if(RMU2-pow(RMU1,2) < 1.0E-12 || A < 1.0E-9)
  {
    B = 1.0;
    RNDC = 1.0-XS0H/XS0;
    if(RNDC < 1.0E-14){ RNDC = 0.0;}
    XS1S = XS1*RNDC;
    XS2S = XS2*RNDC;
    return;
  }

  double RMU1W, RMU2W;
  RMU1W = A*(A+1.0)*log((1.0+A)/A)-A;
  RMU2W = A*(1.0-2.0*RMU1W);
  B = (RMU2W-RMU2)/(RMU2W-RMU1W*RMU1W);

  //  ****  Case I.
  double RMUC;
  if(B > 0.0)
  {
    RNDC = 1.0-XS0H/XS0;
    if(RNDC < 1.0E-6)
    {
      RNDC = 0.0;
      XS0H = XS0;
      XS1S = 0.0;
      XS2S = 0.0;
      return;
    }

    double A1 = A+1.0;
    double B1 = 1.0-B;
    double RND0 = B1*A1*RMU1/(A+RMU1);
    RNDC = 1.0-XS0H/XS0;
    if(RNDC < RND0)
    {
      RMUC = RNDC*A/(B1*A1-RNDC);
      XS1S = B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
      XS2S = B1*(A*A1*pow(RMUC,2)/(A+RMUC))-2.0*A*XS1S;
    }
    else if(RNDC > RND0+B)
    {
      double RNDMB = RNDC-B;
      RMUC = RNDMB*A/(B1*A1-RNDMB);
      XS1S = B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
      XS2S = B1*(A*A1*pow(RMUC,2)/(A+RMUC))-2.0*A*XS1S;
      XS1S = XS1S+B*RMU1;
      XS2S = XS2S+B*pow(RMU1,2);
    }
    else
    {
      RMUC = RMU1;
      double WB = RNDC-RND0;
      XS1S = B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
      XS2S = B1*(A*A1*pow(RMUC,2)/(A+RMUC))-2.0*A*XS1S;
      XS1S = XS1S+WB*RMU1;
      XS2S = XS2S+WB*pow(RMU1,2);
    }
    XS2S = 6.0*XS0*(XS1S-XS2S);
    XS1S = 2.0*XS0*XS1S;
    return;
  }
  if(B > -1.0E-12)
  {
    B = 0.0;
    RNDC = 1.0-XS0H/XS0;
    double A1 = A+1.0;
    RMUC = RNDC*A/(A1-RNDC);
    XS1S = A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
    XS2S = (A*A1*pow(RMUC,2)/(A+RMUC))-2.0*A*XS1S;
    XS2S = 6.0*XS0*(XS1S-XS2S);
    XS1S = 2.0*XS0*XS1S;
    return;
  }

  //  ****  Case II. 

  double D1, D2;
  
  double C1_Aux = 8.333333333333333E-1;
  double C2_Aux = 7.083333333333333E-1;
  D1 = C2_Aux-RMU2;
  D2 = C1_Aux-RMU1;
  double D3 = C2_Aux*RMU1-C1_Aux*RMU2;
  AL = 1.0E-24;
  AU = A;
  Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    A = 0.5*(AL+AU);
    RMU1W = A*(A+1.0)*log((1.0+A)/A)-A;
    RMU2W = A*(1.0-2.0*RMU1W);
    double F = D1*RMU1W-D2*RMU2W-D3;
    if(F < 0.0)
    {
      AL = A;
    }
    else
    {
      AU = A;
    }
    if(AU-AL > 1.0E-14*A){ Eixir = false; continue;}
  }
  B = (RMU1W-RMU1)/(C1_Aux-RMU1W);
  
  RNDC = 1.0-XS0H/XS0;
  if(RNDC < 1.0E-10)
  {
    RNDC = 0.0;
    XS0H = XS0;
    XS1S = 0.0;
    XS2S = 0.0;
    return;
  }
  double A1 = A+1.0;
  double B1 = 1.0+B;
  double RNDCM = B1*A1*0.5/(A+0.5);
  if(RNDC > RNDCM)
  {
    RNDC = RNDCM;
    XS0H = XS0*(1.0-RNDC);
  }
  RMUC = RNDC*A/(B1*A1-RNDC);
  XS1S = B1*A*A1*(log((A+RMUC)/A)-(RMUC/(A+RMUC)));
  XS2S = B1*(A*A1*pow(RMUC,2)/(A+RMUC))-2.0*A*XS1S;
  XS2S = 6.0*XS0*(XS1S-XS2S);
  XS1S = 2.0*XS0*XS1S;
}

//  *********************************************************************
//                       SUBROUTINE DCSEL0
//  *********************************************************************
void DCSEL0(const double E, int IELEC, CDCSEP& dcsep)
{
  //     This subroutine computes a table of the molecular elastic DCSs for
  //  electrons (IELEC=-1) or positrons (IELEC=+1) with kinetic energy  E
  //  (in eV) by log-log cubic spline interpolation in E of the DCS table
  //  prepared by subroutine ELINIT.


  const unsigned int NE = dcsep.NE;
  const unsigned int NA = dcsep.NA;
  
  double Y[NE], A[NE], B[NE], C[NE], D[NE];

  if(E < 49.999 || E > 1.0001E8)
  {
    printf(" Error in DCSEL0: Energy out of range.\n");
    penError(ERR_DCSEL_ENERGY_RANGE); return;
  }
  
  double EL = log(E);
  unsigned int JE;
  FINDI(dcsep.ETL,EL,NE,JE);

  for(unsigned int IA = 0; IA < NA; IA++)
  {
    for(unsigned int IE = 0; IE < NE; IE++)
    {
      if(IELEC == -1)
    {
      Y[IE] = log(dcsep.EDCS[IE][IA]);
    }
      else
    {
      Y[IE] = log(dcsep.PDCS[IE][IA]);
    }
    }
    SPLINE(dcsep.ETL,Y,A,B,C,D,0.0,0.0,NE);
    if(penGetError() != PEN_SUCCESS){ return;}
    dcsep.DCSIL[IA] = A[JE]+EL*(B[JE]+EL*(C[JE]+EL*D[JE]));
    dcsep.DCSI[IA] = exp(dcsep.DCSIL[IA]);
  }

  for(unsigned int IE = 0; IE < NE; IE++)
  {
    if(IELEC == -1)
    {
      Y[IE] = log(dcsep.ECS[IE]);
  }
    else
  {
      Y[IE] = log(dcsep.PCS[IE]);
    }
  }
  SPLINE(dcsep.ETL,Y,A,B,C,D,0.0,0.0,NE);
  if(penGetError() != PEN_SUCCESS){ return;}
  dcsep.CSI = exp(A[JE]+EL*(B[JE]+EL*(C[JE]+EL*D[JE])));
  
  for(unsigned int IE = 0; IE < NE; IE++)
  {
    if(IELEC == -1)
  {
      Y[IE] = log(dcsep.ETCS1[IE]);
  }
    else
  {
      Y[IE] = log(dcsep.PTCS1[IE]);
    }
  }
  SPLINE(dcsep.ETL,Y,A,B,C,D,0.0,0.0,NE);
  if(penGetError() != PEN_SUCCESS){ return;}
  dcsep.TCS1I = exp(A[JE]+EL*(B[JE]+EL*(C[JE]+EL*D[JE])));

  for(unsigned int IE = 0; IE < NE; IE++)
  {
    if(IELEC == -1)
  {
      Y[IE] = log(dcsep.ETCS2[IE]);
  }
    else
  {
      Y[IE] = log(dcsep.PTCS2[IE]);
    }
  }
  SPLINE(dcsep.ETL,Y,A,B,C,D,0.0,0.0,NE);
  if(penGetError() != PEN_SUCCESS){ return;}
  dcsep.TCS2I = exp(A[JE]+EL*(B[JE]+EL*(C[JE]+EL*D[JE])));
}

//  *********************************************************************
//                       FUNCTION DCSEL
//  *********************************************************************
double DCSEL(double RMU, void* arg)
{
  //  This function computes the DCS in (cm**2/sr) by linear log-log inter-
  //  polation in RMU=(1-cos(theta))/2. It is initialised by subroutine
  //  DCSEL0, which must be called before using DCSEL.


  CDCSEP* dcsep = (CDCSEP*)(arg);

  double* XMUL  = dcsep->XMUL;
  double* DCSIL = dcsep->DCSIL;

  const unsigned int NA = dcsep->NA;
  
  double RMUL;
  if(RMU < 1.0E-35){ RMUL = 1.0E-35;}
  else{ RMUL = RMU;}

  if(RMUL > 0.999999999999){ RMUL = 0.999999999999;}
  
  RMUL = log(RMUL);
  unsigned int I;
  FINDI(XMUL,RMUL,NA,I);
  double DCSEL_RETURN = exp(DCSIL[I]+(DCSIL[I+1]-DCSIL[I])*((RMUL-XMUL[I])/(XMUL[I+1]-XMUL[I]))); //Variable que torna la funcio
  //printf("I=%d DCSIL=%12.5E XMUL=%12.5E RMUL=%12.5E DCSEL=%12.5E\n",I,DCSIL[I-1],XMUL[I-1],RMUL,DCSEL_RETURN);
  return DCSEL_RETURN;
}

//  *********************************************************************
//                       FUNCTION GRAaF2
//  *********************************************************************
double GRAaF2( double Q2, void* arg)
{
  //  Squared molecular form factor, as a function of (Q*SL/REV)**2.


  pen_material& mat = *(pen_material*)(arg);
  
  double GRAaF2_RETURN; //Substitueix a la variable GRAaF2 en fortran. Es el valor que torna la funcio
  if(Q2 < 1.0E-9)
  {
    GRAaF2_RETURN = mat.FF0;
  }
  else if(Q2 > mat.QQM)
  {
    GRAaF2_RETURN = 0.0;
  }
  else
  {
    double QL = log(Q2);
    unsigned int I;
    FINDI(mat.QQ,QL,constants::NQ,I);

    double F2 = mat.AR[I]+
      QL*(mat.BR[I]+
	  QL*(mat.CR[I]+
	      QL*mat.DR[I]));
    
    GRAaF2_RETURN = exp(F2);
  }
  return GRAaF2_RETURN;
}

