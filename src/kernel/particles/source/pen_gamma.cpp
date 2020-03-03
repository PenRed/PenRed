
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

#include "pen_gamma.hh" 

// Gamma particle class

pen_gamma::pen_gamma(const pen_context& contextIn,
	    pen_particleStack<pen_particleState>& stackEin,
	    pen_particleStack<pen_particleState>& stackPin,
	    pen_particleStack<pen_state_gPol>& stackGin)

  : abc_particle(contextIn,PEN_PHOTON,4,0.0E0),
    compton(stackEin,stackGin),
    photoelectric(stackEin,stackGin),
    pairProduct(stackEin,stackPin,stackGin)
    
    
{

  //Check ids of interactions
  int auxID = rayleigh.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Gamma rayleigh ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = compton.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Gamma compton ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = photoelectric.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Gamma photelectric ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = pairProduct.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Gamma pair production ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  
  rayleigh.init(context);
  compton.init(context);
  photoelectric.init(context);
  pairProduct.init(context);
  
}

void pen_gamma::JUMP(double &DS, pen_rand& penRand, const double /*DSMAX*/)
{
  //
  //  ************  Photons.
  //
  
  if(state.E < ELAST1)
    {
      //Get material
      const pen_material& mat = *pmat;    
      
      XEL = log(state.E);
      XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
      KE = (int)XE;
      XEK = XE-(double)KE;
      IMFP(mat);      
      ELAST1 = state.E;

      ST = P[0];
      ST += P[1];
      ST += P[2];
      ST += P[3];
    }

  DS = -log(penRand.rand())/ST;  
}

void pen_gamma::JUMPF(double &DS, pen_rand& penRand, const double /*DSMAX*/)
{
  //  ************  Photons.

  //Check energy lost since last call
  if(state.E < ELAST1)
    {
      //Get material
      const pen_material& mat = *pmat;
            
      XEL = log(state.E);
      XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
      KE = (int)XE;
      XEK = XE-(double)KE;
      IMFP(mat);
      ELAST1 = state.E;
      //  ****  Interaction forcing.
      for(unsigned KCOL = 0; KCOL < interactions; KCOL++)
	{
	  if(context.FORCE[state.IBODY][PEN_PHOTON][KCOL] > 1.0 && P[KCOL] > 1.0E-16)
	    {
	      P0[KCOL] = P[KCOL];
	      P[KCOL] = P[KCOL]*context.FORCE[state.IBODY][PEN_PHOTON][KCOL];
	      LFORC[KCOL] = true;
	    }
	  else
	    {
	      LFORC[KCOL] = false;
	    }
	}
      
      ST = P[0];
      ST += P[1];
      ST += P[2];
      ST += P[3];
    }

  DS = -log(penRand.rand())/ST;
    
}

void pen_gamma::KNOCK(double &DE, int &ICOL, pen_rand& penRand)
{

  //  ************  Photons (KPAR=PEN_PHOTON).

  //Get material
  const pen_material& mat = *pmat;
  
  double STS = ST*penRand.rand();
  double SS = 0.0;

  SS += P[GAMMA_RAYLEIGH];
  if(SS > STS)
    {
      ICOL = rayleigh.interact(context,mat,(*this),DE,penRand);
      return;
    }
  
  SS += P[GAMMA_COMPTON];
  if(SS > STS)
    {
      ICOL = compton.interact(context,mat,(*this),DE,penRand);
      return;      
    }

  SS += P[GAMMA_PHOTOELECTRIC];
  if(SS > STS)
    {
      ICOL = photoelectric.interact(context,mat,(*this),DE,penRand);
      return;      
    }

  SS += P[GAMMA_PAIR_PRODUCTION];
  if(SS > STS)
    {
      ICOL = pairProduct.interact(context,mat,(*this),DE,penRand);
      return;      
    }
  
}

void pen_gamma::KNOCKF(double &DE, int &ICOL, pen_rand& penRand)
{
      
  //  ************  Photons (KPAR=PEN_PHOTON).

  //Get material
  const pen_material& mat = *pmat;

  double STS = ST*penRand.rand();
  double SS = 0.0;

  SS += P[GAMMA_RAYLEIGH];
  if(SS > STS)
    {
      ICOL = rayleigh.interactF(context,mat,(*this),DE,penRand);
      return;
    }
  
  SS += P[GAMMA_COMPTON];
  if(SS > STS)
    {
      ICOL = compton.interactF(context,mat,(*this),DE,penRand);
      return;      
    }

  SS += P[GAMMA_PHOTOELECTRIC];
  if(SS > STS)
    {
      ICOL = photoelectric.interactF(context,mat,(*this),DE,penRand);
      return;      
    }

  SS += P[GAMMA_PAIR_PRODUCTION];
  if(SS > STS)
    {
      ICOL = pairProduct.interactF(context,mat,(*this),DE,penRand);
      return;      
    }
}

void pen_gamma::softEloss(double& /*X*/,
			  double& /*Y*/,
			  double& /*Z*/,
			  double& DE,
			  pen_rand& /*penRand*/){
  DE = 0.0;
  return;  
}

void pen_gamma::dpage()
{
  //  The input parameters dsef and dstot are, respectively, the length of
  //  the step in the 'original' material and the total length of the step,
  //  including segments in void volumes. When using PENGEOM, dstot is
  //  given by subroutine STEP.

  //  Usage:
  //  1) CALL PAGE0 just after starting each primary particle.
  //  2) CALL DPAGE(dsef,dstot) at the end of each free flight.

  const double SLCM=2.99792458E10;  // Speed of light (cm/s)
  const double RSLCM=1.0/SLCM;  // Reciprocal of SLCM (1/cm).

  //  ****  Updating the particle's age after each free flight.

      if(state.MAT == 0)
	{
	  state.PAGE = state.PAGE+dsef*RSLCM;
	}
      else
	{
	  state.PAGE = state.PAGE+dstot*RSLCM;
	}  
}

// Gamma Interactions

// photons coherent scattering (Rayleigh)

int pen_gRCS::interact(const pen_context& /*context*/, const pen_material& mat, pen_gamma& gamma, double& DE, pen_rand& penRand) const
{
  //  ****  Rayleigh scattering (ICOL=GAMMA_RAYLEIGH).
  double CDT;
  int IEFF; 
  DE = 0.0;
  GRAa(mat,gamma.state.E,gamma.KE,gamma.XEL,CDT,IEFF,penRand);
  //  ****  Delta interaction. Introduced to correct for the use of an
  //        upper bound of the Rayleigh attenuation coefficient.
  if(IEFF == 0)
    {
      return GAMMA_DELTA;
    }
  double DF;
  if(gamma.state.IPOL == 1)
    {
      DIRPOL(CDT,DF,0.0,gamma.state.SP1,gamma.state.SP2,gamma.state.SP3,
	     gamma.state.U,gamma.state.V,gamma.state.W,penRand);
    }
  else
    {
      DF = constants::TWOPI*penRand.rand();
      DIRECT(CDT,DF,gamma.state.U,gamma.state.V,gamma.state.W);
    }
  return GAMMA_RAYLEIGH;
}

int pen_gRCS::interactF(const pen_context& /*context*/, const pen_material& mat, pen_gamma& gamma, double& DE, pen_rand& penRand) const
{

  //  ****  Rayleigh scattering (ICOL=GAMMA_RAYLEIGH).

  double CDT;
  int IEFF;
  int ICOL;
  gamma.DEA = 0.0;
  DE = 0.0;
  GRAa(mat,gamma.state.E,gamma.KE,gamma.XEL,CDT,IEFF,penRand);
  //  ****  Delta interaction. Introduced to correct for the use of an
  //        upper bound of the Rayleigh attenuation coefficient.
  if(IEFF == 0)
    {
      ICOL = GAMMA_DELTA;
      return ICOL;
    }
  ICOL = GAMMA_RAYLEIGH;
  double DF;
  if(gamma.state.IPOL == 1)
    {
      DIRPOL(CDT,DF,0.0E0,gamma.state.SP1,gamma.state.SP2,gamma.state.SP3,
	     gamma.state.U,gamma.state.V,gamma.state.W,penRand);
    }
  else
    {
      DF = constants::TWOPI*penRand.rand();
      DIRECT(CDT,DF,gamma.state.U,gamma.state.V,gamma.state.W);
    }
  return ICOL;
}
  
// photons incoherent scattering (Compton)

int pen_gCIS::interact(const pen_context& context, const pen_material& mat, pen_gamma& gamma, double& DE, pen_rand& penRand) const
{
  //  ****  Compton scattering (ICOL=GAMMA_COMPTON).

  const double RREV = 1.0/constants::REV;
  
  double EP,CDT,ES,CDTS;
  int IZA,ISA;
  GCOa(mat,gamma.state.E,DE,EP,CDT,ES,CDTS,IZA,ISA,penRand);
  double US = gamma.state.U;
  double VS = gamma.state.V;
  double WS = gamma.state.W;
  double DF = -1.0;
  if(IZA > 0 && ISA < 17)
    {
      int ks;
      RELAX(context.elements,mat,gamma.state,GAMMA_COMPTON,1,
	    IZA,ISA,ks,stackE,stackG,penRand);
      if(penGetError() != PEN_SUCCESS){ return GAMMA_COMPTON;}
    }
  //  ****  New direction and energy.
  if(EP > mat.EABS[PEN_PHOTON])
    {
      if(gamma.state.IPOL == 1)
	{
	  double ECDT = gamma.state.E*RREV*(1.0-CDT);
	  double CONS = ECDT*ECDT/(1.0+ECDT);
	  DIRPOL(CDT,DF,CONS,gamma.state.SP1,gamma.state.SP2,gamma.state.SP3,
		 gamma.state.U,gamma.state.V,gamma.state.W,penRand);
	}
      else
	{
	  DF = constants::TWOPI*penRand.rand();
	  DIRECT(CDT,DF,gamma.state.U,gamma.state.V,gamma.state.W);
	}
      gamma.state.E = EP;
    }
  else
    {
      DE = gamma.state.E;
      gamma.state.E = 0.0;
    }
  //  ****  Compton electron.
  if(ES > mat.EABS[PEN_ELECTRON])
    {
      if(DF < -0.5){ DF = constants::TWOPI*penRand.rand();}
      double DFS = DF+constants::PI;

      DIRECT(CDTS,DFS,US,VS,WS);

      pen_particleState stateStore;

      stateStore.E = ES;

      stateStore.X = gamma.state.X;
      stateStore.Y = gamma.state.Y;
      stateStore.Z = gamma.state.Z;

      stateStore.U = US;
      stateStore.V = VS;
      stateStore.W = WS;

      stateStore.WGHT = gamma.state.WGHT;
      stateStore.IBODY = gamma.state.IBODY;
      stateStore.MAT = gamma.state.MAT;
      
      stateStore.ILB[0] = gamma.state.ILB[0]+1;
      stateStore.ILB[1] = gamma.getKpar();
      stateStore.ILB[2] = GAMMA_COMPTON;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = gamma.state.ILB[4];

      stateStore.LAGE = gamma.state.LAGE;
      stateStore.PAGE = gamma.state.PAGE;

      stackE.store(stateStore);
      
      if(penGetError() != PEN_SUCCESS){ return GAMMA_COMPTON;}
    }
  return GAMMA_COMPTON;
}

int pen_gCIS::interactF(const pen_context& context, const pen_material& mat, pen_gamma& gamma, double& DE, pen_rand& penRand) const
{
  
  //  ****  Compton scattering (ICOL=GAMMA_COMPTON).

  const double RREV = 1.0/constants::REV;
  
  double WFORCE = 1.0;
  bool LCOL;
  int ICOL = GAMMA_COMPTON;
  if(gamma.LFORC[GAMMA_COMPTON])  // Forced interaction.
    {
      WFORCE = gamma.P0[GAMMA_COMPTON]/gamma.P[GAMMA_COMPTON];
      if(penRand.rand() < WFORCE)
	{
	  LCOL = true;
	}
      else
	{
	  LCOL = false;
	}
    }
  else // Unforced interaction.
    {
      LCOL = true;
    }
  double EP,CDT,ES,CDTS;
  int IZA,ISA;
  GCOa(mat,gamma.state.E,DE,EP,CDT,ES,CDTS,IZA,ISA,penRand);
  
  double US = gamma.state.U;
  double VS = gamma.state.V;
  double WS = gamma.state.W;
  double DF = -1.0;
  if(IZA > 0 && ISA < 17)
    {
      int ks;
      double WGHTA = gamma.state.WGHT;
      gamma.state.WGHT = gamma.state.WGHT*WFORCE;
      RELAX(context.elements,mat,gamma.state,ICOL,1,
	    IZA,ISA,ks,stackE,stackG,penRand);      
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
      gamma.state.WGHT = WGHTA;
    }
  //  ****  Compton electron.
  if(ES > mat.EABS[PEN_ELECTRON])
    {
      if(DF < -0.5){ DF = constants::TWOPI*penRand.rand();}
      double DFS = DF+constants::PI;

      DIRECT(CDTS,DFS,US,VS,WS);
      
      pen_particleState stateStore;

      stateStore.E = ES;

      stateStore.X = gamma.state.X;
      stateStore.Y = gamma.state.Y;
      stateStore.Z = gamma.state.Z;

      stateStore.U = US;
      stateStore.V = VS;
      stateStore.W = WS;

      stateStore.WGHT = gamma.state.WGHT*WFORCE;
      stateStore.IBODY = gamma.state.IBODY;
      stateStore.MAT = gamma.state.MAT;
      
      stateStore.ILB[0] = gamma.state.ILB[0]+1;
      stateStore.ILB[1] = gamma.getKpar();
      stateStore.ILB[2] = ICOL;
      stateStore.ILB[3] = IZA*1000000+ISA;
      stateStore.ILB[4] = gamma.state.ILB[4];

      stateStore.LAGE = gamma.state.LAGE;
      stateStore.PAGE = gamma.state.PAGE;
            
      stackE.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
    }
  gamma.DEA = DE;
  DE = gamma.DEA*WFORCE;
  //  ****  New direction and energy.
  if(LCOL)
    {
      if(EP > mat.EABS[PEN_PHOTON])
	{
	  if(gamma.state.IPOL == 1)
	    {
	      double ECDT = gamma.state.E*RREV*(1.0-CDT);
	      double CONS = ECDT*ECDT/(1.0+ECDT);
	      DIRPOL(CDT,DF,CONS,
		     gamma.state.SP1,gamma.state.SP2,gamma.state.SP3,
		     gamma.state.U,gamma.state.V,gamma.state.W,penRand);
	    }
	  else
	    {
	      DF = constants::TWOPI*penRand.rand();
	      DIRECT(CDT,DF,gamma.state.U,gamma.state.V,gamma.state.W);
	    }
	  gamma.state.E = EP;
	}
      else
	{
	  gamma.DEA = gamma.DEA+EP;
	  DE = DE+EP;
	  gamma.state.E = 0.0;
	}
    }
  return ICOL;
}
  
// photons photoelectric interaction

int pen_gPHE::interact(const pen_context& context, const pen_material& mat, pen_gamma& gamma, double& DE, pen_rand& penRand) const
{

  //  ****  Photoelectric absorption (ICOL=GAMMA_PHOTOELECTRIC).

  double ES;
  int IZA, ISA;
  GPHa(mat,context.elements,gamma.state.E,gamma.KE,gamma.XEL,ES,IZA,ISA,penRand);
  //  ****  Delta interaction. Introduced to correct for the use of an
  //        upper bound of the photoelectric attenuation coefficient.
  if(IZA == 0)
    {
      DE = 0.0;
      return GAMMA_DELTA;
    }
    
  if(ES > mat.EABS[PEN_ELECTRON])
    {
      double CDTS;
      double DFS;
      SAUTER(ES,CDTS,DFS,penRand,gamma.state);

      pen_particleState stateStore;

      stateStore.E = ES;

      stateStore.X = gamma.state.X;
      stateStore.Y = gamma.state.Y;
      stateStore.Z = gamma.state.Z;

      stateStore.U = gamma.state.U;
      stateStore.V = gamma.state.V;
      stateStore.W = gamma.state.W;

      stateStore.WGHT = gamma.state.WGHT;
      stateStore.IBODY = gamma.state.IBODY;
      stateStore.MAT = gamma.state.MAT;
      
      stateStore.ILB[0] = gamma.state.ILB[0]+1;
      stateStore.ILB[1] = gamma.getKpar();
      stateStore.ILB[2] = GAMMA_PHOTOELECTRIC;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = gamma.state.ILB[4];

      stateStore.LAGE = gamma.state.LAGE;
      stateStore.PAGE = gamma.state.PAGE;
      
      DIRECT(CDTS,DFS,stateStore.U,stateStore.V,stateStore.W);

      stackE.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return GAMMA_PHOTOELECTRIC;}
    }
  if(ISA < 17)
    {
      int ks;
      RELAX(context.elements,mat,gamma.state,GAMMA_PHOTOELECTRIC,1,
	    IZA,ISA,ks,stackE,stackG,penRand);
      if(penGetError() != PEN_SUCCESS){ return GAMMA_PHOTOELECTRIC;}
    }
  DE = gamma.state.E;
  gamma.state.E = 0.0;
  return GAMMA_PHOTOELECTRIC;
}

int pen_gPHE::interactF(const pen_context& context, const pen_material& mat, pen_gamma& gamma, double& DE, pen_rand& penRand) const
{

  //  ****  Photoelectric absorption (ICOL=GAMMA_PHOTOELECTRIC).

  int ICOL = GAMMA_PHOTOELECTRIC;
  double WFORCE = 1.0;
  bool LCOL;
  double ES;
  int IZA, ISA;
  if(gamma.LFORC[GAMMA_PHOTOELECTRIC])  // Forced interaction.
    {
      WFORCE = gamma.P0[GAMMA_PHOTOELECTRIC]/gamma.P[GAMMA_PHOTOELECTRIC];
      if(penRand.rand() < WFORCE)
	{
	  LCOL = true;
	}
      else
	{
	  LCOL = false;
	}
    }
  else  // Unforced interaction.
    {
      LCOL = true;
    }

  GPHa(mat,context.elements,gamma.state.E,gamma.KE,gamma.XEL,ES,IZA,ISA,penRand);
  //  ****  Delta interaction. Introduced to correct for the use of an
  //        upper bound of the photoelectric attenuation coefficient.
  if(IZA == 0)
    {
      ICOL = GAMMA_DELTA;
      gamma.DEA = 0.0;
      DE = 0.0;
      return ICOL;
    }

  if(ES > mat.EABS[PEN_ELECTRON])
    {
      double CDTS;
      double DFS;
      SAUTER(ES,CDTS,DFS,penRand,gamma.state);

      pen_particleState stateStore;

      stateStore.E = ES;

      stateStore.X = gamma.state.X;
      stateStore.Y = gamma.state.Y;
      stateStore.Z = gamma.state.Z;
      
      stateStore.U = gamma.state.U;
      stateStore.V = gamma.state.V;
      stateStore.W = gamma.state.W;

      stateStore.WGHT = gamma.state.WGHT*WFORCE;
      stateStore.IBODY = gamma.state.IBODY;
      stateStore.MAT = gamma.state.MAT;
      
      stateStore.ILB[0] = gamma.state.ILB[0]+1;
      stateStore.ILB[1] = gamma.getKpar();
      stateStore.ILB[2] = ICOL;
      stateStore.ILB[3] = IZA*1000000+ISA;
      stateStore.ILB[4] = gamma.state.ILB[4];

      stateStore.LAGE = gamma.state.LAGE;
      stateStore.PAGE = gamma.state.PAGE;
      
      DIRECT(CDTS,DFS,stateStore.U,stateStore.V,stateStore.W);

      stackE.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
    }
  if(ISA < 17)
    {
      int ks;
      double WGHTA = gamma.state.WGHT;
      gamma.state.WGHT = gamma.state.WGHT*WFORCE;
      RELAX(context.elements,mat,gamma.state,ICOL,1,
	    IZA,ISA,ks,stackE,stackG,penRand);      
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
      gamma.state.WGHT = WGHTA;
    }
  gamma.DEA = gamma.state.E;
  DE = gamma.state.E*WFORCE;
  if(LCOL){ gamma.state.E = 0.0;}
  return ICOL;
}
  
// photons pair production interaction

int pen_gBBP::interact(const pen_context& context, const pen_material& mat, pen_gamma& gamma, double& DE, pen_rand& penRand) const
{

  //  ****  Electron-positron pair production (ICOL=GAMMA_PAIR_PRODUCTION).

  const double TREV = 2.0*constants::REV;
  
  double EE,CDTE,EP,CDTP;
  int IZA,ISA;
  GPPa(mat,gamma.state.E,gamma.KE,gamma.XEK,EE,CDTE,EP,CDTP,IZA,ISA,penRand);
  DE = gamma.state.E;
  //  ****  Electron.
  if(EE > mat.EABS[PEN_ELECTRON])
    {
      double DF = constants::TWOPI*penRand.rand();

      pen_particleState stateStore;

      stateStore.E = EE;

      stateStore.X = gamma.state.X;
      stateStore.Y = gamma.state.Y;
      stateStore.Z = gamma.state.Z;
      
      stateStore.U = gamma.state.U;
      stateStore.V = gamma.state.V;
      stateStore.W = gamma.state.W;
      
      stateStore.WGHT = gamma.state.WGHT;
      stateStore.IBODY = gamma.state.IBODY;
      stateStore.MAT = gamma.state.MAT;
      
      stateStore.ILB[0] = gamma.state.ILB[0]+1;
      stateStore.ILB[1] = gamma.getKpar();
      stateStore.ILB[2] = GAMMA_PAIR_PRODUCTION;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = gamma.state.ILB[4];

      stateStore.LAGE = gamma.state.LAGE;
      stateStore.PAGE = gamma.state.PAGE;
      
      DIRECT(CDTE,DF,stateStore.U,stateStore.V,stateStore.W);

      stackE.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return GAMMA_PAIR_PRODUCTION;}
    }
  //  ****  Positron.
  if(EP > mat.EABS[PEN_POSITRON])
    {
      double DF = constants::TWOPI*penRand.rand();

      pen_particleState stateStore;

      stateStore.E = EP;
      
      stateStore.X = gamma.state.X;
      stateStore.Y = gamma.state.Y;
      stateStore.Z = gamma.state.Z;
      
      stateStore.U = gamma.state.U;
      stateStore.V = gamma.state.V;
      stateStore.W = gamma.state.W;

      stateStore.WGHT = gamma.state.WGHT;
      stateStore.IBODY = gamma.state.IBODY;
      stateStore.MAT = gamma.state.MAT;
      
      stateStore.ILB[0] = gamma.state.ILB[0]+1;
      stateStore.ILB[1] = gamma.getKpar();
      stateStore.ILB[2] = GAMMA_PAIR_PRODUCTION;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = gamma.state.ILB[4];

      stateStore.LAGE = gamma.state.LAGE;
      stateStore.PAGE = gamma.state.PAGE;
      
      DIRECT(CDTP,DF,stateStore.U,stateStore.V,stateStore.W);

      stackP.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return GAMMA_PAIR_PRODUCTION;}
      //  ****  The positron carries a 'latent' energy of 1022 keV.
      DE = DE-TREV;
    }
  else
    {
      PANaR(gamma.state,mat.EABS[PEN_PHOTON],stackG,penRand);
      if(penGetError() != PEN_SUCCESS){ return GAMMA_PAIR_PRODUCTION;}
    }
  gamma.state.E = 0.0;
  //  ****  Atomic relaxation after triplet production.
  if(ISA < 17)
    {
      int ks;
      RELAX(context.elements,mat,gamma.state,GAMMA_PAIR_PRODUCTION,1,
	    IZA,ISA,ks,stackE,stackG,penRand);
      if(penGetError() != PEN_SUCCESS){ return GAMMA_PAIR_PRODUCTION;}
    }
  return GAMMA_PAIR_PRODUCTION;
}

int pen_gBBP::interactF(const pen_context& context, const pen_material& mat, pen_gamma& gamma, double& DE, pen_rand& penRand) const
{
	  
  //  ****  Electron-positron pair production (ICOL=GAMMA_PAIR_PRODUCTION).

  const double TREV = 2.0*constants::REV;
  
  int ICOL = GAMMA_PAIR_PRODUCTION;
  double WFORCE = 1.0;
  double EE,CDTE,EP,CDTP;
  int IZA,ISA;
  bool LCOL;
  if(gamma.LFORC[GAMMA_PAIR_PRODUCTION])  // Forced interaction.
    {
      WFORCE = gamma.P0[GAMMA_PAIR_PRODUCTION]/gamma.P[GAMMA_PAIR_PRODUCTION];
      if(penRand.rand() < WFORCE)
	{
	  LCOL = true;
	}
      else
	{
	  LCOL = false;
	}
    }
  else  // Unforced interaction.
    {
      LCOL = true;
    }

  GPPa(mat,gamma.state.E,gamma.KE,gamma.XEK,EE,CDTE,EP,CDTP,IZA,ISA,penRand);
  DE = gamma.state.E;
  //  ****  Electron.
  if(EE > mat.EABS[PEN_ELECTRON])
    {
      pen_particleState stateStore;
      
      double DF = constants::TWOPI*penRand.rand();

      stateStore.E = EE;
      
      stateStore.X = gamma.state.X;
      stateStore.Y = gamma.state.Y;
      stateStore.Z = gamma.state.Z;
      
      stateStore.U = gamma.state.U;
      stateStore.V = gamma.state.V;
      stateStore.W = gamma.state.W;

      stateStore.WGHT = gamma.state.WGHT*WFORCE;
      stateStore.IBODY = gamma.state.IBODY;
      stateStore.MAT = gamma.state.MAT;
      
      stateStore.ILB[0] = gamma.state.ILB[0]+1;
      stateStore.ILB[1] = gamma.getKpar();
      stateStore.ILB[2] = ICOL;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = gamma.state.ILB[4];

      stateStore.LAGE = gamma.state.LAGE;
      stateStore.PAGE = gamma.state.PAGE;
      
      DIRECT(CDTE,DF,stateStore.U,stateStore.V,stateStore.W);

      stackE.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
    }
  //  ****  Positron.
  if(EP > mat.EABS[PEN_POSITRON])
    {

      pen_particleState stateStore;

      double DF = constants::TWOPI*penRand.rand();

      stateStore.E = EP;
      
      stateStore.X = gamma.state.X;
      stateStore.Y = gamma.state.Y;
      stateStore.Z = gamma.state.Z;
      
      stateStore.U = gamma.state.U;
      stateStore.V = gamma.state.V;
      stateStore.W = gamma.state.W;

      stateStore.WGHT = gamma.state.WGHT*WFORCE;
      stateStore.IBODY = gamma.state.IBODY;
      stateStore.MAT = gamma.state.MAT;
      
      stateStore.ILB[0] = gamma.state.ILB[0]+1;
      stateStore.ILB[1] = gamma.getKpar();
      stateStore.ILB[2] = ICOL;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = gamma.state.ILB[4];

      stateStore.LAGE = gamma.state.LAGE;
      stateStore.PAGE = gamma.state.PAGE;
      
      DIRECT(CDTP,DF,stateStore.U,stateStore.V,stateStore.W);

      stackP.store(stateStore);
      
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
      //  ****  The positron carries a 'latent' energy of 1022 keV.
      DE = DE-TREV;
    }
  else
    {
      double WGHTA = gamma.state.WGHT;
      gamma.state.WGHT = gamma.state.WGHT*WFORCE;
      PANaR(gamma.state,mat.EABS[PEN_PHOTON],stackG,penRand);
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
      gamma.state.WGHT = WGHTA;
    }
  //  ****  Atomic relaxation after triplet production.
  if(ISA < 17)
    {
      int ks;
      double WGHTA = gamma.state.WGHT;
      gamma.state.WGHT = gamma.state.WGHT*WFORCE;
      RELAX(context.elements,mat,gamma.state,ICOL,1,
	    IZA,ISA,ks,stackE,stackG,penRand);      
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
      gamma.state.WGHT = WGHTA;
    }

  gamma.DEA = DE;
  DE = DE*WFORCE;
  if(LCOL){ gamma.state.E = 0.0;}
  return ICOL;
}


// Auxiliar functions

//  *********************************************************************
//                       SUBROUTINE GPHa
//  *********************************************************************
void GPHa(const pen_material& mat,
	  const pen_elementDataBase& elemDB,
	  const double E, const unsigned KE,
	  const double XEL,
	  double &ES,
	  int &IZZ,
	  int &ISH,
	  pen_rand& penRand)
{
//  Simulation of photoelectric absorption in material M.

//  Output arguments:
//    ES .... kinetic energy of the photoelectron.
//    IZZ ... atomic number of the atom where absorption has occurred.
//    ISH ... atomic electron shell that has been ionised.

//  NOTE: JUMP uses a photoelectric cross section that is slightly larger
//  than its 'true' value. To correct for this, the photon is allowed to
//  'survive' a photoelectric event. Survival of the photon is flagged by
//  setting IZZ=0, ISH=0, ES=0.0D0 (the energy E of the photon is kept
//  unaltered.

  
  double ACP[35], IP[35];

  //  ****  Partial attenuation coefficients.

  bool Eixir;
  int I, IU, IT;
  double DEE, PCSL;
  double PTOT = 0.0;
  for(int IEL = 0; IEL < mat.NELEM; IEL++)
  {
    IZZ = mat.IZ[IEL];
    //  ****  Binary search.
    I = elemDB.IPHF[IZZ-1]-1;
    IU = elemDB.IPHL[IZZ-1]-1;
    do
    {
      IT=(I+IU)/2;
      if(XEL > elemDB.EPH[IT])
      {
        I = IT;
      }
      else
      {
        IU = IT;
      }
    }while(IU-I > 1);
    
    IP[IEL] = I;
    DEE = elemDB.EPH[I+1]-elemDB.EPH[I];
    if(DEE > 1.0E-15)
    {
      PCSL = elemDB.XPH[I][0]+(elemDB.XPH[I+1][0]-elemDB.XPH[I][0])
	*(XEL-elemDB.EPH[I])/DEE;
    }
    else
    {
      PCSL = elemDB.XPH[I][0];
    }
    PTOT = PTOT+mat.STF[IEL]*exp(PCSL);
    ACP[IEL] = PTOT;
  }
  if(PTOT*mat.VMOL > mat.SGPH[KE])
  {
    printf("WARNING: SGPH is less than the actual mac.\n");
  }

  //  ****  Sample the active element.

  int IELAC;
  double TST = penRand.rand()*mat.SGPH[KE]/mat.VMOL;
  bool Return = true;
  for(int IEL = 0; IEL < mat.NELEM; IEL++)
  {
    if(ACP[IEL] > TST)
    {
      IELAC = IEL;
      IZZ = mat.IZ[IEL];
      Return = false;
      break;
    }
  }

  if(Return)
  {
    //  ****  Delta interaction. Introduced to correct for the use of an
    //        upper bound of the photoelectric attenuation coefficient.
    IZZ = 0;  // Flags delta interactions.
    ISH = 0;
    ES = 0.0;
    return;
  }
  
  //  ****  Selection of the active shell.
  I = IP[IELAC];
  DEE = elemDB.EPH[I+1]-elemDB.EPH[I];
  double PIS = 0.0;
  Eixir = false;
  int J;
  if(DEE > 1.0E-15)
  {
    PTOT = exp(elemDB.XPH[I][0]+(elemDB.XPH[I+1][0]-elemDB.XPH[I][0])
	       *(XEL-elemDB.EPH[I])/DEE);
    TST = penRand.rand()*PTOT;
    for(int IS = 0; IS < elemDB.NPHS[IZZ-1]; IS++)
    {
      J = (IS+1);
      PCSL = elemDB.XPH[I][J]+(elemDB.XPH[I+1][J]-elemDB.XPH[I][J])
	*(XEL-elemDB.EPH[I])/DEE;
      PIS = PIS+exp(PCSL);
      if(PIS > TST)
      {
        ISH = IS+1;
        Eixir = true;
        break;
      }
    }
  }
  else
  {
    PTOT = exp(elemDB.XPH[I][0]);
    TST = penRand.rand()*PTOT;
    for(int IS = 0; IS < elemDB.NPHS[IZZ-1]; IS++)
    {
      PIS = PIS+exp(elemDB.XPH[I][IS+1]);
      if(PIS > TST)
      {
        ISH = IS+1;
        Eixir = true;
        break;
      }
    }
  }
  
  if(!Eixir){ ISH = 17;}

  //  ****  Photoelectron emission.

  double EBB;
  if(ISH < 17)
  {
    EBB = elemDB.EB[IZZ-1][ISH-1];
    if(EBB > mat.ECUTR)
    {
      ES = E-EBB;
    }
    else
    {
      ES = E;
      ISH = 17;
    }
  }
  else
  {
    ES = E;
  }
}

//  *********************************************************************
//                       SUBROUTINE GPHa0
//  *********************************************************************
void GPHa0(pen_elementDataBase& elemDB)
{
  //  This subroutine sets all variables in common /CGPH00/ to zero.
  //  It has to be invoked before reading the first material definition
  //  file.

  
  for(unsigned int I = 0; I < elemDB.nElements; I++)
  {
    elemDB.NPHS[I] = 0;
    elemDB.IPHF[I] = constants::NTP;
    elemDB.IPHL[I] = constants::NTP;
  }

  for(unsigned int I = 0; I < constants::NTP; I++)
  {
    elemDB.EPH[I] = 0.0;
    for(unsigned int J = 0; J < 17; J++)
    {
      elemDB.XPH[I][J] = 1.0E-35;
    }
  }
  elemDB.NCUR = 0;
}

//  *********************************************************************
//                       SUBROUTINE GPHaT
//  *********************************************************************
void GPHaT(double &E, double &XS, const pen_material& mat, const pen_elementDataBase& elemDB)
{
  //  Delivers the photoelectric cross section XS (in cm**2) for photons of
  //  energy E in material M.


  double XEL = log(E);
  XS = 0.0;
  double DEE, PCSL;
  int IZZ, I, IU, IT;
  for(int IEL = 0; IEL < mat.NELEM; IEL++)
  {
    IZZ = mat.IZ[IEL];
    I = elemDB.IPHF[IZZ-1]-1;
    IU = elemDB.IPHL[IZZ-1]-1;

    while(IU-I > 1)
    {
      IT = (I+IU)/2;
      if(XEL > elemDB.EPH[IT])
      {
        I = IT;
      }
      else
      {
        IU = IT;
      }
    }
    DEE = elemDB.EPH[I+1]-elemDB.EPH[I];
    if(DEE > 1.0E-15)
    {
      PCSL = elemDB.XPH[I][0]+(elemDB.XPH[I+1][0]-elemDB.XPH[I][0])*(XEL-elemDB.EPH[I])/DEE;
    }
    else
    {
      PCSL = elemDB.XPH[I][0];
    }
    XS = XS+mat.STF[IEL+1]*exp(PCSL);
  }
}

//  *********************************************************************
//                       SUBROUTINE GPPa
//  *********************************************************************
void GPPa(const pen_material& mat,
	  const double E,
	  const unsigned KE,
	  const double XEK,
	  double &EE,
	  double &CDTE,
	  double &EP,
	  double &CDTP,
	  int &IZZ,
	  int &ISH,
	  pen_rand& penRand)
{
  //  Random sampling of electron-positron pair and triplet production by
  //  photons. Bethe-Heitler differential cross section.

  //  Output values:
  //    EE .....  kinetic energy of the electron.
  //    CDTE ...  polar direction cosine of the electron.
  //    EP .....  kinetic energy of the positron.
  //    CDTP ...  polar direction cosine of the positron.
  //    IZZ  ... atomic number of the atom where absorption has occurred.
  //    ISH .... atomic electron shell that has been ionised.


  const double TREV = 2.0*constants::REV;
  const double REV = constants::REV;

  double EKI = REV/E;
  double EPS, ALZ, T, F00, G0, BMIN, G1, G2, G1MIN, G2MIN, A1, P1, XR, RU2M1;
  if(E < 1.1E6)
  {
    EPS = EKI+(1.0-2.0*EKI)*penRand.rand();
  }
  else
  {
    //  ****  Low-energy and Coulomb corrections.
    ALZ = mat.ZEQPP/constants::SL;
    T = sqrt(2.0*EKI);
    F00 = (-1.774-1.210E1*ALZ+1.118E1*ALZ*ALZ)*T
      +(8.523+7.326E1*ALZ-4.441E1*ALZ*ALZ)*pow(T,2)
      -(1.352E1+1.211E2*ALZ-9.641E1*ALZ*ALZ)*pow(T,3)
      +(8.946+6.205E1*ALZ-6.341E1*ALZ*ALZ)*pow(T,4);
    G0 = mat.F0[1]+F00;
    BMIN = 4.0*EKI/mat.BCB;
    SCHIFF(BMIN,G1,G2);
    G1MIN = G1+G0;
    G2MIN = G2+G0;
    XR = 0.5-EKI;
    A1 = 6.666666666666666E-1*G1MIN*pow(XR,2);
    P1 = A1/(A1+G2MIN);
    //  ****  Random sampling of EPS.
    while(true)
    {
      if(penRand.rand() > P1)
      {
	EPS = EKI+2.0*XR*penRand.rand();
	double B = EKI/(mat.BCB*EPS*(1.0-EPS));
	SCHIFF(B,G1,G2);
	G2 += G0;
	if(G2 < 0.0){ G2 = 0.0;}
	if(penRand.rand()*G2MIN > G2){ continue;}
      }
      else{
        RU2M1 = 2.0*penRand.rand()-1.0;
        if(RU2M1 < 0.0)
        {
          EPS = 0.5-XR*pow(fabs(RU2M1),3.333333333333333E-1);
        }
        else
        {
          EPS = 0.5+XR*pow(RU2M1,3.333333333333333E-1);
        }
        double B = EKI/(mat.BCB*EPS*(1.0-EPS));
        SCHIFF(B,G1,G2);
        G1 += G0;
        if(G1 < 0.0){ G1 = 0.0;}
        if(penRand.rand()*G1MIN > G1){ continue;}	
      }
      break;
    }
  }
  //  ****  Electron.
  EE = EPS*E-REV;
  CDTE = 2.0*penRand.rand()-1.0;
  A1 = EE+REV;
  double A2 = sqrt(EE*(EE+TREV));
  CDTE = (CDTE*A1+A2)/(A1+CDTE*A2);
  //  ****  Positron.
  EP = (1.0-EPS)*E-REV;
  CDTP = 2.0*penRand.rand()-1.0;
  A1 = EP+REV;
  A2 = sqrt(EP*(EP+TREV));
  CDTP = (CDTP*A1+A2)/(A1+CDTP*A2);

  //  ****  Triplet production.

  int ISHELL;
  double TRIPL = mat.TRIP[KE]+(mat.TRIP[KE+1]-mat.TRIP[KE])*XEK;
  IZZ = 0;
  ISH = 30;
  if(TRIPL < 1.0E-5){ return;}
  if(penRand.rand() > TRIPL){ return;}
  double TST = penRand.rand();
  //  ****  Binary search.
  if(TST < mat.PTRSH[0])
  {
    ISHELL = 0;
  }
  else
  {
    ISHELL = 0;
    int JO = mat.NOSCCO;
    do
    {
      int I = (ISHELL+JO)/2;
      if(TST > mat.PTRSH[I])
      {
        ISHELL = I;
      }
      else
      {
        JO = I;
      }
    }while(JO-ISHELL > 1);
    ISHELL = ISHELL+1;
  }
  IZZ = mat.KZCO[ISHELL];
  ISH = mat.KSCO[ISHELL];
}

//  *********************************************************************
//                       SUBROUTINE GRAa
//  *********************************************************************
void GRAa(const pen_material& mat,
	  const double E,
	  const unsigned KE,
	  const double XEL,
	  double &CDT,
	  int &IEFF,
	  pen_rand& penRand)
{
  //  Random sampling of coherent (Rayleigh) scattering.


  const double RREV = 1.0/constants::REV;

  //  ****  Binary search.

  int II = mat.IED[KE]-1;
  int IU = mat.IEU[KE]-1;
  int IT;
  while(IU-II > 1)
  {
    IT = (II+IU)/2;
    if(XEL > mat.ERA[IT])
    {
      II = IT;
    }
    else
    {
      IU = IT;
    }
  }
  double XSE = exp(mat.XSRA[II]+(mat.XSRA[II+1]-mat.XSRA[II])*(XEL-mat.ERA[II])/(mat.ERA[II+1]-mat.ERA[II]));
  if(penRand.rand()*mat.SGRA[KE] > XSE)
  {
    IEFF = 0;
    CDT = 1.0;
    return;
  }

  IEFF = 1;
  double QMAX = 2.0*E*RREV;
  if(QMAX < 1.0E-10)
  {
    double G;
    do
    {
      CDT = 1.0-2.0*penRand.rand();
      G = 0.5*(1.0+CDT*CDT);
    }while(penRand.rand() > G);
    return;
  }
  double Q2MAX = QMAX*QMAX;
  if(Q2MAX > mat.QRA[mat.NP_RSC-1]){ Q2MAX = mat.QRA[mat.NP_RSC-1];}

  bool Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    double RU = penRand.rand()*mat.PMAX[KE+1];
  
      //  ****  Selection of the interval
      //        (binary search within pre-calculated limits).

    int ITN = (int)(RU*mat.NPM1_RSC);
    int I = mat.ITLRA[ITN]-1;
    int J = mat.ITURA[ITN]-1;
    if(J-I < 2){}
    else
    {
      while(J-I > 1)
      {
        int K = (I+J)/2;
        if(RU > mat.PRA[K])
        {
          I = K;
        }
        else
        {
          J = K;
        }
      }
    
      //  ****  Sampling from the rational inverse cumulative distribution.
    }

    double XX;
    double RR = RU-mat.PRA[I];
    if(RR > 1.0E-16)
    {
      double D = mat.DPRA[I];
      XX = mat.QRA[I]+((1.0+mat.ARA[I]+mat.BRA[I])*D*RR/(D*D+(mat.ARA[I]*D+mat.BRA[I]*RR)*RR))*(mat.QRA[I+1]-mat.QRA[I]);
    }
    else
    {
      XX = mat.QRA[I];
    }
    if(XX > Q2MAX){ Eixir = false; continue;}
    CDT = 1.0-2.0*XX/Q2MAX;
    //  ****  Rejection.
    double G = 0.5*(1.0+CDT*CDT);
    if(penRand.rand() > G){ Eixir = false; continue;}
  }
}

//  *********************************************************************
//                       SUBROUTINE GCOa
//  *********************************************************************
void GCOa(const pen_material& mat,
	  double E,
	  double &DE,
	  double &EP,
	  double &CDT,
	  double &ES,
	  double &CDTS,
	  int &IZZ,
	  int &ISH,
	  pen_rand& penRand)
{
  //  Random sampling of incoherent (Compton) scattering of photons. Relat-
  //  ivistic impulse approximation with analytical one-electron Compton
  //  profiles.

  //  Input arguments:
  //    E ..... incident photon energy (eV).
  //    M ..... material where photons propagate.
  //  Output argument:
  //    DE .... energy loss (eV).
  //    EP .... energy of the scattered photon (eV).
  //    CDT ... cosine of the polar scattering angle.
  //    ES .... energy of the emitted electron (eV).
  //    CDTS .. polar cosine of direction of the electron.
  //    IZZ ... atomic number of the atom where scattering has occurred.
  //    ISH ... atomic electron shell that has been ionised.


  const double RREV = 1.0/constants::REV;
  const double D2 = 1.4142135623731;
  const double D1 = 1.0/D2;
  const double D12 = 0.5;
  
  double RN[constants::NOCO], PAC[constants::NOCO];
  double TAU, TST;
  int    ISHELL;
  
  double EK = E*RREV;
  double EK2 = EK+EK+1.0;
  double EKS = EK*EK;
  double EK1 = EKS-EK2-1.0;
  double TAUMIN = 1.0/EK2;
  double TAUM2 = TAUMIN*TAUMIN;
  double A1 = log(EK2);
  double A2 = A1+2.0*EK*(1.0+EK)*TAUM2;
  bool Energia_Major_5MeV = false;
  if(E > 5.0E6){ Energia_Major_5MeV = true;}
  else
  {

    //  ****  Incoherent scattering function for theta=PI.
      
    double S0 = 0.0;
    double AUX, PZOMC, RNI;
    for(int I = 0; I < mat.NOSCCO; I++)
    {
      if(mat.UICO[I] < E)
      {
        AUX = E*(E-mat.UICO[I])*2.0;
        PZOMC = mat.FJ0[I]*(AUX-constants::REV*mat.UICO[I])/(constants::REV*sqrt(AUX+AUX+pow(mat.UICO[I],2)));
        if(PZOMC > 0.0)
        {
          RNI = 1.0-0.5*exp(D12-pow(D1+D2*PZOMC,2));
        }
        else
        {
          RNI = 0.5*exp(D12-pow(D1-D2*PZOMC,2));
        }
        S0 = S0+mat.FCO[I]*RNI;
      }
    }
      
    //  ****  Sampling tau.
      
    double CDT1, S, A;
    bool Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      if(penRand.rand()*A2 < A1)
      {
        TAU = pow(TAUMIN,penRand.rand());
      }
      else
      {
        TAU = sqrt(1.0+penRand.rand()*(TAUM2-1.0));
      }
      CDT1 = (1.0-TAU)/(EK*TAU);
      //  ****  Incoherent scattering function.
      S = 0.0;
      for(int I = 0; I < mat.NOSCCO; I++)
      {
        if(mat.UICO[I] < E)
        {
          AUX = E*(E-mat.UICO[I])*CDT1;
          PZOMC = mat.FJ0[I]*(AUX-constants::REV*mat.UICO[I])/(constants::REV*sqrt(AUX+AUX+pow(mat.UICO[I],2)));
          if(PZOMC > 0.0)
          {
            RN[I] = 1.0-0.5*exp(D12-pow(D1+D2*PZOMC,2));
          }
          else
          {
            RN[I] = 0.5*exp(D12-pow(D1-D2*PZOMC,2));
          }
          S = S+mat.FCO[I]*RN[I];
          PAC[I] = S;
        }
        else
        {
          PAC[I] = S;
        }
      }
      //  ****  Rejection function.
      TST = S*(1.0+TAU*(EK1+TAU*(EK2+TAU*EKS)))/(EKS*TAU*(1.0+TAU*TAU));
      if(penRand.rand()*S0 > TST){ Eixir = false; continue;}
    }
    CDT = 1.0-CDT1;
      
    //  ****  Target electron shell.
      
    int JO;
    double XQC, FPZMAX, AF, FPZ, T, B1, B2;
    Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      TST = S*penRand.rand();
      //  ****  Binary search.
      if(TST < PAC[0])
      {
        ISHELL = 0;
      }
      else
      {
        ISHELL = 0;
        JO = mat.NOSCCO;
        int I;
        while(JO-ISHELL > 1)
        {
          I = (ISHELL+JO)/2;
          if(TST > PAC[I])
          {
            ISHELL = I;
          }
          else
          {
            JO = I;
          }
        }
        ISHELL = ISHELL+1;
      }
    
      //  ****  Projected momentum of the target electron.
    
      A = penRand.rand()*RN[ISHELL];
      if(A < 0.5)
      {
        PZOMC = (D1-sqrt(D12-log(A+A)))/(D2*mat.FJ0[ISHELL]);
      }
      else
      {
        PZOMC = (sqrt(D12-log(2.0-A-A))-D1)/(D2*mat.FJ0[ISHELL]);
      }
      if(PZOMC < -1.0){ Eixir = false; continue;}
    
      //  ****  F(EP) rejection.
    
      XQC = 1.0+TAU*(TAU-2.0*CDT);
      AF = sqrt(XQC)*(1.0+TAU*(TAU-CDT)/XQC);
      if(AF > 0.0)
      {
        FPZMAX = 1.0+AF*0.2;
      }
      else
      {
        FPZMAX = 1.0-AF*0.2;
      }
    
      if(PZOMC > 0.2){ FPZ = 0.2;}
      else{ FPZ = PZOMC;}
    
      if(FPZ < -0.2){ FPZ = -0.2;}
    
      FPZ = 1.0+AF*FPZ;
      if(penRand.rand()*FPZMAX > FPZ){ Eixir = false; continue;}
    }
    //  ****  Energy of the scattered photon.
      
    T = pow(PZOMC,2);
    B1 = 1.0-T*TAU*TAU;
    B2 = 1.0-T*TAU*CDT;
    if(PZOMC > 0.0)
    {
      EP = E*(TAU/B1)*(B2+sqrt(fabs(B2*B2-B1*(1.0-T))));
    }
    else
    {
      EP = E*(TAU/B1)*(B2-sqrt(fabs(B2*B2-B1*(1.0-T))));
    }
  }

  if(Energia_Major_5MeV)
  {
    //  ****  No Doppler broadening for E greater than 5 MeV.
    bool Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      if(penRand.rand()*A2 < A1)
      {
        TAU = pow(TAUMIN,penRand.rand());
      }
      else
      {
        TAU = sqrt(1.0+penRand.rand()*(TAUM2-1.0));
      }
      //  ****  Rejection function.
      TST = (1.0+TAU*(EK1+TAU*(EK2+TAU*EKS)))/(EKS*TAU*(1.0+TAU*TAU));
      if(penRand.rand() > TST){ Eixir = false; continue;}
      EP = TAU*E;
      CDT = 1.0-(1.0-TAU)/(EK*TAU);

      //  ****  Target electron shell.

      TST = penRand.rand();
      //  ****  Binary search.
      if(TST < mat.PTRSH[0])
      {
        ISHELL = 0;
      }
      else
      {
        ISHELL = 0;
        int JO = mat.NOSCCO;
        while(JO-ISHELL > 1)
        {
          int I = (ISHELL+JO)/2;
          if(TST > mat.PTRSH[I])
          {
            ISHELL = I;
          }
          else
          {
            JO = I;
          }
        }
        ISHELL = ISHELL+1;
      }
      if(EP > E-mat.UICO[ISHELL]){ Eixir = false; continue;}
    }
  }
  DE = E-EP;
  if(mat.KSCO[ISHELL] < 17)
  {
    if(mat.UICO[ISHELL] > mat.ECUTR)
    {
      ES = DE-mat.UICO[ISHELL];
    }
    else
    {
      ES = DE;
    }
  }
  else
  {
    ES = DE;
  }
  
  double Q2 = E*E+EP*(EP-2.0*E*CDT);
  if(Q2 > 1.0E-12)
  {
    CDTS = (E-EP*CDT)/sqrt(Q2);
  }
  else
  {
    CDTS = 1.0;
  }

  IZZ = mat.KZCO[ISHELL];
  ISH = mat.KSCO[ISHELL];
}

//  *********************************************************************
//                        FUNCTION GCOaD
//  *********************************************************************
double GCOaD(double CDT, void* arg)
{
  //  Single differential cross section for photon Compton scattering by
  //  electrons in the IO-th shell, differential in the direction of the
  //  scattered photon only. Evaluated from the incoherent scattering
  //  function.

  //  The energy E of the primary photon is entered through common CGCO00.
  //  The output value GCOaD is the DCS per electron in units of PIELR2.


  //Extract argument data
  GCO00* co00 = (GCO00*)(arg);

  double EE = co00->EE; 
  int IOSC = co00->IOSC; 
  const pen_material& mat = *(co00->pmat);
  
  const double RREV = 1.0/constants::REV;
  const double D2 = 1.4142135623731;
  const double D1 = 1.0/D2;
  const double D12 = 0.5;

  double GCOaD_RETURN;
  if(EE < mat.UICO[IOSC-1])
  {
    GCOaD_RETURN = 0.0;
    return GCOaD_RETURN;
  }
  //  ****  Energy of the Compton line.
  double CDT1 = 1.0-CDT;
  double EOEC = 1.0+(EE*RREV)*CDT1;
  double ECOE = 1.0/EOEC;
  //  ****  Klein-Nishina X-factor.
  double XKN = EOEC+ECOE-1.0+CDT*CDT;
  //  ****  Incoherent scattering function (analytical profile).
  double AUX = EE*(EE-mat.UICO[IOSC-1])*CDT1;
  double PIMAX = (AUX-constants::REV*mat.UICO[IOSC-1])/(constants::REV*sqrt(AUX+AUX+pow(mat.UICO[IOSC-1],2)));
  double SIA;
  if(PIMAX > 0.0)
  {
    SIA = 1.0-0.5*exp(D12-pow(D1+D2*mat.FJ0[IOSC-1]*PIMAX,2));
  }
  else
  {
    SIA = 0.5*exp(D12-pow(D1-D2*mat.FJ0[IOSC-1]*PIMAX,2));
  }
  //  ****  1st order correction, integral of Pz times the Compton profile.
  //        Calculated approximately using a free-electron gas profile.
  double PF = 3.0/(4.0*mat.FJ0[IOSC-1]);
  double QCOE2, P2, DSPZ;
  if(fabs(PIMAX) < PF)
  {
    QCOE2 = 1.0+pow(ECOE,2)-2.0*ECOE*CDT;
    P2 = pow(PIMAX,2);
    DSPZ = sqrt(QCOE2)*(1.0+ECOE*(ECOE-CDT)/QCOE2)*mat.FJ0[IOSC-1]*0.25*(2*P2-pow(P2,2)/pow(PF,2)-pow(PF,2));
    if(DSPZ > -SIA){ SIA = SIA+DSPZ;}
    else{ SIA = 0.0;}      
  }
  //  ****  Differential cross section (per electron, in units of PIELR2).
  GCOaD_RETURN = pow(ECOE,2)*XKN*SIA;
  return GCOaD_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE GCOaT
//  *********************************************************************
void GCOaT(const pen_material& mat, const double E, double &CS)
{
  //  Total cross section for incoherent (Compton) scattering. Relativistic
  //  Impulse approximation with analytical Compton profiles.

  //  Input arguments:
  //    E ........ photon energy (eV).
  //    M ........ material where photons propagate.
  //  Output argument:
  //    CS ....... incoherent total cross section (cm**2/molecule).


  //double GCOaD(double);

  //Declare the structure to pass as SUMGA argument
  GCO00 gco00;
  gco00.pmat = &mat;
  gco00.EE = E;
  
  const double RREV = 1.0/constants::REV;
  const double PIELR2 = constants::PI*constants::ELRAD*constants::ELRAD;

  double SXCO[constants::NOCO];

  double EK, EKS, EK2, EK1, T0, CSL, TAU, CSKN, CSU;
  CS = 0.0;
  if(E < 5.0E6)
  {
    for(int IO = 0; IO < mat.NOSCCO; IO++)
    {
      gco00.IOSC = IO+1;
      SXCO[IO] = mat.FCO[IO]*PIELR2*SUMGA(GCOaD,&gco00,-1.0,1.0,1.0E-6);
      CS = CS+SXCO[IO];
    }
  }
  else
  {
    //  ****  Klein-Nishina total cross section.
    EK = E*RREV;
    EKS = EK*EK;
    EK2 = 1.0+EK+EK;
    EK1 = EKS-EK2-1.0;
    T0 = 1.0/(1.0+EK+EK);
    CSL = 0.5*EKS*T0*T0+EK2*T0+EK1*log(T0)-1.0/T0;
    for(int IO = 0; IO < mat.NOSCCO; IO++)
    {
      TAU = (E-mat.UICO[IO])/E;
      if(TAU < T0)
      {
        CSKN = 0.0;
      }
      else
      {
        CSU = 0.5*EKS*TAU*TAU+EK2*TAU+EK1*log(TAU)-1.0/TAU;
        CSKN = PIELR2*(CSU-CSL)/(EK*EKS);
      }
      SXCO[IO] = mat.FCO[IO]*CSKN;
      CS = CS+SXCO[IO];        
    }
  }
}

//  *********************************************************************
//                       SUBROUTINE SAUTER
//  *********************************************************************
void SAUTER(const double ES,
      double &CDTS,
      double &DFS,
      pen_rand& penRand,
      pen_state_gPol& gPolstate)
{
  //  Random sampling of the initial direction of photoelectrons from the
  //  Sauter distribution.

  if(ES > 1.0E9)
  {
    CDTS = 1.0;
    DFS = 0.0;
    return;
  }
  double GAM = 1.0+ES/constants::REV;
  double GAM2 = GAM*GAM;
  double BETA = sqrt((GAM2-1.0)/GAM2);
  double A0 = 0.5 * GAM * (GAM - 1.0);
//
  double AC = 1.0/BETA-1.0;
  double A1 = A0 * (GAM - 2.0) * BETA;
  double A2 = AC+2.0;
  double GTMAX = 2.0*(A1+1.0/AC);
  double TSAM;

  double GTR;
  do
  {
    double RU = penRand.rand();
    TSAM = 2.0*AC*(2.0*RU+A2*sqrt(RU))/(A2*A2-4.0*RU);
    GTR = (2.0-TSAM)*(A1+1.0/(AC+TSAM));
  }while(penRand.rand()*GTMAX > GTR);

  CDTS = 1.0-TSAM;
//  ****  Azimuthal angle.
  if (gPolstate.IPOL == 0)    // Unpolarised photons
    {
      DFS = constants::TWOPI * penRand.rand();
    }
  else        // Polarised photons.
    {
      double B0 = A0 * (1.0E0 - BETA * CDTS);
      double B1 = 1.0E0 + B0 * (GAM - 2.0E0);
      double B2 = 1.0E0 - B0;
      double TF = atan2 (gPolstate.SP1, gPolstate.SP3);
      double TERM = B2 * (gPolstate.SP3 * cos (TF) + gPolstate.SP1 * sin (TF));
      double FM = std::max(B1 + TERM,B1 - TERM);

      double F;
      do{
	DFS = constants::TWOPI * penRand.rand();
	double TDFS = DFS + DFS;
	F = B1 + B2 * (gPolstate.SP3 * cos (TDFS) + gPolstate.SP1 * sin (TDFS));
      }while(penRand.rand() * FM > F);
    }
}

//  *********************************************************************
//                       SUBROUTINE SCHIFF
//  *********************************************************************
void SCHIFF(double &B, double &G1, double &G2)
{
  //  Screening functions F1(B) and F2(B) in the Bethe-Heitler differential
  //  cross section for pair production.

  double B2 = B*B;
  double F1 = 2.0-2.0*log(1.0+B2);
  double F2 = F1-6.666666666666666E-1;
  if(B < 1.0E-10)
  {
    F1 = F1-constants::TWOPI*B;
  }
  else
  {
    double A0 = 4.0*B*atan2(1.0,B);
    F1 = F1-A0;
    F2 = F2+2.0*B2*(4.0-A0-3.0*log((1.0+B2)/B2));
  }
  G1 = 0.5*(3.0*F1-F2);
  G2 = 0.25*(3.0*F1+F2);
}
