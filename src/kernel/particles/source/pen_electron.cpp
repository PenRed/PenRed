
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


#include "pen_electron.hh"

// Electrons class functions

const double pen_betaE::REV = 5.10998928e5;  // Electron rest energy (eV)
pen_betaE::pen_betaE(const pen_context& contextIn,
		     pen_particleStack<pen_particleState>& stackEin,
		     pen_particleStack<pen_state_gPol>& stackGin)
  
  : abc_particle(contextIn,PEN_ELECTRON,4,0.0E0,stackEin),
    HinelastCol(stackEin),
    Hbremsstrahlung(stackGin),
    HinnerShell(stackEin,stackGin),
    MHINGE(0),
    DESOFT(0.0),
    SSOFT(0.0)
    
{

  //Check ids of interactions
  int auxID = HelasticCol.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Electron elastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = HinelastCol.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Electron inelastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = Hbremsstrahlung.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Electron inelastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = HinnerShell.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Electron inelastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  
  HelasticCol.init(context);
  HinelastCol.init(context);
  Hbremsstrahlung.init(context);
  HinnerShell.init(context);  
}

void pen_betaE::JUMP(double &DS,
		     pen_rand& penRand,
		     const double DSMAX)
{
  //  ************  Electrons.

  //Get material
  const pen_material& mat = *pmat;
  
  if(MHINGE == 1) // Jump from hinge to hard event
    {
      if(EENDSTEP < ELAST1)
	{
	  XEL = log(EENDSTEP);
	  XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
	  KE = (int)XE;
	  XEK = XE-(double)KE;
	  IMFP(mat);
	  ELAST1 = EENDSTEP;
	}
      DS = DSR;
      return;
    }

  // Jump from hard event to hinge
  
  if(state.E < ELAST2)
    {
      XEL = log(state.E);
      XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
      KE = (int)XE;
      XEK = XE-(double)KE;
      IMFP(mat);
      hingeUpdate(mat);      
      ELAST2 = state.E;
      ELAST1 = state.E;
    }

  //  ****  Inverse hard mean free path (interaction probability per unit
  //        path length).

  ST = P[0];
  ST += P[1];
  ST += P[2];
  ST += P[3];
  
  double DSMAXP = DSMAX;

  //  ****  Soft stopping interactions.
  //        KSOFTI=1, soft stopping is active,
  //        KSOFTI=0, soft stopping is not active.
  if(W1 > 1.0E-20)
    {
      KSOFTI = 1;
      //  ****  The maximum step length, DSMAXP, is determined in terms of the
      //        input DSMAX value (which is specified by the user) and the mean
      //        free path for hard interactions (1/ST).
      double DSMC = 4.0/ST;
      if(DSMAXP > DSMC)
	{
	  DSMAXP = DSMC;
	}
      else if(DSMAXP < 1.0E-8)
	{
	  DSMAXP = DSMC;
	}

      //  ****  Upper bound for the interaction probability along the step
      //        (including soft energy straggling).

      double EDE0 = W1*DSMAXP;
      double VDE0 = W2*DSMAXP;

      double FSEDE = 1.0-mat.DW1EL[KE]*EDE0;
      if(FSEDE < 0.75){ FSEDE = 0.75;}

      double FSVDE = 1.0-mat.DW2EL[KE]*EDE0;
      if(FSVDE < 0.75){ FSVDE = 0.75;}
    
      double EDEM = EDE0*FSEDE;
      double VDEM = VDE0*FSVDE;
      double W21 = VDEM/EDEM;
      double ELOWER;
      if(EDEM > 9.0*W21)
	{
	  ELOWER = state.E-(EDEM+3.0*sqrt(VDEM));
	  if(ELOWER < context.grid.EMIN){ ELOWER = context.grid.EMIN;}
	}
      else if(EDEM > 3.0*W21)
	{
	  ELOWER = state.E-(EDEM+sqrt(3.0*VDEM));
	  if(ELOWER < context.grid.EMIN){ ELOWER = context.grid.EMIN;}
	}
      else
	{
	  ELOWER = state.E-1.5*(EDEM+W21);
	  if(ELOWER < context.grid.EMIN){ ELOWER = context.grid.EMIN;}
	}
      double XE1 = (log(ELOWER)-context.grid.DLEMP1)*context.grid.DLFC;
      int KE1 = (int)XE1;
      double XEK1 = XE1-(double)KE1;
      double STLWR = exp(mat.SETOT[KE1]+mat.DSETOT[KE1]*XEK1);
      if(ST < STLWR){ ST = STLWR;}
    }
  else
    {
      KSOFTI = 0;
      DESOFT = 0.0;
      SSOFT = 0.0;
    }

  //  ****  Soft elastic scattering.
  //        KSOFTE=1, soft scattering is active,
  //        KSOFTE=0, soft scattering is not active.

  if(T1 > 1.0E-20)
    {
      KSOFTE = 1;
    }
  else
    {
      KSOFTE = 0;
    }

  //  ****  Delta interactions.
  //        KDELTA=0, a hard interaction follows,
  //        KDELTA=1, a delta interaction follows.

  DST = -log(penRand.rand())/ST;
  
  if(DST < DSMAXP)
    {
      KDELTA = 0;
    }
  else
    {
      DST = DSMAXP;
      KDELTA = 1;
    }

  if(KSOFTI == 0 && KSOFTE == 0)
    {
      MHINGE = 1;
      DS = DST;
    }
  else
    {
      DS = DST*penRand.rand();
      
      DSR = DST-DS;
      if(KSOFTI == 1)
	{
	  if(DST < 1.0E-8)
	    {
	      SSOFT = W1;
	      DESOFT = SSOFT*DST;
	    }
	  else
	    {
	      double EDE0 = W1*DST;
	      double VDE0 = W2*DST;
	      double FSEDE = 1.0-mat.DW1EL[KE]*EDE0;
	      if(FSEDE < 0.75){ FSEDE = 0.75;}

	      double FSVDE = 1.0-mat.DW2EL[KE]*EDE0;
	      if(FSVDE < 0.75){ FSVDE = 0.75;}
	      double EDE = EDE0*FSEDE;
	      double VDE = VDE0*FSVDE;
	      //  ****  Generation of random values DE with mean EDE and variance VDE.
	      double SIGMA = sqrt(VDE);
	      
	      if(SIGMA < 0.333333333*EDE)
		{
		  //  ****  Truncated Gaussian distribution.
		  DESOFT = EDE+RNDG3(context.rndg3,penRand)*SIGMA;
		}
	      else
		{
		  double RU = penRand.rand();
		  double EDE2 = EDE*EDE;
		  double VDE3 = 3.0*VDE;
		  if(EDE2 < VDE3)
		    {
		      double PNULL = (VDE3-EDE2)/(VDE3+3.0*EDE2);
		      if(RU < PNULL)
			{
			  DESOFT = 0.0;
			  SSOFT = 0.0;

			  if(KSOFTE == 0)
			    {
			      MHINGE = 1;
			      DS = DST;
			    }
			  else
			    {
			      KSOFTI = 0;
			    }
			  return;
			}
		      else
			{
			  //  ****  Uniform distribution.
			  DESOFT = 1.5*(EDE+VDE/EDE)*(RU-PNULL)/(1.0-PNULL);
			}
		    }
		  else
		    {
		      DESOFT = EDE+(2.0*RU-1.0)*sqrt(VDE3);
		    }
		}
	      SSOFT = DESOFT/DST;
	    }
	}
    }
}

void pen_betaE::JUMPF(double &DS,
		      pen_rand& penRand,
		      const double DSMAX)
{
  //  ************  Electrons

  //Get material
  const pen_material& mat = *pmat;
  
  if(MHINGE == 1)
    {
      if(EENDSTEP < ELAST1)
	{
	  XEL = log(EENDSTEP);
	  XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
	  KE = (int)XE;
	  XEK = XE-(double)KE;
	  IMFP(mat);
	  for(unsigned KCOL = 0; KCOL < interactions; KCOL++)
	    {
	      POR[KCOL] = P[KCOL];  // Int. forcing modifies the IMFPs.
	    }
	  ELAST1 = EENDSTEP;
	}
      DS = DSR;
      return;
    }

  if(state.E < ELAST2)
    {
      XEL = log(state.E);
      XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
      KE = (int)XE;
      XEK = XE-(double)KE;
      IMFP(mat);
      hingeUpdate(mat);
      for(unsigned KCOL = 0; KCOL < interactions; KCOL++)
	{
	  POR[KCOL] = P[KCOL];  // Int. forcing modifies the IMFPs.
	}
      ELAST2 = state.E;
      ELAST1 = state.E;
    }
  //  ****  Interaction forcing.
  for(unsigned KCOL = 0; KCOL < interactions; KCOL++)
    {
      if(context.FORCE[state.IBODY][PEN_ELECTRON][KCOL] > 1.0 && P[KCOL] > 1.0E-16)
	{
	  P0[KCOL] = POR[KCOL];
	  P[KCOL] = POR[KCOL]*context.FORCE[state.IBODY][PEN_ELECTRON][KCOL];
	  LFORC[KCOL] = true;
	}
      else
	{
	  LFORC[KCOL] = false;
	}
    }

  //  ****  Inverse hard mean free path (interaction probability per unit
  //        path length).

  ST = P[0];
  ST += P[1];
  ST += P[2];
  ST += P[3];
  
  double DSMAXP = DSMAX;

  //  ****  Soft stopping interactions.
  //        KSOFTI=1, soft stopping is active,
  //        KSOFTI=0, soft stopping is not active.

  if(W1 > 1.0E-20)
    {
      KSOFTI = 1;
      //  ****  The maximum step length, DSMAXP, is determined in terms of the
      //        input DSMAX value (which is specified by the user) and the mean
      //        free path for hard interactions (1/ST).
      double DSMC = 4.0/ST;
      if(DSMAXP > DSMC)
	{
	  DSMAXP = DSMC;
	}
      else if(DSMAXP < 1.0E-8)
	{
	  DSMAXP = DSMC;
	}

      //  ****  Upper bound for the interaction probability along the step
      //        (including soft energy straggling).

      double EDE0 = W1*DSMAXP;
      double VDE0 = W2*DSMAXP;

      double FSEDE = 1.0-mat.DW1EL[KE]*EDE0;
      if(FSEDE < 0.75){ FSEDE = 0.75;}

      double FSVDE = 1.0-mat.DW2EL[KE]*EDE0;
      if(FSVDE < 0.75){ FSVDE = 0.75;}
	  
      double EDEM = EDE0*FSEDE;
      double VDEM = VDE0*FSVDE;
      double W21 = VDEM/EDEM;
      double ELOWER;
      if(EDEM > 9.0*W21)
	{
	  ELOWER = state.E-(EDEM+3.0*sqrt(VDEM));
	  if(ELOWER < context.grid.EMIN){ ELOWER = context.grid.EMIN;}
	}
      else if(EDEM > 3.0*W21)
	{
	  ELOWER = state.E-(EDEM+sqrt(3.0*VDEM));
	  if(ELOWER < context.grid.EMIN){ ELOWER = context.grid.EMIN;}
	}
      else
	{
	  ELOWER = state.E-1.5*(EDEM+W21);
	  if(ELOWER < context.grid.EMIN){ ELOWER = context.grid.EMIN;}
	}
	  
      double XE1 = (log(ELOWER)-context.grid.DLEMP1)*context.grid.DLFC;
      int KE1 = (int)XE1;
      double XEK1 = XE1-(double)KE1;
      double STLWR = exp(mat.SETOT[KE1]+mat.DSETOT[KE1]*XEK1);
      if(ST < STLWR){ ST = STLWR;}
    }
  else
    {
      KSOFTI = 0;
      DESOFT = 0.0;
      SSOFT = 0.0;
    }

  //  ****  Soft elastic scattering.
  //        KSOFTE=1, soft scattering is active,
  //        KSOFTE=0, soft scattering is not active.

  if(T1 > 1.0E-20)
    {
      KSOFTE = 1;
    }
  else
    {
      KSOFTE = 0;
    }

  //  ****  Delta interactions.
  //        KDELTA=0, a hard interaction follows,
  //        KDELTA=1, a delta interaction follows.

  DST = -log(penRand.rand())/ST;
  if(DST < DSMAXP)
    {
      KDELTA = 0;
    }
  else
    {
      DST = DSMAXP;
      KDELTA = 1;
    }

  if(KSOFTI == 0 && KSOFTE == 0)
    {
      MHINGE = 1;
      DS = DST;
    }
  else
    {
      DS = DST*penRand.rand();
      DSR = DST-DS;
      if(KSOFTI == 1)
	{
	  if(DST < 1.0E-8)
	    {
	      SSOFT = W1;
	      DESOFT = SSOFT*DST;
	    }
	  else
	    {
	      double EDE0 = W1*DST;
	      double VDE0 = W2*DST;
	      double FSEDE = 1.0-mat.DW1EL[KE]*EDE0;
	      if(FSEDE < 0.75){ FSEDE = 0.75;}

	      double FSVDE = 1.0-mat.DW2EL[KE]*EDE0;
	      if(FSVDE < 0.75){ FSVDE = 0.75;}     
	      double EDE = EDE0*FSEDE;
	      double VDE = VDE0*FSVDE;
	      //  ****  Generation of random values DE with mean EDE and variance VDE.
	      double SIGMA = sqrt(VDE);
	      if(SIGMA < 0.333333333*EDE)
		{
		  //  ****  Truncated Gaussian distribution.
		  DESOFT = EDE+RNDG3(context.rndg3,penRand)*SIGMA;
		}
	      else
		{
		  double RU = penRand.rand();
		  double EDE2 = EDE*EDE;
		  double VDE3 = 3.0*VDE;
		  if(EDE2 < VDE3)
		    {
		      double PNULL = (VDE3-EDE2)/(VDE3+3.0*EDE2);
		      if(RU < PNULL)
			{
			  DESOFT = 0.0;
			  SSOFT = 0.0;
			  if(KSOFTE == 0)
			    {
			      MHINGE = 1;
			      DS = DST;
			    }
			  else
			    {
			      KSOFTI = 0;
			    }
			  return;
			}
		      else
			{
			  //  ****  Uniform distribution.
			  DESOFT = 1.5*(EDE+VDE/EDE)*(RU-PNULL)/(1.0-PNULL);
			}
		    }
		  else
		    {
		      DESOFT = EDE+(2.0*RU-1.0)*sqrt(VDE3);
		    }
		}
	      SSOFT = DESOFT/DST;
	    }
	}
    }
}


void pen_betaE::KNOCK(double &DE,
		      int &ICOL,
		      pen_rand& penRand)
{
  //  ************  Electrons.

  //Get material
  const pen_material& mat = *pmat;
    
  if(MHINGE == 1){}
  else
    {
      //  ****  Hinge, artificial soft event (ICOL=BETAe_SOFT_INTERACTION).
      
      ICOL = BETAe_SOFT_INTERACTION;
      MHINGE = 1;
      
      //  ****  Energy loss.
      DE = 0.0;

      // Knock must calculate the energy at the end of step (EENDSTEP)
      // to allow JUMP to compute the mean free paths for that energy.
      // Notice that the particle energy will reache EENDSTEP at
      // softEloss function if no interface is crossed.
      
      if(KSOFTI == 1){

	//Check if particle has sufficient energy to
	//end the whole step
	EENDSTEP = state.E-SSOFT*DSR;      
	if(EENDSTEP < mat.EABS[kpar]){
	  DE = state.E;
	  state.E = 0.0;
	  return;
	}
	if(KSOFTE == 0){return;}
	  
	//Get energy grid values at the hinge to
	//calculate angular deflection
	XEL = log(state.E);
	XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
	KE = (int)XE;
	XEK = XE-(double)KE;	  
	  
      }
      else{
	//Energy at the end of step is the actual one
	EENDSTEP = state.E;
      }

      //  ****  Angular deflection.

      if(mat.T1E[KE+1] > -78.3)
	{
	  T1 = exp(mat.T1E[KE]+(mat.T1E[KE+1]-mat.T1E[KE])*XEK);
	  T2 = exp(mat.T2E[KE]+(mat.T2E[KE+1]-mat.T2E[KE])*XEK);
	}
      else
	{
	  T1 = 0.0;
	  T2 = 0.0;
	}
      if(T1 < 1.0E-20){ return;}
      //  ****  1st and 2nd moments of the angular distribution.
      double EMU1 = 0.5*(1.0-exp(-DST*T1));
      double EMU2 = EMU1-(1.0-exp(-DST*T2))/6.0;
      //  ****  Sampling from a two-bar histogram with these moments.
      double PNUM = 2.0*EMU1-3.0*EMU2;
      double PDEN = 1.0-2.0*EMU1;
      double PMU0 = PNUM/PDEN;
      double PA = PDEN+PMU0;
      double RND = penRand.rand();
      double CDT;
      if(RND < PA)
	{
	  CDT = 1.0-2.0*PMU0*(RND/PA);
	}
      else
	{
	  CDT = 1.0-2.0*(PMU0+(1.0-PMU0)*((RND-PA)/(1.0-PA)));
	}
      double DF = constants::TWOPI*penRand.rand();
      DIRECT(CDT,DF,state.U,state.V,state.W);
      return;
    }

  //  ************  Hard event.

  MHINGE = 0;
  //  ****  A delta interaction (ICOL=7) occurs when the maximum
  //        allowed step length is exceeded.
  if(KDELTA == 1)
    {
      ICOL = BETAe_DELTA;
      DE = 0.0;
      return;
    }
  //  ****  Random sampling of the interaction type.
  
  double STNOW = P[0];
  STNOW += P[1];
  STNOW += P[2];
  STNOW += P[3];
  
  
  double STS = (STNOW > ST ? STNOW : ST)*penRand.rand();

  double SS = P[BETAe_HARD_ELASTIC];
  if(SS > STS)
    {
      ICOL = HelasticCol.interact(context,mat,(*this),DE,penRand);
      return;
    }
  
  SS += P[BETAe_HARD_INELASTIC];
  if(SS > STS)
    {
      ICOL = HinelastCol.interact(context,mat,(*this),DE,penRand);
      return;
    }

  SS += P[BETAe_HARD_BREMSSTRAHLUNG];
  if(SS > STS)
    {
      ICOL = Hbremsstrahlung.interact(context,mat,(*this),DE,penRand);
      return;
    }  

  SS += P[BETAe_HARD_INNER_SHELL];
  if(SS > STS)
    {
      ICOL = HinnerShell.interact(context,mat,(*this),DE,penRand);
      return;
    }
  
  //  ****  A delta interaction (ICOL=BETAe_DELTA) may occur when the total
  //        interaction probability per unit path length, ST, is
  //        larger than STNOW.
  ICOL = BETAe_DELTA;
  DE = 0.0;
}

void pen_betaE::KNOCKF(double &DE,
		       int &ICOL,
		       pen_rand& penRand)
{
  //  ************  Electrons (KPAR=PEN_ELECTRON).

  //Get material
  const pen_material& mat = *pmat;  
  
  if(MHINGE == 1){}
  else
    {
      //  ****  Hinge, artificial soft event (ICOL=1).

      ICOL = BETAe_SOFT_INTERACTION;
      MHINGE = 1;

      //  ****  Energy loss.	  
      DE = 0.0;

      // Knock must calculate the energy at the end of step (EENDSTEP)
      // to allow JUMP to compute the mean free paths for that energy.
      // Notice that the particle energy will reache EENDSTEP at
      // softEloss function if no interface is crossed.

      if(KSOFTI == 1)
	{

	  //Check if particle has sufficient energy to
	  //end the whole step
	  EENDSTEP = state.E-SSOFT*DSR;      
	  if(EENDSTEP < mat.EABS[kpar]){
	    DE = state.E;
	    DEA = DE;
	    state.E = 0.0;
	    return;
	  }
	  if(KSOFTE == 0){ return;}
	  XEL = log(state.E);
	  XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
	  KE = (int)XE;
	  XEK = XE-(double)KE;
	}
      else
	{
	  //Energy at the end of step is the actual one
	  EENDSTEP = state.E;
	  DEA = 0.0;
	}

      //  ****  Angular deflection.

      if(mat.T1E[KE+1] > -78.3)
	{
	  T1 = exp(mat.T1E[KE]+(mat.T1E[KE+1]-mat.T1E[KE])*XEK);
	  T2 = exp(mat.T2E[KE]+(mat.T2E[KE+1]-mat.T2E[KE])*XEK);
	}
      else
	{
	  T1 = 0.0;
	  T2 = 0.0;
	}
      if(T1 < 1.0E-20){ return;}
      //  ****  1st and 2nd moments of the angular distribution.
      double EMU1 = 0.5*(1.0-exp(-DST*T1));
      double EMU2 = EMU1-(1.0-exp(-DST*T2))/6.0;
      //  ****  Sampling from a two-bar histogram with these moments.
      double PNUM = 2.0*EMU1-3.0*EMU2;
      double PDEN = 1.0-2.0*EMU1;
      double PMU0 = PNUM/PDEN;
      double PA = PDEN+PMU0;
      double RND = penRand.rand();
      double CDT;
      if(RND < PA)
	{
	  CDT = 1.0-2.0*PMU0*(RND/PA);
	}
      else
	{
	  CDT = 1.0-2.0*(PMU0+(1.0-PMU0)*((RND-PA)/(1.0-PA)));
	}
      double DF = constants::TWOPI*penRand.rand();
      DIRECT(CDT,DF,state.U,state.V,state.W);
      return;
    }

  //  ************  Hard event.
	  
  MHINGE = 0;
  //  ****  A delta interaction (ICOL=7) occurs when the maximum
  //        allowed step length is exceeded.
  if(KDELTA == 1)
    {
      ICOL = BETAe_DELTA;
      DEA = 0.0;
      DE = 0.0;
      return;
    }
  
  //  ****  Interaction forcing.
  for(unsigned int KCOL = 0; KCOL < interactions; KCOL++)
    {
      if(context.FORCE[state.IBODY][PEN_ELECTRON][KCOL] > 1.0 && POR[KCOL] > 1.0E-16)
	{
	  P0[KCOL] = POR[KCOL];
	  P[KCOL] = POR[KCOL]*context.FORCE[state.IBODY][PEN_ELECTRON][KCOL];
	  LFORC[KCOL] = true;
	}
      else
	{
	  LFORC[KCOL] = false;
	}
    }
  IBR = (context.IBRSPL[state.IBODY] < 1 ? 1 : context.IBRSPL[state.IBODY]);
  //  ****  Random sampling of the interaction type.

  double STNOW = P[0];
  STNOW += P[1];
  STNOW += P[2];
  STNOW += P[3];
  
  
  double STS = (STNOW > ST ? STNOW : ST)*penRand.rand();

  double SS = P[BETAe_HARD_ELASTIC];
  if(SS > STS)
    {
      ICOL = HelasticCol.interactF(context,mat,(*this),DE,penRand);
      return;
    }
  
  SS += P[BETAe_HARD_INELASTIC];
  if(SS > STS)
    {
      ICOL = HinelastCol.interactF(context,mat,(*this),DE,penRand);
      return;
    }

  SS += P[BETAe_HARD_BREMSSTRAHLUNG];
  if(SS > STS)
    {
      ICOL = Hbremsstrahlung.interactF(context,mat,(*this),DE,penRand);
      return;
    }  

  SS += P[BETAe_HARD_INNER_SHELL];
  if(SS > STS)
    {
      ICOL = HinnerShell.interactF(context,mat,(*this),DE,penRand);
      return;
    }
  
  //  ****  A delta interaction (ICOL=7) may occur when the total
  //        interaction probability per unit path length, ST, is
  //        larger than STNOW.
  ICOL = BETAe_DELTA;
  DEA = 0.0;
  DE = 0.0;
  return;
}

void pen_betaE::softEloss(double& X,
			  double& Y,
			  double& Z,
			  double& DE,
			  pen_rand& penRand){

  if(DESOFT > 0.0){
    //Calculate energy loss
    DE = SSOFT*dsef;
    //Sample random position
    double dsRand = penRand.rand()*dsef;
    X = XL+dsRand*state.U;
    Y = YL+dsRand*state.V;
    Z = ZL+dsRand*state.W;
  
    if(DE >= state.E){
      DE = state.E;
      state.E = 0.0;
    }else{
      state.E -= DE;
    }
  }else{
    DE = 0.0;
  }
}

void pen_betaE::dpage()
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
  const double TREV=2.0E0*constants::REV;

  //  ****  Updating the particle's age after each free flight.

  if(SSOFT > 1.0E-16)
    {
      double E0=state.E;
      double SQE=sqrt(E0*(E0+TREV));
      double E1=E0-SSOFT*dsef;  // Energy at the end of the step, S=dsef.
      if(E0 < 1.000001*E1)
	{
	  state.PAGE = state.PAGE+dstot*RSLCM*(E0+constants::REV)/SQE;
	}
      else
	{
	  if(E1 < 50.0){ E1=50.0;}
	  double SQE1=sqrt(E1*(E1+TREV));
	  state.PAGE = state.PAGE+RSLCM*(SQE-SQE1)/SSOFT;
	  if(dstot > dsef && state.MAT != 0){
	    state.PAGE = state.PAGE+(dstot-dsef)*RSLCM*(E1+constants::REV)/SQE1;
	  }
	}
    }
  else
    {
      double SQE=sqrt(state.E*(state.E+TREV));
      if(state.MAT == 0)
	{
	  state.PAGE = state.PAGE+dsef*RSLCM*(state.E+constants::REV)/SQE;
	}
      else
	{
	  state.PAGE = state.PAGE+dstot*RSLCM*(state.E+constants::REV)/SQE;
	}
    }
}

// Electrons hard elastic collision

int pen_eHEC::interact(const pen_context& context,
		       const pen_material& mat,
		       pen_betaE& beta,
		       double& DE,
		       pen_rand& penRand) const
{
  //  ****  Hard elastic collision (ICOL=BETAe_HARD_ELASTIC).
  DE = 0.0;
  double RMU;
  if(beta.state.E >= mat.EELMAX)
    {
      double TRNDC = mat.RNDCE[beta.KE]+(mat.RNDCE[beta.KE+1]-mat.RNDCE[beta.KE])*beta.XEK;
      double TA = exp(mat.AE[beta.KE]+(mat.AE[beta.KE+1]-mat.AE[beta.KE])*beta.XEK);
      double TB = mat.BE[beta.KE]+(mat.BE[beta.KE+1]-mat.BE[beta.KE])*beta.XEK;
      EELa(TA,TB,TRNDC,RMU,penRand);
    }
  else
    {

      double TRNDC = mat.RNDCEd[beta.KE]+(mat.RNDCEd[beta.KE+1]-mat.RNDCEd[beta.KE])*beta.XEK;
      EELd(mat,beta.KE,beta.XEL,TRNDC,RMU,context.grid.DLEMP,context.grid.DLFC,penRand);  // Uses the ELSEPA database.
    }
  double CDT = 1.0-(RMU+RMU);
  double DF = constants::TWOPI*penRand.rand();
  DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
  return BETAe_HARD_ELASTIC;
}

int pen_eHEC::interactF(const pen_context& context,
			const pen_material& mat,
			pen_betaE& beta,
			double& DE,
			pen_rand& penRand) const
{

  //  ****  Hard elastic collision (ICOL=BETAe_HARD_ELASTIC).
  DE = 0.0;
  double RMU;
  if(beta.state.E >= mat.EELMAX)
    {
      double TRNDC = mat.RNDCE[beta.KE]+(mat.RNDCE[beta.KE+1]-mat.RNDCE[beta.KE])*beta.XEK;
      double TA = exp(mat.AE[beta.KE]+(mat.AE[beta.KE+1]-mat.AE[beta.KE])*beta.XEK);
      double TB = mat.BE[beta.KE]+(mat.BE[beta.KE+1]-mat.BE[beta.KE])*beta.XEK;
      EELa(TA,TB,TRNDC,RMU,penRand);
    }
  else
    {
      double TRNDC = mat.RNDCEd[beta.KE]+(mat.RNDCEd[beta.KE+1]-mat.RNDCEd[beta.KE])*beta.XEK;
      EELd(mat,beta.KE,beta.XEL,TRNDC,RMU,context.grid.DLEMP,context.grid.DLFC,penRand);  // Uses the ELSEPA database.
    }
  double CDT = 1.0-(RMU+RMU);
  double DF = constants::TWOPI*penRand.rand();
  DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
  beta.DEA = 0.0;
  return BETAe_HARD_ELASTIC;
}

// Electrons hard inelastic collision

int pen_eHIC::interact(const pen_context& context,
		       const pen_material& mat,
		       pen_betaE& beta,
		       double& DE,
		       pen_rand& penRand) const
{

  //  ****  Hard inelastic collision (ICOL=BETAe_HARD_INELASTIC).
  double EP,CDT,ES,CDTS;
  int IOSC;
  double DELTA = mat.DEL[beta.KE]+(mat.DEL[beta.KE+1]-mat.DEL[beta.KE])*beta.XEK;
  EINa(mat,beta.state.E,beta.KE,beta.XEL,DE,EP,CDT,ES,CDTS,IOSC,DELTA,context.grid.DLEMP,context.grid.DLFC,penRand);
  //  ****  Scattering angles (primary electron).
  double DF = constants::TWOPI*penRand.rand();
  //  ****  Delta ray.
  if(ES > mat.EABS[PEN_ELECTRON])
    {
      double DFS = DF+constants::PI;
      pen_particleState storeState;

      storeState.E = ES;
      
      storeState.X = beta.state.X;
      storeState.Y = beta.state.Y;
      storeState.Z = beta.state.Z;

      storeState.U = beta.state.U;
      storeState.V = beta.state.V;
      storeState.W = beta.state.W;

      storeState.WGHT = beta.state.WGHT;
      storeState.IBODY = beta.state.IBODY;
      storeState.MAT = beta.state.MAT;
      
      storeState.ILB[0] = beta.state.ILB[0]+1;
      storeState.ILB[1] = beta.getKpar();
      storeState.ILB[2] = BETAe_HARD_INELASTIC;
      storeState.ILB[3] = 0;
      storeState.ILB[4] = beta.state.ILB[4];

      storeState.LAGE = beta.state.LAGE;
      storeState.PAGE = beta.state.PAGE;
      
      DIRECT(CDTS,DFS,storeState.U,storeState.V,storeState.W);

      stack.store(storeState);
      if(penGetError() != PEN_SUCCESS){ return BETAe_HARD_INELASTIC;}        
    }
  //  ****  New energy and direction.
  if(EP > mat.EABS[PEN_ELECTRON])
    {
      beta.state.E = EP;
      DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
    }
  else
    {
      DE = beta.state.E;
      beta.state.E = 0.0;
    }
  return BETAe_HARD_INELASTIC;     
}

int pen_eHIC::interactF(const pen_context& context,
			const pen_material& mat,
			pen_betaE& beta,
			double& DE,
			pen_rand& penRand) const
{

  //  ****  Hard inelastic collision (ICOL=BETAe_HARD_INELASTIC).
  double EP,CDT,ES,CDTS;
  double WFORCE = 1.0;
  bool LCOL;
  int IOSC;
  int ICOL = BETAe_HARD_INELASTIC;
  if(beta.LFORC[BETAe_HARD_INELASTIC])  // Forced interaction.
    {
      WFORCE = beta.P0[BETAe_HARD_INELASTIC]/beta.P[BETAe_HARD_INELASTIC];
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

  double DELTA = mat.DEL[beta.KE]+(mat.DEL[beta.KE+1]-mat.DEL[beta.KE])*beta.XEK;
  EINa(mat,beta.state.E,beta.KE,beta.XEL,DE,EP,CDT,ES,CDTS,IOSC,DELTA,context.grid.DLEMP,context.grid.DLFC,penRand);
  //  ****  Scattering angles (primary electron).
  double DF = constants::TWOPI*penRand.rand();
  //  ****  Delta ray.
  if(ES > mat.EABS[PEN_ELECTRON])
    {
      double DFS = DF+constants::PI;

      pen_particleState stateStore;

      stateStore.E = ES;
      
      stateStore.X = beta.state.X;
      stateStore.Y = beta.state.Y;
      stateStore.Z = beta.state.Z;      
      
      stateStore.U = beta.state.U;
      stateStore.V = beta.state.V;
      stateStore.W = beta.state.W;

      stateStore.WGHT = beta.state.WGHT*WFORCE;
      stateStore.IBODY = beta.state.IBODY;
      stateStore.MAT = beta.state.MAT;
      
      stateStore.ILB[0] = beta.state.ILB[0]+1;
      stateStore.ILB[1] = beta.getKpar();
      stateStore.ILB[2] = ICOL;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = beta.state.ILB[4];

      stateStore.LAGE = beta.state.LAGE;
      stateStore.PAGE = beta.state.PAGE;      

      DIRECT(CDTS,DFS,stateStore.U,stateStore.V,stateStore.W);

      stack.store(stateStore);      
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
    }
  //  ****  New energy and direction.
  beta.DEA = DE;
  DE = beta.DEA*WFORCE;
  if(LCOL)
    {
      if(EP > mat.EABS[PEN_ELECTRON])
	{
	  beta.state.E = EP;
	  DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
	}
      else
	{
	  beta.DEA = beta.DEA+EP;
	  DE = DE+EP;
	  beta.state.E = 0.0;
	}
    }
  return ICOL;
}
  
// Electrons hard bremsstrahlung emission

pen_eHBE::pen_eHBE(pen_particleStack<pen_state_gPol>& stackIn)
  : abc_interaction(BETAe_HARD_BREMSSTRAHLUNG),
    stack(stackIn)
{
  const double TREV = 2.0*constants::REV;
  
  const int NE = 7;
  double E[NE];

  E[0] = 1.0E3;
  E[1] = 5.0E3;
  E[2] = 1.0E4;
  E[3] = 5.0E4;
  E[4] = 1.0E5;
  E[5] = 5.0E5;
  E[6] = 1.0E6;

  for(int IE = 0; IE < NET; IE++)
    {
      BET[IE] = sqrt(E[IE]*(E[IE]+TREV))/(E[IE]+constants::REV);
    }
}

int pen_eHBE::interact(const pen_context& /*context*/,
		       const pen_material& mat,
		       pen_betaE& beta,
		       double& DE,
		       pen_rand& penRand) const
{   
  //  ****  Hard bremsstrahlung emission (ICOL=BETAe_HARD_BREMSSTRAHLUNG).
    
  EBRa(mat,beta.state.E,beta.KE,beta.XEK,DE,constants::WB,penRand);
  //  ****  Bremsstrahlung photon.
  if(DE > mat.EABS[PEN_PHOTON])
    {
      double CDTS;
      EBRaA(mat,beta.state.E,DE,CDTS,BET,penRand);
      double DFS = constants::TWOPI*penRand.rand();

      pen_state_gPol storeState;

      storeState.E = DE;
      
      storeState.X = beta.state.X;
      storeState.Y = beta.state.Y;
      storeState.Z = beta.state.Z;

      storeState.U = beta.state.U;
      storeState.V = beta.state.V;
      storeState.W = beta.state.W;
      
      storeState.WGHT = beta.state.WGHT;
      storeState.IBODY = beta.state.IBODY;
      storeState.MAT = beta.state.MAT;

      storeState.ILB[0] = beta.state.ILB[0]+1;
      storeState.ILB[1] = beta.getKpar();
      storeState.ILB[2] = BETAe_HARD_BREMSSTRAHLUNG;
      storeState.ILB[3] = 0;
      storeState.ILB[4] = beta.state.ILB[4];

      storeState.LAGE = beta.state.LAGE;
      storeState.PAGE = beta.state.PAGE;

      DIRECT(CDTS,DFS,storeState.U,storeState.V,storeState.W);

      stack.store(storeState);
      
      if(penGetError() != PEN_SUCCESS){ return BETAe_HARD_BREMSSTRAHLUNG;}
    }
  //  ****  New energy.
  beta.state.E = beta.state.E-DE;
  if(beta.state.E < mat.EABS[PEN_ELECTRON])
    {
      DE = beta.state.E+DE;
      beta.state.E = 0.0;
    }
  return BETAe_HARD_BREMSSTRAHLUNG;
}

int pen_eHBE::interactF(const pen_context& /*context*/,
			const pen_material& mat,
			pen_betaE& beta,
			double& DE,
			pen_rand& penRand) const
{
	  
  //  ****  Hard bremsstrahlung emission (ICOL=BETAe_HARD_BREMSSTRAHLUNG).
  double WFORCE = 1.0;
  bool LCOL;
  unsigned int ICOL = BETAe_HARD_BREMSSTRAHLUNG;
  if(beta.LFORC[BETAe_HARD_BREMSSTRAHLUNG])  // Forced interaction.
    {
      WFORCE = beta.P0[BETAe_HARD_BREMSSTRAHLUNG]/beta.P[BETAe_HARD_BREMSSTRAHLUNG];
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
	  
  EBRa(mat,beta.state.E,beta.KE,beta.XEK,DE,constants::WB,penRand);
  //  ****  Bremsstrahlung photon.
  if(DE > mat.EABS[PEN_PHOTON])
    {
      double CDTS;
      EBRaA(mat,beta.state.E,DE,CDTS,BET,penRand);
      double WSPLIT = beta.state.WGHT*WFORCE/beta.IBR;
      for(int I = 0; I < beta.IBR; I++)
	{
	  double DFS = constants::TWOPI*penRand.rand();

	  pen_state_gPol stateStore;

	  stateStore.E = DE;

	  stateStore.X = beta.state.X;
	  stateStore.Y = beta.state.Y;
	  stateStore.Z = beta.state.Z;
	  
	  stateStore.U = beta.state.U;
	  stateStore.V = beta.state.V;
	  stateStore.W = beta.state.W;

	  stateStore.WGHT = WSPLIT;
	  stateStore.IBODY = beta.state.IBODY;
	  stateStore.MAT = beta.state.MAT;
	  
	  stateStore.ILB[0] = beta.state.ILB[0]+1;
	  stateStore.ILB[1] = beta.getKpar();
	  stateStore.ILB[2] = ICOL;
	  stateStore.ILB[3] = 0;
	  stateStore.ILB[4] = beta.state.ILB[4];

	  stateStore.LAGE = beta.state.LAGE;
	  stateStore.PAGE = beta.state.PAGE;
	  
	  DIRECT(CDTS,DFS,stateStore.U,stateStore.V,stateStore.W);

	  stack.store(stateStore);
	  if(penGetError() != PEN_SUCCESS){ return ICOL;}
	}
    }
  //  ****  New energy.
  beta.DEA = DE;
  DE = DE*WFORCE;
  if(LCOL)
    {
      beta.state.E = beta.state.E-beta.DEA;
      if(beta.state.E < mat.EABS[PEN_ELECTRON])
	{
	  beta.DEA = beta.DEA+beta.state.E;
	  DE = DE+beta.state.E;
	  beta.state.E = 0.0;
	}
    }
  return ICOL;	  
}
  
// Electrons inner-shell impact ionisation

int pen_eSII::interact(const pen_context& context,
		       const pen_material& mat,
		       pen_betaE& beta,
		       double& DE,
		       pen_rand& penRand) const
{
       
  //  ****  Ionisation of an inner shell (ICOL=BETAe_HARD_INNER_SHELL).
  double EP,CDT,ES,CDTS;
  int IZA,ISA;
  double DELTA = mat.DEL[beta.KE]+(mat.DEL[beta.KE+1]-mat.DEL[beta.KE])*beta.XEK;
  ESIa(mat,beta.state.E,beta.KE,beta.XEL,DE,EP,CDT,ES,CDTS,IZA,ISA,context.grid.DLEMP,context.grid.DLFC,DELTA,penRand);
  
  //  ****  Atomic relaxation.
  if(IZA > 2)
    {
      int ks;
      RELAX(context.elements,mat,beta.state,BETAe_HARD_INNER_SHELL,
	    1,IZA,ISA,ks,stackE,stackG,penRand);
      
      if(penGetError() != PEN_SUCCESS){ return BETAe_HARD_INNER_SHELL;}
    }
  //  ****  Scattering angles (primary electron).
  double DF = constants::TWOPI*penRand.rand();
  //  ****  Delta ray.
  if(ES > mat.EABS[PEN_ELECTRON])
    {
      double DFS = DF+constants::PI;

      pen_particleState storeState;

      storeState.E = ES;

      storeState.X = beta.state.X;
      storeState.Y = beta.state.Y;
      storeState.Z = beta.state.Z;
      
      storeState.U = beta.state.U;
      storeState.V = beta.state.V;
      storeState.W = beta.state.W;
      
      storeState.WGHT = beta.state.WGHT;
      storeState.IBODY = beta.state.IBODY;
      storeState.MAT = beta.state.MAT;

      storeState.ILB[0] = beta.state.ILB[0]+1;
      storeState.ILB[1] = beta.getKpar();
      storeState.ILB[2] = BETAe_HARD_INNER_SHELL;
      storeState.ILB[3] = 0;
      storeState.ILB[4] = beta.state.ILB[4];

      storeState.LAGE = beta.state.LAGE;
      storeState.PAGE = beta.state.PAGE;
	  
      DIRECT(CDTS,DFS,storeState.U,storeState.V,storeState.W);

      stackE.store(storeState);

      if(penGetError() != PEN_SUCCESS){ return BETAe_HARD_INNER_SHELL;}
    }
  //  ****  New energy and direction.
  if(EP > mat.EABS[PEN_ELECTRON])
    {
      beta.state.E = EP;
      DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
    }
  else
    {
      DE = beta.state.E;
      beta.state.E = 0.0;
    }
  
  return BETAe_HARD_INNER_SHELL;
}

int pen_eSII::interactF(const pen_context& context,
			const pen_material& mat,
			pen_betaE& beta,
			double& DE,
			pen_rand& penRand) const
{

  //  ****  Ionization of an inner shell (ICOL=BETAe_HARD_INNER_SHELL).
  double EP,CDT,ES,CDTS;
  int IZA,ISA;
  bool LCOL;
  double WFORCE = 1.0;
  int ICOL = BETAe_HARD_INNER_SHELL;
  if(beta.LFORC[BETAe_HARD_INNER_SHELL])  // Forced interaction.
    {
      WFORCE = beta.P0[BETAe_HARD_INNER_SHELL]/beta.P[BETAe_HARD_INNER_SHELL];
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

  double DELTA = mat.DEL[beta.KE]+(mat.DEL[beta.KE+1]-mat.DEL[beta.KE])*beta.XEK;
  ESIa(mat,beta.state.E,beta.KE,beta.XEL,DE,EP,CDT,ES,CDTS,IZA,ISA,context.grid.DLEMP,context.grid.DLFC,DELTA,penRand);
  //  ****  Scattering angles (primary electron).
  double DF = constants::TWOPI*penRand.rand();
  //  ****  Delta ray.
  if(ES > mat.EABS[PEN_ELECTRON])
    {
      double DFS = DF+constants::PI;

      pen_particleState stateStore;

      stateStore.E = ES;
      
      stateStore.X = beta.state.X;
      stateStore.Y = beta.state.Y;
      stateStore.Z = beta.state.Z;
      
      stateStore.U = beta.state.U;
      stateStore.V = beta.state.V;
      stateStore.W = beta.state.W;

      stateStore.WGHT = beta.state.WGHT*WFORCE;
      stateStore.IBODY = beta.state.IBODY;
      stateStore.MAT = beta.state.MAT;
      
      stateStore.ILB[0] = beta.state.ILB[0]+1;
      stateStore.ILB[1] = beta.getKpar();
      stateStore.ILB[2] = ICOL;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = beta.state.ILB[4];

      stateStore.LAGE = beta.state.LAGE;
      stateStore.PAGE = beta.state.PAGE;
      
      DIRECT(CDTS,DFS,stateStore.U,stateStore.V,stateStore.W);

      stackE.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
    }
  //  ****  Atomic relaxation.
  if(IZA > 2)
    {
      int ks;
      double WGHTA = beta.state.WGHT;
      beta.state.WGHT = beta.state.WGHT*WFORCE;
      RELAX(context.elements,mat,beta.state,ICOL,
	    1,IZA,ISA,ks,stackE,stackG,penRand);
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
      beta.state.WGHT = WGHTA;
    }
  //  ****  New energy and direction.
  beta.DEA = DE;
  DE = beta.DEA*WFORCE;
  if(LCOL)
    {
      if(EP > mat.EABS[PEN_ELECTRON])
	{
	  beta.state.E = EP;
	  DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
	}
      else
	{
	  beta.DEA = beta.DEA+EP;
	  DE = DE+EP;
	  beta.state.E = 0.0;
	}
    }
  return ICOL;
}


//  *********************************************************************
//                       SUBROUTINE EBRa
//  *********************************************************************
void EBRa(const pen_material& mat,
	  const double E,
	  const unsigned KE,
	  const double XEK,
	  double &W,
	  const double* WB,
	  pen_rand& penRand)
{
  //  Simulation of bremsstrahlung emission by electrons or positrons in
  //  material M.

  
  if(mat.WCR > E)
  {
    W = 0.0;
    return;
  }

  //  ****  Selection of the energy grid point.
  int IE;
  if(penRand.rand() < XEK)
  {
    IE = KE+1;
  }
  else
  {
    IE = KE;
  }
  //  ****  Pointer.
  double PT, W1, W2, DW, B, A, PMAX;
  int I, J, K;
  bool Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    PT = mat.PBCUT[IE]+penRand.rand()*(mat.PACB[IE][constants::NBW-1]-mat.PBCUT[IE]);
      //  ****  Binary search of the W-interval.
    I = 0;
    J = constants::NBW-1;
    while(J-I > 1)
    {
      K = (I+J)/2;
      if(PT > mat.PACB[IE][K])
      {
        I = K;
      }
      else
      {
        J = K;
      }
    }
      //  ****  Sampling the photon energy (rejection method).
    W1 = WB[I];
    W2 = WB[I+1];
    DW = W2-W1;
    B = mat.DPDFB[IE][I]/DW;
    A = mat.PDFB[IE][I]-B*W1;
    if(W1 < mat.WBCUT[IE]){ W1 = mat.WBCUT[IE];}
    if(W2 < W1)
    {
      printf(" **** WARNING: EBR. Conflicting end-point values.\n");
      W = W1;
      return;
    }

    PMAX = A+B*W1;
    if(PMAX < A+B*W2){ PMAX = A+B*W2;}

    do
    {
      W = W1*pow(W2/W1,penRand.rand());
    }while(penRand.rand()*PMAX > A+B*W);
      
    W = W*E;
    if(W < mat.WCR){Eixir = false; continue;}
  }
     
}

//  *********************************************************************
//                       SUBROUTINE EELd
//  *********************************************************************
void EELd(const pen_material& mat,
	  const unsigned KE,
	  const double XEL,
	  double &RNDC,
	  double &RMU,
	  const double* DLEMP,
	  const double DLFC,
	  pen_rand& penRand)
{
  //     Simulation of electron hard elastic events. Cross sections from
  //  the ELSEPA numerical database.

  //  Argument value:
  //    RNDC ... cutoff value of the uniform random number
  //             (only hard events are simulated).
  //    RMU .... sampled angular deflection, =(1-CDT)/2.

  
  const unsigned int NPM1 = mat.NPM1_P;
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
  //  ****  Pointer.
  double RU = RNDC+penRand.rand()*(1.0-RNDC);

  //  ****  Selection of the interval (binary search in a restricted
  //        interval).
  int ITN = (int)(RU*NPM1);
  int I = mat.ITLE[ITN][JE]-1;
  int J = mat.ITUE[ITN][JE]-1;
  if(J-I < 2){}
  else
    {
      do
	{
	  int K = (I+J)/2;
	  if(RU > mat.PSE[K][JE])
	    {
	      I = K;
	    }
	  else
	    {
	      J = K;
	    }
	}while(J-I > 1);
    }
  //  ****  Sampling from the rational inverse cumulative distribution.
  double PP = mat.PSE[I][JE];
  double RR = RU-PP;
  if(RR > 1.0E-16)
    {
      double XX = mat.XSE[I][JE];
      double AA = mat.ASE[I][JE];
      double BB = mat.BSE[I][JE];
      double D = mat.PSE[I+1][JE]-PP;
      RMU = XX+((1.0+AA+BB)*D*RR/(D*D+(AA*D+BB*RR)*RR))
	*(mat.XSE[I+1][JE]-XX);
    }
  else
    {
      RMU = mat.XSE[I][JE];
    }
}

//  *********************************************************************
//                       SUBROUTINE EINa
//  *********************************************************************
void EINa(const pen_material& mat,
	  const double E,
	  const int KE,
	  const double XEL,
	  double &DE,
	  double &EP,
	  double &CDT,
	  double &ES,
	  double &CDTS,
	  int &IOSC,
	  const double DELTA,
	  const double* DLEMP,
	  const double DLFC,
	  pen_rand& penRand)
{
  //  Random sampling of hard inelastic collisions of electrons.
  //
  //  Sternheimer-Liljequist GOS model.
  //
  //  Input arguments:
  //    DELTA ... Fermi's density effect correction.
  //  Output arguments:
  //    DE ...... energy loss (eV).
  //    EP ...... energy of the scattered electron (eV).
  //    CDT ..... cosine of the polar scattering angle.
  //    ES ...... energy of the emitted secondary electron (eV).
  //    CDTS .... polar cosine of direction of the secondary electron.
  //    IOSC .... index of the oscillator that has been 'ionised'.

  
  bool LDIST;
  const double RREV = 1.0/constants::REV;
  const double TREV = 2.0*constants::REV;
  const double RTREV = 1.0/TREV;
  
  double WCCM = mat.WCC;
  if(WCCM > E)
  {
    DE = 0.0;
    EP = E;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    IOSC = constants::NO;
    return;
  }
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
  int IT;
  int IO = 0;
  int JO = mat.NEIN;
  while(JO-IO > 1)
  {
    IT = (IO+JO)/2;
    if(TST > mat.EINAC[JE][IT])
    {
      IO = IT;
    }
    else
    {
      JO = IT;
    }
  }
  IOSC = mat.IEIN[IO];
  double UK = mat.UI[IOSC-1];
  double WK = mat.WRI[IOSC-1];
  double WTHR;
  if(UK > 1.0E-3)
  {
    WTHR = WCCM;
    if(WTHR < UK){ WTHR = UK;}
  }
  else
  {
    WTHR = WCCM;
    if(WTHR < WK){ WTHR = WK;}
  }

  if(E < WTHR+1.0E-6)
  {
    DE = 0.0;
    EP = E;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    IOSC = constants::NO;
    return;
  }

  //  ****  Trick: The resonance energy and the cutoff recoil energy of
  //        inner shells are varied to yield a smooth threshold.

  LDIST = true;
  double WM, WKP, QKP;
  double EE, WCMAX, WDMAX;
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
    if(WCCM > WM){ LDIST = false;}
    EE = E+UK;
    WCMAX = 0.5*EE;
    WDMAX = WM;
    if(WDMAX > WCMAX){ WDMAX = WCMAX;}
    if(WTHR > WDMAX){ LDIST = false;}
  }
  else
  {
    if(WCCM > WK){ LDIST = false;}
    WKP = WK;
    QKP = WK;
    WM = E;
    EE = E;
    WCMAX = 0.5*EE;
    WDMAX = WKP+1.0;
  }

  //  ****  Constants.

  double RB = E+TREV;
  double GAM = 1.0+E*RREV;
  double GAM2 = GAM*GAM;
  double BETA2 = (GAM2-1.0)/GAM2;
  double AMOL = pow((GAM-1.0E0)/GAM,2);
  double CPS = E*RB;
  double CP = sqrt(CPS);

  //  ************  Partial cross sections of the active oscillator.

  //  ****  Distant excitations.
  double CPP, CPPS;
  double QM, XHDL, XHDT;
  if(LDIST)
  {
      
    CPPS = (E-WKP)*(E-WKP+TREV);
    CPP = sqrt(CPPS);
    if(WKP > 1.0E-6*E)
    {
      QM = sqrt(pow(CP-CPP,2)+constants::REV*constants::REV)-constants::REV;
    }
    else
    {
      QM = pow(WKP,2)/(BETA2*TREV);
      QM = QM*(1.0-QM*RTREV);
    }
    if(QM < QKP)
    {
      double RWKP = 1.0/WKP;
      XHDL = log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;
      XHDT = log(GAM2)-BETA2-DELTA;
      if(XHDT < 0.0){ XHDT = 0.0;}
      XHDT = XHDT*RWKP;
      if(UK > 1.0E-3)
      {
      double F0 = (WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/pow(WM-UK,2);
      XHDL = F0*XHDL;
      XHDT = F0*XHDT;
      }
    }
    else
    {
      XHDL = 0.0;
      XHDT = 0.0;
    }
  }
  else
  {
    QM = 0.0;    // Defined to avoid compilation warnings.
    CPP = 0.0;   // same
    CPPS = 0.0;  // same
    XHDL = 0.0;
    XHDT = 0.0;
  }
  //  ****  Close collisions.
  double RCL = WTHR/EE;
  double RL1, RRL1, XHC;
  if(RCL < 0.5)
  {
    RL1 = 1.0-RCL;
    RRL1 = 1.0/RL1;
    XHC = (AMOL*(0.5-RCL)+1.0/RCL-RRL1+(1.0-AMOL)*log(RCL*RRL1))/EE;
  }
  else
  {
    XHC = 0.0;
  }

  double XHTOT = XHC+XHDL+XHDT;
  if(XHTOT < 1.0E-35)
  {
    DE = 0.0;
    EP = E;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    IOSC = constants::NO;
    return;
  }

  //  ************  Sampling of final-state variables.

  TST = penRand.rand()*XHTOT;
  
  //  ****  Hard close collision.

  double TS1 = XHC;
  double A, ARCL, RK, RK2, RKF, PHI;
  if(TST < TS1)
  {
    A = 5.0*AMOL;
    ARCL = A*0.5*RCL;

    do
    {
      double FB = (1.0+ARCL)*penRand.rand();
      if(FB < 1.0)
      {
        RK = RCL/(1.0E0-FB*(1.0E0-(RCL+RCL)));
      }
      else
      {
        RK = RCL+(FB-1.0E0)*(0.5E0-RCL)/ARCL;
      }
      RK2 = RK*RK;
      RKF = RK/(1.0-RK);
      PHI = 1.0+pow(RKF,2)-RKF+AMOL*(RK2+RKF);
    }while(penRand.rand()*(1.0+A*RK2) > PHI);
      //  ****  Energy and scattering angle (primary electron).
    DE = RK*EE;
    EP = E-DE;
    CDT = sqrt(EP*RB/(E*(RB-DE)));
    //  ****  Energy and emission angle of the delta ray.
    if(mat.KS[IOSC-1] < 17)
    {
      if(UK > mat.ECUTR)
      {
        ES = DE-UK;  // Inner shells only.
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
    CDTS = sqrt(DE*RB/(E*(DE+TREV)));
    return;
  }

  //  ****  Hard distant longitudinal interaction.

  TS1 = TS1+XHDL;
  if(UK > 1.0E-3)
  {
    DE = WM-sqrt(pow(WM-WTHR,2)-penRand.rand()*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
  }
  else
  {
    DE = WKP;
  }
  EP = E-DE;
  double QS, Q, QTREV;
  if(TST < TS1)
  {
    QS = QM/(1.0+QM*RTREV);
    Q = QS/(pow((QS/QKP)*(1.0+QKP*RTREV),penRand.rand())-(QS*RTREV));
    QTREV = Q*(Q+TREV);
    CDT = (CPPS+CPS-QTREV)/(2.0*CP*CPP);
    if(CDT > 1.0){ CDT = 1.0;}
    //  ****  Energy and emission angle of the delta ray.
    if(mat.KS[IOSC-1] < 17)
    {
      ES = DE-UK;  // Inner shells only.
    }
    else
    {
      ES = DE;
    }
    CDTS = 0.5*(WKP*(E+RB-WKP)+QTREV)/sqrt(CPS*QTREV);
    if(CDTS > 1.0){ CDTS = 1.0;}
    return;
  }

  //  ****  Hard distant transverse interaction.

  CDT = 1.0;
  //  ****  Energy and emission angle of the delta ray.
  if(mat.KS[IOSC-1] < 17)
  {
    if(UK > mat.ECUTR)
    {
      ES = DE-UK;  // Inner shells only.
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
  CDTS = 1.0;

}

//  *********************************************************************
//                       SUBROUTINE ESIa
//  *********************************************************************
void ESIa(const pen_material& mat,
	  const double E,
	  const int KE,
	  const double XEL,
	  double &DE,
	  double &EP,
	  double &CDT,
	  double &ES,
	  double &CDTS,
	  int &IZZ,
	  int &ISH,
	  const double* DLEMP,
	  const double DLFC,
	  const double DELTA,
	  pen_rand& penRand)
{
  //  Random sampling of inner-shell ionisation by electron impact.
  //
  //  Sternheimer-Liljequist GOS model.
  //
  //  Input arguments:
  //    E ....... electron energy (eV).
  //    M ....... material where electrons propagate.
  //    DELTA ... Fermi's density effect correction.
  //  Output arguments:
  //    DE ...... energy loss (eV).
  //    EP ...... energy of the scattered electron (eV).
  //    CDT ..... cosine of the polar scattering angle.
  //    ES ...... energy of the emitted secondary electron (eV).
  //    CDTS .... polar cosine of direction of the secondary electron.
  //    IZZ ..... atomic number of the element where ionisation has
  //             occurred.
  //    ISH ..... atomic electron shell that has been ionised.


  const double RREV = 1.0/constants::REV;
  const double TREV = 2.0*constants::REV;
  const double RTREV = 1.0/TREV;

  //  ****  Energy grid point.
  int JE;
  double PK = (XEL-DLEMP[KE])*DLFC;
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
  int IO, JO, IT, IOSC;
  //  ****  Binary search.
  IO = 0;
  JO = mat.NESI;
  while(JO-IO > 1)
  {
    IT = (IO+JO)/2;
    if(TST > mat.ESIAC[JE][IT])
    {
      IO = IT;
    }
    else
    {
      JO = IT;
    }
  }
  IOSC = mat.IESI[IO];

  IZZ = mat.KZ[IOSC-1];
  ISH = mat.KS[IOSC-1];

  double UK = mat.UI[IOSC-1];
  double WK = mat.WRI[IOSC-1];
  double WTHR;
  if(UK > 1.0E-3)
  {
    WTHR = UK;
  }
  else
  {
    WTHR = WK;
  }

  if(E < WTHR+1.0E-6)
  {
    DE = UK;
    EP = E-DE;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    return;
  }

  //  ****  Trick: The resonance energy and the cutoff recoil energy of
  //        inner shells are varied to yield a smooth threshold.

  double WM = 3.0*WK-2.0*UK;
  double WKP, QKP;
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
  double EE = E+UK;
  double WCMAX = 0.5*EE;
  double WDMAX;
  if(WM < WCMAX){ WDMAX = WM;}else{ WDMAX = WCMAX;}
      

  //  ****  Constants.

  const double RB =E+TREV;
  const double GAM = 1.0+E*RREV;
  const double GAM2 = GAM*GAM;
  const double BETA2 = (GAM2-1.0)/GAM2;
  const double AMOL = pow((GAM-1.0)/GAM,2);
  const double CPS = E*RB;
  const double CP = sqrt(CPS);

  //  ************  Partial cross sections of the active oscillator.

  //  ****  Distant excitations.
      
  double CPPS =(E-WKP)*(E-WKP+TREV);
  double CPP = sqrt(CPPS);
  double QM;
  double RWKP, XHDL, XHDT, F0;
  if(WKP > 1.0E-6*E)
  {
    QM = sqrt(pow(CP-CPP,2)+constants::REV*constants::REV)-constants::REV;
  }
  else
  {
    QM = pow(WKP,2)/(BETA2*TREV);
    QM = QM*(1.0-QM*RTREV);
  }
  if(QM < QKP)
  {
    RWKP = 1.0/WKP;
    XHDL = log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;

    XHDT = log(GAM2)-BETA2-DELTA;
    if(XHDT < 0.0){ XHDT = 0.0;}
    
    XHDT = XHDT*RWKP;
    F0 = (WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/pow(WM-UK,2);
    XHDL = F0*XHDL;
    XHDT = F0*XHDT;
  }
  else
  {
    XHDL = 0.0;
    XHDT = 0.0;
  }
      //  ****  Close collisions.
  double RCL = WTHR/EE;
  double RL1 = 1.0-RCL;
  double RRL1 = 1.0/RL1;
  double XHC = (AMOL*(0.5-RCL)+1.0/RCL-RRL1+(1.0-AMOL)*log(RCL*RRL1))/EE;

  double XHTOT = XHC+XHDL+XHDT;
  if(XHTOT < 1.0E-35)
  {
    DE = UK;
    EP = E-DE;
    CDT = 1.0;
    ES = 0.0;
    CDTS = 0.0;
    return;
  }

      //  ************  Sampling of final-state variables.

  TST=penRand.rand()*XHTOT;

      //  ****  Hard close collision.

  double TS1 = XHC;
  double A, ARCL, FB, RK, RK2, RKF, PHI;
  if(TST < TS1)
  {
    A = 5.0*AMOL;
    ARCL = A*0.5*RCL;
    do
    {
      FB = (1.0+ARCL)*penRand.rand();
      if(FB < 1.0)
      {
        RK = RCL/(1.0-FB*(1.0-(RCL+RCL)));
      }
      else
      {
        RK = RCL+(FB-1.0)*(0.5-RCL)/ARCL;
      }
      RK2 = RK*RK;
      RKF = RK/(1.0-RK);
      PHI = 1.0+pow(RKF,2)-RKF+AMOL*(RK2+RKF);
    }while(penRand.rand()*(1.0+A*RK2) > PHI);
    //  ****  Energy and scattering angle (primary electron).
    DE = RK*EE;
    EP = E-DE;
    CDT = sqrt(EP*RB/(E*(RB-DE)));
    //  ****  Energy and emission angle of the delta ray.
    ES = DE-UK;
    CDTS = sqrt(DE*RB/(E*(DE+TREV)));
    return;
  }

  //  ****  Hard distant longitudinal interaction.

  double QS, Q, QTREV;
  TS1 = TS1+XHDL;
  DE = WM-sqrt(pow(WM-WTHR,2)-penRand.rand()*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
  EP = E-DE;
  if(TST < TS1)
  {
    QS = QM/(1.0+QM*RTREV);
    Q = QS/(pow((QS/QKP)*(1.0+QKP*RTREV),penRand.rand())-(QS*RTREV));
    QTREV = Q*(Q+TREV);
    CDT = (CPPS+CPS-QTREV)/(2.0*CP*CPP);
    if(CDT > 1.0){ CDT = 1.0;}
    //  ****  Energy and emission angle of the delta ray.
    ES = DE-UK;  // Inner shells only.
    CDTS = 0.5*(WKP*(E+RB-WKP)+QTREV)/sqrt(CPS*QTREV);
    if(CDTS > 1.0){ CDTS = 1.0;}
    return;
  }

  //  ****  Hard distant transverse interaction.

  CDT = 1.0;
  //  ****  Energy and emission angle of the delta ray.
  ES = DE-UK;  // Inner shells only.
  CDTS = 1.0;
      
}

//  *********************************************************************
//                       SUBROUTINE EELa
//  *********************************************************************
void EELa(double A,
	  double B,
	  double RNDC,
	  double &RMU,
	  pen_rand& penRand)
{
  //  Simulation of hard elastic events. Modified-Wentzel model.
  //
  //  Input arguments:
  //    A, B ... angular distribution parameters.
  //    RNDC ... cutoff probability.
  //  Output values:
  //    RMU .... angular deflection, =(1-CDT)/2.

  double A1 = A+1.0;
  double B1, BB;
  double RND0, RND, RNDMB, RNDRC;
  if(B >= 0.0)
  {
      
    //  ****  Case I.
      
    double RMUAV = A*A1*log(A1/A)-A;
    B1 = 1.0-B;
    RND0 = B1*A1*RMUAV/(A+RMUAV);
    RND = RNDC+penRand.rand()*(1.0-RNDC);
    if(RND < RND0)
    {
      RMU = RND*A/(B1*A1-RND);
    }
    else if(RND > RND0+B)
    {
      RNDMB = RND-B;
      RMU = RNDMB*A/(B1*A1-RNDMB);
    }
    else
    {
      RMU = RMUAV;
    }
  }
  else
  {

    //  ****  Case II.

    BB = -B;
    B1 = 1.0-BB;
    double RMUC = RNDC*A/(B1*A1-RNDC);
    double PW = B1*A*(1.0-RMUC)/(A+RMUC);
    if(penRand.rand()*(BB+PW) < BB)
    {
      RMU = 0.5*(1.0+sqrt(penRand.rand()));
    }
    else
    {
      RNDRC = penRand.rand()*(1.0-RMUC);
      RMU = (A*RNDRC+A1*RMUC)/(A1-RNDRC);
    }
  }
  
}

//  *********************************************************************
//                       SUBROUTINE EBRaA
//  *********************************************************************
void EBRaA(const pen_material& mat,
	   double &E,
	   double &DE,
	   double &CDT,
	   const double* BET,
	   pen_rand& penRand)
{
  //  Random sampling of the initial direction of bremss photons, relative
  //  to the direction of the projectile.
  //  Numerical fit/interpolation of partial-wave shape functions generated
  //  by the program BREMS of A. Poskus, Comp. Phys. Commun. (2018).

  //  Input parameters:
  //    M ..... material where the projectile moves.
  //    E ..... kinetic energy of the projectile.
  //    DE .... energy of the emitted photon.
  //  Output parameter:
  //    CDT ... cosine of the polar emission angle.

  
  const double TREV = 2.0*constants::REV;
  
  //  ****  Distribution parameters.
  
  double BETA = sqrt(E*(E+TREV))/(E+constants::REV);

  //  A pure dipole distribution is used for E>500 keV.
  if(E > 1.0E6)
  {
    CDT = 2.0*penRand.rand()-1.0;
    if(penRand.rand() > 0.75)
    {
      if(CDT < 0.0)
      {
        CDT = -pow(-CDT,0.333333333333333);
      }
      else
      {
      CDT = pow(CDT,0.333333333333333);
      }
    }
    CDT = (CDT+BETA)/(1.0+BETA*CDT);

  }
  else
  {
    int IE, IET, IE1;
    if(BETA > BET[6])
    {
      BETA = BET[6];
      IE = 6;
    }
    else if(BETA < BET[0])
    {
      BETA = BET[0];
      IE = 1;
    }
    else
    {
      IE = 1;
      IE1 = 7;
      do
	{
	  IET = (IE+IE1)/2;
	  if(BETA > BET[IET - 1])
	    {
	      IE = IET;
	    }
	  else
	    {
	      IE1 = IET;
	    }
	}while(IE1-IE > 1);
    }
    IE--;

    double RK = 1.0 + 40.0 * DE/E;
    int IK = std::min(int(RK),40);
    //int IK = int(RK);
    //if(IK > 40){ IK = 40;}
    IK--;

    double P10 = mat.BP1[IE][IK][0]+
    BETA*(mat.BP1[IE][IK][1]+
    BETA*(mat.BP1[IE][IK][2]+
    BETA*mat.BP1[IE][IK][3]));
  
    double P11 = mat.BP1[IE][IK+1][0]+
    BETA*(mat.BP1[IE][IK+1][1]+
    BETA*(mat.BP1[IE][IK+1][2]+
    BETA*mat.BP1[IE][IK+1][3]));
  
    double P1 = exp(P10+(RK-(IK+1))*(P11-P10));

    double P20 = mat.BP2[IE][IK][0]+
    BETA*(mat.BP2[IE][IK][1]+
    BETA*(mat.BP2[IE][IK][2]+
    BETA*mat.BP2[IE][IK][3]));
  
    double P21 = mat.BP2[IE][IK+1][0]+
    BETA*(mat.BP2[IE][IK+1][1]+
    BETA*(mat.BP2[IE][IK+1][2]+
    BETA*mat.BP2[IE][IK+1][3]));
  
    double P2 = P20+(RK-(IK+1))*(P21-P20);
  
  //  ****  Sampling from the Lorentz-transformed dipole distributions.

    if (penRand.rand() < P1)
      {
        bool Eixir = false;
        while (!Eixir)
          {
            Eixir = true;
            CDT = 2.0 * penRand.rand() - 1.0;
            if (2.0 * penRand.rand() > 1.0 + CDT * CDT)
              {
                Eixir = false;
                continue;
              }
          }
      }
    else
      {
        bool Eixir = false;
        while (!Eixir)
          {
            Eixir = true;
            CDT = 2.0 * penRand.rand() - 1.0;
            if (penRand.rand() > 1.0 - CDT * CDT)
              {
                Eixir = false;
                continue;
              }
          }
      }
    CDT = (CDT+P2)/(1.0+P2*CDT);    
  }  
}
/*
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
*/
