
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


#include "pen_positron.hh"

const double pen_betaP::REV = 5.10998928e5;  // Electron rest energy (eV)
const double pen_betaP::mc2 = 5.10998918E5;
const double pen_betaP::twomc2 = 2.0*mc2;    

pen_betaP::pen_betaP(const pen_context& contextIn,
		     pen_particleStack<pen_particleState>& stackEin,
		     pen_particleStack<pen_state_gPol>& stackGin)
  : abc_particle(contextIn,PEN_POSITRON,5,twomc2), 
    HinelastCol(stackEin),
    Hbremsstrahlung(stackGin),
    HinnerShell(stackEin,stackGin),
    annihilation(stackGin),
    MHINGE(0),
    DESOFT(0.0),
    SSOFT(0.0),
    stackG(stackGin)
    
{
  //Check ids of interactions
  int auxID = HelasticCol.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Positron elastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = HinelastCol.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Positron inelastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = Hbremsstrahlung.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Positron inelastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = HinnerShell.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Positron inelastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  auxID = annihilation.getID();
  if(auxID < 0 || auxID >= (int)interactions)
    {
      char error[300];
      sprintf(error,"Positron inelastic colision ID (%d) out of range [0,%d)\n  Check enums in 'interactions.hh'",auxID,interactions);
      throw std::out_of_range(error);
    }
  
  HelasticCol.init(context);
  HinelastCol.init(context);
  Hbremsstrahlung.init(context);
  HinnerShell.init(context);
  annihilation.init(context);
}

void pen_betaP::JUMP(double &DS,
		     pen_rand& penRand,
		     const double DSMAX)
{
  //  ************  Positrons.
  
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
  ST += P[4];
  
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

      //
      //  ****  Upper bound for the interaction probability along the step
      //        (including soft energy straggling).

      double EDE0 = W1*DSMAXP;
      double VDE0 = W2*DSMAXP;

      double FSEDE = 1.0-mat.DW1PL[KE]*EDE0;
      if(FSEDE < 0.75){ FSEDE = 0.75;}

      double FSVDE = 1.0-mat.DW2PL[KE]*EDE0;
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
	  if( ELOWER < context.grid.EMIN){ ELOWER = context.grid.EMIN;}
	}
      double XE1 = (log(ELOWER)-context.grid.DLEMP1)*context.grid.DLFC;
      int KE1 = (int)XE1;
      double XEK1 = XE1-(double)KE1;
      double STLWR = exp(mat.SPTOT[KE1]+(mat.SPTOT[KE1+1]-mat.SPTOT[KE1])*XEK1);
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
	      double FSEDE = 1.0-mat.DW1PL[KE]*EDE0;
	      if(FSEDE < 0.75){ FSEDE = 0.75;}

	      double FSVDE = 1.0-mat.DW2PL[KE]*EDE0;
	      if(FSVDE < 0.75){ FSVDE = 0.75;}      
	      double EDE = EDE0*FSEDE;
	      double VDE = VDE0*FSVDE;
	      //  ****  Generation of random values DE with mean EDE and variance VDE.
	      double SIGMA = sqrt(VDE);
	      if(SIGMA < 0.333333333E0*EDE)
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

void pen_betaP::JUMPF(double &DS,
		      pen_rand& penRand,
		      const double DSMAX)
{
  //  ************  Positrons.

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
      if(context.FORCE[state.IBODY][PEN_POSITRON][KCOL] > 1.0 && P[KCOL] > 1.0E-16)
	{
	  P0[KCOL] = POR[KCOL];
	  P[KCOL] = POR[KCOL]*context.FORCE[state.IBODY][PEN_POSITRON][KCOL];
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
  ST += P[4];
  
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

      double FSEDE = 1.0-mat.DW1PL[KE]*EDE0;
      if(FSEDE < 0.75){ FSEDE = 0.75;}

      double FSVDE = 1.0-mat.DW2PL[KE]*EDE0;
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
      double STLWR = exp(mat.SPTOT[KE1]+(mat.SPTOT[KE1+1]-mat.SPTOT[KE1])*XEK1);
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
	      double FSEDE = 1.0-mat.DW1PL[KE]*EDE0;
	      if(FSEDE < 0.75){ FSEDE = 0.75;}
		    
	      double FSVDE = 1.0-mat.DW2PL[KE]*EDE0;
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

void pen_betaP::KNOCK(double &DE,
		      int &ICOL,
		      pen_rand& penRand)
{
  //  ************  Positrons (KPAR=PEN_POSITRON).

  const double TREV=2.0E0*constants::REV;
  
  //Get material
  const pen_material& mat = *pmat;
    
  if(MHINGE == 1){}
  else
    {
  
      //  ****  Hinge, artificial soft event (ICOL=BETAp_SOFT_INTERACTION).
  
      ICOL = BETAp_SOFT_INTERACTION;
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
	    PANaR(this->state,mat.EABS[PEN_PHOTON],stackG,penRand);  // Annihilation at rest.
	    if(penGetError() != PEN_SUCCESS){ return;}	  
	    DE = state.E+TREV;
	    state.E = 0.0;
	    return;
	  }
	
	  if(KSOFTE == 0){ return;}

	  //Get energy grid values at the hinge to
	  //calculate angular deflection
	  XEL = log(state.E);
	  XE = (XEL-context.grid.DLEMP1)*context.grid.DLFC;
	  KE = (int)XE;
	  XEK = XE-(double)KE;	  	  
	}
      else
	{
	  //Energy at the end of step is the actual one
	  EENDSTEP = state.E;
	}

      //  ****  Angular deflection.

      if(mat.T1P[KE+1] > -78.3)
	{
	  T1 = exp(mat.T1P[KE]+(mat.T1P[KE+1]-mat.T1P[KE])*XEK);
	  T2 = exp(mat.T2P[KE]+(mat.T2P[KE+1]-mat.T2P[KE])*XEK);
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
  //  ****  A delta interaction (ICOL=BETAp_DELTA) occurs when the maximum
  //        allowed step length is exceeded.
  if(KDELTA == 1)
    {
      ICOL = BETAp_DELTA;
      DE = 0.0;
      return;
    }
  //  ****  Random sampling of the interaction type.

  double STNOW = P[0];
  STNOW += P[1];
  STNOW += P[2];
  STNOW += P[3];
  STNOW += P[4];
  
  
  double STS = (STNOW > ST ? STNOW : ST)*penRand.rand();
  
  double SS = P[BETAp_HARD_ELASTIC];
  if(SS > STS)
    {
      ICOL = HelasticCol.interact(context,mat,(*this),DE,penRand);
      return;
    }
  
  SS += P[BETAp_HARD_INELASTIC];
  if(SS > STS)
    {
      ICOL = HinelastCol.interact(context,mat,(*this),DE,penRand);
      return;
    }

  SS += P[BETAp_HARD_BREMSSTRAHLUNG];
  if(SS > STS)
    {
      ICOL = Hbremsstrahlung.interact(context,mat,(*this),DE,penRand);
      return;
    }  

  SS += P[BETAp_HARD_INNER_SHELL];
  if(SS > STS)
    {
      ICOL = HinnerShell.interact(context,mat,(*this),DE,penRand);
      return;
    }

  SS += P[BETAp_ANNIHILATION];
  if(SS > STS)
    {
      ICOL = annihilation.interact(context,mat,(*this),DE,penRand);
      return;
    }
  
  //  ****  A delta interaction (ICOL=BETAp_DELTA) may occur when the total
  //        interaction probability per unit path length, ST, is
  //        larger than STNOW.
  ICOL = BETAp_DELTA;
  DE = 0.0;
}

void pen_betaP::KNOCKF(double &DE,
		       int &ICOL,
		       pen_rand& penRand)
{

  //  ************  Positrons (KPAR=PEN_POSITRON).

  const double TREV=2.0E0*constants::REV;
  
  //Get material
  const pen_material& mat = *pmat;
  
  if(MHINGE == 1){}
  else
    {

      //  ****  Hinge, artificial soft event (ICOL=1).
	  
      ICOL = BETAp_SOFT_INTERACTION;
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
	    // Annihilation at rest.
	    PANaR(this->state,mat.EABS[PEN_PHOTON],stackG,penRand); 
	    if(penGetError() != PEN_SUCCESS){ return;}	  
	    DE = state.E+TREV;
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
	  T1 = exp(mat.T1P[KE]+(mat.T1P[KE+1]-mat.T1P[KE])*XEK);
	  T2 = exp(mat.T2P[KE]+(mat.T2P[KE+1]-mat.T2P[KE])*XEK);
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
      ICOL = BETAp_DELTA;
      DEA = 0.0;
      DE = 0.0;
      return;
    }
  //  ****  Interaction forcing.

  for(unsigned int i = 0; i < interactions; i++)
    {
      if(context.FORCE[state.IBODY][getKpar()][i] > 1.0 && POR[i] > 1.0E-16)
	{
	  P0[i] = POR[i];
	  P[i] = POR[i]*context.FORCE[state.IBODY][getKpar()][i];
	  LFORC[i] = true;
	}
      else
	{
	  LFORC[i] = false;
	}      
    }
  
  IBR = (context.IBRSPL[state.IBODY] > 1 ? context.IBRSPL[state.IBODY] : 1);

  //  ****  Random sampling of the interaction type.

  double STNOW = P[0];
  STNOW += P[1];
  STNOW += P[2];
  STNOW += P[3];
  STNOW += P[4];
  
  
  double STS = (STNOW > ST ? STNOW : ST)*penRand.rand();

  double SS = P[BETAp_HARD_ELASTIC];
  if(SS > STS)
    {
      ICOL = HelasticCol.interactF(context,mat,(*this),DE,penRand);
      return;
    }
  
  SS += P[BETAp_HARD_INELASTIC];
  if(SS > STS)
    {
      ICOL = HinelastCol.interactF(context,mat,(*this),DE,penRand);
      return;
    }

  SS += P[BETAp_HARD_BREMSSTRAHLUNG];
  if(SS > STS)
    {
      ICOL = Hbremsstrahlung.interactF(context,mat,(*this),DE,penRand);
      return;
    }  

  SS += P[BETAp_HARD_INNER_SHELL];
  if(SS > STS)
    {
      ICOL = HinnerShell.interactF(context,mat,(*this),DE,penRand);
      return;
    }

  SS += P[BETAp_ANNIHILATION];
  if(SS > STS)
    {
      ICOL = annihilation.interactF(context,mat,(*this),DE,penRand);
      return;
    }

  
  //  ****  A delta interaction (ICOL=7) may occur when the total
  //        interaction probability per unit path length, ST, is
  //        larger than STNOW.
  ICOL = BETAp_DELTA;
  DEA = 0.0;
  DE = 0.0;
  return;
}

void pen_betaP::softEloss(double& X,
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

void pen_betaP::dpage()
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

void pen_betaP::annihilate(pen_rand& penRand){

  //Get material
  PANaR(state,pmat->EABS[PEN_PHOTON],stackG,penRand);
}


// Interactions

// Positron elastic collision

int pen_pHEC::interact(const pen_context& context, const pen_material& mat, pen_betaP& beta, double& DE, pen_rand& penRand) const
{

  //  ****  Hard elastic collision (ICOL=BETAp_HARD_ELASTIC).
  DE = 0.0;
  double RMU;
  if(beta.state.E >= mat.PELMAX)
    {
      double TRNDC = mat.RNDCP[beta.KE]+(mat.RNDCP[beta.KE+1]-mat.RNDCP[beta.KE])*beta.XEK;
      double TA = exp(mat.AP[beta.KE]+(mat.AP[beta.KE+1]-mat.AP[beta.KE])*beta.XEK);
      double TB = mat.BP[beta.KE]+(mat.BP[beta.KE+1]-mat.BP[beta.KE])*beta.XEK;
      PELa(TA,TB,TRNDC,RMU,penRand);
    }
  else
    {
      double TRNDC = mat.RNDCPd[beta.KE]+(mat.RNDCPd[beta.KE+1]-mat.RNDCPd[beta.KE])*beta.XEK;
      PELd(mat,beta.KE,beta.XEL,TRNDC,RMU,context.grid.DLEMP,context.grid.DLFC,penRand);  // Uses the ELSEPA database.
    }
  double CDT = 1.0-(RMU+RMU);
  double DF = constants::TWOPI*penRand.rand();
  DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
  DE = 0.0;
  return BETAp_HARD_ELASTIC;
}

int pen_pHEC::interactF(const pen_context& context,
			const pen_material& mat,
			pen_betaP& beta,
			double& DE,
			pen_rand& penRand) const
{

  //  ****  Hard elastic collision (ICOL=BETAp_HARD_ELASTIC).
  DE = 0.0;
  double RMU;
  int ICOL = BETAp_HARD_ELASTIC;
  if(beta.state.E >= mat.PELMAX)
    {
      double TRNDC = mat.RNDCP[beta.KE]+(mat.RNDCP[beta.KE+1]-mat.RNDCP[beta.KE])*beta.XEK;
      double TA = exp(mat.AP[beta.KE]+(mat.AP[beta.KE+1]-mat.AP[beta.KE])*beta.XEK);
      double TB = mat.BP[beta.KE]+(mat.BP[beta.KE+1]-mat.BP[beta.KE])*beta.XEK;
      PELa(TA,TB,TRNDC,RMU,penRand);
    }
  else
    {
      double TRNDC = mat.RNDCPd[beta.KE]+(mat.RNDCPd[beta.KE+1]-mat.RNDCPd[beta.KE])*beta.XEK;
      PELd(mat,beta.KE,beta.XEL,TRNDC,RMU,context.grid.DLEMP,context.grid.DLFC,penRand);  // Uses the ELSEPA database.
    }
  double CDT = 1.0-(RMU+RMU);
  double DF = constants::TWOPI*penRand.rand();
  DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
  beta.DEA = 0.0;
  return ICOL;
}

// Positron inelastic collision

int pen_pHIC::interact(const pen_context& context,
		       const pen_material& mat,
		       pen_betaP& beta,
		       double& DE,
		       pen_rand& penRand) const
{
  //  ****  Hard inelastic collision (ICOL=BETAp_HARD_INELASTIC).
  double EP, CDT, ES, CDTS;
  int IOSC;
  const double TREV = 2.0*constants::REV;
  const double DELTA = mat.DEL[beta.KE]+(mat.DEL[beta.KE+1]-mat.DEL[beta.KE])*beta.XEK;
  PINa(context,mat,beta.state.E,beta.KE,beta.XEL,DELTA,DE,EP,CDT,ES,CDTS,IOSC,penRand);
  //  ****  Scattering angles (primary positron).
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
      storeState.ILB[2] = BETAp_HARD_INELASTIC;
      storeState.ILB[3] = 0;
      storeState.ILB[4] = beta.state.ILB[4];

      storeState.LAGE = beta.state.LAGE;
      storeState.PAGE = beta.state.PAGE;
      
      DIRECT(CDTS,DFS,storeState.U,storeState.V,storeState.W);

      stack.store(storeState);

      if(penGetError() != PEN_SUCCESS){ return BETAp_HARD_INELASTIC;}
    }
  //  ****  New energy and direction.
  if(EP > mat.EABS[PEN_POSITRON])
    {
      beta.state.E = EP;
      DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
    }
  else
    {
      PANaR(beta.state,mat.EABS[PEN_PHOTON],beta.stackG,penRand);  // Annihilation at rest.
      if(penGetError() != PEN_SUCCESS){ return BETAp_HARD_INELASTIC;}
      DE = beta.state.E+TREV;
      beta.state.E = 0.0;
    }
  return BETAp_HARD_INELASTIC;
}

int pen_pHIC::interactF(const pen_context& context,
			const pen_material& mat,
			pen_betaP& beta,
			double& DE,
			pen_rand& penRand) const
{

  //  ****  Hard inelastic collision (ICOL=BETAp_HARD_INELASTIC).
  double EP, CDT, ES, CDTS;
  double WFORCE = 1.0;
  int IOSC;
  bool LCOL;
  int ICOL = BETAp_HARD_INELASTIC;
  if(beta.LFORC[BETAp_HARD_INELASTIC])  // Forced interaction.
    {
      WFORCE = beta.P0[BETAp_HARD_INELASTIC]/beta.P[BETAp_HARD_INELASTIC];
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

  const double TREV = 2.0*constants::REV;
  const double DELTA = mat.DEL[beta.KE]+(mat.DEL[beta.KE+1]-mat.DEL[beta.KE])*beta.XEK;
  PINa(context,mat,beta.state.E,beta.KE,beta.XEL,DELTA,DE,EP,CDT,ES,CDTS,IOSC,penRand);

  //  ****  Scattering angles (primary positron).
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
      if(EP > mat.EABS[PEN_POSITRON])
	{
	  beta.state.E = EP;
	  DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
	}
      else
	{
	  PANaR(beta.state,mat.EABS[PEN_PHOTON],beta.stackG,penRand);  // Annihilation at rest.
	  if(penGetError() != PEN_SUCCESS){ return ICOL;}
	  beta.DEA = beta.DEA+EP+TREV;
	  DE = DE+EP+TREV;
	  beta.state.E = 0.0;
	}
    }
  return ICOL;
}
  
// Positron bremsstrahlung emission
pen_pHBE::pen_pHBE(pen_particleStack<pen_state_gPol>& stackIn) :
  abc_interaction(BETAp_HARD_BREMSSTRAHLUNG),
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

int pen_pHBE::interact(const pen_context& /*context*/,
		       const pen_material& mat,
		       pen_betaP& beta,
		       double& DE,
		       pen_rand& penRand) const
{
    
  //  ****  Hard bremsstrahlung emission (ICOL=BETAp_HARD_BREMSSTRAHLUNG).
  const double TREV = 2.0*constants::REV;
  PBRa(mat,beta.state.E,beta.KE,beta.XEK,DE,constants::WB,penRand);
  
  //  ****  Bremsstrahlung photon.
  if(DE > mat.EABS[PEN_PHOTON])
    {      
      double CDTS;
      PBRaA(mat,beta.state.E,DE,CDTS,BET,penRand);      
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
      storeState.ILB[2] = BETAp_HARD_BREMSSTRAHLUNG;
      storeState.ILB[3] = 0;
      storeState.ILB[4] = beta.state.ILB[4];

      storeState.LAGE = beta.state.LAGE;
      storeState.PAGE = beta.state.PAGE;
      
      DIRECT(CDTS,DFS,storeState.U,storeState.V,storeState.W);
      
      stack.store(storeState);
      
      if(penGetError() != PEN_SUCCESS){ return BETAp_HARD_BREMSSTRAHLUNG;}
    }
  //  ****  New energy.
  beta.state.E = beta.state.E-DE;
  if(beta.state.E < mat.EABS[PEN_POSITRON])
    {
      PANaR(beta.state,mat.EABS[PEN_PHOTON],beta.stackG,penRand);  // Annihilation at rest.
      if(penGetError() != PEN_SUCCESS){ return BETAp_HARD_BREMSSTRAHLUNG;}
      DE = beta.state.E+DE+TREV;
      beta.state.E = 0.0;
    }  
  return BETAp_HARD_BREMSSTRAHLUNG;
}

int pen_pHBE::interactF(const pen_context& /*context*/,
			const pen_material& mat,
			pen_betaP& beta,
			double& DE,
			pen_rand& penRand) const
{

  //  ****  Hard bremsstrahlung emission (ICOL=BETAp_HARD_BREMSSTRAHLUNG).

  const double TREV = 2.0*constants::REV;

  double WFORCE = 1.0;
  bool LCOL;
  int ICOL = BETAp_HARD_BREMSSTRAHLUNG;
  if(beta.LFORC[BETAp_HARD_BREMSSTRAHLUNG])  // Forced interaction.
    {
      WFORCE = beta.P0[BETAp_HARD_BREMSSTRAHLUNG]/beta.P[BETAp_HARD_BREMSSTRAHLUNG];
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
	  
  PBRa(mat,beta.state.E,beta.KE,beta.XEK,DE,constants::WB,penRand);
  //  ****  Bremsstrahlung photon.
  if(DE > mat.EABS[PEN_PHOTON])
    {
      double CDTS;
      PBRaA(mat,beta.state.E,DE,CDTS,BET,penRand);
      double WSPLIT = beta.state.WGHT*WFORCE/beta.IBR;
      for(int I = 0; I < beta.IBR; I++)
	{
	  double DFS = constants::TWOPI*penRand.rand();

	  pen_state_gPol storeState;

	  storeState.E = DE;
	  
	  storeState.X = beta.state.X;
	  storeState.Y = beta.state.Y;
	  storeState.Z = beta.state.Z;

	  storeState.U = beta.state.U;
	  storeState.V = beta.state.V;
	  storeState.W = beta.state.W;

	  storeState.WGHT = WSPLIT;
	  storeState.IBODY = beta.state.IBODY;
	  storeState.MAT = beta.state.MAT;
	  
	  storeState.ILB[0] = beta.state.ILB[0]+1;
	  storeState.ILB[1] = beta.getKpar();
	  storeState.ILB[2] = ICOL;
	  storeState.ILB[3] = 0;
	  storeState.ILB[4] = beta.state.ILB[4];

	  storeState.LAGE = beta.state.LAGE;
	  storeState.PAGE = beta.state.PAGE;
	  
	  DIRECT(CDTS,DFS,storeState.U,storeState.V,storeState.W);

	  stack.store(storeState);	  
	  if(penGetError() != PEN_SUCCESS){ return ICOL;}
	}
    }
  //  ****  New energy.
  beta.DEA = DE;
  DE = DE*WFORCE;
  if(LCOL)
    {
      beta.state.E = beta.state.E-beta.DEA;
      if(beta.state.E < mat.EABS[beta.getKpar()])
	{
	  PANaR(beta.state,mat.EABS[PEN_PHOTON],beta.stackG,penRand);  // Annihilation at rest.
	  if(penGetError() != PEN_SUCCESS){ return ICOL;}
	  beta.DEA = beta.DEA+beta.state.E+TREV;
	  DE = DE+beta.state.E+TREV;
	  beta.state.E = 0.0;
	}
    }
  return ICOL;	  
}
  
// Positron inner-shell impact ionisation

int pen_pSII::interact(const pen_context& context,
		       const pen_material& mat,
		       pen_betaP& beta,
		       double& DE,
		       pen_rand& penRand) const
{

  //  ****  Ionisation of an inner shell (ICOL=BETAp_HARD_INNER_SHELL).

  const double TREV = 2.0*constants::REV;
  
  double EP,CDT,ES,CDTS;
  int IZA,ISA;
  double const DELTA = mat.DEL[beta.KE]+(mat.DEL[beta.KE+1]-mat.DEL[beta.KE])*beta.XEK;
  PSIa(context,mat,beta.state.E,beta.KE,beta.XEL,DELTA,DE,EP,CDT,ES,CDTS,IZA,ISA,penRand);
  //  ****  Atomic relaxation.
  if(IZA > 2)
    {
      int ks;
      RELAX(context.elements,mat,beta.state,BETAp_HARD_INNER_SHELL,1,
	    IZA,ISA,ks,stackE,stackG,penRand);
      if(penGetError() != PEN_SUCCESS){ return BETAp_HARD_INNER_SHELL;}
    }
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

      stateStore.WGHT = beta.state.WGHT;
      stateStore.IBODY = beta.state.IBODY;
      stateStore.MAT = beta.state.MAT;
      
      stateStore.ILB[0] = beta.state.ILB[0]+1;
      stateStore.ILB[1] = beta.getKpar();
      stateStore.ILB[2] = BETAp_HARD_INNER_SHELL;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = beta.state.ILB[4];

      stateStore.LAGE = beta.state.LAGE;
      stateStore.PAGE = beta.state.PAGE;
      
      DIRECT(CDTS,DFS,stateStore.U,stateStore.V,stateStore.W);

      stackE.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return BETAp_HARD_INNER_SHELL;}
    }
  //  ****  New energy and direction.
  if(EP > mat.EABS[PEN_POSITRON])
    {
      beta.state.E = EP;
      DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
    }
  else
    {
      PANaR(beta.state,mat.EABS[PEN_PHOTON],beta.stackG,penRand);  // Annihilation at rest.
      if(penGetError() != PEN_SUCCESS){ return BETAp_HARD_INNER_SHELL;}
      DE = beta.state.E+TREV;
      beta.state.E = 0.0;
    }
  return BETAp_HARD_INNER_SHELL;
}

int pen_pSII::interactF(const pen_context& context,
			const pen_material& mat,
			pen_betaP& beta,
			double& DE,
			pen_rand& penRand) const
{

  //  ****  Ionization of an inner shell (ICOL=BETAp_HARD_INNER_SHELL).

  const double TREV = 2.0*constants::REV;
  
  double EP,CDT,ES,CDTS;
  double WFORCE = 1.0;
  int IZA,ISA;
  int ICOL = BETAp_HARD_INNER_SHELL;
  bool LCOL;
  if(beta.LFORC[BETAp_HARD_INNER_SHELL])  // Forced interaction.
    {
      WFORCE = beta.P0[BETAp_HARD_INNER_SHELL]/beta.P[BETAp_HARD_INNER_SHELL];
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

  const double DELTA = mat.DEL[beta.KE]+(mat.DEL[beta.KE+1]-mat.DEL[beta.KE])*beta.XEK;
  PSIa(context,mat,beta.state.E,beta.KE,beta.XEL,DELTA,DE,EP,CDT,ES,CDTS,IZA,ISA,penRand);
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
      RELAX(context.elements,mat,beta.state,ICOL,1,
	    IZA,ISA,ks,stackE,stackG,penRand);      
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
      beta.state.WGHT = WGHTA;
    }
  //  ****  New energy and direction.
  beta.DEA = DE;
  DE = beta.DEA*WFORCE;
  if(LCOL)
    {
      if(EP > mat.EABS[PEN_POSITRON])
	{
	  beta.state.E = EP;
	  DIRECT(CDT,DF,beta.state.U,beta.state.V,beta.state.W);
	}
      else
	{
	  PANaR(beta.state,mat.EABS[PEN_PHOTON],beta.stackG,penRand);  // Annihilation at rest.
	  if(penGetError() != PEN_SUCCESS){ return ICOL;}
	  beta.DEA = beta.DEA+EP+TREV;
	  DE = DE+EP+TREV;
	  beta.state.E = 0.0;
	}
    }
  return ICOL;  
}
  
// Positron annihilation

int pen_pHAN::interact(const pen_context& /*context*/, const pen_material& mat, pen_betaP& beta, double& DE, pen_rand& penRand) const
{
  //  ****  Positron annihilation in flight (ICOL=BETAp_ANNIHILATION).
  const double TREV = 2.0*constants::REV;
  
  double E1,CDT1,E2,CDT2;
  PANa(mat.EABS[beta.getKpar()],beta.state.E,E1,CDT1,E2,CDT2,penRand);
  double DF = constants::TWOPI*penRand.rand();
  if(E1 > mat.EABS[PEN_PHOTON])
    {
      pen_state_gPol stateStore;

      stateStore.E = E1;
      
      stateStore.X = beta.state.X;
      stateStore.Y = beta.state.Y;
      stateStore.Z = beta.state.Z;      
      
      stateStore.U = beta.state.U;
      stateStore.V = beta.state.V;
      stateStore.W = beta.state.W;

      stateStore.WGHT = beta.state.WGHT;
      stateStore.IBODY = beta.state.IBODY;
      stateStore.MAT = beta.state.MAT;
     
      stateStore.ILB[0] = beta.state.ILB[0]+1;
      stateStore.ILB[1] = beta.getKpar();
      stateStore.ILB[2] = BETAp_ANNIHILATION;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = beta.state.ILB[4];

      stateStore.LAGE = beta.state.LAGE;
      stateStore.PAGE = beta.state.PAGE;

      DIRECT(CDT1,DF,stateStore.U,stateStore.V,stateStore.W);

      stack.store(stateStore);
      
      if(penGetError() != PEN_SUCCESS){ return BETAp_ANNIHILATION;}
    }
  if(E2 > mat.EABS[PEN_PHOTON])
    {

      pen_state_gPol stateStore;

      DF = DF+constants::PI;
      stateStore.E = E2;
      
      stateStore.X = beta.state.X;
      stateStore.Y = beta.state.Y;
      stateStore.Z = beta.state.Z;      
      
      stateStore.U = beta.state.U;
      stateStore.V = beta.state.V;
      stateStore.W = beta.state.W;

      stateStore.WGHT = beta.state.WGHT;
      stateStore.IBODY = beta.state.IBODY;
      stateStore.MAT = beta.state.MAT;
     
      stateStore.ILB[0] = beta.state.ILB[0]+1;
      stateStore.ILB[1] = beta.getKpar();
      stateStore.ILB[2] = BETAp_ANNIHILATION;
      stateStore.ILB[3] = 0;
      stateStore.ILB[4] = beta.state.ILB[4];

      stateStore.LAGE = beta.state.LAGE;
      stateStore.PAGE = beta.state.PAGE;

      DIRECT(CDT2,DF,stateStore.U,stateStore.V,stateStore.W);

      stack.store(stateStore);
      
      if(penGetError() != PEN_SUCCESS){ return BETAp_ANNIHILATION;}
    }
  DE = beta.state.E+TREV;
  beta.state.E = 0.0;
  return BETAp_ANNIHILATION;
}

int pen_pHAN::interactF(const pen_context& /*context*/, const pen_material& mat, pen_betaP& beta, double& DE, pen_rand& penRand) const
{

  //  ****  Positron annihilation in flight (ICOL=BETAp_ANNIHILATION).

  const double TREV = 2.0*constants::REV;
  
  int ICOL = BETAp_ANNIHILATION;
  double WFORCE = 1.0;
  bool LCOL;
  if(beta.LFORC[BETAp_ANNIHILATION])  // Forced interaction.
    {
      WFORCE = beta.P0[BETAp_ANNIHILATION]/beta.P[BETAp_ANNIHILATION];
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
	  
  double E1,CDT1,E2,CDT2;  
  PANa(mat.EABS[beta.getKpar()],beta.state.E,E1,CDT1,E2,CDT2,penRand);
  double DF = constants::TWOPI*penRand.rand();
  if(E1 > mat.EABS[PEN_PHOTON])
    {
      pen_state_gPol stateStore;

      stateStore.E = E1;
      
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
      
      DIRECT(CDT1,DF,stateStore.U,stateStore.V,stateStore.W);

      stack.store(stateStore);
      
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
    }
  if(E2 > mat.EABS[PEN_PHOTON])
    {
      pen_state_gPol stateStore;

      DF = DF+constants::PI;

      stateStore.E = E2;

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
      
      DIRECT(CDT2,DF,stateStore.U,stateStore.V,stateStore.W);
     
      stack.store(stateStore);
      if(penGetError() != PEN_SUCCESS){ return ICOL;}
    }
  beta.DEA = beta.state.E+TREV;
  DE = beta.DEA*WFORCE;
  if(LCOL){ beta.state.E = 0.0;}
  return ICOL;	  
}


// Auxiliar functions

//  *********************************************************************
//                       SUBROUTINE PANaR 
//  *********************************************************************
void PANaR(const pen_particleState& betaState,
	   const double ECUT,
	   pen_particleStack<pen_state_gPol>& stackG,
	   pen_rand& penRand)
{
  //  Simulation of annihilation of positrons at rest. Annihilation quanta
  //  are stored in the secondary stack only when ECUT is less than REV.


  if(constants::REV < ECUT){ return;}

  pen_state_gPol storeState;

  storeState.E = constants::REV;
  
  storeState.X = betaState.X;
  storeState.Y = betaState.Y;
  storeState.Z = betaState.Z;
  
  storeState.U = betaState.U;
  storeState.V = betaState.V;
  storeState.W = betaState.W;

  storeState.WGHT = betaState.WGHT;
  storeState.IBODY = betaState.IBODY;
  storeState.MAT = betaState.MAT;
  
  double CDT1 = -1.0+2.0*penRand.rand();
  double DF = constants::TWOPI*penRand.rand();
  
  storeState.ILB[0] = betaState.ILB[0]+1;
  storeState.ILB[1] = PEN_POSITRON;
  storeState.ILB[2] = BETAp_ANNIHILATION;
  storeState.ILB[3] = 0;
  storeState.ILB[4] = betaState.ILB[4];

  storeState.LAGE = betaState.LAGE;
  storeState.PAGE = betaState.PAGE;

  DIRECT(CDT1,DF,storeState.U,storeState.V,storeState.W);
  
  stackG.store(storeState);
  if(penGetError() != PEN_SUCCESS){ return;}

  storeState.U = -storeState.U;
  storeState.V = -storeState.V;
  storeState.W = -storeState.W;
  
  stackG.store(storeState);
  if(penGetError() != PEN_SUCCESS){ return;}

}

//  *********************************************************************
//                       SUBROUTINE PELd
//  *********************************************************************
void PELd(const pen_material& mat,
	  const unsigned KE,
	  const double XEL,
	  double &RNDC,
	  double &RMU,
	  const double* DLEMP,
	  const double DLFC,
	  pen_rand& penRand)
{
  //     Simulation of positron hard elastic events. Cross sections from
  //  the ELSEPA numerical database.

  //  Argument value:
  //    RNDC ... cutoff value of the uniform random number
  //             (only hard events are simulated).
  //    RMU .... sampled angular deflection, =(1-CDT)/2.


  const unsigned int NPM1= mat.NPM1_P;
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
  int I = mat.ITLP[ITN][JE]-1;
  int J = mat.ITUP[ITN][JE]-1;
  if(J-I < 2){}
  else
  {
    while(J-I > 1)
    {
      int K=(I+J)/2;
      if(RU > mat.PSP[K][JE])
      {
        I = K;
      }
      else
      {
        J = K;
      }
    }
  }
  //  ****  Sampling from the rational inverse cumulative distribution.
  double PP = mat.PSP[I][JE];
  double RR = RU-PP;
  if(RR > 1.0E-16)
  {
    double XX = mat.XSP[I][JE];
    double AA = mat.ASP[I][JE];
    double BB = mat.BSP[I][JE];
    double D = mat.PSP[I+1][JE]-PP;
    RMU = XX+((1.0+AA+BB)*D*RR/(D*D+(AA*D+BB*RR)*RR))
      *(mat.XSP[I+1][JE]-XX);
  }
  else
  {
    RMU = mat.XSP[I][JE];
  }
}

//  *********************************************************************
//                       SUBROUTINE PINa
//  *********************************************************************
void PINa(const pen_context& context,
	  const pen_material& mat,
	  const double E,
	  const unsigned KE,
	  const double XEL,
	  const double DELTA,
	  double &DE,
	  double &EP,
	  double &CDT,
	  double &ES,
	  double &CDTS,
	  int &IOSC,
	  pen_rand& penRand)
{
  //  Random sampling of hard inelastic collisions of positrons.
  //
  //  Sternheimer-Liljequist GOS model.
  //
  //  Input arguments:
  //    E ....... positron energy (eV).
  //    M ....... material where positrons propagate.
  //    DELTA ... Fermi's density effect correction.
  //  Output arguments:
  //    DE ...... energy loss (eV).
  //    EP ...... energy of the scattered positron (eV).
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
  double PK;
  int JE;
  PK = (XEL-context.grid.DLEMP[KE])*context.grid.DLFC;
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
  int JO = mat.NPIN;
  while(JO-IO > 1)
  {
    int IT=(IO+JO)/2;
    if(TST > mat.PINAC[JE][IT])
    {
      IO = IT;
    }
    else
    {
      JO = IT;
    }
  }
  IOSC = mat.IPIN[IO];

  double UK = mat.UI[IOSC-1];
  double WK = mat.WRI[IOSC-1];
  double WTHR;
  if(UK > 1.0E-3)
  {
    if(WCCM > UK){ WTHR = WCCM;}
    else{ WTHR = UK;}
  }
  else
  {
    if(WCCM > WK){ WTHR = WCCM;}
    else{ WTHR = WK;}
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
  //
  //  ****  Trick: The resonance energy and the cutoff recoil energy of
  //        inner shells are varied to yield a smooth threshold.

  LDIST = true;
  double WM, WKP, QKP;
  double WCMAX, WDMAX;
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
    WCMAX = E;
    if(WM < WCMAX){ WDMAX = WM;}else{ WDMAX = WCMAX;}
    if(WTHR > WDMAX){ LDIST = false;}
  }
  else
  {
    if(WCCM > WK){ LDIST = false;}
    WKP = WK;
    QKP = WK;
    WM = E;
    WCMAX = E;
    WDMAX = WKP+1.0;
  }

  //  ****  Constants.

  const double RB = E+TREV;
  const double GAM = 1.0+E*RREV;
  const double GAM2 = GAM*GAM;
  const double BETA2 = (GAM2-1.0)/GAM2;
  const double G12 = pow(GAM+1.0,2);
  const double AMOL = pow((GAM-1.0)/GAM,2);
  const double BHA1 = AMOL*(2.0*G12-1.0)/(GAM2-1.0);
  const double BHA2 = AMOL*(3.0+1.0/G12);
  const double BHA3 = AMOL*2.0*GAM*(GAM-1.0)/G12;
  const double BHA4 = AMOL*pow(GAM-1.0,2)/G12;
  const double CPS = E*RB;
  const double CP = sqrt(CPS);

  //  ************  Partial cross sections of the active oscillator.

  //  ****  Distant excitations.
  double CPPS, CPP, QM, RWKP, XHDL, XHDT, F0;
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
      RWKP = 1.0/WKP;
      XHDL = log(QKP*(QM+TREV)/(QM*(QKP+TREV)))*RWKP;

      XHDT = log(GAM2)-BETA2-DELTA;
      if(XHDT < 0.0){ XHDT = 0.0;}
    
      XHDT = XHDT*RWKP;
      if(UK > 1.0E-3)
      {
        F0 = (WDMAX-WTHR)*(WM+WM-WDMAX-WTHR)/pow(WM-UK,2);
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
    CPP = 0.0;   // Defined to avoid compilation warnings.
    CPPS = 0.0;  // Defined to avoid compilation warnings.
    XHDL = 0.0;
    XHDT = 0.0;
  }
  //  ****  Close collisions.
  double RCL = WTHR/E;
  double RL1 = 1.0-RCL;
  double XHC =((1.0/RCL-1.0)+BHA1*log(RCL)+BHA2*RL1+(BHA3/2.0)*(pow(RCL,2)-1.0)+(BHA4/3.0)*(1.0-pow(RCL,3)))/E;
  
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
  double RK, PHI;
  if(TST < TS1)
  {
    bool Eixir = false;
    while(!Eixir)
    {
      Eixir = true;
      RK = RCL/(1.0-penRand.rand()*RL1);
      PHI = 1.0-RK*(BHA1-RK*(BHA2-RK*(BHA3-BHA4*RK)));
      if(penRand.rand() > PHI){ Eixir = false; continue;}
    }
    //  ****  Energy and scattering angle (primary positron).
    DE = RK*E;
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
//                       SUBROUTINE PSIa
//  *********************************************************************
void PSIa(const pen_context& context,
	  const pen_material& mat,
	  const double E,
	  const unsigned KE,
	  const double XEL,
	  const double DELTA,
	  double &DE,
	  double &EP,
	  double &CDT,
	  double &ES,
	  double &CDTS,
	  int &IZZ,
	  int &ISH,
	  pen_rand& penRand)
{
  //  Random sampling of inner-shell ionisation by positron impact.

  //  Sternheimer-Liljequist GOS model.

  //  Input arguments:
  //    E ....... positron energy (eV).
  //    M ....... material where positrons propagate.
  //    DELTA ... Fermi's density effect correction.
  //  Output arguments:
  //    DE ...... energy loss (eV).
  //    EP ...... energy of the scattered positron (eV).
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
  double PK =(XEL-context.grid.DLEMP[KE])*context.grid.DLFC;
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
  int IT;
  while(JO-IO > 1)
    {
      IT = (IO+JO)/2;
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
  //
  double WM = 3.0*WK-2.0*UK;
  double WKP, QKP, WCMAX, WDMAX;
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
  
  //  ****  Constants.

  const double RB = E+TREV;
  const double GAM = 1.0+E*RREV;
  const double GAM2 = GAM*GAM;
  const double BETA2 = (GAM2-1.0)/GAM2;
  const double G12 = pow(GAM+1.0,2);
  const double AMOL = pow((GAM-1.0)/GAM,2);
  const double BHA1 = AMOL*(2.0*G12-1.0)/(GAM2-1.0);
  const double BHA2 = AMOL*(3.0+1.0/G12);
  const double BHA3 = AMOL*2.0*GAM*(GAM-1.0)/G12;
  const double BHA4 = AMOL*pow(GAM-1.0,2)/G12;
  const double CPS = E*RB;
  const double CP = sqrt(CPS);

  //  ************  Partial cross sections of the active oscillator.

  //  ****  Distant excitations.
  double CPPS = (E-WKP)*(E-WKP+TREV);
  double CPP = sqrt(CPPS);
  double RWKP, XHDL, XHDT, F0, QM;
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
  double RCL = WTHR/E;
  double RL1 = 1.0-RCL;
  double XHC = ((1.0/RCL-1.0)+BHA1*log(RCL)+BHA2*RL1+(BHA3/2.0)*(pow(RCL,2)-1.0)+(BHA4/3.0)*(1.0-pow(RCL,3)))/E;
  
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

  TST = penRand.rand()*XHTOT;

  //  ****  Hard close collision.

  double TS1 = XHC;
  double RK, PHI;
  if(TST < TS1)
  {
    bool Eixir = false;
    while(!Eixir)
      {
	Eixir = true;
	RK = RCL/(1.0-penRand.rand()*RL1);
	PHI = 1.0-RK*(BHA1-RK*(BHA2-RK*(BHA3-BHA4*RK)));
	if(penRand.rand() > PHI){ Eixir = false; continue;}
      }
    //  ****  Energy and scattering angle (primary positron).
    DE = RK*E;
    EP = E-DE;
    CDT = sqrt(EP*RB/(E*(RB-DE)));
    //  ****  Energy and emission angle of the delta ray.
    ES = DE-UK;
    CDTS = sqrt(DE*RB/(E*(DE+TREV)));
    return;
  }

  //  ****  Hard distant longitudinal interaction.

  TS1 = TS1+XHDL;
  DE = WM-sqrt(pow(WM-WTHR,2)-penRand.rand()*(WDMAX-WTHR)*(WM+WM-WDMAX-WTHR));
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
//                       SUBROUTINE PANa
//  *********************************************************************
void PANa(const double EABS,
	  const double E,
	  double &E1,
	  double &CDT1,
	  double &E2,
	  double &CDT2,
	  pen_rand& penRand)
{
  //  Simulation of positron annihilation (either at rest or in flight) in
  //  material M. Ei and CDTi are the energies and polar direction cosines
  //  of the two annihilation photons.


  const double TREV = 2.0*constants::REV;

  //  ****  Slow positrons (assumed at rest).

  if(E < EABS)
  {
    E1 = 0.5*(E+TREV);
    E2 = E1;
    CDT1 = -1.0+2.0*penRand.rand();
    CDT2 = -CDT1;
  }
  else
  {
    //  ****  Annihilation in flight (two photons with energy and directions
    //        determined from the dcs and energy-momentum conservation).
    double GAM;
    if(E < 1.0){ GAM = 1.0+1.0/constants::REV;}
    else{ GAM = 1.0+E/constants::REV;}
      
    double GAM21 = sqrt(GAM*GAM-1.0);
    double ANI = 1.0+GAM;
    double CHIMIN = 1.0/(ANI+GAM21);
    double RCHI = (1.0-CHIMIN)/CHIMIN;
    double GT0 = ANI*ANI-2.0;
    bool Eixir = false;
    double CHI, GREJ;
    while(!Eixir)
    {
      Eixir = true;
      CHI = CHIMIN*pow(RCHI,penRand.rand());
      GREJ = ANI*ANI*(1.0-CHI)+GAM+GAM-1.0/CHI;
      if(penRand.rand()*GT0 > GREJ){ Eixir = false; continue;}
    }
      
    double DET = E+TREV;
    E1 = CHI*DET;
    CDT1 = (ANI-1.0/CHI)/GAM21;
    double CHIP = 1.0-CHI;
    E2 = DET-E1;
    CDT2 = (ANI-1.0/CHIP)/GAM21;
  }

}

//  *********************************************************************
//                       SUBROUTINE EELa
//  *********************************************************************
void PELa(double A,
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
//                       SUBROUTINE EBRa
//  *********************************************************************
void PBRa(const pen_material& mat,
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
//                       SUBROUTINE EBRaA
//  *********************************************************************
void PBRaA(const pen_material& mat,
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
