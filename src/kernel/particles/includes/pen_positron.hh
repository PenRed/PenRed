
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

 
#ifndef __PEN_POSITRON_
#define __PEN_POSITRON_

class pen_betaP; 

// Interactions

class pen_pHEC final : public abc_interaction<pen_betaP, pen_context, pen_material>{
 private:
 public:
 pen_pHEC() : abc_interaction(BETAp_HARD_ELASTIC){}
  double iMeanFreePath(const pen_material&,
		       const pen_betaP&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaP&,
	       double& DE,
	       pen_rand& penRand) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaP&,
		double&,
		pen_rand& penRand) const final override;  
};

class pen_pHIC final : public abc_interaction<pen_betaP, pen_context, pen_material>{
 private:
  pen_particleStack<pen_particleState>& stack;  //b- stack
 public:

 pen_pHIC(pen_particleStack<pen_particleState>& stackIn)
   : abc_interaction(BETAp_HARD_INELASTIC),
    stack(stackIn){}
  
  double iMeanFreePath(const pen_material&,
		       const pen_betaP&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaP&,
	       double&,
	       pen_rand& penRand) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaP&,
		double&,
		pen_rand& penRand) const final override;  
};

class pen_pHBE final : public abc_interaction<pen_betaP, pen_context, pen_material>{
 private:
  pen_particleStack<pen_state_gPol>& stack; //Gamma stack

  // **** Bremsstrahlung angular distribution parameters
  static const int NET = 7;
  double BET[NET];      
 public:

  //  ****  Bremsstrahlung emission.
  static const double WB[constants::NBW];
    
  pen_pHBE(pen_particleStack<pen_state_gPol>& stackIn);
  
  double iMeanFreePath(const pen_material&,
		       const pen_betaP&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaP&,
	       double&,
	       pen_rand& penRand) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaP&,
		double&,
		pen_rand& penRand) const final override;  
};

class pen_pSII final : public abc_interaction<pen_betaP, pen_context, pen_material>{
 private:
  pen_particleStack<pen_particleState>& stackE; //b- stack
  pen_particleStack<pen_state_gPol>&    stackG; //Gamma stack
 public:

 pen_pSII(pen_particleStack<pen_particleState>& stackEin,
	  pen_particleStack<pen_state_gPol>& stackGin)
   : abc_interaction(BETAp_HARD_INNER_SHELL),
    stackE(stackEin),
    stackG(stackGin){}
  
  double iMeanFreePath(const pen_material&,
		       const pen_betaP&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaP&,
	       double&,
	       pen_rand& penRand) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaP&,
		double&,
		pen_rand& penRand) const final override;  
};

class pen_pHAN final : public abc_interaction<pen_betaP, pen_context, pen_material>{
 private:
  pen_particleStack<pen_state_gPol>& stack; //Gamma stack
 public:
  
 pen_pHAN(pen_particleStack<pen_state_gPol>& stackIn)
   : abc_interaction(BETAp_ANNIHILATION),
    stack(stackIn)
    {}
  
  double iMeanFreePath(const pen_material&,
		       const pen_betaP&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaP&,
	       double&,
	       pen_rand& penRand) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaP&,
		double&,
		pen_rand& penRand) const final override;  
};

class pen_betaP final : public abc_particle<pen_particleState, pen_context, pen_material>{

  friend class pen_pHEC;
  friend class pen_pHIC;
  friend class pen_pHBE;
  friend class pen_pSII;
  friend class pen_pHAN;
  
private:

  pen_pHEC HelasticCol;
  pen_pHIC HinelastCol;
  pen_pHBE Hbremsstrahlung;
  pen_pSII HinnerShell;
  pen_pHAN annihilation;

  //  ****  Random-hinge slowing-down parameters (output from subroutines
  //        JUMP and KNOCK).
  //
  //            <----------------- step ----------------->
  //           hard event -------- hinge ------- hard event
  //            <---- segment 0 ----><---- segment 1 ---->
  //
  //     MHINGE ... labels the two segments (jumps) of a step between
  //                two hard events;
  //                = 0 (1) for the segment before (after) the hinge,
  //     DESOFT ... energy loss along the whole step,
  //     SSOFT .... effective stopping power, = DESOFT/step_length.
  int MHINGE;
  double DESOFT, SSOFT;
  
  //Clase II simulation variables
  double ELAST2, DST, DSR, W1, W2, T1, T2;
  int KSOFTE, KDELTA;
  double EENDSTEP;

  //Variance reduction variables
  double POR[constants::MAXINTERACTIONS];
  int    IBR;
    
  // Stack for annihilation produced photons
  pen_particleStack<pen_state_gPol>& stackG;

  //  ----  Energy deposited in the last event (analogue simulation).
  double DEA;
  
public:

  static const double REV;  // Electron rest energy (eV)
  static const double mc2;
  static const double twomc2;    
  
  pen_betaP(const pen_context& context,
	    pen_particleStack<pen_particleState>& stackEin,
	    pen_particleStack<pen_state_gPol>& stackGin,
        pen_particleStack<pen_particleState>& stackPin);

  void JUMP(double &DS,
	    pen_rand& penRand,
	    const double DSMAX) final override;
  void JUMPF(double &DS,
	     pen_rand& penRand,
	     const double DSMAX) final override;
  void KNOCK(double &DE,
	     int &ICOL,
	     pen_rand& penRand) final override;
  void KNOCKF(double &DE,
	      int &ICOL,
	      pen_rand& penRand) final override;

  void softEloss(double& X,
		 double& Y,
		 double& Z,
		 double& DE,
		 pen_rand& penRand) final override;  
  
  void dpage() final override;
  inline void page0() final override{
    SSOFT = 0.0;
  }
  
  void IMFP(const pen_material& mat);  
  void hingeUpdate(const pen_material& mat);

  void annihilate(pen_rand& penRand) final override;
  
  inline void START() final override{
    MHINGE = 0;
    ELAST1 = state.E+1.0E30;
    ELAST2 = ELAST1;
  }  
};

inline double pen_pHEC::iMeanFreePath(const pen_material& mat,
				      const pen_betaP& particle) const{
  return exp(mat.SPHEL[particle.KE]+mat.DSPHEL[particle.KE]*particle.XEK);
}

inline double pen_pHIC::iMeanFreePath(const pen_material& mat,
				      const pen_betaP& particle) const{
  return exp(mat.SPHIN[particle.KE]+mat.DSPHIN[particle.KE]*particle.XEK);
}

inline double pen_pHBE::iMeanFreePath(const pen_material& mat,
				      const pen_betaP& particle) const{
  return exp(mat.SPHBR[particle.KE]+mat.DSPHBR[particle.KE]*particle.XEK);
}

inline double pen_pSII::iMeanFreePath(const pen_material& mat,
				      const pen_betaP& particle) const{
  return exp(mat.SPISI[particle.KE]+mat.DSPISI[particle.KE]*particle.XEK);
}

inline double pen_pHAN::iMeanFreePath(const pen_material& mat,
				      const pen_betaP& particle) const{
  return exp(mat.SPAN[particle.KE]+mat.DSPAN[particle.KE]*particle.XEK);
}

inline void pen_betaP::IMFP(const pen_material& mat)
{
  P[BETAp_HARD_ELASTIC]        = HelasticCol.iMeanFreePath(mat,(*this));
  P[BETAp_HARD_INELASTIC]      = HinelastCol.iMeanFreePath(mat,(*this));
  P[BETAp_HARD_BREMSSTRAHLUNG] = Hbremsstrahlung.iMeanFreePath(mat,(*this));
  P[BETAp_HARD_INNER_SHELL]    = HinnerShell.iMeanFreePath(mat,(*this));  
  P[BETAp_ANNIHILATION]        = annihilation.iMeanFreePath(mat,(*this));
}

inline void pen_betaP::hingeUpdate(const pen_material& mat)
{
  if(mat.W1P[KE+1] > -78.3)
    {
      W1 = exp(mat.W1P[KE]+(mat.W1P[KE+1]-mat.W1P[KE])*XEK);
      W2 = exp(mat.W2P[KE]+(mat.W2P[KE+1]-mat.W2P[KE])*XEK);
    }
  else
    {
      W1 = 0.0;
      W2 = 0.0;
    }
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
}

// Auxiliar functions

void PELd(const pen_material& mat,
	  const unsigned KE,
	  const double XEL,
	  double &RNDC,
	  double &RMU,
	  const double* DLEMP,
	  const double DLFC,
	  pen_rand& penRand);


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
	  pen_rand& penRand);

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
	  pen_rand& penRand);

void PANa(const double EABS,
	  const double E,
	  double &E1,
	  double &CDT1,
	  double &E2,
	  double &CDT2,
	  pen_rand& penRand);

void PELa(double A,
	  double B,
	  double RNDC,
	  double &RMU,
	  pen_rand& penRand);

void PBRa(const pen_material& mat,
	  const double E,
	  const unsigned KE,
	  const double XEK,
	  double &W,
	  const double* WB,
	  pen_rand& penRand);

void PBRaA(const pen_material& mat,
	   double &E,
	   double &DE,
	   double &CDT,
	   const double* BET,
	   pen_rand& penRand);

#endif
