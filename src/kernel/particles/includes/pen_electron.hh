
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

 
#ifndef __PEN_ELECTRON_
#define __PEN_ELECTRON_

class pen_betaE;

class pen_eHEC final : public abc_interaction<pen_betaE, pen_context, pen_material>{
private:
public:
  pen_eHEC() : abc_interaction(BETAe_HARD_ELASTIC){}
  double iMeanFreePath(const pen_material& mat,
		       const pen_betaE& particle) const final override;
  
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaE&,
	       double&,
	       pen_rand&) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaE&,
		double&,
		pen_rand&) const final override;
};

class pen_eHIC final : public abc_interaction<pen_betaE, pen_context, pen_material>{
 private:
  pen_particleStack<pen_particleState>& stack;  //b- stack
 public:

 pen_eHIC(pen_particleStack<pen_particleState>& stackIn)
   : abc_interaction(BETAe_HARD_INELASTIC),
    stack(stackIn){}
  
  double iMeanFreePath(const pen_material& mat,
		       const pen_betaE& particle) const final override;
  
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaE&,
	       double&,
	       pen_rand&) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaE&,
		double&,
		pen_rand&) const final override;
};

class pen_eHBE final : public abc_interaction<pen_betaE, pen_context, pen_material>{
 private:
  pen_particleStack<pen_state_gPol>& stack; //Gamma stack

  // **** Bremsstrahlung angular distribution parameters
  static const int NET = 7;
  double BET[NET];      
 public:

  //  ****  Bremsstrahlung emission.
  static const double WB[constants::NBW];
    
  pen_eHBE(pen_particleStack<pen_state_gPol>& stackIn);
  
  double iMeanFreePath(const pen_material& mat,
		       const pen_betaE& particle) const final override;
  
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaE&,
	       double&,
	       pen_rand&) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaE&,
		double&,
		pen_rand&) const final override;
};

class pen_eSII final : public abc_interaction<pen_betaE, pen_context, pen_material>{
 private:
  pen_particleStack<pen_particleState>& stackE; //b- stack
  pen_particleStack<pen_state_gPol>&    stackG; //Gamma stack
 public:

 pen_eSII(pen_particleStack<pen_particleState>& stackEin,
	  pen_particleStack<pen_state_gPol>& stackGin)
   : abc_interaction(BETAe_HARD_INNER_SHELL),
    stackE(stackEin),
    stackG(stackGin){}
  
  double iMeanFreePath(const pen_material& mat,
		       const pen_betaE& particle) const final override;
  
  int interact(const pen_context&,
	       const pen_material&,
	       pen_betaE&,
	       double&,
	       pen_rand&) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_betaE&,
		double&,
		pen_rand&) const final override;
};

class pen_betaE final : public abc_particle<pen_particleState, pen_context, pen_material>{

  friend class pen_eHEC;
  friend class pen_eHIC;
  friend class pen_eHBE;
  friend class pen_eSII;
  
private:
  pen_eHEC HelasticCol;
  pen_eHIC HinelastCol;
  pen_eHBE Hbremsstrahlung;
  pen_eSII HinnerShell;

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

  //  ----  Energy deposited in the last event (analogue simulation).
  double DEA;
  
public:

  static const double REV;  // Electron rest energy (eV)

      
  pen_betaE(const pen_context& contextIn,
	    pen_particleStack<pen_particleState>& stackEin,
	    pen_particleStack<pen_state_gPol>& stackGin);

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

    
  inline void START() final override{
    MHINGE = 0;
    ELAST1 = state.E+1.0E30;
    ELAST2 = ELAST1;
  }
};

inline double pen_eHEC::iMeanFreePath(const pen_material& mat,
				      const pen_betaE& particle) const{
  //Compute inverse mean free pathz
  return exp(mat.SEHEL[particle.KE]+mat.DSEHEL[particle.KE]*particle.XEK);
}

inline double pen_eHIC::iMeanFreePath(const pen_material& mat,
				      const pen_betaE& particle) const{
  //Compute inverse mean free path
  return exp(mat.SEHIN[particle.KE]+mat.DSEHIN[particle.KE]*particle.XEK);
}

inline double pen_eHBE::iMeanFreePath(const pen_material& mat,
				      const pen_betaE& particle) const{
  //Compute inverse mean free path
  return exp(mat.SEHBR[particle.KE]+mat.DSEHBR[particle.KE]*particle.XEK);
}

inline double pen_eSII::iMeanFreePath(const pen_material& mat,
				      const pen_betaE& particle) const{
  //Compute inverse mean free path
  return exp(mat.SEISI[particle.KE]+mat.DSEISI[particle.KE]*particle.XEK);
}

inline void pen_betaE::IMFP(const pen_material& mat){
  P[BETAe_HARD_ELASTIC]        = HelasticCol.iMeanFreePath(mat,(*this));
  P[BETAe_HARD_INELASTIC]      = HinelastCol.iMeanFreePath(mat,(*this));
  P[BETAe_HARD_BREMSSTRAHLUNG] = Hbremsstrahlung.iMeanFreePath(mat,(*this));
  P[BETAe_HARD_INNER_SHELL]    = HinnerShell.iMeanFreePath(mat,(*this));
}

inline void pen_betaE::hingeUpdate(const pen_material& mat){
  if(mat.W1E[KE+1] > -78.3)
    {
      W1 = exp(mat.W1E[KE]+(mat.W1E[KE+1]-mat.W1E[KE])*XEK);
      W2 = exp(mat.W2E[KE]+(mat.W2E[KE+1]-mat.W2E[KE])*XEK);
    }
  else
    {
      W1 = 0.0;
      W2 = 0.0;
    }
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
}

void EBRa(const pen_material& mat,
	  const double E,
	  const unsigned KE,
	  const double XEK,	
	  double &W,
	  const double* WB,
	  pen_rand& penRand);

void EELd(const pen_material& mat,
	  const unsigned KE,
	  const double XEL,
	  double &RNDC,
	  double &RMU,
	  const double* DLEMP,
	  const double DLFC,
	  pen_rand& penRand);

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
	  pen_rand& penRand);

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
	  pen_rand& penRand);

void EELa(double A,
	  double B,
	  double RNDC,
	  double &RMU,
	  pen_rand& penRand);

void EBRaA(const pen_material& mat,
	   double &E,
	   double &DE,
	   double &CDT,
	   const double* BET,
	   pen_rand& penRand);

void ESIb(const pen_material& mat,
	  const int KE,
	  const double XEL,
	  const double* DLEMP,
	  const double DLFC,
	  int &IZZ,
	  int &ISH,
	  pen_rand& penRand);

#endif
