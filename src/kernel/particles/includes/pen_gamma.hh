
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

 
#ifndef __PEN_GAMMA_
#define __PEN_GAMMA_

class pen_gamma;

struct GCO00{
  double EE;
  int IOSC;
  const pen_material* pmat;
};

class pen_gRCS final : public abc_interaction<pen_gamma, pen_context, pen_material>{
 private:
  
 public:
 pen_gRCS() : abc_interaction(GAMMA_RAYLEIGH){}
  double iMeanFreePath(const pen_material&, const pen_gamma&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_gamma&,
	       double&,
	       pen_rand&) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_gamma&,
		double&,
		pen_rand&) const final override;
};

class pen_gCIS final : public abc_interaction<pen_gamma, pen_context, pen_material>{
 private:
  pen_particleStack<pen_particleState>& stackE; //b- stack
  pen_particleStack<pen_state_gPol>& stackG; //Gamma stack  
 public:

 pen_gCIS(pen_particleStack<pen_particleState>& stackInE,
	  pen_particleStack<pen_state_gPol>& stackInG)
   : abc_interaction(GAMMA_COMPTON),
    stackE(stackInE),
    stackG(stackInG){}
  
  double iMeanFreePath(const pen_material&, const pen_gamma&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_gamma&,
	       double&,
	       pen_rand&) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_gamma&,
		double&,
		pen_rand&) const final override;
};

class pen_gPHE final : public abc_interaction<pen_gamma, pen_context, pen_material>{
 private:
  pen_particleStack<pen_particleState>& stackE; //b- stack
  pen_particleStack<pen_state_gPol>& stackG; //Gamma stack  
 public:
  
 pen_gPHE(pen_particleStack<pen_particleState>& stackInE,
	  pen_particleStack<pen_state_gPol>& stackInG)
   : abc_interaction(GAMMA_PHOTOELECTRIC),
    stackE(stackInE),
    stackG(stackInG){}
  
  double iMeanFreePath(const pen_material&, const pen_gamma&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_gamma&,
	       double&,
	       pen_rand&) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_gamma&,
		double&,
		pen_rand&) const final override;
};

class pen_gBBP final : public abc_interaction<pen_gamma, pen_context, pen_material>{
 private:
  pen_particleStack<pen_particleState>& stackE; //b- stack
  pen_particleStack<pen_particleState>& stackP; //b+ stack
  pen_particleStack<pen_state_gPol>& stackG; //gamma stack
 public:
  
 pen_gBBP(pen_particleStack<pen_particleState>& stackEin,
	  pen_particleStack<pen_particleState>& stackPin,
	  pen_particleStack<pen_state_gPol>& stackGin)
   : abc_interaction(GAMMA_PAIR_PRODUCTION),
    stackE(stackEin),
    stackP(stackPin),
    stackG(stackGin){}
  
  double iMeanFreePath(const pen_material&, const pen_gamma&) const final override;
  int interact(const pen_context&,
	       const pen_material&,
	       pen_gamma&,
	       double&,
	       pen_rand&) const final override;
  int interactF(const pen_context&,
		const pen_material&,
		pen_gamma&,
		double&,
		pen_rand&) const final override;
};


class pen_gamma final : public abc_particle<pen_state_gPol, pen_context, pen_material>{

  friend class pen_gRCS;
  friend class pen_gCIS;
  friend class pen_gPHE;
  friend class pen_gBBP;
  
private:
  
  pen_gRCS rayleigh;
  pen_gCIS compton;
  pen_gPHE photoelectric;
  pen_gBBP pairProduct;

  //  ----  Energy deposited in the last event (analogue simulation).
  double DEA;
  
public:
    
  pen_gamma(const pen_context& contextIn,
	    pen_particleStack<pen_particleState>& stackEin,
	    pen_particleStack<pen_particleState>& stackPin,
	    pen_particleStack<pen_state_gPol>& stackGin);
  
  void JUMP(double &DS, pen_rand& penRand, const double /*DSMAX*/) final override;
  void JUMPF(double &DS, pen_rand& penRand, const double /*DSMAX*/) final override;
  void KNOCK(double &DE, int &ICOL, pen_rand& penRand) final override;
  void KNOCKF(double &DE, int &ICOL, pen_rand& penRand) final override;

  void softEloss(double& /*X*/,
		 double& /*Y*/,
		 double& /*Z*/,
		 double& DE,
		 pen_rand& /*penRand*/) final override;  
  
  void dpage() final override;
  inline void page0() final override{}
  void IMFP(const pen_material& mat);
  inline void START() final override{ELAST1 = state.E+1.0E30;}
};

inline double pen_gRCS::iMeanFreePath(const pen_material& mat,
				      const pen_gamma& particle) const{
  //Compute inverse mean free path
  return mat.SGRA[particle.KE];
}

inline double pen_gCIS::iMeanFreePath(const pen_material& mat,
				      const pen_gamma& particle) const{
  //Compute inverse mean free path
  return exp(mat.SGCO[particle.KE]+mat.DSGCO[particle.KE]*particle.XEK);
}

inline double pen_gPHE::iMeanFreePath(const pen_material& mat,
				      const pen_gamma& particle) const{
  //Compute inverse mean free path
  return mat.SGPH[particle.KE];
}

inline double pen_gBBP::iMeanFreePath(const pen_material& mat,
				      const pen_gamma& particle) const{
  if(particle.state.E < 1.023e6)
    return 0.0;
  
  //Compute inverse mean free path
  return exp(mat.SGPP[particle.KE]+mat.DSGPP[particle.KE]*particle.XEK);
}

inline void pen_gamma::IMFP(const pen_material& mat){
  P[GAMMA_RAYLEIGH]        = rayleigh.iMeanFreePath(mat,(*this));
  P[GAMMA_COMPTON]         = compton.iMeanFreePath(mat,(*this));
  P[GAMMA_PHOTOELECTRIC]   = photoelectric.iMeanFreePath(mat,(*this));
  P[GAMMA_PAIR_PRODUCTION] = pairProduct.iMeanFreePath(mat,(*this));
}

void GPHa(const pen_material& mat,
	  const pen_elementDataBase& elemDB,
	  const double E,
	  const unsigned KE,
	  const double XEL,
	  double &ES,
	  int &IZZ,
	  int &ISH,
	  pen_rand& penRand);

void GPHa0(pen_elementDataBase& elemDB);

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
	  pen_rand& penRand);


void GRAa(const pen_material& mat,
	  const double E,
	  const unsigned KE,
	  const double XEL,
	  double &CDT,
	  int &IEFF,
	  pen_rand& penRand);


void GCOa(const pen_material& mat,
	  double E,
	  double &DE,
	  double &EP,
	  double &CDT,
	  double &ES,
	  double &CDTS,
	  int &IZZ,
	  int &ISH,
	  pen_rand& penRand);

double GCOaD(double CDT,
	     void* arg);

void GCOaT(const pen_material& mat,
	   const double E,
	   double &CS);

void SAUTER(const double ES,
	    double &CDTS,
	    double &DFS,
	    pen_rand& penRand,
	    pen_state_gPol& gPolstate);

void SCHIFF(double &B,
	    double &G1,
	    double &G2);

#endif
