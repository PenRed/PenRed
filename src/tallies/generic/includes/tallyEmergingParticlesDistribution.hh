
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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


#ifndef __PEN_EMERGIN_PART_DISTRIB_TALLY__
#define __PEN_EMERGIN_PART_DISTRIB_TALLY__

#include "pen_constants.hh" 

class pen_EmergingPartDistrib : public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_EmergingPartDistrib,pen_particleState,EMERGING_PART_DISTRIB,
		std::pair<double, penred::tally::Dim<1>>, //Upbound electrons
		std::pair<double, penred::tally::Dim<1>>, //Upbound gamma
		std::pair<double, penred::tally::Dim<1>>, //Upbound positrons
		std::pair<double, penred::tally::Dim<1>>, //Downbound electrons
		std::pair<double, penred::tally::Dim<1>>, //Downbound gamma
		std::pair<double, penred::tally::Dim<1>>, //Downbound positrons
		std::pair<double, penred::tally::Dim<2>>, //Angular electrons
		std::pair<double, penred::tally::Dim<2>>, //Angular gamma
		std::pair<double, penred::tally::Dim<2>>  //Angular positrons
		)
  
private:
  bool energyLogscale, angularLogscale;
  static const int nbinmax=32000;
  int nBinsE, nBinsTheta, nBinsPhi, nAngBins; //Number of bins
  double ebin, tbin, pbin;  //Bins size
  double iebin, itbin, ipbin;  //Bins inverse size
  double emin, emax;
  double configEmin;
  double tmin, tmax;
  const double pmin;
  const double pmax;
  
  double eHistUp[constants::nParTypes][nbinmax],
    eHistTmpUp[constants::nParTypes][nbinmax],
    eHist2Up[constants::nParTypes][nbinmax];
  
  double eHistDown[constants::nParTypes][nbinmax],
    eHistTmpDown[constants::nParTypes][nbinmax],
    eHist2Down[constants::nParTypes][nbinmax];
  
  unsigned long long enlastUp[constants::nParTypes][nbinmax];
  unsigned long long enlastDown[constants::nParTypes][nbinmax];
  
  double* angHist[constants::nParTypes];
  double* angHistTmp[constants::nParTypes];
  double* angHist2[constants::nParTypes];
  unsigned long long* angnlast[constants::nParTypes];
    
    
public:
    
  pen_EmergingPartDistrib();

  penred::measurements::results<double, 1> generateUpbound(const pen_KPAR kpar,
							   const unsigned long long nhist);

  penred::measurements::results<double, 1> generateDownbound(const pen_KPAR kpar,
							     const unsigned long long nhist);

  penred::measurements::results<double, 2> generateAngular(const pen_KPAR kpar,
							   const unsigned long long nhist);  
  
  void scapedParticle(const unsigned long long nhist,
				 const pen_KPAR /*kpar*/,
				 const pen_particleState& /*state*/);
  
  void tally_move2geo(const unsigned long long nhist,
		      const unsigned /*kdet*/,
		      const pen_KPAR /*kpar*/,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);
  
  void tally_matChange(const unsigned long long nhist,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const unsigned /*prevMat*/);
  
  void tally_endSim(const unsigned long long /*nhist*/);
  
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],     
		const pen_parserSection& config, const unsigned verbose);
  void flush();
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_EmergingPartDistrib& tally);

  ~pen_EmergingPartDistrib(){
    for(unsigned int i = 0; i < constants::nParTypes; i++){

      if(angHist[i] != nullptr){
	free(angHist[i]);
	angHist[i] = nullptr;
      }

      if(angHistTmp[i] != nullptr){
	free(angHistTmp[i]);
	angHistTmp[i] = nullptr;
      }
      
      if(angHist2[i] != nullptr){
	free(angHist2[i]);
	angHist2[i] = nullptr;
      }
      
      if(angnlast[i] != nullptr){
	free(angnlast[i]);
	angnlast[i] = nullptr;
      }      
    }
  }
};


#endif
