
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

 
#ifndef __PEN_FLN_TRACK_LENGTH_TALLY__
#define __PEN_FLN_TRACK_LENGTH_TALLY__

#include "pen_constants.hh"

class pen_ImpactDetector: public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_ImpactDetector,pen_particleState)
    
  private:
  
  bool isLinScaleFlu, isLinScaleSpc,
    isLinScaleAge, isLinScaleDep,
    inside, spc, ageActive, fln,
    eDepActive;
  int  nbin, nbinAge;
  unsigned idet;
  static const int nbinmax=32000; 

  double flutmp[constants::nParTypes][nbinmax],
         flu[constants::nParTypes][nbinmax],
         flu2[constants::nParTypes][nbinmax],
         flutmptotal[nbinmax],
         flutotal[nbinmax],
         flu2total[nbinmax];
    
  double spectrumTmp[constants::nParTypes][nbinmax],
         spectrum[constants::nParTypes][nbinmax],
         spectrum2[constants::nParTypes][nbinmax],
         spectrumTotalTmp[nbinmax],
         spectrumTotal[nbinmax],
         spectrum2Total[nbinmax];
         
  double edepCounterTmp,
         edepCounter, 
         edep2Counter;
        
  long unsigned int detected[nbinmax];
         
  double agetmp[nbinmax], 
         age[nbinmax],
         age2[nbinmax];
         

  double egrid[nbinmax],ebingrd[nbinmax];
  double ebin,iebin,eratio,unclimit,
         emin,emax,ageMin,ageMax, ageBinW, iageBinW;
  double configEmin, configAgeMin;
  double configEmax;
    
    
public:
    
  pen_ImpactDetector() : pen_genericTally( USE_STEP |
					   USE_LOCALEDEP |
					   USE_INTERFCROSS |
					   USE_MOVE2GEO |
					   USE_BEGINPART |
					   USE_ENDHIST),
			 isLinScaleFlu(true), isLinScaleSpc(true),
			 isLinScaleAge(true), isLinScaleDep(true),
			 inside(false), spc(false), ageActive(false),
			 fln(false), eDepActive(false),			 
			 nbin(0),
			 nbinAge(0)			 
    
    
  {}
    
  //void trackl(double ds, double energy, double de, const pen_KPAR kpar);
  void discreteTrackL(const pen_KPAR kpar,
		      const double ds,
		      const double energy);
    
  void continuousTrackL(const pen_KPAR kpar,
                        const pen_particleState& state,
                        const double ds,
                        const double de,
                        const double energy);
  
  void countSpectrum(const pen_KPAR kpar,
                     const pen_particleState& state);
  
  void energyDep(const double eDep);
  
  void ageSpectrum(const pen_particleState& state);

    
  void tally_step(const double nhist,
		  const pen_KPAR kpar,
		  const pen_particleState& state,
		  const tally_StepData& stepData);
  
  void tally_localEdep(const double /*nhist*/,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state,
		       const double dE);
    
  void tally_interfCross(const double nhist,
			 const unsigned kdet,
			 const pen_KPAR kpar,
			 const  pen_particleState& state);
			  
  void tally_move2geo(const double nhist,
		      const unsigned kdet,
		      const pen_KPAR kpar,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);
  
  
  void tally_beginPart(const double /*nhist*/,
		       const unsigned kdet,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state);
  
  
  void tally_endHist(const double /*nhist*/);
    
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  void flush();
  void saveData(const double nhist) const;
  int sumTally(const pen_ImpactDetector& tally);
};

#endif
