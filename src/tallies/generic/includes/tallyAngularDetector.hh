
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

 
#ifndef __PEN_ANGULAR_DET_TALLY__
#define __PEN_ANGULAR_DET_TALLY__

#include "pen_constants.hh" 


class pen_AngularDet : public
pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_AngularDet,pen_particleState)
  
private:
    bool isLinScale;
    static const int nbinmax = 32000;
    double phi1, phi2,
           theta1, theta2;
    double emin, emax;
    int nBinsE;
    unsigned idet;
    
    double ebin,iebin;
    double angDetSpcTmp[constants::nParTypes][nbinmax],
           angDetSpc[constants::nParTypes][nbinmax], 
           angDetSpc2[constants::nParTypes][nbinmax];
           
           
public:
    
  pen_AngularDet() : pen_genericTally( USE_ENDHIST |
				       USE_MOVE2GEO |
				       USE_MATCHANGE),
		     nBinsE(0),
		     ebin(0.0),
		     iebin(1.0e35)		     
    {}
    
      
  void scapedParticle(const pen_KPAR kpar,
				 const pen_particleState& state);
  
  void tally_move2geo(const unsigned long long /*nhist*/,
		      const unsigned /*kdet*/,
		      const pen_KPAR kpar,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);
  
  void tally_matChange(const unsigned long long /*nhist*/,
		       const pen_KPAR kpar,
		       const pen_particleState& state,
		       const unsigned /*prevMat*/);
  
  void tally_endHist(const unsigned long long /*nhist*/);
  
  
  int configure(const wrapper_geometry& /*geometry*/,
		const abc_material* const /*materials*/[constants::MAXMAT],     
		const pen_parserSection& config, const unsigned verbose);
  void flush();
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_AngularDet& tally);
  
};


#endif
