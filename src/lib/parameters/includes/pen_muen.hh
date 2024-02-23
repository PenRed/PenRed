 
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

#ifndef __PEN_MUEN__
#define __PEN_MUEN__

#include <thread>
#include <atomic>
#include "PenRed.hh"
#include "pen_geometries.hh"
#include "pen_samplers.hh"

namespace pen_muen{

  struct muData{
    double muRho, muenRho;
    double E0, f, fg;
    double err;

    inline double muen() const {return fg*muRho/E0;}

    static inline std::string header(unsigned shift = 0) {
      char aux[400];
      sprintf(aux,"#%*s    Mean E       mu/rho      mu_en/rho     "
	      "(1-g)*f        f            1-g      uncert.      \n"
	      "#%*s     (eV)       (cm^2/g)     (cm^2/g)%45.45s(%%)       \n",
	      shift, " ", shift, " "," ");
      return std::string(aux);
    }
    inline std::string stringify() const {
      char aux[150];
      sprintf(aux,"%.5E  %.5E  %.5E  %.5E  %.5E  %.5E  %.1E",
	      E0,muRho,muen(),fg/E0,f/E0,fg/f,err);
      return std::string(aux);
    }
  };

  struct muenSimTally{
    double EDPT, EDPT2;
    double ETRT, ETRT2;
    double RANGET, RANGET2;
    double ET, ET2;
    unsigned long long nhist;
  };

  int calculate(const double Emin,
		const double Emax,
		const unsigned nBins,
		const double tolerance,
		const double simTime,
		const char* matFilename,
		std::vector<double>& EData,
		std::vector<double>& muenData);

  int calculate(const char** energySpectrums,
		const unsigned nSpectrums,
		const double tolerance,
		const double simTime,
		const char* matFilename,
		std::vector<pen_muen::muData>& muenData,
		unsigned verbose = 1);
  
  double simulate(const pen_context& context,
		  const double E0,
		  const double simTime,
		  const double tolerance,
		  int& seed1, int& seed2,
		  const unsigned ithread,
		  const unsigned long long minHists = 100000);
  
  double simulate(const pen_context& context,
		  const fileSpectrum_energySampling& spectrum,
		  const double simTime,
		  const double tolerance,
		  int& seed1, int& seed2,
		  pen_muen::muenSimTally &tally,
		  const unsigned long long minHists = 100000);
  
};

#endif
