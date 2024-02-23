
//
//
//    Copyright (C) 2021 Universitat de València - UV
//    Copyright (C) 2021 Universitat Politècnica de València - UPV
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


#ifndef __FILE_ENERGY_SPECTRUM_SAMPLING__
#define __FILE_ENERGY_SPECTRUM_SAMPLING__

#include "pen_auxiliar.hh"

class fileSpectrum_energySampling : public abc_energySampler{
  DECLARE_SAMPLER(fileSpectrum_energySampling)
private:

  std::vector<double> energies;
  std::vector<double> dE;
  std::vector<double> cummulative;
  unsigned nEBins;
  
  
public:
  void energySampling(double& energy, pen_rand& random) const;
  int configure(double& Emax, const pen_parserSection& config, const unsigned verbose);
  inline double minE(){return energies.front();}
  inline const std::vector<double>& readEnergy() const {return energies;}
  inline const std::vector<double>& readDE() const {return dE;}
  inline const std::vector<double>& readCummulative() const {return cummulative;}
  inline unsigned nBins() const {return nEBins;}
};

#endif
