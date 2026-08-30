
//
//
//    Copyright (C) 2026 Vicent Giménez Alventosa
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
//    
//


#ifndef __MEASURE_ENERGY_SPECTRUM_SAMPLING__
#define __MEASURE_ENERGY_SPECTRUM_SAMPLING__

class measure_energySampling : public abc_energySampler{
  DECLARE_SAMPLER(measure_energySampling)
private:

  penred::sampling::aliasing<1> sampler;
  
public:
  void energySampling(double& energy, pen_rand& random) const;
  int configure(double& Emax,
                const pen_parserSection& config,
                const unsigned verbose) final override;
};

#endif
