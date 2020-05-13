
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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



#ifndef __RANDOM_SPECIFIC_SAMPLER__
#define __RANDOM_SPECIFIC_SAMPLER__

class random_specificSampler : public abc_specificSampler<pen_particleState>{
  DECLARE_SAMPLER(random_specificSampler)
  private:
  public:

  random_specificSampler() : abc_specificSampler<pen_particleState>(USE_NONE)
  {}
  
  void sample(pen_particleState& state,
	      pen_KPAR& genKpar,
	      unsigned long long& dhist,
	      pen_rand& random);
  
  int configure(double& Emax,
		const abc_spatialSampler* /*pSpatial*/,
		const abc_directionSampler* /*pDirection*/,
		const abc_energySampler* /*pEnergy*/,
		const abc_timeSampler* /*pTime*/,
		const pen_parserSection& /*config*/,
		const unsigned /*verbose*/);
};

#endif
