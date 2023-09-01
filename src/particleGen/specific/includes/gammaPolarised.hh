
//
//
//    Copyright (C) 2019-2023 Universitat de València - UV
//    Copyright (C) 2019-2023 Universitat Politècnica de València - UPV
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



#ifndef __GAMMA_POLARISED_SPECIFIC_SAMPLER__
#define __GAMMA_POLARISED_SPECIFIC_SAMPLER__

class gammaPolarised_specificSampler : public abc_specificSampler<pen_state_gPol>{
  DECLARE_SAMPLER(gammaPolarised_specificSampler)
  private:

  double SP10,SP20,SP30;
  int IPOL0;
  
  public:

  gammaPolarised_specificSampler() : abc_specificSampler<pen_state_gPol>(USE_GENERIC),
				     SP10(0.0),
				     SP20(0.0),
				     SP30(0.0),
				     IPOL0(0)
  {}
  
  void sample(pen_state_gPol& state,
	      pen_KPAR& /*genKpar*/,
	      unsigned long long& /*dhist*/,
	      pen_rand& /*random*/);
  
  int configure(double& /*Emax*/,
		const pen_parserSection& config,
		const unsigned /*nthreads*/,
		const unsigned verbose);
};

#endif
