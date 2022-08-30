
//
//
//    Copyright (C) 2019-2022 Universitat de València - UV
//    Copyright (C) 2019-2022 Universitat Politècnica de València - UPV
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

#include "randomState.hh" 

void random_specificSampler::sample(pen_particleState& state,
				    pen_KPAR& genKpar,
				    unsigned long long& dhist,
				    pen_rand& random){
    
  state.X = random.rand();
  state.Y = random.rand();
  state.Z = random.rand();
    
  state.U = random.rand()+0.1;
  state.V = random.rand()+0.1;
  state.W = random.rand()+0.1;

  double modul = sqrt(state.U*state.U+
		      state.U*state.V+
		      state.W*state.W);

  state.U /= modul;
  state.V /= modul;
  state.W /= modul;
    
  state.E = random.rand()*1.0e4+1.0e3;
  state.PAGE = random.rand();
    
  genKpar = PEN_ELECTRON;
  dhist = 1;
}

int random_specificSampler::configure(double& Emax,
				      const pen_parserSection& /*config*/,
				      const unsigned /*verbose*/){
  Emax = 1.0e4+2.0e3;
  return 0;
}

REGISTER_SPECIFIC_SAMPLER(random_specificSampler,pen_particleState, RANDOM)
