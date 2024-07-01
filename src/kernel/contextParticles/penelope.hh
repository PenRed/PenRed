
//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
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

 
#ifndef __PENRED_CONTEXT_PARTICLES_PENELOPE_INCLUDE__
#define __PENRED_CONTEXT_PARTICLES_PENELOPE_INCLUDE__

//Implement the 'createContextParticles' partial specialization for this context
namespace penred{
  namespace context{

    //Define particle types for this context
    template<>
    struct particleTypes<pen_context> {
      typedef std::tuple<pen_gamma,pen_betaE,pen_betaP> type;
    };

    //Define particles constructor for this context
    template<> inline particles<pen_context>::particles(const pen_context& c) :
      particleTuple({c, getStack<1>(), getStack<2>(), getStack<0>()},
		    {c, getStack<1>(), getStack<0>()               },
		    {c, getStack<1>(), getStack<0>(), getStack<2>()}) {}
  } // namespace context
} // namespace penred

#endif
