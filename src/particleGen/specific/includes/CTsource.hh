
//
//
//    Copyright (C) 2020-2021 Universitat de València - UV
//    Copyright (C) 2020-2021 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#ifndef __CT_SPECIFIC_SAMPLER__
#define __CT_SPECIFIC_SAMPLER__

#include "PSFsource.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace CTsource{
  struct trans{
    double x,y,z;

    inline void translate(pen_particleState& state){
      state.X += x;
      state.Y += y;
      state.Z += z;
    }
  };
  struct rotz{
    double c,s;
    inline void rotate(pen_particleState& state) const {
      const double xcopy = state.X;
      state.X =  c*xcopy - s*state.Y;
      state.Y =  s*xcopy + c*state.Y;

      const double ucopy = state.U;
      state.U =  c*ucopy - s*state.V;
      state.V =  s*ucopy + c*state.V;      
    }
  };
}

class ct_specificSampler : public abc_specificSampler<pen_particleState>{
  DECLARE_SAMPLER(ct_specificSampler)

  psf_specificSampler psf;
  double  tmin, dt;
  unsigned long nphi;
  double polarRot, azimRot;
  double azim, polar, r;

  //CTsource::trans psf2originRing;
  CTsource::trans part2psfOrigin;
  CTsource::trans origin2CT;
  //CTsource::rotz particleRotPolar;
  double particleRot[9];
  
  std::vector<CTsource::rotz> rotations;
  std::vector<std::pair<pen_particleState,pen_KPAR>> histStates;
  std::pair<pen_particleState,pen_KPAR> nextHistFirstState;
  unsigned long savedDhist;
  unsigned long actualCTpos;
  size_t actualState;
  
public:
    
    
ct_specificSampler() :abc_specificSampler<pen_particleState>(USE_NONE){}    
    
void skip(const unsigned long long dhists);
    
int configure(double& Emax,
	      const abc_spatialSampler* /*pSpatial*/,
	      const abc_directionSampler* /*pDirection*/,
	      const abc_energySampler* /*pEnergy*/,
	      const abc_timeSampler* /*pTime*/,
	      const pen_parserSection& config,
	      const unsigned verbose);  

void sample(pen_particleState& state,
	      pen_KPAR& genKpar,
	      unsigned long long& dhist,
	      pen_rand& random);  
};

#endif
