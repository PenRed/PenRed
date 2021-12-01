
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

#ifndef __PEN_SPATIAL_DOSE_DISTRIB_TALLY__
#define __PEN_SPATIAL_DOSE_DISTRIB_TALLY__


#include "pen_constants.hh"

class pen_SpatialDoseDistrib: public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_SpatialDoseDistrib,pen_particleState)
  
  private:
  int  nx,ny,nz;
  long int nxy,nbin;
  double xmin,ymin,zmin,dx,dy,dz,idx,idy,idz;
  double voxVol;
  unsigned long long* nlast;
  double* edptmp;
  double* edep;
  double* edep2;
  double* ivoxMass;

  unsigned long long* nlastdepth;
  double* edepthtmp;
  double* edepth;
  double* edepth2;
  
  double imatDens[constants::MAXMAT];
public:
    
      pen_SpatialDoseDistrib() : pen_genericTally( USE_LOCALEDEP |
						   USE_BEGINPART |
						   USE_SAMPLEDPART |
						   USE_STEP |
						   USE_ENDSIM |
						   USE_MOVE2GEO)
                         

  {
      nx = ny = nz = nxy = nbin = 0;
      dx = dy = dz = 0.0;
      idx = idy = idz = 1.0e35;
      xmin = ymin = zmin = 0.0;
      
      nlast = nullptr;
      edptmp = nullptr;
      edep = nullptr;
      edep2 = nullptr;
      ivoxMass = nullptr;

      nlastdepth = nullptr;
      edepthtmp = nullptr;
      edepth = nullptr;
      edepth2 = nullptr;      
}
  
  void updateEdepCounters(const double dE,
                          const unsigned long long nhist,
                          const double X,
			  const double Y,
			  const double Z,
			  const double WGHT,
			  const unsigned MAT);
    
  void clear(void);
  
  void tally_localEdep(const unsigned long long nhist,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state,
		       const double dE);
  
  void tally_beginPart(const unsigned long long nhist,
		       const unsigned /*kdet*/,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state);
  
  void tally_sampledPart(const unsigned long long nhist,
			 const unsigned long long /*dhist*/,
			 const unsigned /*kdet*/,
			 const pen_KPAR /*kpar*/,
			 const pen_particleState& state);

  void tally_step(const unsigned long long nhist,
		  const pen_KPAR /*kpar*/,
		  const pen_particleState& state,
		  const tally_StepData& stepData);
  
  void tally_move2geo(const unsigned long long nhist,
		      const unsigned /*kdet*/,
		      const pen_KPAR /*kpar*/,
		      const pen_particleState& state,
		      const double /*dsef*/,
		      const double /*dstot*/);
  
  
  void tally_endSim(const unsigned long long /*nhist*/);
    
  int configure(const wrapper_geometry& geometry,
		const abc_material* const materials[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  void saveData(const unsigned long long nhist) const;
  void flush(void);
  int sumTally(const pen_SpatialDoseDistrib& tally);
  
  ~pen_SpatialDoseDistrib(){clear();}
  
};


#endif
