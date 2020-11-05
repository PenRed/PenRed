
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

#ifndef __PEN_DICOM_DOSE_DISTRIB_TALLY__
#define __PEN_DICOM_DOSE_DISTRIB_TALLY__

#ifdef _PEN_USE_DICOM_

class pen_DICOMDoseDistrib: public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_DICOMDoseDistrib,pen_particleState)
  
  private:
  long int  nx,ny,nz;
  long int nxy,nbin;
  double xmin,ymin,zmin,dx,dy,dz,idx,idy,idz;
  double voxVol;

  //Spatial distribution
  unsigned long long* nlast;
  double* edptmp;
  double* edep;
  double* edep2;

  double* ivoxMass;

  const int* contourVox;
  
  //Contour dose tallies
  long int ncontours;
  double* contEdptmp;
  double* contEdep;
  double* contEdep2;

  double* contMass;
  double* contVol;
  
  double matDens[constants::MAXMAT];

  std::vector<std::string> contNames;
public:
    
      pen_DICOMDoseDistrib() : pen_genericTally( USE_LOCALEDEP |
						   USE_BEGINPART |
						   USE_BEGINHIST |
						   USE_STEP |
						   USE_ENDSIM |
						   USE_MOVE2GEO |
						   USE_ENDHIST )
                         

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

      contEdptmp = nullptr;
      contEdep = nullptr;
      contEdep2 = nullptr;
 
      contMass = nullptr;
      contVol  = nullptr;

      contourVox = nullptr;
      ncontours = 0;
      
}
  
  void updateEdepCounters(const double dE,
                          const unsigned long long nhist,
                          const double X,
			  const double Y,
			  const double Z,
			  const double WGHT);
    
  void clear(void);
  
  void tally_localEdep(const unsigned long long nhist,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state,
		       const double dE);
  
  void tally_beginPart(const unsigned long long nhist,
		       const unsigned /*kdet*/,
		       const pen_KPAR /*kpar*/,
		       const pen_particleState& state);
  
  void tally_beginHist(const unsigned long long nhist,
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
  
  
  void tally_endHist(const unsigned long long /*nhist*/);
  
  void tally_endSim(const unsigned long long /*nhist*/);
  
  int configure(const wrapper_geometry& geometry,
		const abc_material* const materials[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);
  void saveData(const unsigned long long nhist) const;
  void flush(void);
  int sumTally(const pen_DICOMDoseDistrib& tally);
  
  ~pen_DICOMDoseDistrib(){clear();}
  
};

#endif

#endif
