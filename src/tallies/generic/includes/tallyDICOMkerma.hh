
//
//
//    Copyright (C) 2021-2022 Universitat de València - UV
//    Copyright (C) 2021-2022 Universitat Politècnica de València - UPV
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

#ifndef __PEN_DICOM_KERMA_TRACK_LENGTH_TALLY__
#define __PEN_DICOM_KERMA_TRACK_LENGTH_TALLY__

#ifdef _PEN_USE_DICOM_

#include "tallyKermaTrackLength.hh"

class pen_tallyDICOMkerma: public pen_genericTally<pen_particleState> {
  DECLARE_TALLY(pen_tallyDICOMkerma,pen_particleState)
  
  private: 
  
  long int nx,ny,nz;
  long int nxy,nbin;
  double dx,dy,dz;
  double xmin,ymin,zmin;
  double voxVol;
  
  double* contVol;
  
  const int* contourVox;
  long int ncontours;

  std::vector<std::string> contNames;
  
  //DVH
  double prescribedDose;
  double DVHfactor;
  unsigned long DVHnbins;
  double DVHbinWidth, DVHmaxDose;
  

  pen_tallyKermaTrackLength tallyKerma;

  
public:
    
  pen_tallyDICOMkerma() : pen_genericTally(USE_JUMP | USE_STEP),
			  contourVox(nullptr),
			  ncontours(0)
  {}

  void tally_step(const unsigned long long nhist,
		  const pen_KPAR kpar,
		  const pen_particleState& state,
		  const tally_StepData& stepData);  
    
  void tally_jump(const unsigned long long nhist,
		  const pen_KPAR kpar,
		  const pen_particleState& state,
		  const double ds);
    
  int configure(const wrapper_geometry& geometry,
		const abc_material* const materials[constants::MAXMAT],
		const pen_parserSection& config,
		const unsigned verbose);

  inline int sharedConfig(const pen_tallyDICOMkerma& tally){
    return tallyKerma.sharedConfig(tally.tallyKerma);
  }

  void flush();
  void saveData(const unsigned long long nhist) const;
  int sumTally(const pen_tallyDICOMkerma& tally);
    
};

#endif

#endif
