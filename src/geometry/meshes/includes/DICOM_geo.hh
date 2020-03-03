
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
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

 
#ifndef __PENRED_DICOM_GEOMETRY__
#define __PENRED_DICOM_GEOMETRY__

#ifdef _PEN_USE_DICOM_

#include "pen_dicom.hh"
#include "voxel_geo.hh"

class pen_dicomGeo : public pen_voxelGeo{
  DECLARE_GEOMETRY(pen_dicomGeo)
  
  private:

  //Store dicom information
  pen_dicom dicom;

  //Store reference material densities
  double densities[constants::MAXMAT];
  
  public:

  int configure(const pen_parserSection& config, const unsigned verbose);

  virtual int printImage(const char* filename) const;
  
  inline const pen_dicom& readDicom() const {return dicom;}
  
};

#endif
#endif
