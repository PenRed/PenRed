
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

 
#ifndef __PENRED_DICOM_GEOMETRY__
#define __PENRED_DICOM_GEOMETRY__

#ifdef _PEN_USE_DICOM_

#include <limits>
#include "pen_dicom.hh"
#include "voxel_geo.hh"
#include "pen_image.hh"

struct intensityRange{
  unsigned mat;
  double low, top, dens;

  intensityRange() = default;
  intensityRange(const unsigned inMat,
		 const double inLow,
		 const double inTop,
		 const double inDens) : mat(inMat),
					low(inLow),
					top(inTop),
					dens(inDens) {}

  inline bool inner(const double value) const{
    return value >= low && value < top;
  }
};

struct densityRange{
  unsigned mat;
  double low, top;

  densityRange() = default;
  densityRange(const unsigned inMat,
	       const double inLow,
	       const double inTop) : mat(inMat),
				     low(inLow),
				     top(inTop) {}

  inline bool inner(const double value) const{
    return value >= low && value < top;
  }  
};

struct contourAssign{

  unsigned defaultMat;
  double defaultDens;
  double priority;

  std::vector<intensityRange> intensityRanges;
  std::vector<densityRange> densityRanges;
  
};

class pen_dicomGeo : public pen_voxelGeo{
  DECLARE_GEOMETRY(pen_dicomGeo)
  
  private:

  //Store dicom information
  pen_dicom dicom;

  //Store reference material densities
  double densities[constants::MAXMAT];
  
  public:

  int configure(const pen_parserSection& config, const unsigned verbose) override;

  virtual int printImage(const char* filename) const override;

  int printContourMasks(const char* filename) const;

  int printContourMaskSummary(const char* filename) const;
  
  inline const pen_dicom& readDicom() const {return dicom;}

  inline void getOffset(double* offset) const override{
    offset[0] = dicom.getOriginX();
    offset[1] = dicom.getOriginY();
    offset[2] = dicom.getOriginZ();
  }
};

int readIntensityRanges(const pen_parserSection& config,
			std::vector<intensityRange>& data,
			const unsigned verbose);

int readDensityRanges(const pen_parserSection& config,
		      std::vector<densityRange>& data,
		      const unsigned verbose);

#endif
#endif
