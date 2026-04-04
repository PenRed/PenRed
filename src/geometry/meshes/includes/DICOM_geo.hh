
//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2026 Vicent Giménez Alventosa
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
#include <algorithm>
#include <set>
#include <unordered_set>
#include <iterator>
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

struct voxelInfo{
  unsigned long index;
  unsigned mat;
  double dens;

  constexpr voxelInfo(const unsigned long i,
		      const unsigned matIn,
		      const double densIn) : index(i),
					     mat(matIn),
					     dens(densIn){}
};

struct segmentConstraints{

  unsigned mat;
  double minVolume; //cm**3
  double maxVolume; //cm**3
  unsigned maxClusters;

  constexpr segmentConstraints(const unsigned matIn,
			       const double minV,
			       const double maxV,
			       const unsigned maxC) : mat(matIn),
						      minVolume(minV),
						      maxVolume(maxV),
						      maxClusters(maxC){}
  inline std::string stringify() const {
    std::string out;
    out = "Material: " + std::to_string(mat);
    double auxMin = std::max(minVolume,0.0);
    out += ", Volume: [" + std::to_string(auxMin);
    if(maxVolume <= 0.0){
      out += ", Inf]";
    }
    else{
      out += ", " + std::to_string(maxVolume) + "]";
    }
    if(maxClusters > 0)
      out += ", Clusters: " + std::to_string(maxClusters);

    return out;
  }
};

class pen_dicomGeo : public pen_voxelGeo{
  DECLARE_GEOMETRY(pen_dicomGeo)
  
  private:

  //Store dicom information
  pen_dicom dicom;

  //Store reference material densities
  double densities[constants::MAXMAT];
  
  public:

  enum errors{
    SUCCESS = 0,
    MISSING_CONFIG_PARAMETER,
    CORRUPTED_FILE,
    BAD_VALUE,
    UNKNOWN_BODY,
    INTENSITY_RANGE_CONFIG_ERROR,
    DENSITY_RANGE_CONFIG_ERROR,
    DICOM_LOAD_ERROR,
    DICOM_CONTOUR_ERROR,
    LIMIT_EXCEEDED,
    VOXEL_ASSIGN_ERROR,
    NULL_FILENAME,
    UNABLE_TO_CREATE_FILE,
  };

  static constexpr const char* errorMessage(const int val) noexcept {
    switch(val){
    case SUCCESS: return "Success";
    case MISSING_CONFIG_PARAMETER: return "Missing configuration parameter";
    case CORRUPTED_FILE: return "Corrupted file";
    case BAD_VALUE: return "Bad parameter value";
    case UNKNOWN_BODY: return "Unknown body";
    case INTENSITY_RANGE_CONFIG_ERROR: return "Error configuring intensity range";
    case DENSITY_RANGE_CONFIG_ERROR: return "Error configuring density range";
    case DICOM_LOAD_ERROR: return "Error loading DICOMs";
    case DICOM_CONTOUR_ERROR: return "Error on DICOM contours";
    case LIMIT_EXCEEDED: return "Limit exceeded";
    case VOXEL_ASSIGN_ERROR: return "Voxel assign error";
    case NULL_FILENAME: return "Null filename provided";
    case UNABLE_TO_CREATE_FILE: return "Unable to create file";
    default: return "Unknown error";
    }
  }

  penred::errors::Error specificConfigure(const pen_parserSection& config,
					  const unsigned verbose) override;

  virtual penred::errors::Error printImage(const char* filename) const override;

  penred::errors::Error printContourMasks(const char* filename) const;

  penred::errors::Error printContourMaskSummary(const char* filename) const;
  
  inline const pen_dicom& readDicom() const {return dicom;}

  inline void getOffset(double* offset) const override{
    offset[0] = dicom.getOriginX();
    offset[1] = dicom.getOriginY();
    offset[2] = dicom.getOriginZ();
  }

  static penred::errors::Error readIntensityRanges(const pen_parserSection& config,
						   std::vector<intensityRange>& data,
						   const unsigned verbose);

  static penred::errors::Error readDensityRanges(const pen_parserSection& config,
						 std::vector<densityRange>& data,
						 const unsigned verbose);

  static penred::errors::Error readSegmentConstraints(const pen_parserSection& config,
						      std::vector<segmentConstraints>& data,
						      const unsigned verbose);  
};

//Define pen_dicomGeo error message function
template<>
constexpr const char* penred::errors::errorMessage<pen_dicomGeo>(const int val) noexcept {
  return pen_dicomGeo::errorMessage(val);
}

#endif
#endif
