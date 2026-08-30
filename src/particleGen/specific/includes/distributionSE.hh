
//
//
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//    
//

#ifndef __SE_DISTRIB_SPECIFIC_SAMPLER__
#define __SE_DISTRIB_SPECIFIC_SAMPLER__

class distribSE_specificSampler : public abc_specificSampler<pen_particleState>{
  DECLARE_SPECIFIC_SAMPLER(distribSE_specificSampler, pen_particleState)

  std::shared_ptr<abc_spatialSampler> shiftSampler;
  std::shared_ptr<penred::measurements::results<double, 4>> distribSE;

  vector3D<double> distribOrigin;
  
  size_t ebins;
  double eBinWidth;
  
  std::pair<double, double> elimits;
  std::pair<double, double> xlimits;
  std::pair<double, double> ylimits;
  std::pair<double, double> zlimits;

  double xbins, ybins, zbins;
  double xBinWidth, yBinWidth, zBinWidth;
  
public:
    
  distribSE_specificSampler() : abc_specificSampler<pen_particleState>(USE_SPATIAL){}    
    
  inline void skip(const unsigned long long /*dhists*/){}
    
  int configure(double& Emax,
                const pen_parserSection& config,
                const unsigned nthreads,
                const unsigned verbose);  

  void sample(pen_particleState& state,
              pen_KPAR& genKpar,
              unsigned long long& dhist,
              pen_rand& random);

  int sharedConfig(const distribSE_specificSampler& o);
  void updateGeometry(const wrapper_geometry* geometryIn){
    if(shiftSampler){
      shiftSampler->updateGeometry(geometryIn);
    }
  }

};

//Tally configuration reader
class tallyReader_DistributionSE : public pen_configReader<tallyReader_DistributionSE>{

public:
  
  enum errors{
    SUCCESS = 0,
    UNHANDLED = 1,
  };

  vector3D<double> distribOrigin;
  std::string distribFilename;
  
  int storeElement(const std::string& pathInSection,
                   const pen_parserData& element,
                   const unsigned verbose);
  
  int storeString(const std::string& pathInSection,
                  const std::string& element,
                  const unsigned verbose);  
};

template<>
struct pen_format<tallyReader_DistributionSE>{
  static constexpr const char* format = R"===(

# Tally "DetectionEDep" reader configuration

reader-description "This sampler uses a 4D distribution to assign particles' energy depending on the spatial position. Those distributions can be generated using the 'DETECTION_SPATIAL_DISTRIB' tally, which dimensions are, in order, (E,X,Y,Z). The spatial position is sampled using a generic spatial sampler, and it should fall within the energetic 4D distribution. Finally, the direction is set using an 'origin' point, being the vector from the origin and the sampled particle position the final direction. The particle position can be then shifted using a spatial sampler. Notice that energy and direction are assigned using the initial particle's sampled position."

# Distribution

## Filename
distribution/filename/reader-description "Spatial-energy distribution filename"
distribution/filename/reader-value ""
distribution/filename/reader-required/type "required"

## X
distribution/origin/x/reader-description "Distribution's origin along the X axis, in cm"
distribution/origin/x/reader-value 0.0
distribution/origin/x/reader-required/type "required"

## Y
distribution/origin/y/reader-description "Distribution's origin along the Y axis, in cm"
distribution/origin/y/reader-value 0.0
distribution/origin/y/reader-required/type "required"

## Z
distribution/origin/z/reader-description "Distribution's origin along the Z axis, in cm"
distribution/origin/z/reader-value 0.0
distribution/origin/z/reader-required/type "required"

# Origin Sampler
position/shift/config/${subsection}/reader-description "Section to configure a post-sampling random translation. If used, must be a valid spatial sampler."
position/shift/config/${subsection}/reader-required/type "optional"
)===";
};

#endif
