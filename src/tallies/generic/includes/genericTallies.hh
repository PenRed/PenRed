
//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
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


#ifndef __PEN_GENERIC_TALLIES__ 
#define __PEN_GENERIC_TALLIES__ 

#include "tallyEnergyDepositionMat.hh"
#include "tallyEnergyDepositionBody.hh"
#include "tallySphericalDoseDistrib.hh"
#include "tallyTracking.hh"
#include "tallyEmergingParticlesDistribution.hh"
#include "tallyCylindricalDoseDistrib.hh"
#include "tallyImpactDetector.hh"
#include "tallySecondaryGen.hh"
#include "tallySpatialDoseDistrib.hh"
#include "tallyAngularDetector.hh"
#include "tallyPhaseSpaceFile.hh"
#include "tallyKermaTrackLength.hh"
#include "tallyDICOMDoseDistrib.hh"
#include "tallyCTsinogram.hh"
#include "tallyDICOMkerma.hh"
#include "tallyPSS.hh"
#include "tallyDetectionSpatialDistrib.hh"
#include "tallyEmergingSphericalDistribution.hh"
#include "tallySingles.hh"

namespace penred{
  namespace tally{

#ifdef _PEN_USE_DICOM_
    using typesGenericTallies = std::tuple<pen_EdepMat,
					   pen_EdepBody,
					   pen_SphericalDoseDistrib,
					   pen_tallyTracking,
					   pen_EmergingPartDistrib,
					   pen_CylindricalDoseDistrib,
					   pen_ImpactDetector,
					   pen_tallySecondary,
					   pen_SpatialDoseDistrib,
					   pen_AngularDet,
					   pen_tallyPhaseSpaceFile,
					   pen_tallyKermaTrackLength,
					   pen_DICOMDoseDistrib,
					   pen_CTsinogram,
					   pen_tallyDICOMkerma,
					   pen_PSS,
					   pen_DetectionSpatialDistrib,
					   pen_EmergingSphericalDistrib,
					   pen_Singles>;
#else
    using typesGenericTallies = std::tuple<pen_EdepMat,
					   pen_EdepBody,
					   pen_SphericalDoseDistrib,
					   pen_tallyTracking,
					   pen_EmergingPartDistrib,
					   pen_CylindricalDoseDistrib,
					   pen_ImpactDetector,
					   pen_tallySecondary,
					   pen_SpatialDoseDistrib,
					   pen_AngularDet,
					   pen_tallyPhaseSpaceFile,
					   pen_tallyKermaTrackLength,
					   pen_CTsinogram,
					   pen_PSS,
					   pen_DetectionSpatialDistrib,
					   pen_EmergingSphericalDistrib,
					   pen_Singles>;    
#endif

    template<size_t T>
    typename std::enable_if<T >= std::tuple_size<typesGenericTallies>::value, bool>::type
    checkRegistersGeneric(const unsigned){ return true; }
    
    template<size_t T>
    typename std::enable_if<T < std::tuple_size<typesGenericTallies>::value, bool>::type
    checkRegistersGeneric(const unsigned verbose){
      using tupleType = typename std::tuple_element<T, typesGenericTallies>::type;
      int val = tupleType::registerStatus();
      if(val != 0){
	if(verbose > 0){
	  printf("Warning: generic tally type '%s' register failed."
		 " Error code: %d\n",
		 tupleType::___ID, val);
	}
	checkRegistersGeneric<T+1>(verbose);
	return false;
      }else{
	return checkRegistersGeneric<T+1>(verbose);
      }
    }
    
    template<>
    bool checkRegistered<pen_particleState>(const unsigned verbose);
  }
}

#endif
