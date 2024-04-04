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
//        sanolgi@upvnet.upv.es
//    
//

 
#ifndef __PEN_DETECTOR_SIMULATION_FUNCTIONS__
#define __PEN_DETECTOR_SIMULATION_FUNCTIONS__

#include "base_simulation_functions.hh"
#include <cstdio>
#include <vector>

namespace penred{

  namespace simulation{

    //Auxiliary structure to store a detected particle
    struct detectedPart{
      pen_particleState state; //Particle state
      unsigned long long dh;        //History increment
      pen_KPAR kpar;           //Particle type

      detectedPart() = default;
      inline detectedPart(const pen_particleState& stateIn,
			  const unsigned long long dhIn,
			  const pen_KPAR kparIn) : state(stateIn),
						   dh(dhIn),
						   kpar(kparIn){}
      
      static int printSpectrum(const std::vector<detectedPart>& data,
			       const std::string& filename,
			       const double eMin,
			       const double eMax,
			       const unsigned nBins,
			       unsigned long long nHists = 0);

      
      static inline int printSpectrum(const std::vector<detectedPart>& data,
				      const std::string& filename,
				      const unsigned nBins,
				      unsigned long long nHists = 0){
	
	//Get maximum and minimum energy
	double maxE = 0.0;
	double minE = 1.0e35;
	for(const detectedPart& part : data){
	  if(maxE < part.state.E)
	    maxE = part.state.E;
	  if(minE > part.state.E)
	    minE = part.state.E;
	}

	return printSpectrum(data, filename, minE, maxE, nBins, nHists);	
      }


      static int printSpatialXY(const std::vector<detectedPart>& data,
				const std::string& filename,
				const unsigned nBins,
				const double xMin, const double xMax,
				const double yMin, const double yMax,
				const double zMin = -1.0e35,
				const double zMax =  1.0e35,
				unsigned long long nHists = 0);

      static inline int printSpatialXY(const std::vector<detectedPart>& data,
				       const std::string& filename,
				       const unsigned nBins,
				       unsigned long long nHists = 0){
	
	//Get maximum and minimum energy
	double xMin =  1.0e35;
	double xMax = -1.0e35;
	double yMin =  1.0e35;
	double yMax = -1.0e35;
	
	for(const detectedPart& part : data){

	  //X
	  if(xMax < part.state.X)
	    xMax = part.state.X;
	  if(xMin > part.state.X)
	    xMin = part.state.X;

	  //Y
	  if(yMax < part.state.Y)
	    yMax = part.state.Y;
	  if(yMin > part.state.Y)
	    yMin = part.state.Y;	  
	}

	return printSpatialXY(data, filename, nBins,
			      xMin, xMax, yMin, yMax, nHists);	
      }
      
    };

    template<class stateType, class contextType, class... vrTypes>
    inline int simulateAndDetect(simConfig& config,
				 const contextType& context,
				 std::vector<detectedPart>& result,
				 const std::string& sourceName,
				 sampleFuncType<stateType>& fsource,
				 const unsigned long long nhists,
				 const unsigned detector,
				 const vrTypes&... vr){

      //This function simulates the sampled particles across the geometry
      //and stores all the particles impinging the specified detector
      //
      //  result: Vector with the registered particles
      //  fsource: Sampling function
      //  nhists: Histories to be simulated
      //  detector: Detector index where the particles are registered
      //  particles: List of the particle instances involved in the simulation
      //
      
      //Create the tally function
      unsigned long long lastTalliedHist = 0;
      tallyFuncType ftallyFuncs = [&result, &lastTalliedHist]
	(const pen_particleState& state,
	 const pen_KPAR kpar,
	 const unsigned long long nhist,
	 const int ret){
	if(ret == 1){
	  result.emplace_back(state,nhist-lastTalliedHist,kpar);
	  lastTalliedHist = nhist; 
	}
      };

      return sampleAndSimulateCondContext(config, context, nhists,
					  sourceName, fsource,
					  finishTypes::DETECTOR_REACHED, //Detector finish type
					  detector, //Detector index as finish value
					  ftallyFuncs,
					  vr...);
      
    }


    template<class contextType, class... vrTypes>
    inline int simulateVectorAndDetect(simConfig& config,
				       const contextType& context,
				       std::vector<detectedPart>& result,
				       const std::vector<detectedPart>& initParticles,
				       const unsigned long long initHist,
				       const unsigned long long nhists,
				       const std::string& sourceName,
				       const unsigned detector,
				       const vrTypes&... vr){

      //Skip particles until desired history      
      size_t nextPart = 0;
      unsigned long long actualHist = 0;
      for(; nextPart < initParticles.size(); ++nextPart){
	actualHist += initParticles[nextPart].dh;
	if(actualHist > initHist){
	  if(nextPart > 0)
	    --nextPart;
	  break;
	}
      }

      sampleFuncType<pen_particleState> fsource =
	[&nextPart, &initParticles]
	(pen_particleState& state,    //Generated state
	 pen_KPAR& kpar,          //Generated kpar
	 unsigned long long& dh,  //History increment
	 const unsigned,          //Thread number
	 pen_rand&){              //Random generator
	  
	  if(nextPart < initParticles.size()){
	    const detectedPart& p = initParticles[nextPart++];
	    state = p.state;
	    kpar = p.kpar;
	    dh = p.dh;
	  }
	  else{
	    kpar = ALWAYS_AT_END;
	  }
	};

      return simulateAndDetect(config, context, result, sourceName, fsource,
			       nhists, detector, vr...);
    }

    template<class contextType, class... vrTypes>
    inline int simulateVectorAndDetect(simConfig& config,
				       const contextType& context,
				       std::vector<detectedPart>& result,
				       const std::vector<detectedPart>& initParticles,
				       const unsigned long long nhists,
				       const std::string& sourceName,
				       const unsigned detector,
				       vrTypes&... vr){

      return simulateSubVectorAndDetect(config, context, result, initParticles,
					0, nhists, sourceName, detector,
					vr...);
    }

  };
};

#endif
