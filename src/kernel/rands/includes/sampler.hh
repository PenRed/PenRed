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


#ifndef __PEN_SAMPLING_CLASSES__
#define __PEN_SAMPLING_CLASSES__

#include <numeric>
#include <array>
#include <vector>

namespace penred{

  namespace sampling{

    template<size_t dim = 1>
    class aliasing{
      //Implements Walker's aliasing algorithm
      //for sampling discrete probability distributions

    public:
      enum errors{    
	SUCCESS = 0,
	INVALID_NUMBER_OF_BINS,
	NUMBER_OF_BINS_MISMATCH,
	INVALID_LIMITS,
      };            
    private:

      std::array<unsigned long, dim> nBins;
      unsigned long totalBins;
      std::array<unsigned long, dim> binsPerIncrement;
      std::array<double, dim> binWidths;
      std::array<std::pair<double, double>, dim> limits;

      std::vector<double> cutoff;
      std::vector<unsigned> alias;
      
    public:

      static constexpr size_t dimensions = dim;      

      int init(const std::vector<double>& data,
	       const std::array<unsigned long, dim>& nBinsIn,
	       const std::array<std::pair<double, double>, dim>& limitsIn){

	//Init aliasing algorithm. Notice that provided probabilities must be positive or zero.
	//
	//
	// Input:
	//
	// data    : Probabilities of each bin (not necessarily normalized to unity)
	// nBinsIn : Number of bins in each dimension. Total number of bins must match with "data" size 
	// limitsIn: Pairs with limits for each dimension
	//
	// Return error code or "SUCCESS" on success execution

	//Calculate total number of bins
	totalBins = std::accumulate(nBinsIn.begin(),
				    nBinsIn.end(), 1,
				    std::multiplies<unsigned long>());
	if(totalBins == 0)
	  return errors::INVALID_NUMBER_OF_BINS;

	if(totalBins != data.size()){
	  return errors::NUMBER_OF_BINS_MISMATCH;
	}
    
	//Check limits
	for(size_t i = 0; i < dim; ++i){
	  if(limitsIn[i].first >= limitsIn[i].second){
	    return errors::INVALID_LIMITS;
	  }
	}

	//Save bins
	nBins = nBinsIn;

	//Calculate bins per increment in each dimension
	binsPerIncrement[0] = 1;
	for(size_t i = 1; i < dim; ++i){
	  binsPerIncrement[i] = binsPerIncrement[i-1]*nBins[i-1];
	}

	//Save limits
	limits = limitsIn;

	//Calculate bin widths
	for(size_t i = 0; i < dim; ++i){
	  binWidths[i] = (limits[i].second - limits[i].first)/static_cast<double>(nBins[i]);
	}

	//Save initial cutoffs
	cutoff = data;
	
	//Resize alias vector
	alias.resize(totalBins);
	//Fill alias with the next bin index
	std::iota(alias.begin(), alias.end(), 1);

	//Get normalization factor
	const double sum = std::accumulate(data.cbegin(), data.cend(), 0.0);
	const double normFact = static_cast<double>(totalBins)/sum;

	//Normalize initial cutoffs
	for(double& c : cutoff){
	  c *= normFact;
	}

	//Compute alias and cutoffs
	for(unsigned long i = 0; i < totalBins-1; ++i){

	  double lowVal = 1.0;
	  double highVal = 1.0;
	  unsigned ilow = totalBins;
	  unsigned ihigh = totalBins;
	  
	  //Get maximum and minimum value
	  for(unsigned long j = 0; j < totalBins; ++j){

	    if(alias[j] == j+1){
	      
	      if(cutoff[j] < lowVal){
		lowVal = cutoff[j];
		ilow = j;
	      }
	      else if(cutoff[j] > highVal){
		highVal = cutoff[j];
		ihigh = j;		
	      }		
	    }
	  }
	  if(ilow == totalBins || ihigh == totalBins){ return errors::SUCCESS; }
	  alias[ilow] = ihigh;
	  cutoff[ihigh] = highVal + lowVal - 1.0;
	}

	return errors::SUCCESS;    
      }
      
      unsigned long sample(pen_rand& random) const {

	double r = random.rand()*totalBins;
	unsigned long rInt = static_cast<unsigned long>(r);
	double tst = r - static_cast<double>(rInt);
	if(tst > cutoff[rInt]){
	  return alias[rInt];
	}
	else{
	  return rInt;
	}	
      }

      inline std::array<unsigned long, dim> sampleByDim(pen_rand& random) const {

	//Sample global bin
	unsigned long globBin = sample(random);

	std::array<unsigned long, dim> result;
	for(long int i = dim-1; i >= 0; --i){
	  result[i] = globBin/binsPerIncrement[i];
	  globBin -= result[i]*binsPerIncrement[i];
	}

	return result;
      }

      inline std::array<double, dim> samplePositions(pen_rand& random) const {

	//Sample index in each dimension
	const std::array<unsigned long, dim> indexes = sampleByDim(random);

	//Sample uniform position inside each dimension
	std::array<double, dim> result;
	for(size_t i = 0; i < dim; ++i){
	  result[i] = limits[i].first + (static_cast<double>(indexes[i]) + random.rand())*binWidths[i];
	}

	return result;
      }      
      
    };
    
  }; // namespace sampling
  
}; // namespace penred


#endif
