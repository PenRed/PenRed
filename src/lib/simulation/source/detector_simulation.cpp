 

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
//        vicent.gimenez.alventosa@gmail.com
//        sanolgi@upvnet.upv.es
//        vicente.gimenez@uv.es
//    
//

#include "detector_simulation.hh"

namespace penred{

  namespace simulation{

    int detectedPart::printSpectrum(const std::vector<detectedPart>& data,
				    const std::string& filename,
				    const double eMin,
				    const double eMax,
				    const unsigned nBins,
				    unsigned long long nHists){

      if(eMin >= eMax)
	return -1;

      if(data.size() == 0)
	return -2;
	
      FILE* fspec = nullptr;
      fspec = fopen(filename.c_str(), "w");
      if(fspec != nullptr){

	//Get total number of histories
	if(nHists == 0){
	  for(const detectedPart& part : data){
	    nHists += part.dh;
	  }
	}
	const double dE =
	  (eMax - eMin)/static_cast<double>(nBins);
	std::vector<double> spectrum(nBins,0.0);
	for(const detectedPart& p : data){
	  if(p.state.E >= eMin && p.state.E < eMax){
	    unsigned i = static_cast<unsigned>((p.state.E-eMin)/dE);
	    spectrum[i] += p.state.WGHT;
	  }
	}

	//Normalize
	for(double& b : spectrum){
	  b /= static_cast<double>(nHists);
	}
	  
	fprintf(fspec,"# Particles registered: %llu\n"
		"# Histories registered: %llu\n"
		"# Elow (eV) | weight (1/eV) \n",
		static_cast<unsigned long long>(data.size()),
		nHists);
	for(size_t i = 0; i < spectrum.size(); ++i){
	  fprintf(fspec, " %15.5E    %15.5E\n",
		  eMin + dE*static_cast<double>(i),
		  spectrum[i]);
	}
	fclose(fspec);
      }
      else{
	printf("penred:simulation:detectedPart:printSpectrum:Error: "
	       "Unable to open file '%s'\n",filename.c_str());
	return -3;
      }

      return 0;
    }

    
    int detectedPart::printSpatialXY(const std::vector<detectedPart>& data,
				     const std::string& filename,
				     const unsigned nBins,
				     const double xMin, const double xMax,
				     const double yMin, const double yMax,
				     const double zMin, const double zMax,
				     unsigned long long nHists){
      
      if(xMin >= xMax ||
	 yMin >= yMax ||
	 zMin >= zMax)
	return -1;
      
      if(data.size() == 0)
	return -2;
	
      FILE* fmap = nullptr;
      fmap = fopen(filename.c_str(), "w");
      if(fmap != nullptr){

	if(nHists == 0){
	  for(const detectedPart& part : data){
	    nHists += part.dh;
	  }
	}
	
	const double dx =
	  (xMax - xMin)/static_cast<double>(nBins);

	const double dy =
	  (yMax - yMin)/static_cast<double>(nBins);
	
	std::vector<double> map(nBins*nBins,0.0);
	for(const detectedPart& p : data){
	  if(p.state.X >= xMin && p.state.X < xMax &&
	     p.state.Y >= yMin && p.state.Y < yMax &&
	     p.state.Z >= zMin && p.state.Z < zMax){
	    unsigned ix = static_cast<unsigned>((p.state.X-xMin)/dx);
	    unsigned iy = static_cast<unsigned>((p.state.Y-yMin)/dy);
	    map[iy*nBins + ix] += p.state.WGHT;
	  }
	}

	//Normalize
	for(double& b : map){
	  b /= static_cast<double>(nHists);
	}
	  
	fprintf(fmap,"# Particles registered: %llu\n"
		"# Histories registered: %llu\n"
		"# x range: [%E, %E) cm\n"
		"# y range: [%E, %E) cm\n"
		"# z range: [%E, %E) cm\n"
		"# Units are 1/hist*cm**3\n",
		static_cast<unsigned long long>(data.size()),
		nHists,
		xMin, xMax,
		yMin, yMax,
		zMin, zMax);

	//Print first line with bin X positions
	fprintf(fmap, " %u ", nBins);
	for(unsigned i = 0; i < nBins; ++i){
	  fprintf(fmap, " %.5f", xMin + (static_cast<double>(i) + 0.5)*dx);
	}
	fprintf(fmap, "\n");
	
	
	for(unsigned iy = 0; iy < nBins; ++iy){
	  //Print Y bin position
	  fprintf(fmap, " %.5f", yMin + (static_cast<double>(iy) + 0.5)*dy);
	  for(unsigned ix = 0; ix < nBins; ++ix){
	    fprintf(fmap, " %15.5E", map[iy*nBins + ix]);
	  }
	  fprintf(fmap, "\n");
	}
	fclose(fmap);
      }
      else{
	printf("penred:simulation:detectedPart:printSpatialXY:Error: "
	       "Unable to open file '%s'\n", filename.c_str());
	return -3;
      }

      return 0;
    }
    
  };
  
};