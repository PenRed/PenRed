 

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

    int detectedPart::getSpectrum(std::vector<double>& spectrum,
				  const std::vector<detectedPart>& data,
				  const double eMin,
				  const double eMax,
				  const unsigned nBins,
				  unsigned long long& nHists){

      if(eMin >= eMax)
	return -1;

      if(data.size() == 0)
	return -2;
	     
      //Get total number of histories
      if(nHists == 0){
	for(const detectedPart& part : data){
	  nHists += part.dh;
	}
      }
      
      const double dE =
	(eMax - eMin)/static_cast<double>(nBins);
      
      spectrum.resize(nBins);
      std::fill(spectrum.begin(), spectrum.end(), 0.0);
      for(const detectedPart& p : data){
	if(p.state.E >= eMin && p.state.E < eMax){
	  unsigned i = static_cast<unsigned>((p.state.E-eMin)/dE);
	  spectrum[i] += p.state.WGHT;
	}
      }

      //Normalize
      for(double& b : spectrum){
	b /= static_cast<double>(nHists)*dE;
      }
	  
      return 0;
    }

    int detectedPart::printSpectrum(const std::vector<detectedPart>& data,
				    const std::string& filename,
				    const double eMin,
				    const double eMax,
				    const unsigned nBins,
				    unsigned long long& nHists){

      if(data.size() == 0)
	return -2;
	
      FILE* fspec = nullptr;
      fspec = fopen(filename.c_str(), "w");
      if(fspec != nullptr){

	//Generate spectrum
	std::vector<double> spectrum;
	int err = getSpectrum(spectrum, data, eMin, eMax, nBins, nHists);
	if(err != 0)
	  return err;

	const double dE =
	  (eMax - eMin)/static_cast<double>(nBins);
	
	//Print spectrum
	fprintf(fspec,"# Particles registered: %llu\n"
		"# Histories registered: %llu\n"
		"# Energy range: [%E, %E) eV\n"
		"# Elow (eV) | weight (1/eV) \n",
		static_cast<unsigned long long>(data.size()),
		nHists, eMin, eMax);
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

    int detectedPart::getProfileX(std::vector<double>& profile,
				  const std::vector<detectedPart>& data,
				  const unsigned nBins,
				  const double xMin, const double xMax,
				  const double yMin, const double yMax,
				  const double zMin, const double zMax,
				  unsigned long long& nHists){

      if(xMin >= xMax ||
	 yMin >= yMax ||
	 zMin >= zMax)
	return -1;
      
      if(data.size() == 0)
	return -2;

      const double dx =
	(xMax - xMin)/static_cast<double>(nBins);

      profile.resize(nBins);
      std::fill(profile.begin(), profile.end(), 0.0);
      for(const detectedPart& p : data){
	if(p.state.X >= xMin && p.state.X < xMax &&
	   p.state.Y >= yMin && p.state.Y < yMax &&
	   p.state.Z >= zMin && p.state.Z < zMax){
	  unsigned ix = static_cast<unsigned>((p.state.X-xMin)/dx);
	  profile[ix] += p.state.WGHT;
	}
      }

      //Normalize
      for(double& b : profile){
	b /= static_cast<double>(nHists)*dx;
      }      
      
      return 0;
    }

    int detectedPart::printProfileX(const std::vector<detectedPart>& data,
				    const std::string& filename,
				    const unsigned nBins,
				    const double xMin, const double xMax,
				    const double yMin, const double yMax,
				    const double zMin, const double zMax,
				    unsigned long long& nHists){

      if(data.size() == 0)
	return -2;
	
      FILE* fspec = nullptr;
      fspec = fopen(filename.c_str(), "w");
      if(fspec != nullptr){

	//Generate profile
	std::vector<double> profile;
	int err = getProfileX(profile, data, nBins,
			      xMin, xMax,
			      yMin, yMax,
			      zMin, zMax,
			      nHists);
	if(err != 0)
	  return err;

	const double dx =
	  (xMax - xMin)/static_cast<double>(nBins);
	
	//Print profile
	fprintf(fspec,"# Particles registered: %llu\n"
		"# Histories registered: %llu\n"
		"# x range: [%E, %E) cm\n"
		"# Xlow (cm) | weight (1/cm) \n",
		static_cast<unsigned long long>(data.size()),
		nHists, xMin, xMax);
	for(size_t i = 0; i < profile.size(); ++i){
	  fprintf(fspec, " %15.5E    %15.5E\n",
		  xMin + dx*static_cast<double>(i),
		  profile[i]);
	}
	fclose(fspec);
      }
      else{
	printf("penred:simulation:detectedPart:printProfileX:Error: "
	       "Unable to open file '%s'\n",filename.c_str());
	return -3;
      }

      return 0;
    }    

    int detectedPart::getProfileY(std::vector<double>& profile,
				  const std::vector<detectedPart>& data,
				  const unsigned nBins,
				  const double xMin, const double xMax,
				  const double yMin, const double yMax,
				  const double zMin,
				  const double zMax,
				  unsigned long long& nHists){

      if(xMin >= xMax ||
	 yMin >= yMax ||
	 zMin >= zMax)
	return -1;
      
      if(data.size() == 0)
	return -2;

      const double dy =
	(yMax - yMin)/static_cast<double>(nBins);

      profile.resize(nBins);
      std::fill(profile.begin(), profile.end(), 0.0);
      for(const detectedPart& p : data){
	if(p.state.X >= xMin && p.state.X < xMax &&
	   p.state.Y >= yMin && p.state.Y < yMax &&
	   p.state.Z >= zMin && p.state.Z < zMax){
	  unsigned iy = static_cast<unsigned>((p.state.Y-yMin)/dy);
	  profile[iy] += p.state.WGHT;
	}
      }

      //Normalize
      for(double& b : profile){
	b /= static_cast<double>(nHists)*dy;
      }      
      
      return 0;
    }

    int detectedPart::printProfileY(const std::vector<detectedPart>& data,
				    const std::string& filename,
				    const unsigned nBins,
				    const double xMin, const double xMax,
				    const double yMin, const double yMax,
				    const double zMin, const double zMax,
				    unsigned long long& nHists){

      if(data.size() == 0)
	return -2;
	
      FILE* fspec = nullptr;
      fspec = fopen(filename.c_str(), "w");
      if(fspec != nullptr){

	//Generate profile
	std::vector<double> profile;
	int err = getProfileY(profile, data, nBins,
			      xMin, xMax,
			      yMin, yMax,
			      zMin, zMax,
			      nHists);
	if(err != 0)
	  return err;

	const double dy =
	  (yMax - yMin)/static_cast<double>(nBins);
	
	//Print profile
	fprintf(fspec,"# Particles registered: %llu\n"
		"# Histories registered: %llu\n"
		"# y range: [%E, %E) cm\n"
		"# Ylow (cm) | weight (1/cm) \n",
		static_cast<unsigned long long>(data.size()),
		nHists, yMin, yMax);
	for(size_t i = 0; i < profile.size(); ++i){
	  fprintf(fspec, " %15.5E    %15.5E\n",
		  yMin + dy*static_cast<double>(i),
		  profile[i]);
	}
	fclose(fspec);
      }
      else{
	printf("penred:simulation:detectedPart:printProfileX:Error: "
	       "Unable to open file '%s'\n",filename.c_str());
	return -3;
      }

      return 0;
    }        
    
    int detectedPart::printSpatialXY(const std::vector<detectedPart>& data,
				     const std::string& filename,
				     const unsigned nBinsX, const unsigned nBinsY,
				     const double xMin, const double xMax,
				     const double yMin, const double yMax,
				     const double zMin, const double zMax,
				     unsigned long long& nHists){
      
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
	  (xMax - xMin)/static_cast<double>(nBinsX);

	const double dy =
	  (yMax - yMin)/static_cast<double>(nBinsY);
	
	std::vector<double> map(nBinsX*nBinsY,0.0);
	for(const detectedPart& p : data){
	  if(p.state.X >= xMin && p.state.X < xMax &&
	     p.state.Y >= yMin && p.state.Y < yMax &&
	     p.state.Z >= zMin && p.state.Z < zMax){
	    unsigned ix = static_cast<unsigned>((p.state.X-xMin)/dx);
	    unsigned iy = static_cast<unsigned>((p.state.Y-yMin)/dy);
	    map[iy*nBinsX + ix] += p.state.WGHT;
	  }
	}

	//Normalize
	for(double& b : map){
	  b /= static_cast<double>(nHists)*dx*dy;
	}
	  
	fprintf(fmap,"# Particles registered: %llu\n"
		"# Histories registered: %llu\n"
		"# x range: [%E, %E) cm\n"
		"# y range: [%E, %E) cm\n"
		"# z range: [%E, %E) cm\n"
		"# Units are 1/hist*cm**2\n",
		static_cast<unsigned long long>(data.size()),
		nHists,
		xMin, xMax,
		yMin, yMax,
		zMin, zMax);

	//Print first line with bin X positions
	fprintf(fmap, " %u ", nBinsX);
	for(unsigned i = 0; i < nBinsX; ++i){
	  fprintf(fmap, " %.5f", xMin + (static_cast<double>(i) + 0.5)*dx);
	}
	fprintf(fmap, "\n");
	
	
	for(unsigned iy = 0; iy < nBinsY; ++iy){
	  //Print Y bin position
	  fprintf(fmap, " %.5f", yMin + (static_cast<double>(iy) + 0.5)*dy);
	  for(unsigned ix = 0; ix < nBinsX; ++ix){
	    fprintf(fmap, " %15.5E", map[iy*nBinsX + ix]);
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
    
  } // namespace simulation
  
} // namespace penred
