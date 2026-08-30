
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


#include "distributionSE.hh"

//Reader functions
int tallyReader_DistributionSE::storeElement(const std::string& pathInSection,
                                             const pen_parserData& element,
                                             const unsigned){

  if(pathInSection.compare("distribution/origin/x") == 0){
    distribOrigin.x = element;
  }
  else if(pathInSection.compare("distribution/origin/y") == 0){
    distribOrigin.y = element;
  }
  else if(pathInSection.compare("distribution/origin/z") == 0){
    distribOrigin.z = element;
  }  
  else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
  
}

int tallyReader_DistributionSE::storeString(const std::string& pathInSection,
                                            const std::string& element,
                                            const unsigned /*verbose*/){

  if(pathInSection.compare("distribution/filename") == 0){
    distribFilename = element;
  }
  else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
}

int distribSE_specificSampler::configure(double& Emax,
                                         const pen_parserSection& config,
                                         const unsigned /*nthreads*/,
                                         const unsigned verbose){

  
  //Only the thread 0 will initialize the shared samplers and distributions
  if(getThread() == 0){

    //Try to configure the post-sampling shift spatial sampler
    pen_parserSection shiftConf;
    if(config.readSubsection("position/shift/config",shiftConf) == INTDATA_SUCCESS){
      //Read sampler type
      std::string shiftType;
      if(shiftConf.read("type", shiftType) != INTDATA_SUCCESS){
        
      }

      //Instantiate the sampler
      shiftSampler = pen_genericStateGen::instanceSpatial(shiftType.c_str());
      if(!shiftSampler){
        printf("distribSE:config:Error: Unable to instantiate spatial sampler with type '%s'",
               shiftType.c_str());
        return -1;
      }

      //Configure the sampler
      int errConfig = shiftSampler->configure(shiftConf, verbose);
      if(errConfig != 0){
        if(verbose > 0){
          printf("distribSE:config:Error: Unable to "
                 "configurate spatial sampler '%s'\n",shiftType.c_str());
        }
        return -2;
      }

      //Update sampler geometry
      if(geo() != nullptr){
        shiftSampler->updateGeometry(geo());
      }
    }
    
    //Read information from config section
    tallyReader_DistributionSE reader;
    int err = reader.read(config,verbose);
    if(err != tallyReader_DistributionSE::SUCCESS){
      return err;
    }

    //Save distribution origin
    distribOrigin = reader.distribOrigin;

    //Read distribution
    distribSE = std::make_shared<penred::measurements::results<double, 4>>();
    std::ifstream fin(reader.distribFilename, std::ifstream::in);
    if(!fin){
      if(verbose > 0){
        printf("distribSE_specificSampler: Error: Unable to "
               "open distribution file '%s'",
               reader.distribFilename.c_str());
      }
      return -1;
    }
  
    err = distribSE->read(fin);
    if(err != 0){
      if(verbose > 0){
        printf("distribSE_specificSampler: Error: Unable to "
               "read distribution file '%s'.\n  Error code: %d\n",
               reader.distribFilename.c_str(), err);
      }
      return -2;
    }

    //Get energy bins
    ebins = distribSE->getNBins(0);
    eBinWidth = distribSE->readBinWidth(0);

    //Get limits
    elimits = distribSE->readLimits()[0];
    xlimits = distribSE->readLimits()[1];
    ylimits = distribSE->readLimits()[2];
    zlimits = distribSE->readLimits()[3];

    //Get spatial bins
    xbins = distribSE->getNBins(1);
    xBinWidth = distribSE->readBinWidth(1);

    ybins = distribSE->getNBins(2);
    yBinWidth = distribSE->readBinWidth(2);

    zbins = distribSE->getNBins(3);
    zBinWidth = distribSE->readBinWidth(3);
    
    //Set Emax
    Emax = elimits.second;
  }
  
  return 0;
}

void distribSE_specificSampler::sample(pen_particleState& state,
                                       pen_KPAR& /*genKpar*/,
                                       unsigned long long& dhist,
                                       pen_rand& random){

  //Sample position
  spatial()->sample(state,random);

  //Calculate direction
  vector3D<double> dir(state.X - distribOrigin.x,
                       state.Y - distribOrigin.y,
                       state.Z - distribOrigin.z);
  if(dir.mod2() < 1.0e-8){
    //Sampled position at the distribution origin. Set direction randomly
    const double phi = random.rand()*2.0*constants::PI;
    const double cosp = 1.0 - random.rand()*2.0;

    const double sinp = sqrt(1.0 - cosp*cosp);
    const double cosa = cos(phi);
    const double sina = sin(phi);

    dir.x = sinp * cosa;
    dir.y = sinp * sina;
    dir.z = cosp;
  }
    
  //Set direction
  dir.normalize();
  state.U = dir.x;
  state.V = dir.y;
  state.W = dir.z;

  //Obtain the energy distribution for the sampled position
  penred::measurements::results<double, 1> spectrum;
  distribSE->extractSpectrum1D(0, state.X, state.Y, state.Z, spectrum);

  //Sample and apply post-sampling shift
  if(shiftSampler){
    pen_particleState shift;
    shiftSampler->sample(shift,random);

    state.X += shift.X;
    state.Y += shift.Y;
    state.Z += shift.Z;
  }
  
  //Set ILB as primary
  state.ILB[0] = 1;
  
  //Increase history counter
  dhist = 1;

  // * Set the energy

  //Get the spectrum data
  const std::vector<double>& sData = spectrum.readData();

  //Find global maximum
  const double maxVal = *std::max_element(sData.begin(), sData.end());

  // Scale them up for better precision
  double scaleFactor;
  if(maxVal > 0.0) {
    scaleFactor = 1.0 / maxVal;
  }
  else{
    //Handle zero distribution
    state.E = elimits.first + random.rand() * (elimits.second - elimits.first);
    return;
  }

  //Build CDF using scaled values
  std::vector<double> cdf(sData.size() + 1);
  cdf[0] = 0.0;
  for(size_t i = 0; i < sData.size(); ++i) {
    cdf[i+1] = cdf[i] + sData[i] * scaleFactor;
  }

  //Normalize CDF
  const double totalSum = cdf.back();
  const double invTotalSum = 1.0 / totalSum;
  for(auto& val : cdf) {
    val *= invTotalSum;
  }

  //Sample
  double rval = random.rand();
  auto it = std::lower_bound(cdf.begin(), cdf.end(), rval);
  size_t bin = std::distance(cdf.begin(), it) - 1;
  if(bin >= sData.size()) bin = sData.size() - 1;  
  state.E = elimits.first + eBinWidth * (static_cast<double>(bin) + random.rand());
}

int distribSE_specificSampler::sharedConfig(const distribSE_specificSampler& o){
  //Share distribution and origin sampler
  shiftSampler = o.shiftSampler;
  distribSE = o.distribSE;

  //Set energy bins and limits
  ebins = o.ebins;
  eBinWidth = o.eBinWidth;
  elimits = o.elimits;
  xlimits = o.xlimits;
  ylimits = o.ylimits;
  zlimits = o.zlimits;

  xbins = o.xbins;
  xBinWidth = o.xBinWidth;

  ybins = o.ybins;
  yBinWidth = o.yBinWidth;

  zbins = o.zbins;
  zBinWidth = o.zBinWidth;
  
  return 0;
}

REGISTER_SPECIFIC_SAMPLER(distribSE_specificSampler,pen_particleState, COMBINED_DISTRIBUTIONS)
