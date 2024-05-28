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
//    
//
 
#include "device.hh"


namespace penred{

  namespace xray{

    int constructDevice(std::ostream& out,
			const double focalSpot,
			const double source2det,
			const double source2filter,
			const double source2bowtie,
			const double detectorDx,
			const double detectorDy,
			const double inherentFilterSize,
			const std::vector<double>& filters,
			std::vector<double> bowtieDz,
			const vector3D<double> sourcePos,
			const bool constructAnode,
			const double anodeAngle,
			const unsigned verbose){

      //This function constructs a mesh based geometry of a x-ray device
      //following the specifications provided by the function parameters

      // ** Distances
      
      //Calculate and check distances

      constexpr double collHeight = 1.0;
      constexpr double elementSpacing = 0.5;
      constexpr double source2inherentCollTop = 1.0;
      constexpr double source2inherentCollBot = source2inherentCollTop + collHeight;
      constexpr double source2inherentFilter = source2inherentCollBot + elementSpacing;

      const double inherentFilter2filters =
	source2filter - (source2inherentFilter + inherentFilterSize);
	
      if(inherentFilter2filters < 4.0*elementSpacing+2.0*collHeight){
	if(verbose > 1){
	  printf("constructDevice: Error: The minimum distance between "
		 "inherent filter and first added filter must be %f cm\n"
		 "    Inherent filter to source: %f\n"
		 "    Added filter to source   : %f\n",
		 4.0*elementSpacing+2.0*collHeight,
		 source2inherentFilter,
		 source2filter);
	}
	return errors::INVALID_DISTANCE;
      }

      const double source2filterCollBot = source2filter-elementSpacing;
      const double source2filterCollTop = source2filterCollBot-collHeight;

      //Calculate total filters width
      double filtersWidth = 0.0;
      for(const double& f : filters){
	filtersWidth += f;
      }

      double source2filtersEnd = source2filter + filtersWidth;
      if(source2bowtie > 0.0 && bowtieDz.size() > 0){
	const double maxDz = *std::max_element(bowtieDz.cbegin(), bowtieDz.cend());
	const double source2BowtieBot = source2bowtie + maxDz;
	
	source2filtersEnd = source2BowtieBot;
      }
      const double source2detCollTop = source2filtersEnd + elementSpacing;
      const double source2detCollBot = source2detCollTop + collHeight;

      if(verbose > 1){
	printf("\n + X-ray characteristics:\n"
	       "     Source to detector   : %.4f cm\n"
	       "     Source to inherent f.: %.4f cm\n"
	       "     Source to filters    : %.4f cm\n"
	       "     Detector X size      : %.4f cm\n"
	       "     Detector Y size      : %.4f cm\n",
	       source2det,
	       source2inherentFilter,
	       source2filter,
	       detectorDx,
	       detectorDy);
      }
      
      // ** Number of objects
      unsigned nBodies = 2; //World and detector
      if(constructAnode)
	++nBodies; //Anode
      
      ++nBodies; //First collimator

      if(inherentFilterSize > 0.0)
	nBodies += 2; //Inherent filter and collimator

      if(filters.size() > 0)
	nBodies += filters.size() + 1; //Filters and collimator

      if(source2bowtie > 0.0)
	nBodies += 1; //Bowtie

      out << "# Number of bodies" << std::endl;
      out << " " << nBodies << std::endl;
      out << "#" << std::endl;

      // ** Material index counter
      unsigned collMat = 1; //Collimator material
      unsigned nextMat = 2;
      
      // ** Anode
      //-----------------

      //Check if the anode must be simulated
      if(constructAnode){

	if(verbose > 1){
	  printf("\n + Creating anode: \n"
		 "     Angle          : %.3f deg\n",
		 anodeAngle);
	}

	createAnode(out, anodeAngle, nextMat++,
		    1.0,2.0,1.0,
		    "anode", "world", false,
		    sourcePos);
      }

      // ** Anode collimator
      //

      // Add a collimator before the innerent filter to speed-up the simulation

      //Calculate field size at top
      const std::pair<double, double> inherentCollTopSizes =
	fieldSize(focalSpot,
		  detectorDx,
		  detectorDy,
		  source2det,
		  source2inherentCollTop);

      //Calculate field size at bot
      const std::pair<double, double> inherentCollBotSizes =
	fieldSize(focalSpot,
		  detectorDx,
		  detectorDy,
		  source2det,
		  source2inherentCollBot);

      //Create collimator
      vector3D<double> inherentCollCenter = sourcePos;
      inherentCollCenter.z -= source2inherentCollTop + collHeight/2.0;
      createBaseCollimator(std::min(2.0*inherentCollBotSizes.first, detectorDx),
			   inherentCollTopSizes.first,
			   inherentCollBotSizes.first,
			   std::min(2.0*inherentCollBotSizes.second, detectorDy),
			   inherentCollTopSizes.second,
			   inherentCollBotSizes.second,
			   collHeight,
			   out,
			   collMat,
			   "inherent-collimator",
			   "world",
			   false,
			   inherentCollCenter);

      // ** Inherent filter

      if(inherentFilterSize > 0.0){
      
	//Create an inherent filter

	//Calculate field size at filter bottom
	const std::pair<double, double> inherentFilterFieldSizes =
	  fieldSize(focalSpot,
		    detectorDx,
		    detectorDy,
		    source2det,
		    source2inherentFilter + inherentFilterSize);

	//Create filter
	vector3D<double> inherentFilterCenter = sourcePos;
	inherentFilterCenter.z -= source2inherentFilter + inherentFilterSize/2.0;
	createBaseFilter(1.2*inherentFilterFieldSizes.first,
			 1.2*inherentFilterFieldSizes.second,
			 inherentFilterSize,
			 1,
			 out,
			 nextMat++,
			 "inherent-filter",
			 "world",
			 false,
			 inherentFilterCenter);

	// ** Inherent filter collimator

	// Create a collimator before the filters

	//Calculate field size at top
	const std::pair<double, double> filterCollTopSizes =
	  fieldSize(focalSpot,
		    detectorDx,
		    detectorDy,
		    source2det,
		    source2filterCollTop);

	//Calculate field size at bot
	const std::pair<double, double> filterCollBotSizes =
	  fieldSize(focalSpot,
		    detectorDx,
		    detectorDy,
		    source2det,
		    source2filterCollBot);

	//Create collimator
	vector3D<double> filterCollCenter = sourcePos;
	filterCollCenter.z -= source2filterCollTop + collHeight/2.0;
	createBaseCollimator(std::min(2.0*filterCollBotSizes.first, detectorDx),
			     filterCollTopSizes.first,
			     filterCollBotSizes.first,
			     std::min(2.0*filterCollBotSizes.second, detectorDy),
			     filterCollTopSizes.second,
			     filterCollBotSizes.second,
			     collHeight,
			     out,
			     collMat,
			     "filter-collimator",
			     "world",
			     false,
			     filterCollCenter);

      }

      // ** Filters
      //

      //Calculate field size at filters bottom
      const double filtersEnd = source2filter + filtersWidth;
      const std::pair<double, double> filterFieldSizes =
	fieldSize(focalSpot,
		  detectorDx,
		  detectorDy,
		  source2det,
		  filtersEnd);

      if(filters.size() > 0){
      
	size_t ifilter = 0;
	double zorigin = source2filter;
	for(const double& f : filters){

	  //Create ith-filter
	  std::string filterName("filter-");
	  filterName += std::to_string(ifilter);
	  vector3D<double> filterCenter = sourcePos;
	  filterCenter.z -= zorigin + f/2.0;
	  createBaseFilter(1.2*filterFieldSizes.first,
			   1.2*filterFieldSizes.second,
			   f,
			   1,
			   out,
			   nextMat++,
			   filterName,
			   "world",
			   false,
			   filterCenter);	

	  zorigin += f;
	  ++ifilter;
	}

	// ** Filters collimation

	// Create a collimator before the irradiated object and detector

	//Calculate field size at top
	const std::pair<double, double> detCollTopSizes =
	  fieldSize(focalSpot,
		    detectorDx,
		    detectorDy,
		    source2det,
		    source2detCollTop);

	//Calculate field size at bot
	const std::pair<double, double> detCollBotSizes =
	  fieldSize(focalSpot,
		    detectorDx,
		    detectorDy,
		    source2det,
		    source2detCollBot);

	//Create collimator
	vector3D<double> detCollCenter = sourcePos;
	detCollCenter.z -= source2detCollTop + collHeight/2.0;
	createBaseCollimator(std::min(2.0*detCollBotSizes.first, detectorDx),
			     detCollTopSizes.first,
			     detCollBotSizes.first,
			     std::min(2.0*detCollBotSizes.second, detectorDy),
			     detCollTopSizes.second,
			     detCollBotSizes.second,
			     collHeight,
			     out,
			     collMat,
			     "detector-collimator",
			     "world",
			     false,
			     detCollCenter);

      }

      // ** Bowtie
      //

      if(source2bowtie > 0.0 && bowtieDz.size() > 0){

	if(source2bowtie <= filtersEnd){
	  if(verbose > 1){
	    printf("constructDevice: Error: The distance between source"
		   "and bowtie filter is lesser than the distance between the "
		   "source and the collimator after flat filters\n"
		   "    Distance source to bowtie         : %f cm\n"
		   "    Distance source to last collimator: %f cm\n",
		   source2bowtie,
		   filtersEnd);
	  }
	  return errors::INVALID_DISTANCE;
	}
	
	//Get the maximum bowtie height
	const double maxDz = *std::max_element(bowtieDz.cbegin(), bowtieDz.cend());
	const double source2BowtieBot = source2bowtie + maxDz;

	//Check if the maximum dz is positive
	if(maxDz <= 0.0){
	  if(verbose > 1){
	    printf("constructDevice: Error: The maximum bowtie height must be greater than 0.\n"
		   "                        Bowtie maximum height: %f cm\n",
		   maxDz);
	  }
	  return errors::INVALID_DISTANCE;
	}

	//Calculate field size at bowtie bot face
	const std::pair<double, double> bowtieCollBotSizes =
	  fieldSize(focalSpot,
		    detectorDx,
		    detectorDy,
		    source2det,
		    source2BowtieBot);

	vector3D<double> bowtieCenter = sourcePos;
	bowtieCenter.z -= source2bowtie + maxDz/2.0;
	createTopFaceIrregularFilter(bowtieCollBotSizes.first,
				     bowtieCollBotSizes.second,
				     bowtieDz,
				     out,
				     nextMat++,
				     "bowtie",
				     "world",
				     false,
				     bowtieCenter);	
      }

      // ** Detector
      //

      //Create final detector
      vector3D<double> detCenter = sourcePos;
      detCenter.z -= source2det + 0.5;
      createBaseFilter(detectorDx,
		       detectorDy,
		       1.0,
		       1,
		       out,
		       nextMat++,
		       "detector",
		       "world",
		       false,
		       detCenter);

      // ** World
      //

      //Create world body
      vector3D<double> worldCenter = sourcePos;
      worldCenter.z -= source2det/2.0;
      createBaseFilter(2.0*detectorDx,
		       2.0*detectorDy,
		       1.2*source2det,
		       1,
		       out,
		       0,
		       "world",
		       "void",
		       false,
		       worldCenter);
      
      return 0;
    }
    
    int simDevice(const pen_parserSection& config,
		  const unsigned verbose){
      
      // ** Parse configuration
      
      //Read information from config section
      readerXRayDeviceSimulate reader;
      int err = reader.read(config,verbose);
      if(err != readerXRayDeviceSimulate::SUCCESS){
	return err;
      }

      //Save added geometry configuration
      config.readSubsection("geometry/config", reader.addedGeoConf);

      // ** Comon simulation configuration
      penred::simulation::simConfig baseSimConfig;
      err = baseSimConfig.configure("simulation",config);
      if(err != penred::simulation::errors::SUCCESS){
	if(verbose > 0)
	  printf("simDevice: Error: Unable to parse 'simulation' "
		 "section: Invalid seeds\n");
	return -1;
      }
      //Set verbose level
      baseSimConfig.verbose = verbose;      
      if(verbose > 1){
	printf("%s\n", baseSimConfig.stringifyConfig().c_str());
      }

      // ** Sampling function
      
      measurements::results<double, 2> spatialDistrib;
      measurements::results<double, 1> energyDistrib;

      sampling::aliasing<2> spatialSampler;
      sampling::aliasing<1> energySampler;

      simulation::sampleFuncType<pen_particleState> fsample;
      double maxE;

      const vector3D<double> sourcePos = reader.sourcePosition;      

      //Calculate beam radius
      constexpr double pi = 3.141592653589793;
      constexpr double pi05 = pi/2.0;
      constexpr double pi2 = 2.0*pi;
      constexpr double deg2rad = pi/180.0;
      
      const double anodeAngleRad = deg2rad*reader.anodeAngle;
      const double beamRad = reader.focalSpot*tan(pi05-anodeAngleRad);
      const double tanAnodeAngle = tan(anodeAngleRad);
      
      if(reader.simAnode){
	
	const double beamE = reader.kvp*1.0e3;
	maxE = beamE;

	//It is a primary (source) particle, set ILB[0] = 1
	pen_particleState baseState;
	baseState.X = sourcePos.x;
	baseState.Y = sourcePos.y;
	baseState.Z = sourcePos.z;

	baseState.U = 0.0;
	baseState.V = -1.0;
	baseState.W = 0.0;

	baseState.E = beamE;
	
	//Flag it as primary particle
	baseState.ILB[0] = 1;
	
	fsample = [beamRad, baseState]
	  (pen_particleState& state,    //Generated state
	   pen_KPAR& kpar,          //Generated kpar
	   unsigned long long& dh,  //History increment
	   const unsigned,          //Thread number
	   pen_rand& random) -> void{

	  dh = 1;
	  kpar = PEN_ELECTRON;
	  
	  //Sample particle position in the beam radius
	  double r = beamRad * sqrt(random.rand());
	  double theta = random.rand() * pi2; 

	  state = baseState;
	  
	  state.X += r * cos(theta);
	  state.Z += r * sin(theta);
	};
      }
      else{

	//Read energy distribution
	std::ifstream fin(reader.energyDistribFile, std::ifstream::in);
	if(!fin){
	  printf("Unable to open energy distribution file '%s'\n",
		 reader.energyDistribFile.c_str());
	  return -1;
	}
  
	err = energyDistrib.read(fin);
	if(err != 0){
	  if(verbose > 0){
	    printf("Error reading energy distribution data file '%s'.\n"
		   "  Error code: %d\n"
		   "  Error message: %s\n",
		   reader.energyDistribFile.c_str(),
		   err,
		   penred::measurements::errorToString(err));
	  }
	  fin.close();
	  return -2;
	}
	fin.close();

	if(verbose > 1){
	  printf("Configuring samplers...\n");
	}

	//Configure energy sampler
	err = energySampler.init(energyDistrib.readData(),
				 energyDistrib.readDimBins(),
				 energyDistrib.readLimits());
	if(err != 0){
	  if(verbose > 0){
	    printf("Error on energy sampler initialization from file '%s'.\n"
		   "  Error code: %d\n"
		   "  Please, report this error\n",
		   reader.energyDistribFile.c_str(),
		   err);
	  }
	  return -3;
	}

	//Read spatial distribution
	fin.open(reader.spatialDistribFile, std::ifstream::in);
	if(!fin){
	  printf("Unable to open spatial distribution file '%s'\n",
		 reader.spatialDistribFile.c_str());
	  return -1;
	}
  
	err = spatialDistrib.read(fin);
	if(err != 0){
	  printf("Error reading spatial distribution data file '%s'.\n"
		 "  Error code: %d\n"
		 "  Error message: %s\n",
		 reader.spatialDistribFile.c_str(),
		 err,
		 penred::measurements::errorToString(err));
	  fin.close();
	  return -2;
	}
	fin.close();

	//Configure spatial sampler
	err = spatialSampler.init(spatialDistrib.readData(),
				  spatialDistrib.readDimBins(),
				  spatialDistrib.readLimits());
	if(err != 0){
	  if(verbose > 0){
	    printf("Error on spatial sampler initialization from file '%s'.\n"
		   "  Error code: %d\n"
		   "  Please, report this error\n",
		   reader.energyDistribFile.c_str(),
		   err);
	  }
	  return -3;
	}

	if(verbose > 1){
	  printf("Samplers configured!\n");
	}
	
	const double distrib2source = reader.distrib2source;
	maxE = energyDistrib.readLimits()[0].second+50;
	
	fsample = [tanAnodeAngle, pi2, beamRad, sourcePos, distrib2source, &spatialSampler, &energySampler]
	  (pen_particleState& state,    //Generated state
	   pen_KPAR& kpar,          //Generated kpar
	   unsigned long long& dh,  //History increment
	   const unsigned,          //Thread number
	   pen_rand& random) -> void{

	  dh = 1;
	  kpar = PEN_PHOTON;

	  //Sample a position in the beam radius
	  double r = beamRad * sqrt(random.rand());
	  double theta = random.rand() * pi2;

	  const double dx = r * cos(theta);
	  const double dz = r * sin(theta);

	  //Calculate the displacement in y axis
	  //_________                        n +Z
	  //   | |  /                        |
	  //   | | /____________beam center  |--> +Y
	  //   | |/_|dz
	  //   |a/dy
	  //___|/
	  //
	  //      tan(a) = dy/dz -> dy = tan(a)*dz
	  //
	  const double dy = tanAnodeAngle*dz; //Notice the dz sign

	  //Sample particle position
	  std::array<double, 2> posXY = spatialSampler.samplePositions(random);

	  //shift the distribution
	  posXY[0] += dx;
	  posXY[1] += dy;

	  state.reset();
	  state.ILB[0] = 1;
	  
	  state.X = posXY[0] + sourcePos.x;
	  state.Y = posXY[1] + sourcePos.y;
	  state.Z =  sourcePos.z - distrib2source + dz;

	  //Sample energy
	  const std::array<double, 1> energy = energySampler.samplePositions(random);

	  state.E = energy[0];

	  //Set direction
	  vector3D<double> dir(posXY[0], posXY[1], -distrib2source + dz);
	  dir.normalize();
	  state.U = dir.x;
	  state.V = dir.y;
	  state.W = dir.z;

	};
	
      }

      auto itMinMax = std::minmax_element(reader.bowtieDz.cbegin(), reader.bowtieDz.cend());
      const double bowtieMin = *itMinMax.first;
      const double bowtieMax = *itMinMax.second;
      const unsigned nSpatBinsX = reader.detBinsX;
      const unsigned nSpatBinsY = reader.detBinsY;
      const double tolerance = reader.tolerance;

      if(reader.bowtieAutoDesign){
	if(reader.bowtieDesignBins > reader.bowtieDz.size()){
	  reader.bowtieDz.resize(reader.bowtieDesignBins);
	  std::fill(reader.bowtieDz.begin(), reader.bowtieDz.end(), bowtieMin);
	  reader.bowtieDz[0] = bowtieMax;
	  reader.bowtieDz.back() = bowtieMax;
	}
	reader.detBinsX = reader.bowtieDz.size();
	reader.detBinsY = 1;

	//Set initial tolerance to 10%
	reader.tolerance = 0.1;

	//Set verbose to 1
	baseSimConfig.verbose = 1;
      }

      
      if(reader.bowtieAutoDesign){	
	double midVal;
	const double maxChangeFactor = 1.0;
	const double minChangeFactor = 0.2;
	double changeFactor = minChangeFactor;

	const double maxTolerance = 0.1;
	const double minTolerance = tolerance/10.0;

	bool reached = false;
	for(unsigned i = 0; i < reader.bowtieDesignIterations; ++i){

	  if(verbose > 1){
	    printf("\n+ Bowtie design iteration %u:\n"
		   "   - Tolerance: %.2f %%\n",
		   i, reader.tolerance*100.0);
	  }
	  
	  measurements::measurement<double, 2> detFluence;
	  measurements::measurement<double, 2> detEdep;
	  measurements::measurement<double, 1> detSpec;
	  unsigned long long simHists;
	
	  err = simDevice(reader, maxE, fsample, baseSimConfig,
			  detFluence, detEdep, detSpec, simHists,
			  1);
	  if(err != errors::SUCCESS){
	    return err;
	  }

	  std::string sufix("_");
	  sufix += std::to_string(i);
	
	  //Generate results
	  penred::measurements::results<double, 2> fluenceResults;
	  detFluence.results(simHists, fluenceResults);
	  
	  //Generate fluence profile
	  penred::measurements::results<double, 1> profile;

	  err = fluenceResults.profile1D(0, profile);
	  if(err != 0){
	    printf("Unexpected error profiling data.\n"
		   "  Error code: %d\n"
		   "  Error message: %s\n",
		   err,
		   penred::measurements::errorToString(err));
	    return err;
	  }

	  FILE* fout = nullptr;
	  std::string filename = reader.outputPrefix + "fluenceProfile.dat" + sufix;
	  fout = fopen(filename.c_str(), "w");
	  profile.print(fout, 2, true, false);
	  fclose(fout);

	  fout = nullptr;
	  filename = reader.outputPrefix + "bowtie.dat" + sufix;
	  fout = fopen(filename.c_str(), "w");
	  for(size_t j = 0; j < reader.bowtieDz.size(); ++j){
	    fprintf(fout, "%E\n", reader.bowtieDz[j]);
	  }
	  fclose(fout);
	      
	  //Get fluence relative differences
	  std::vector<double> relDiff = profile.readData();

	  //At the first iteration, calculate the medium value at the detector center
	  if(i == 0){
	    //Get central value
	    if(relDiff.size() % 2 == 0){
	      midVal = (relDiff[relDiff.size()/2] + relDiff[relDiff.size()/2+1])/2.0;
	    }else{
	      midVal = (relDiff[relDiff.size()/2] +
			relDiff[relDiff.size()/2+1] +
			relDiff[relDiff.size()/2-1])/3.0;	    
	    }	  
	  }

	  double maxRelDiff = 0.0;
	  for(double& val : relDiff){
	    val = (val-midVal)/midVal;
	    if(val > maxRelDiff)
	      maxRelDiff = val;
	  }

	  if(verbose > 1){
	    printf("   - Maximum relative difference: %.3f %%\n",
		   maxRelDiff*100.0);
	  }

	  //Update change bowtie reshape factor
	  changeFactor = maxRelDiff/2.0;
	  if(changeFactor > maxChangeFactor)
	    changeFactor = maxChangeFactor;
	  if(changeFactor < minChangeFactor)
	    changeFactor = minChangeFactor;

	  //Update tolerance
	  reader.tolerance = maxRelDiff/10.0;
	  if(reader.tolerance > maxTolerance)
	    reader.tolerance = maxTolerance;
	  if(reader.tolerance < minTolerance)
	    reader.tolerance = minTolerance;

	  fout = nullptr;
	  filename = reader.outputPrefix + "fluenceRelativeDiff.dat" + sufix;
	  fout = fopen(filename.c_str(), "w");
	  for(size_t j = 0; j < relDiff.size(); ++j){
	    fprintf(fout, "%E\n", relDiff[j]);
	  }
	  fclose(fout);

	  if(maxRelDiff < tolerance){
	    //Required tolerance reached

	    if(reached){
	      //Smooth already done
	      break;
	    }
	    
	    if(verbose > 1){
	      printf("   - Tolerance reached, smooth the bowtie and simulate again\n");
	    }

	    //Flag tolerance as reached and smooth the bowtie
	    reached = true;
	    //Smooth the bowtie
	    for(size_t k = 0; k < 10; ++k){
	      std::vector<double> bowtieDz = reader.bowtieDz;
	      for(size_t j = 1; j < relDiff.size()-1; ++j){
		reader.bowtieDz[j] =
		  bowtieDz[j]*0.7 + bowtieDz[j-1]*0.15 + bowtieDz[j+1]*0.15;
	      }
	    }
	    continue;
	  }
	  reached = false;

	  if(verbose > 1){
	    printf("   - Maximum bowtie change factor: %.3f %%\n",
		   changeFactor*100.0);
	  }

	  //Reshape bowtie
	  for(size_t j = 0; j < relDiff.size(); ++j){
	    double factor = relDiff[j];
	    if(std::fabs(factor) > changeFactor){
	      if(std::signbit(factor))
		reader.bowtieDz[j] -= reader.bowtieDz[j]*changeFactor;
	      else
		reader.bowtieDz[j] += reader.bowtieDz[j]*changeFactor;
	    }
	    else{
	      reader.bowtieDz[j] += reader.bowtieDz[j]*factor;
	    }
	      
	      
	    if(reader.bowtieDz[j] > bowtieMax)
	      reader.bowtieDz[j] = bowtieMax;
	    if(reader.bowtieDz[j] < bowtieMin)
	      reader.bowtieDz[j] = bowtieMin;
	  }

	  //Smooth bowtie with a smooth filter
	  for(size_t k = 0; k < 5; ++k){
	    std::vector<double> bowtieDz = reader.bowtieDz;
	    for(size_t j = 1; j < relDiff.size()-1; ++j){
	      reader.bowtieDz[j] =
		bowtieDz[j]*0.7 + bowtieDz[j-1]*0.15 + bowtieDz[j+1]*0.15;
	    }
	  }
	}
      }

      //Perform the final simulation
      reader.detBinsX = nSpatBinsX;
      reader.detBinsY = nSpatBinsY;
      reader.tolerance = tolerance;
      baseSimConfig.verbose = verbose;
      
      measurements::measurement<double, 2> detFluence;
      measurements::measurement<double, 2> detEdep;
      measurements::measurement<double, 1> detSpec;
      unsigned long long simHists;
	
      err = simDevice(reader, maxE, fsample, baseSimConfig,
		      detFluence, detEdep, detSpec, simHists,
		      verbose);
      if(err != errors::SUCCESS){
	return err;
      }

      // ** Print results
      FILE* fout = nullptr;
      std::string filename = reader.outputPrefix + "detectedFluence.dat";
      fout = fopen(filename.c_str(), "w");
      detFluence.print(fout, simHists, 2, true, false);
      fclose(fout);

      fout = nullptr;
      filename = reader.outputPrefix + "detectedEdep.dat";
      fout = fopen(filename.c_str(), "w");
      detEdep.print(fout, simHists, 2, true, false);
      fclose(fout);

      fout = nullptr;
      filename = reader.outputPrefix + "detectedSpectrum.dat";
      fout = fopen(filename.c_str(), "w");
      detSpec.print(fout, simHists, 2, true, false);
      fclose(fout);      

      return errors::SUCCESS;
      
    }

    int simDevice(const readerXRayDeviceSimulate& reader,
		  const double maxE,
		  const simulation::sampleFuncType<pen_particleState>& fsample,
		  const penred::simulation::simConfig& baseSimConfig,
		  measurements::measurement<double, 2>& detFluence,
		  measurements::measurement<double, 2>& detEdep,
		  measurements::measurement<double, 1>& detSpec,
		  unsigned long long& simHistsOut,
		  const unsigned verbose){
      
      // ** Threads number

      //Get the number of threads to use
      unsigned nThreads = reader.nThreads;
#ifdef _PEN_USE_THREADS_
      if(nThreads == 0){
	nThreads = std::max(static_cast<unsigned int>(2),
			    std::thread::hardware_concurrency());
      }
#else
      if(nThreads != 0){
	if(verbose > 1)
	  printf("simDevice: Warning: Number of threads has been specified,"
		 " but the code has been compiled with no multithreading support.\n"
		 " Only one thread will be used.\n");
      }
      nThreads = 1;
#endif

      // ** Tally function

      const vector3D<double> sourcePos = reader.sourcePosition;      
      
      std::vector<measurements::measurement<double, 2>> detectedFluence(nThreads);
      for(size_t i = 0; i < nThreads; ++i){
	detectedFluence[i].
	  init({reader.detBinsX, reader.detBinsY},
	       {std::pair<double,double>(sourcePos.x-reader.detectorDx/2.0,
					 sourcePos.x+reader.detectorDx/2.0),
		std::pair<double,double>(sourcePos.y-reader.detectorDy/2.0,
					 sourcePos.y+reader.detectorDy/2.0)});
	
	detectedFluence[i].setDimHeader(0, "X (cm)");
	detectedFluence[i].setDimHeader(1, "Y (cm)");
	detectedFluence[i].setValueHeader("Value (prob)");	
	
      }

      std::vector<measurements::measurement<double, 2>> detectedEdep(nThreads);
      for(size_t i = 0; i < nThreads; ++i){
	detectedEdep[i].
	  init({reader.detBinsX, reader.detBinsY},
	       {std::pair<double,double>(sourcePos.x-reader.detectorDx/2.0,
					 sourcePos.x+reader.detectorDx/2.0),
		std::pair<double,double>(sourcePos.y-reader.detectorDy/2.0,
					 sourcePos.y+reader.detectorDy/2.0)});

	detectedEdep[i].setDimHeader(0, "X (cm)");
	detectedEdep[i].setDimHeader(1, "Y (cm)");
	detectedEdep[i].setValueHeader("Edep (eV)");	
      }

      std::vector<measurements::measurement<double, 1>> detectedSpectrum(nThreads);
      for(size_t i = 0; i < nThreads; ++i){
	detectedSpectrum[i].
	  init({reader.eBins},
	       {std::pair<double,double>(reader.minEnergy,maxE+50)});
	detectedSpectrum[i].setDimHeader(0, "Energy (eV)");
	detectedSpectrum[i].setValueHeader("Value (prob)");
      }

      std::vector<simulation::tallyFuncType> ftallies(nThreads);
	
      for(size_t i = 0; i < nThreads; ++i){
	ftallies[i] = [&detectedFluence, &detectedEdep, &detectedSpectrum, i]
	  (const pen_particleState& state,  //Resulting state
	   const pen_KPAR kpar,  //State kpar
	   const unsigned long long& hist, //History number
	   const int kdet){        //Returned value by 'simulatePartCond' function

	  if(kdet == 1 && kpar == PEN_PHOTON){
	    detectedFluence[i].add({state.X, state.Y}, state.WGHT, hist);
	    detectedEdep[i].add({state.X, state.Y}, state.E*state.WGHT, hist);
	    detectedSpectrum[i].add({state.E}, state.WGHT, hist);
	  }
	  
	};
      
      }

      // ** Simulation configuration for each thread
      
      std::vector<penred::simulation::simConfig> simConfigs(nThreads);

      //Copy basic configuration
      const double tolerance = reader.tolerance*sqrt(static_cast<double>(nThreads))*0.8;
      for(unsigned i = 0; i < nThreads; i++){
	simConfigs[i].iThread = i;
	simConfigs[i].copyCommonConfig(baseSimConfig);
	//Set, by default, std::cout as output stream
	simConfigs[i].setOutstream(std::cout);
	//Set seed pair
	simConfigs[i].setSeeds(reader.seedPair+i);
	simConfigs[i].fSimFinish =
	  [tolerance, i, &detectedFluence]
	  (const unsigned long long ihist){
	    
	    if(ihist % 10000 == 0){
	      //Get mean relative error
	      const double erel = detectedFluence[i].errorRel(ihist);
	      if(erel < tolerance)
		return false;
	    }
	    return true;
	  };
      }
      
      // ** Device Geometry

      //Create device geometry
      std::stringstream geoStream;
      int err = constructDevice(geoStream,
				reader.focalSpot,
				reader.source2det,
				reader.source2filter,
				reader.source2bowtie,
				reader.detectorDx,
				reader.detectorDy,
				reader.inherentFilterWidth,
				reader.filtersWidth,
				reader.bowtieDz,
				reader.sourcePosition,
				reader.simAnode,
				reader.anodeAngle,
				verbose);

      if(err != 0){
	if(verbose > 0)
	  printf("simDevice: Error: Unable to create device geometry\n");
	return -1;
      }

      // ** Geometry configuration

      pen_parserSection geoConfig;
      std::shared_ptr<wrapper_geometry> geometry;

      if(reader.addedGeoType.compare("-") == 0){
	//Use Mesh geometry

	//Create a mesh geometry instance
	std::shared_ptr<pen_meshBodyGeo> geometryMesh =
	  std::make_shared<pen_meshBodyGeo>();

	//Set geometry file in geometry instance
	geometryMesh->preloadGeo.assign(geoStream.str());

	if(reader.printGeo){
	  //Print geometry file
	  std::ofstream out("device.msh", std::ofstream::out);
	  out << geoStream.str() << std::endl;
	  out.close();
	}
	
	geometry = geometryMesh;
	if(!geometry){
	  printf("Unexpected Error: Unable to create a geometry instance "
		 "of type 'MESH_BODY'\n"
		 "                  Please, report this error");
	  return -2;
	}

	//Set in-memory geometry
	geoConfig.set("memory-file", true);

	//Assign kdet 1 to detector object
	geoConfig.set("kdet/detector", 1);

	//Limit anode dsmax if it is simulated
	if(reader.simAnode)
	  geoConfig.set("dsmax/anode", 2.0e-2);
	
      }
      else{

	//Create combo geometry
	std::shared_ptr<pen_comboGeo> geometryCombo =
	  std::make_shared<pen_comboGeo>();

	//Print geometry file
	std::ofstream out("device.msh", std::ofstream::out);
	out << geoStream.str() << std::endl;
	out.close();	
	  
	geometry = geometryCombo;
	if(!geometry){
	  printf("Unexpected Error: Unable to create a geometry instance "
		 "of type 'COMBO'\n"
		 "                  Please, report this error");
	  return -2;
	}
  
	// * Configure device geometry
	
	geoConfig.set("geometries/device/priority", 0);
	geoConfig.set("geometries/device/config/type", "MESH_BODY");
	geoConfig.set("geometries/device/config/input-file", "device.msh");

	//Assign kdet 1 to detector object
	geoConfig.set("geometries/device/config/kdet/detector", 1);

	//Limit anode dsmax if it is simulated
	if(reader.simAnode)
	  geoConfig.set("geometries/device/config/dsmax/anode", 2.0e-2);

	// * Configure added geometry

	//Read configuration section
	if(reader.addedGeoConf.size() == 0){
	  if(verbose > 0){
	    printf("simDevice: Error: Configuration for added "
		   "geometry not provided.\n"
		   "                  Missing section at: 'geometry/config'\n");
	  }
	  return -2;
	}

	geoConfig.addSubsection("geometries/added/config", reader.addedGeoConf);
	geoConfig.set("geometries/added/priority", 1);
	geoConfig.set("geometries/added/config/type", reader.addedGeoType);	
      }

      //Clear geometry stream
      geoStream.clear();      
      
      // ** Context configuration

      pen_parserSection simConf;
      int nextMat = 1;
      
      // ** Collimators (perfect absorber)

      //Create collimators material file
      std::string errorString;
      err = penred::penMaterialCreator::createMat(82,
						  "collimator.mat",
						  errorString);

      if(err != 0){
	printf ("simDevice: Error: Unable to create "
		"collimator material: %s\n", errorString.c_str());
	printf ("IRETRN =%d\n", err);
	return -2;
      }
	
      simConf.set("materials/collimators/number", nextMat++);
      for(unsigned i = 0; i < constants::nParTypes; ++i){
	std::string path("materials/collimators/eabs/");
	path += particleName(i);
	simConf.set(path, 1.0e35);
      }
      simConf.set("materials/collimators/filename", "collimator.mat");

      // ** Anode
      
      if(reader.simAnode){

	//Create anode material file
	err = penred::penMaterialCreator::createMat(reader.anodeZ,
						    "anode.mat",
						    errorString);
	if(err != 0){
	  printf ("simDevice: Error: Unable to create "
		  "anode material: %s\n", errorString.c_str());
	  printf ("IRETRN =%d\n", err);
	  return -2;
	}
	
	//Configure anode material
	simConf.set("materials/anode/number", nextMat++);
	for(unsigned i = 0; i < constants::nParTypes; ++i){
	  std::string path("materials/anode/eabs/");
	  path += particleName(i);
	  simConf.set(path, reader.minEnergy);
	}
	
	simConf.set("materials/anode/C1", 0.05);
	simConf.set("materials/anode/C2", 0.05);
	simConf.set("materials/anode/WCC", std::min(5e3,reader.minEnergy/100.0));
	simConf.set("materials/anode/WCR", std::min(5e3,reader.minEnergy/100.0));
	simConf.set("materials/anode/filename", "anode.mat");

	//VR
	simConf.set("VR/IForcing/bremss/particle", "electron");
	simConf.set("VR/IForcing/bremss/interaction", BETAe_HARD_BREMSSTRAHLUNG);
	simConf.set("VR/IForcing/bremss/factor", 400);
	simConf.set("VR/IForcing/bremss/min-weight", 0.1);
	simConf.set("VR/IForcing/bremss/max-weight", 2.0);
	simConf.set("VR/IForcing/bremss/bodies/anode", true);

	simConf.set("VR/IForcing/innerShell/particle", "electron");
	simConf.set("VR/IForcing/innerShell/interaction", BETAe_HARD_INNER_SHELL);
	simConf.set("VR/IForcing/innerShell/factor", 400);
	simConf.set("VR/IForcing/innerShell/min-weight", 0.1);
	simConf.set("VR/IForcing/innerShell/max-weight", 2.0);
	simConf.set("VR/IForcing/innerShell/bodies/anode", true);

	simConf.set("VR/bremss/split4/splitting", 4);
	simConf.set("VR/bremss/split4/bodies/anode", true);      	
      }

      // ** Inherent filter

      if(reader.inherentFilterWidth > 0.0){

	//Create inherent filter material
	err = penred::penMaterialCreator::createMat(13,
						    "inherent.mat",
						    errorString);
	if(err != 0){
	  printf ("simDevice: Error: Unable to create "
		  "inherent filter material: %s\n",
		  errorString.c_str());
	  printf ("IRETRN =%d\n", err);
	  return -2;
	}

	//Configure inherent filter material
	simConf.set("materials/inherentFilter/number", nextMat++);
	for(unsigned j = 0; j < constants::nParTypes; ++j){
	  std::string path = "materials/inherentFilter/eabs/";
	  path += particleName(j);
	  if(j == PEN_PHOTON)
	    simConf.set(path, reader.minEnergy);
	  else
	    simConf.set(path, 1.0e35);
	}
	simConf.set("materials/inherentFilter/filename", "inherent.mat");

      }

      // ** Filters
      for(size_t i = 0; i < reader.filtersZ.size(); ++i){

	//Create filter material file
	std::string filterName("filter_");
	filterName += std::to_string(i);
	
	std::string filterMatFile = filterName + ".mat";
	
	if(reader.filtersMatFile[i].compare("-") == 0){
	  err = penred::penMaterialCreator::createMat(reader.filtersZ[i],
						      filterMatFile.c_str(),
						      errorString);
	  if(err != 0){
	    printf ("simDevice: Error: Unable to create "
		    "filter %lu material: %s\n",
		    static_cast<unsigned long>(i),
		    errorString.c_str());
	    printf ("IRETRN =%d\n", err);
	    return -2;
	  }
	}else{
	  filterMatFile = reader.filtersMatFile[i];
	}
	
	//Configure material
	std::string prefix = "materials/" + filterName + "/";
	simConf.set((prefix + "number").c_str(), nextMat++);
	for(unsigned j = 0; j < constants::nParTypes; ++j){
	  std::string path = prefix + "eabs/";
	  path += particleName(j);
	  if(j == PEN_PHOTON)
	    simConf.set(path, reader.minEnergy);
	  else
	    simConf.set(path, 1.0e35);
	}

	//simConf.set((prefix + "C1").c_str(), 0.05);
	//simConf.set((prefix + "C2").c_str(), 0.05);
	//simConf.set((prefix + "WCC").c_str(), std::min(5e3,reader.minEnergy/100.0));
	//simConf.set((prefix + "WCR").c_str(), std::min(5e3,reader.minEnergy/100.0));
	simConf.set((prefix + "filename").c_str(), filterMatFile);
		
      }

      // ** Bowtie filter

      if(reader.source2bowtie > 0.0 && reader.bowtieDz.size() > 0){

	//Create bowtie material
	std::string matFilename;
	if(reader.bowtieMatFile.compare("-") == 0){
	  matFilename.assign("bowtie.mat");
	  err = penred::penMaterialCreator::createMat(reader.bowtieZ,
						      matFilename.c_str(),
						      errorString);
	  if(err != 0){
	    printf ("simDevice: Error: Unable to create "
		    "bowtie filter material: %s\n",
		    errorString.c_str());
	    printf ("IRETRN =%d\n", err);
	    return -2;
	  }
	}else{
	  matFilename = reader.bowtieMatFile;
	}

	//Configure inherent filter material
	simConf.set("materials/bowtie/number", nextMat++);
	for(unsigned j = 0; j < constants::nParTypes; ++j){
	  std::string path = "materials/bowtie/eabs/";
	  path += particleName(j);
	  if(j == PEN_PHOTON)
	    simConf.set(path, reader.minEnergy);
	  else
	    simConf.set(path, 1.0e35);
	}
	simConf.set("materials/bowtie/filename", matFilename);
      }      

      // ** Detector
      simConf.set("materials/detector/number", nextMat++);
      for(unsigned i = 0; i < constants::nParTypes; ++i){
	std::string path("materials/detector/eabs/");
	path += particleName(i);
	simConf.set(path, 1.0e35);
      }
      simConf.set("materials/detector/filename", "collimator.mat");

      // ** Added geometry

      if(reader.addedGeoType.compare("-") != 0){
	//Provided geometry
	for(const readerXRayDeviceSimulate::materialData& mat : reader.addedGeoMats){
	  if(mat.index < nextMat && mat.index != 0){
	    if(verbose > 0){
	      printf("simDevice: Error: Added geometry uses a non void material "
		     "reserved for device geometry (mat %d).\n"
		     "           Number of reserved materials: %d\n",
		     mat.index, nextMat);
	    }
	    return -3;
	  }

	  //Add configuration for this material
	  std::string prefix = "materials/added_" + mat.name + "/";
	  simConf.set((prefix + "number").c_str(), mat.index);
	  for(unsigned j = 0; j < constants::nParTypes; ++j){
	    std::string path = prefix + "eabs/";
	    path += particleName(j);
	    if(j == PEN_PHOTON)
	      simConf.set(path, reader.minEnergy);
	    else
	      simConf.set(path, 1.0e35);
	  }

	  for(const penred::massFraction& mf : mat.composition){
	    std::string prefixElement = prefix + "elements/" + std::to_string(mf.Z);
	    simConf.set(prefixElement, mf.fraction);
	  }

	  simConf.set((prefix + "density").c_str(), mat.density);
	  simConf.set((prefix + "filename").c_str(), "added_" + mat.name + ".mat");
	  simConf.set((prefix + "force-creation").c_str(), true);
	}
      }

      // ** Create context

      pen_context context;

      //Run context configuration step with no geometry
      pen_parserSection matInfo;
      if(context.configure(maxE,
			   simConf,
			   matInfo,
			   verbose > 2 ? verbose : 1) != pen_context::SUCCESS){
	printf("simDevice: Error at simulation context initialization. "
	       "See context report.\n");
	return -3;
      }

      //Configure geometry
      geoConfig.addSubsection("materials", matInfo);
      
      if(geometry->configure(geoConfig,
			     verbose > 2 ? verbose : 1) != PEN_MESHBODY_GEO_SUCCESS){
	if(verbose > 0){
	  printf("Unexpected Error: Unable to construct the geometry. "
		 "Please, report this error\n");
	}
	return -2;
      }      
      
      //Set the geometry to the simulation context
      context.setGeometry(geometry.get());      

      //Run context configuration step with geometry
      if(context.configureWithGeo(simConf,
				  verbose > 2 ? verbose : 1) != pen_context::SUCCESS){
	printf("simDevice: Error at simulation context initialization with geometry. "
	       "Please, report this error.\n");
	return -3;
      }

      if(reader.printGeo){
	//Print configuration also
	geoConfig.remove("materials");
	simConf.addSubsection("geometry", geoConfig);
	std::ofstream out("device.conf", std::ofstream::out);
	out << "## Simulation " << std::endl;
	out << "###############\n" << std::endl;
	out << "simulation/threads " << reader.nThreads << std::endl;
	out << "simulation/max-time " << reader.maxTime << std::endl;
	out << "\n" << std::endl;
	out << "## Source " << std::endl;
	out << "###############\n" << std::endl;
	out << "sources/generic/source1/nhist " << reader.nHists << std::endl;
	out << "sources/generic/source1/kpar \"gamma\"" << std::endl;
	out << "\n" << std::endl;
	out << "## Tallies " << std::endl;
	out << "###############\n" << std::endl;
	out << "tallies/FluenceDetectior/type \"DETECTION_SPATIAL_DISTRIB\" " << std::endl;
	out << "tallies/SpatialDetectior/spatial/nx " << reader.detBinsX << std::endl;
	out << "tallies/SpatialDetectior/spatial/xmin " <<
	  detectedFluence[0].readLimits()[0].first << std::endl;
	out << "tallies/SpatialDetectior/spatial/xmax " <<
	  detectedFluence[0].readLimits()[0].second << std::endl;
	out << "tallies/SpatialDetectior/spatial/ny " << reader.detBinsY << std::endl;
	out << "tallies/SpatialDetectior/spatial/ymin " <<
	  detectedFluence[0].readLimits()[1].first << std::endl;
	out << "tallies/SpatialDetectior/spatial/ymax  " <<
	  detectedFluence[0].readLimits()[1].second << std::endl;
	out << "tallies/SpatialDetectior/detector 1 " << std::endl;
	out << "tallies/SpatialDetectior/particle \"gamma\" " << std::endl;
	out << "\n" << std::endl;
	out << "tallies/SpectrumDetectior/type \"DETECTION_SPATIAL_DISTRIB\" " << std::endl;
	out << "tallies/SpectrumDetectior/spatial/nbins " << reader.eBins << std::endl;
	out << "tallies/SpectrumDetectior/spatial/emin " <<
	  detectedSpectrum[0].readLimits()[0].first << std::endl;
	out << "tallies/SpectrumDetectior/spatial/emax " <<
	  detectedSpectrum[0].readLimits()[0].second << std::endl;
	out << "tallies/SpectrumDetectior/detector 1 " << std::endl;
	out << "tallies/SpectrumDetectior/particle \"gamma\" " << std::endl;	
	out << "\n" << std::endl;
	out << "## Geometry " << std::endl;
	out << "###############\n" << std::endl;
	pen_parserSection geoConfOut;
	geoConfOut.addSubsection("geometry", geoConfig);
	if(reader.addedGeoType.compare("-") == 0){
	  geoConfOut.set("geometry/type", "MESH_BODY");
	}
	else{
	  geoConfOut.set("geometry/type", "COMBO");
	}
	out << geoConfOut.stringify() << std::endl;
	out << "\n" << std::endl;
	out << "## Materials " << std::endl;
	out << "###############\n" << std::endl;
	simConf.remove("geometry");
	out << simConf.stringify() << std::endl;
	out.close();
      }
      
      // ** Simulate

#ifdef _PEN_USE_THREADS_

      std::vector<std::thread> simThreads;

      std::string sourceName = "deviceSource";
      for(size_t i = 0; i < nThreads; ++i){
	unsigned long long hists = reader.nHists/nThreads;
	if(nThreads == 0)
	  hists += reader.nHists % nThreads;
	simThreads.emplace_back([i, hists, &simConfigs, &context, &sourceName,
				 &fsample, &ftallies]() -> void {
	  
	  simulation::sampleAndSimulateCondContext(simConfigs[i],
						   context,
						   hists,
						   sourceName,
						   fsample,
						   simulation::finishTypes::DETECTOR_REACHED,
						   1,
						   ftallies[i]);
	});
      }

      //Join threads
      for(size_t i = 0; i < nThreads; ++i){
	simThreads[i].join();
	if(i > 0){
	  detectedFluence[0].add(detectedFluence[i]);
	  detectedEdep[0].add(detectedEdep[i]);
	  detectedSpectrum[0].add(detectedSpectrum[i]);
	}
      }
#else
      simulation::sampleAndSimulateCondContext(simConfigs[0],
					       context,
					       reader.nHists,
					       sourceName,
					       fsample,
					       simulation::finishTypes::DETECTOR_REACHED,
					       1,
					       ftallies[0]);
#endif

      unsigned long long simHists = 0;
      for(size_t i = 0; i < nThreads; ++i){
	simHists += simConfigs[i].getTotalSimulated();
      }

      //Save results
      detFluence = detectedFluence[0];
      detEdep = detectedEdep[0];
      detSpec = detectedSpectrum[0];
      simHistsOut = simHists;
      
      return errors::SUCCESS;
      
    }
    

    // ** Device creation. Reader functions
    int readerXRayDeviceCreate::beginSectionFamily(const std::string& pathInSection,
						   const size_t,
						   const unsigned){

      if(family == -1){
	if(pathInSection.compare("filters") == 0){
	  family = 0;
	}else{
	  return errors::UNHANDLED;
	}
      }else{
	return errors::UNHANDLED;
      }

      return errors::SUCCESS;
    }

    int readerXRayDeviceCreate::endSectionFamily(const unsigned){
      if(family == 0){
	family = -1;
      }else{
	return errors::UNHANDLED;
      }
      return errors::SUCCESS;
    }

    int readerXRayDeviceCreate::beginSection(const std::string&,
					     const unsigned){
      if(family == 0){ //Filters
	return errors::SUCCESS;
      }
      return errors::UNHANDLED;  
    }

    int readerXRayDeviceCreate::endSection(const unsigned){

      if(family == 0){ //Filters
	return errors::SUCCESS;
      }
  
      return errors::UNHANDLED;  
    }

    int readerXRayDeviceCreate::storeElement(const std::string& pathInSection,
					     const pen_parserData& element,
					     const unsigned){

      if(family == -1){ //Root
	if(pathInSection.compare("anode/create") == 0){
	  createAnode = element;
	}
	else if(pathInSection.compare("anode/angle") == 0){
	  anodeAngle = element;
	}
	else if(pathInSection.compare("focalSpot") == 0){
	  focalSpot = element;
	}    
	else if(pathInSection.compare("source2det") == 0){
	  source2det = element;
	}
	else if(pathInSection.compare("source2bowtie") == 0){
	  source2bowtie = element;
	}	
	else if(pathInSection.compare("detector/dx") == 0){
	  detectorDx = element;
	}
	else if(pathInSection.compare("detector/dy") == 0){
	  detectorDy = element;
	}
	else if(pathInSection.compare("inherentFilter/width") == 0){
	  inherentFilterWidth = element;
	}
	else if(pathInSection.compare("source2filter") == 0){
	  source2filter = element;
	}
	else if(pathInSection.compare("source/pos/x") == 0){
	  sourcePosition.x = element;
	}    
	else if(pathInSection.compare("source/pos/y") == 0){
	  sourcePosition.y = element;
	}    
	else if(pathInSection.compare("source/pos/z") == 0){
	  sourcePosition.z = element;
	}
	else{
	  return errors::UNHANDLED;
	}
      }
      else if(family == 0){ //Filters
	if(pathInSection.compare("width") == 0){
	  double width = element;
	  filters.push_back(width);
	}
	else{
	  return errors::UNHANDLED;
	}
      }
      else{
	return errors::UNHANDLED;
      }

      return errors::SUCCESS;
  
}

    int readerXRayDeviceCreate::beginArray(const std::string& pathInSection,
					   const size_t,
					   const unsigned){
      if(family == -1){
	if(pathInSection.compare("bowtie/dz") == 0){
	  return errors::SUCCESS;
	}
      }
      return errors::UNHANDLED;
    }
      
    int readerXRayDeviceCreate::endArray(const unsigned){

      if(family == -1){
	return errors::SUCCESS;
      }
      return errors::UNHANDLED;
    }

    int readerXRayDeviceCreate::storeArrayElement(const std::string& pathInSection,
						  const pen_parserData& element,
						  const size_t,
						  const unsigned){
      if(family == -1){
	if(pathInSection.compare("bowtie/dz") == 0){
	  bowtieDz.push_back(static_cast<float>(element));
	  return errors::SUCCESS;
	}
      }
      return errors::UNHANDLED;
    }


    // ** Device simulation. Reader functions
    int readerXRayDeviceSimulate::beginSectionFamily(const std::string& pathInSection,
						     const size_t,
						     const unsigned){

      if(family == -1){
	//In root section
	if(pathInSection.compare("x-ray/filters") == 0){
	  //Enters to filters section
	  family = 0;
	}
	else if(pathInSection.compare("geometry/materials") == 0){
	  //Enters to materials section
	  family = 1;
	}
	else{
	  return errors::UNHANDLED;
	}
      }
      else if(family == 1){
	//In materials section
	if(pathInSection.compare("elements") == 0){
	  //Enters to elements section
	  family = 2;
	}
	else{
	  return errors::UNHANDLED;
	}	
      }
      else{
	return errors::UNHANDLED;
      }

      return errors::SUCCESS;
    }

    int readerXRayDeviceSimulate::endSectionFamily(const unsigned){
      if(family == 0 || family == 1){
	//Return from filters or materials to root section
	family = -1;
      }
      else if(family == 2){
	//Return from elements to materials section
	family = 1;
      }
      else{
	return errors::UNHANDLED;
      }
      return errors::SUCCESS;
    }

    int readerXRayDeviceSimulate::beginSection(const std::string& name,
					       const unsigned verbose){
      if(family == 0){ //Filters
	return errors::SUCCESS;
      }
      else if(family == 1){
	//Begin material
	addedGeoMats.emplace_back(name);
	return errors::SUCCESS;
      }
      else if(family == 2){
	//Begin filter
	//Try to convert the section name to an integer
	unsigned Z;
	if(sscanf(name.c_str(), "%u", &Z) != 1){
	  if(verbose > 0){
	    printf("Error: Unknown atomic number: %s\n", name.c_str());
	  }
	  return errors::INVALID_ATOMIC_NUMBER;      
	}

	addedGeoMats.back().composition.emplace_back(Z);	
	return errors::SUCCESS;
      }
      return errors::UNHANDLED;  
    }

    int readerXRayDeviceSimulate::endSection(const unsigned){

      if(family == 0 || family == 1 || family == 2){ //Filters or Materials or Elements
	return errors::SUCCESS;
      }
  
      return errors::UNHANDLED;  
    }

    int readerXRayDeviceSimulate::storeElement(const std::string& pathInSection,
					       const pen_parserData& element,
					       const unsigned){

      if(family == -1){ //Root
	if(pathInSection.compare("simulation/sim-anode") == 0){
	  simAnode = element;
	}
	else if(pathInSection.compare("simulation/histories") == 0){
	  nHists = element;
	}
	else if(pathInSection.compare("simulation/max-time") == 0){
	  maxTime = element;
	}
	else if(pathInSection.compare("simulation/min-energy") == 0){
	  minEnergy = element;
	}
	else if(pathInSection.compare("simulation/print-geometry") == 0){
	  printGeo = element;
	}
	else if(pathInSection.compare("simulation/nthreads") == 0){
	  nThreads = element;
	}
	else if(pathInSection.compare("simulation/seedPair") == 0){
	  seedPair = element;
	}	
	else if(pathInSection.compare("simulation/detBins/nx") == 0){
	  detBinsX = element;
	}
	else if(pathInSection.compare("simulation/detBins/ny") == 0){
	  detBinsY = element;
	}	
	else if(pathInSection.compare("simulation/eBins") == 0){
	  eBins = element;
	}
	else if(pathInSection.compare("simulation/tolerance") == 0){
	  tolerance = element;
	  tolerance /= 100.0;
	}
	else if(pathInSection.compare("x-ray/focal-spot") == 0){
	  focalSpot = element;
	}	
	else if(pathInSection.compare("x-ray/source/distrib2source") == 0){
	  distrib2source = element;
	}
	else if(pathInSection.compare("x-ray/source2det") == 0){
	  source2det = element;
	}
	else if(pathInSection.compare("x-ray/detector/dx") == 0){
	  detectorDx = element;
	}
	else if(pathInSection.compare("x-ray/detector/dy") == 0){
	  detectorDy = element;
	}
	else if(pathInSection.compare("x-ray/inherentFilter/width") == 0){
	  inherentFilterWidth = element;
	}
	else if(pathInSection.compare("x-ray/source2filter") == 0){
	  source2filter = element;
	}
	else if(pathInSection.compare("x-ray/source2bowtie") == 0){
	  source2bowtie = element;
	}
	else if(pathInSection.compare("x-ray/bowtie/z") == 0){
	  bowtieZ = element;
	}
	else if(pathInSection.compare("x-ray/bowtie/auto-design") == 0){
	  bowtieAutoDesign = element;
	}
	else if(pathInSection.compare("x-ray/bowtie/design-bins") == 0){
	  bowtieDesignBins = element;
	}
	else if(pathInSection.compare("x-ray/bowtie/design-iterations") == 0){
	  bowtieDesignIterations = element;
	}	
	else if(pathInSection.compare("anode/angle") == 0){
	  anodeAngle = element;
	}
	else if(pathInSection.compare("anode/z") == 0){
	  anodeZ = element;
	}	
	else if(pathInSection.compare("x-ray/kvp") == 0){
	  kvp = element;
	}
	else{
	  return errors::UNHANDLED;
	}
      }
      else if(family == 0){ //Filters
	if(pathInSection.compare("width") == 0){	  
	  double width = element;
	  filtersWidth.push_back(width);
	}
	else if(pathInSection.compare("z") == 0){	  
	  unsigned z = element;
	  filtersZ.push_back(z);
	}
	else{
	  return errors::UNHANDLED;
	}
      }
      else if(family == 1){ //Materials
	if(pathInSection.compare("density") == 0){
	  addedGeoMats.back().density = element;
	}
	else if(pathInSection.compare("number") == 0){
	  addedGeoMats.back().index = element;
	}	
	else{
	  return errors::UNHANDLED;
	}
      }
      else if(family == 2){ //Elements
	if(pathInSection.empty()){
	  //Read fraction by weight
	  addedGeoMats.back().composition.back().fraction = element;
	}
	else{
	  return errors::UNHANDLED;
	}
      }
      else{
	return errors::UNHANDLED;
      }

      return errors::SUCCESS;
  
    }

    int readerXRayDeviceSimulate::storeString(const std::string& pathInSection,
					      const std::string& element,
					      const unsigned){

      if(family == -1){ //Root
	if(pathInSection.compare("x-ray/source/spatial-distrib-file") == 0){	  
	  spatialDistribFile = element;
	}
	else if(pathInSection.compare("x-ray/source/energy-distrib-file") == 0){	  
	  energyDistribFile = element;
	}
	else if(pathInSection.compare("geometry/type") == 0){
	  addedGeoType = element;
	}
	else if(pathInSection.compare("x-ray/bowtie/mat-file") == 0){
	  bowtieMatFile = element;
	}
	else if(pathInSection.compare("simulation/output-prefix") == 0){
	  outputPrefix = element;
	}
	else{
	  return errors::UNHANDLED;
	}
      }
      else if(family == 0){ //Filters
	if(pathInSection.compare("mat-file") == 0){	  
	  filtersMatFile.push_back(element);
	}
	else{
	  return errors::UNHANDLED;
	}
      }
      else{
	return errors::UNHANDLED;
      }

      return errors::SUCCESS;
  
    }

    int readerXRayDeviceSimulate::beginArray(const std::string& pathInSection,
					     const size_t size,
					     const unsigned verbose){

      if(family == -1){
	if(pathInSection.compare("x-ray/source/position") == 0){
	  if(size != 3){
	    if(verbose > 0){
	      printf("Error: Bad position vector (x,y,z).\n"
		     "      Required coordinates: 3\n"
		     "      provided coordinates: %lu\n",
		     static_cast<unsigned long>(size));
	    }
	    return errors::BAD_DIMENSIONS;
	  }
	}
	else if(pathInSection.find("x-ray/bowtie/dz") == 0){
	  return errors::SUCCESS;
	}	
	else{
	  return errors::UNHANDLED;
	}
      }else{
	return errors::UNHANDLED;
      }
      return errors::SUCCESS;
    }

    int readerXRayDeviceSimulate::endArray(const unsigned){
      return errors::SUCCESS;
    }

    int readerXRayDeviceSimulate::storeArrayElement(const std::string& pathInSection,
						    const pen_parserData& element,
						    const size_t pos,
						    const unsigned){

      if(family == -1){
	if(pathInSection.compare("x-ray/source/position") == 0){
	  if(pos == 0){
	    sourcePosition.x = element;
	  }
	  else if(pos == 1){
	    sourcePosition.y = element;	    
	  }
	  else if(pos == 2){
	    sourcePosition.z = element;	    
	  }
	  else{
	    return errors::UNHANDLED;
	  }
	}
	else if(pathInSection.compare("x-ray/bowtie/dz") == 0){
	  bowtieDz.push_back(static_cast<float>(element));
	}
	else{
	  return errors::UNHANDLED;
	}	
      }
      else{
	return errors::UNHANDLED;
      }
      return errors::SUCCESS;
    }
    
    
  };
};
