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
 
#include "utilities.hh"


namespace penred{

  namespace xray{

    int filterWithCollimation(std::vector<detectedPart>& particlesIn,
			      const double filterZorigin,
			      const double filter2col,
			      const std::vector<std::pair<unsigned, double>>& filters,
			      const double coldz,
			      const double coldx1, const double coldy1,
			      const double coldx2, const double coldy2,
			      const double emin,
			      const unsigned nthreadsIn,
			      const unsigned verbose,
			      const std::string& geoFilename){

      //Ensure all distances are positive
      if(filter2col < 0.0 ||
	 coldz < 0.0 ||
	 coldx1 < 0.0 || coldy1 < 0.0 ||
	 coldx2 < 0.0 || coldy2 < 0.0){

	return NEGATIVE_DISTANCE;
      }

      if(filters.size() == 0){
	return NO_FILTER_PROVIDED;
      }

      if(particlesIn.size() == 0){
	//No particles to simulate
	return SUCCESS;
      }

      //Get maximum energy and total number of histories
      double maxE = 0.0;
      unsigned long long nHists = 1;
      for(const detectedPart& part : particlesIn){
	if(maxE < part.state.E)
	  maxE = part.state.E;
	nHists += part.dh;
      }

      //Check filter data
      double filtersEnd = filterZorigin;
      for(const auto& f : filters){
	if(f.first == 0)
	  return INVALID_Z;
	if(f.second <= 0.0)
	  return INVALID_SIZE;
	filtersEnd -= f.second;
      }
      filtersEnd -= 0.1;

      //Calculate filter size
      double filterSizeX = 2.0*std::max(coldx1, coldx2);
      double filterSizeY = 2.0*std::max(coldy1, coldy2);

      //Create a string to store the created geometries
      std::stringstream ss;

      //Create the void world containing the filtes and detector geometry
      ss << "# Number of objects:\n " << filters.size() + 2 << std::endl;

      // ** Create the world
      
      //Create a void world to contain all filters and the detector
      double worldSize = filterZorigin-filtersEnd + 2.0*coldz + 1.0;
      createBaseFilter(filterSizeX+0.1,  //dx
		       filterSizeY+0.1,  //dy
		       worldSize,  //dz
		       1,  //Number of subdivisions in X axis
		       ss, //Output stream
		       0,  //Material index
		       "world",  //"Filter" name
		       "void",   //Parent name
		       false,    //Do not print number of objects
		       vector3D<double>(0.0,0.0,filterZorigin+0.1-worldSize/2.0)); //Origin

      // ** Create all filters
      double forigin = filterZorigin;
      unsigned fcount = 0;
      for(const auto& f : filters){

	//Create filter material
	std::string filterName = "filter" + std::to_string(++fcount);
	std::string errorString;
	int err = penred::penMaterialCreator::createMat(f.first,
							filterName + ".mat",
							errorString);

	if(err != 0){
	  if(verbose > 0){
	    printf("Unable to create material for filter '%s'.\n"
		   "%s\n"
		   "IRETRN =%d\n",
		   filterName.c_str(),
		   errorString.c_str(),
		   err);
	  }
	  return UNABLE_TO_CREATE_MATERIAL;
	}

	//Create filter geometry
	createBaseFilter(filterSizeX,  //dx
			 filterSizeY,  //dy
			 f.second,  //dz
			 1,  //Number of subdivisions in X axis
			 ss, //Output stream
			 fcount,  //Material index
			 filterName,  //Filter name
			 "world",   //Parent name
			 false,    //Do not print number of objects
			 vector3D<double>(0.0,0.0,forigin - f.second/2.0));  //Origin
	forigin -= f.second;
      }

      // ** Create the detector

      //Create filter geometry
      createBaseFilter(filterSizeX,  //dx
		       filterSizeY,  //dy
		       0.05,  //dz
		       1,  //Number of subdivisions in X axis
		       ss, //Output stream
		       1,  //Material index
		       "detector",  //Filter name
		       "world",   //Parent name
		       false,    //Do not print number of objects
		       vector3D<double>(0.0,0.0,forigin - filter2col - 0.05/2.0));//Origin

      
      // ** Configure simulation

      // * Create simulation context
  
      //Create a context
      std::shared_ptr<pen_context> pcontext(new pen_context()); 
      pen_context& contextSim = *pcontext.get();

      //Create context configuration
      pen_parserSection contextConfig;
      contextConfig.set("material-eabs/electron", 1.0e35);
      contextConfig.set("material-eabs/positron", 1.0e35);
      contextConfig.set("material-eabs/gamma", emin);
      
      //Define material files
      for(int i = 0; i < static_cast<int>(filters.size()); ++i){
	std::string path = "materials/filter" + std::to_string(i+1) + "/number"; 
	contextConfig.set(path, i+1);
	path = "materials/filter" + std::to_string(i+1) + "/filename";
	contextConfig.set(path, "filter" + std::to_string(i+1) + ".mat");
      }

      //Configure context with no geometry
      pen_parserSection mateInfoSection;
      int err = contextSim.configure(maxE, contextConfig, mateInfoSection, verbose);
      if(err != pen_context::SUCCESS){
	return ERROR_ON_CONTEXT_INITIALIZATION;
      }

      // * Configure geometry

      //Create a mesh geometry instance
      pen_meshBodyGeo* geometry = new pen_meshBodyGeo;

      //Set geometry file in geometry instance
      geometry->preloadGeo.assign(ss.str());

      if(!geoFilename.empty()){
	FILE* fout = fopen(geoFilename.c_str(), "w");
	if(fout != nullptr){
	  fprintf(fout,"%s", ss.str().c_str());
	  fclose(fout);
	}
      }

      //Set the geometry configuration
      pen_parserSection geoConfig;

      //Set inmemory geometry
      geoConfig.set("memory-file", true);

      //Set detector object as perfect absorber
      geoConfig.set("eabs/detector/electron", 1.0e35);      
      geoConfig.set("eabs/detector/positron", 1.0e35);      
      geoConfig.set("eabs/detector/gamma"   , 1.0e35);

      //Set detector kdet
      geoConfig.set("kdet/detector", 1);

      //Configure geometry
      if(geometry->configure(geoConfig,verbose) != PEN_MESHBODY_GEO_SUCCESS){
	printf("Unexpected Error: Unable to construct the geometry. "
	       "Please, report this error\n");
	delete geometry;
	printf("      Geometry string: \n%s\n.", ss.str().c_str());
	return ERROR_ON_GEOMETRY_INITIALIZATION;
      }

      //Set geometry in context
      if(contextSim.setGeometry(geometry) != 0){
	printf("Unexpected Error: Unable to set the geometry to context. "
	       "Please, report this error\n");
	delete geometry;
	return ERROR_ON_GEOMETRY_INITIALIZATION;
      }

      //Set the number of threads to use
#ifndef _PEN_USE_THREADS_
      //Without multithreading
      unsigned nThreads = 1;
#else
      //With multithreading
      unsigned nThreads = nthreadsIn;
      if(nThreads == 0){
	nThreads = std::max(std::thread::hardware_concurrency(),1u);
      }
#endif

      //Create a base configuration
      penred::simulation::simConfig basicConf;
      basicConf.verbose = verbose;

      //Generate configurations for each thread
      std::vector<penred::simulation::simConfig> simConf =
	basicConf.generateThreadConfigs(nThreads);

      //Create vectors to store resulting particles in each thread
      std::vector<std::vector<detectedPart>> result(nThreads);

      //Calculate number of histories per thread
      unsigned long long nHistsPerThread = nHists/nThreads;
      
#ifdef _PEN_USE_THREADS_

      std::vector<std::thread> threads;

      for(size_t ith = 0; ith < nThreads; ++ith){

	threads.emplace_back([&, ith, nHists, nThreads, nHistsPerThread]
				  (){	  

	  unsigned long long hists2simulate = nHistsPerThread;
	  if(ith == nThreads-1)
	    hists2simulate += nHists % nThreads;
	  
	  simulateVectorAndDetect(simConf[ith],contextSim,
				  result[ith],
				  particlesIn,
				  ith*nHistsPerThread,
				  hists2simulate,
				  "filter-and-collimate",
				  1); //Detector index
	  
	});
      }

      for(size_t ith = 0; ith < nThreads; ++ith){
	threads[ith].join();	
      }
      
#else

      simulateVectorAndDetect(simConf[0], context, result[0],
			      particlesIn,
			      nHists,
			      "filter-and-collimate",
			      1); //Detector index      
#endif

      //Clear input particles and store the resulting ones
      particlesIn.clear();
      for(size_t ith = 0; ith < nThreads; ++ith){
	particlesIn.insert(particlesIn.end(),
			   result[ith].begin(), result[ith].end());
      }

      //Delete geometry
      delete geometry;

      //Collimate particles after filter

      const double colz1 = forigin - filter2col;
      const double colz2 = colz1 - coldz;

      if(verbose > 1){
	printf("\n + Apply collimation after filters:\n"
	       "      Z up        : %.4f cm\n"
	       "      Z bot       : %.4f cm\n"
	       "      Top aperture: (%.4f x %.4f) cm\n"
	       "      Bot aperture: (%.4f x %.4f) cm\n"
	       ,colz1, colz2,
	       coldx1, coldy1,
	       coldx2, coldy2);
      }
      
      collimateInVoid(particlesIn,
		      colz1, colz2,
		      coldx1, coldy1,
		      coldx2, coldy2);
      
      return SUCCESS;
    }

    int readerHVL::beginSectionFamily(const std::string& pathInSection,
				      const size_t,
				      const unsigned verbose){
      if(family == -1){
	if(pathInSection.compare("filters") == 0){
	  family = 0;
	}
	else{
	  if(verbose > 0)
	    printf("Error: Unhandled section '%s'\n", pathInSection.c_str());
	  return UNHANDLED;
	}
      }
      else{
	if(verbose > 0)
	  printf("Error: Unknown section family '%s'\n", pathInSection.c_str());
	return UNHANDLED;
      }
      return SUCCESS;
    }

    int readerHVL::endSectionFamily(const unsigned){

      if(family == 0){
	family = -1;
      }
      else{
	return UNHANDLED;
      }
      return SUCCESS;
    }

    int readerHVL::beginSection(const std::string&,
				const unsigned){
      
      if(family == 0){
	filters.emplace_back(13,0.1);
      }
      else{
	return UNHANDLED;
      }
      return SUCCESS;      
    }

    int readerHVL::endSection(const unsigned){ return SUCCESS; }
    
    
    int readerHVL::storeElement(const std::string& pathInSection,
				const pen_parserData& element,
				const unsigned verbose){
      
      if(family == -1){
	//Root section
	if(pathInSection.compare("anode/simulation/enabled") == 0){
	  simAnode = element;
	}
	else if(pathInSection.compare("anode/simulation/nhists") == 0){
	  nhists = element;
	}
	else if(pathInSection.compare("anode/simulation/energy") == 0){
	  beamEnergy = element;
	}      
	else if(pathInSection.compare("minE") == 0){
	  minE = element;
	}
	else if(pathInSection.compare("anode/simulation/focalSpot") == 0){
	  focalSpot = element;
	}
	else if(pathInSection.compare("anode/simulation/angle") == 0){
	  anodeAngle = element;
	}
	else if(pathInSection.compare("anode/material/z") == 0){
	  anodeZ = element;
	}
	else if(pathInSection.compare("anode/material/density") == 0){
	  anodeDensity = element;
	}      
	else if(pathInSection.compare("source/zpos") == 0){
	  sourcePositionZ = element;
	}
	else if(pathInSection.compare("distance/source-detector") == 0){
	  source2det = element;
	}
	else if(pathInSection.compare("distance/source-filter") == 0){
	  source2filter = element;
	}
	else if(pathInSection.compare("distance/source-register") == 0){
	  source2register = element;
	}
	else if(pathInSection.compare("detector/radius") == 0){
	  detectorRad = element;
	}
	else if(pathInSection.compare("nThreads") == 0){
	  nThreads = element;
	}
	else if(pathInSection.compare("print-geometries") == 0){
	  printGeometries = element;
	}	
	else{
	  if(verbose > 0){
	    printf("Error: Unexpected path: '%s'\n", pathInSection.c_str());
	  }
	  return UNHANDLED;
	}
      }
      else if(family == 0){
	//Filters section
	if(pathInSection.compare("z") == 0){
	  filters.back().first = element;
	}
	else if(pathInSection.compare("width") == 0){
	  filters.back().second = element;
	}
	else{
	  return UNHANDLED;
	}
      }else{
	return UNHANDLED;
      }

      return SUCCESS;
    }

    int readerHVL::storeString(const std::string& pathInSection,
			       const std::string& element,
			       const unsigned){
      
      if(family == -1){
	if(pathInSection.compare("anode/material/filename") == 0){
	  anodeMatFilename = element;
	}
	else{
	  return UNHANDLED;
	}
      }
      else{
	return UNHANDLED;
      }

      return SUCCESS;
    }

    int HVL(const pen_parserSection config,
	    const std::vector<detectedPart>& particlesIn,
	    double& hvl,
	    const unsigned verbose){
      
      //This function computes the HVL value from a x-ray configuration.
      //The returned value is in cm
      
      //Read material information from config section
      readerHVL reader;
      int err = reader.read(config,verbose);
      if(err != readerHVL::SUCCESS){
	return err;
      }

      //Calculate the beam maximum aperture
      //
      //          fs
      //         |--|        n
      //        /|  |\       |
      //       /_|  | \      |
      //      / a|  |  \     | source to detector
      //     /   |  |   \    |
      //    /    |  |    \   |
      //   |--------------|  u
      //     | detector
      //     u
      // (detector-fs)/2
      //
      //
      double tanBeamAperture;
      if(reader.focalSpot < 2.0*reader.detectorRad)
	tanBeamAperture = (reader.detectorRad - reader.focalSpot/2.0)/reader.source2det;
      else
	tanBeamAperture = 0.0;
	  
      const double beamSemiApertureAngle = atan(tanBeamAperture);
      

      // ** Anode
      //-----------------

      std::vector<detectedPart> actualParticles;
      unsigned long long nHists = reader.nhists;

      //Check if the anode must be simulated
      if(reader.simAnode){

	//Calculate anode reject angle
	double rejectAngle;
	if(reader.focalSpot < 2.0*reader.detectorRad){
	  double rejectAngleTan =
	    (reader.detectorRad + reader.focalSpot/2.0)/reader.source2det;
	  rejectAngle = atan(rejectAngleTan);
	}else{
	  double rejectAngleTan =
	    (reader.focalSpot/2.0)/reader.source2det;
	  rejectAngle = atan(rejectAngleTan);	  
	}
	

	//Simulate anode
	if(verbose > 1){
	  printf(" + Simulating anode: \n"
		 "     Atomic number        : %u\n"
		 "     Beam energy          : %.3f keV\n"
		 "     Detection threshold  : %.3f keV\n"
		 "     Focal spot           : %.3f cm\n"
		 "     Anode angle          : %.3f deg\n"
		 "     Max histories        : %llu\n"
		 "     Reject angle         : %.6f deg\n",
		 reader.anodeZ,
		 reader.beamEnergy/1000.0,
		 reader.minE/1000.0,
		 reader.focalSpot,
		 reader.anodeAngle,
		 reader.nhists,
		 rejectAngle*180.0/constants::PI);
	}

	//Set the source origin to 0.0
	reader.sourcePositionZ = 0.0;

	//Check if the particle filename has been provided and is accessible
	bool fileExists = false;
	if(reader.anodeMatFilename.compare("-") != 0){

	  //Try to open material file
	  FILE* fmat = nullptr;
	  fmat = fopen(reader.anodeMatFilename.c_str(),"r");
	  if(fmat != nullptr){
	    fclose(fmat);
	    fileExists = true;
	  }
	  
	}

	if(!fileExists){
	  //If material does not exist, create it
	  if(reader.anodeMatFilename.compare("-") == 0)
	    reader.anodeMatFilename = "anode.mat";

	  std::string errorString;

	  if(reader.anodeDensity > 0.0){
	    err = penred::penMaterialCreator::createMat("anode",
							reader.anodeDensity,
							reader.anodeZ,
							errorString,
							reader.anodeMatFilename);
	  }else{
	    err = penred::penMaterialCreator::createMat(reader.anodeZ,
							reader.anodeMatFilename,
							errorString);
	  }

	  if(err != 0){
	    if(verbose > 0){
	      printf("Unable to create anode material.\n"
		     "%s\n"
		     "IRETRN =%d\n",
		     errorString.c_str(),
		     err);
	    }
	    return -1;
	  }
	}

	err = simAnode(reader.anodeMatFilename.c_str(),
		       reader.beamEnergy,
		       reader.minE,
		       reader.focalSpot,
		       reader.anodeAngle,
		       reader.nhists,
		       100000000,
		       reader.source2register,
		       actualParticles,
		       rejectAngle*180.0/constants::PI,
		       true,  //Only register photons
		       verbose > 3 ? verbose : 1);
	if(err != 0){
	  return -2;
	}

	if(verbose > 1){
	  printf("     Distance to register : %.5f\n", reader.source2register);
	}
	
      }
      else{

	if(verbose > 1){
	  printf(" + Anode pre-simulated. \n");
	}
	//Save particles from input
	actualParticles = particlesIn;
	//Calculate histories
	nHists = 0;
	for(const detectedPart& part : actualParticles){
	  nHists += part.dh;
	}	
      }

      if(verbose > 1){
	//Print anode spectrum

	printf("   - Printing anode spectrum to 'x-ray-anode-spectrum.dat'.\n");
	
	if(detectedPart::printSpectrum(actualParticles,
				       "x-ray-anode-spectrum.dat",
				       200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-anode-spectrum.dat'\n");
	  return ERROR_UNABLE_TO_OPEN_FILE;
	}

	printf("   - Printing anode spatial distribution to 'x-ray-anode-spatial.dat'.\n");
	
	if(detectedPart::printSpatialXY(actualParticles,
					"x-ray-anode-spatial.dat",
					200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-anode-spatial.dat'\n");
	  return ERROR_UNABLE_TO_OPEN_FILE;
	}	
	
      }
      

      // ** Collimator
      //-----------------
      
      //Now, collimate the resulting particles from the anode simulation

      //Calculate the first collimator apertures
      const double dcol1 = 0.1;
      const double ds2col1 = reader.source2register + dcol1;
      
      double col1z1 = reader.sourcePositionZ - ds2col1;
      double col1z2 = col1z1 - dcol1;
      
      double col1r1 = reader.focalSpot + 2.0*ds2col1*sin(beamSemiApertureAngle);
      double col1r2 = reader.focalSpot + 2.0*(ds2col1 + dcol1)*sin(beamSemiApertureAngle);


      if(verbose > 1){
	printf("\n + Apply anode collimation:\n"
	       "      Z up        : %.4f cm\n"
	       "      Z bot       : %.4f cm\n"
	       "      Top aperture: %.4f cm\n"
	       "      Bot aperture: %.4f cm\n"
	       ,col1z1, col1z2,
	       col1r1, col1r2);
      }
      
      collimateInVoid(actualParticles,
		      col1z1, col1z2,
		      col1r1, col1r1,
		      col1r2, col1r2);

      if(verbose > 1){

	//Print collimated anode spectrum
	printf("   - Printing collimated anode spectrum to"
	       " 'x-ray-collimated-anode-spectrum.dat'.\n");
	
	if(detectedPart::printSpectrum(actualParticles,
				       "x-ray-collimated-anode-spectrum.dat",
				       200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-collimated-anode-spectrum.dat'\n");
	  printf(" Number of particles: %lu\n",
		 static_cast<unsigned long>(actualParticles.size()));
	  return ERROR_UNABLE_TO_OPEN_FILE;
	}

	printf("   - Printing collimated anode spectrum to"
	       " 'x-ray-collimated-anode-spatial.dat'.\n");
	
	if(detectedPart::printSpatialXY(actualParticles,
					"x-ray-collimated-anode-spatial.dat",
					200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-collimated-anode-spatial.dat'\n");
	  printf(" Number of particles: %lu\n",
		 static_cast<unsigned long>(actualParticles.size()));
	  return ERROR_UNABLE_TO_OPEN_FILE;
	}
		
      }

      // ** Filtration and collimation
      //--------------------------------
      double filtersOrigin = col1z2-dcol1;

      //Calculate total filters width
      double filtersWidth = 0.0;
      for(const auto& f : reader.filters){
	filtersWidth += f.second;
      }
      const double filtersEnd = filtersOrigin-filtersWidth;


      const double col2z1 = filtersEnd - dcol1;
      //const double col2z2 = col2z1 - dcol1;
      
      const double ds2col2 = reader.sourcePositionZ - col2z1;
      const double ds2col2Bot = ds2col2 + dcol1;
      
      const double col2r1 = reader.focalSpot + 2.0*ds2col2*sin(beamSemiApertureAngle);
      const double col2r2 = reader.focalSpot + 2.0*ds2col2Bot*sin(beamSemiApertureAngle);

      if(verbose > 1){
	printf("\n + Filter anode beam:\n"
	       "      Filters start   : %.4f cm\n"
	       "      Filters width   : %.4f cm\n"
	       "      Filters end     : %.4f cm\n"
	       "      Collimator start: %.4f cm\n"
	       "      Collimator width: %.4f cm\n"
	       "      Top aperture    : %.4f cm\n"
	       "      Bot aperture    : %.4f cm\n"	       
	       "      Filters         : \n"
	       ,filtersOrigin, filtersWidth, filtersEnd,
	       col2z1, dcol1, col2r1, col2r2);

	for(const auto& f : reader.filters){
	  printf("         - Z: %d; width: %.4f cm\n",
		 f.first, f.second);
	}

	if(reader.printGeometries)
	  printf("      Filters geometry saved in 'filters.geo'\n");
      }

      std::string filterGeoFilename;
      if(reader.printGeometries)
	filterGeoFilename = "filters.geo";
      err = filterWithCollimation(actualParticles,
				  filtersOrigin,
				  dcol1,
				  reader.filters,
				  dcol1,
				  col2r1, col2r1,
				  col2r2, col2r2,
				  reader.minE,
				  reader.nThreads,
				  verbose > 3 ? verbose : 1,
				  filterGeoFilename);

      if(err != SUCCESS){
	return -3;
      }

      if(verbose > 1){

	//Print filtered spectrum
	printf("   - Printing filtered spectrum to"
	       " 'x-ray-filtered-spectrum.dat'.\n");
	
	if(detectedPart::printSpectrum(actualParticles,
				       "x-ray-filtered-spectrum.dat",
				       200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-filtered-spectrum.dat'\n");
	  return ERROR_UNABLE_TO_OPEN_FILE;
	}

	printf("   - Printing filtered spatial distribution to"
	       " 'x-ray-filtered-spatial.dat'.\n");
	
	if(detectedPart::printSpatialXY(actualParticles,
				       "x-ray-filtered-spatial.dat",
					200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-filtered-spatial.dat'\n");
	  return ERROR_UNABLE_TO_OPEN_FILE;
	}	
      }
      

      // ** Iterate to obtain the HVL value
      //-------------------------------------


      hvl = 1.0; //CACA
      return SUCCESS;
    }
    
  };
};
