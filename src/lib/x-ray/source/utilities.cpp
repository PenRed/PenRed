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

	return errors::NEGATIVE_DISTANCE;
      }

      if(filters.size() == 0){
	return errors::NO_FILTER_PROVIDED;
      }

      if(particlesIn.size() == 0){
	//No particles to simulate
	return errors::SUCCESS;
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
	  return errors::INVALID_Z;
	if(f.second <= 0.0)
	  return errors::INVALID_SIZE;
	filtersEnd -= f.second;
      }

      //Calculate filter size
      double filterSizeX = 2.0*std::max(coldx1, coldx2);
      double filterSizeY = 2.0*std::max(coldy1, coldy2);

      if(verbose > 1){
	printf("\n + Simulate filters:\n"
	       "      Minimum energy  : %.4E keV\n"
	       "      Filters start   : %.4f cm\n"
	       "      Filters width   : %.4f cm\n"
	       "      Filters end     : %.4f cm\n"
	       "      Filters         : \n",
	       emin,
	       filterZorigin,
	       filterZorigin-filtersEnd,
	       filtersEnd);

	for(const auto& f : filters){
	  printf("         - Z: %d; width: %.4f cm\n",
		 f.first, f.second);
	}

	if(!geoFilename.empty())
	  printf("      Filters geometry saved in '%s'\n", geoFilename.c_str());
      }     

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
	  return errors::UNABLE_TO_CREATE_MATERIAL;
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
      int err = contextSim.configure(maxE, contextConfig,
				     mateInfoSection, verbose >= 3 ? verbose : 1);
      if(err != pen_context::SUCCESS){
	return errors::ERROR_ON_CONTEXT_INITIALIZATION;
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
      if(geometry->configure(geoConfig,verbose >= 3 ? verbose : 1) != PEN_MESHBODY_GEO_SUCCESS){
	printf("Unexpected Error: Unable to construct the geometry. "
	       "Please, report this error\n");
	delete geometry;
	printf("      Geometry string: \n%s\n.", ss.str().c_str());
	return errors::ERROR_ON_GEOMETRY_INITIALIZATION;
      }

      //Set geometry in context
      if(contextSim.setGeometry(geometry) != 0){
	printf("Unexpected Error: Unable to set the geometry to context. "
	       "Please, report this error\n");
	delete geometry;
	return errors::ERROR_ON_GEOMETRY_INITIALIZATION;
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
      basicConf.verbose = verbose >= 3 ? verbose : 1;

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
      
      return errors::SUCCESS;
    }

    int bowtieWithCollimation(std::vector<detectedPart>& particlesIn,
			      const double filterZbot,
			      const double filter2col,
			      const unsigned matZ,
			      const std::vector<double>& bowtieWidths,
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

	return errors::NEGATIVE_DISTANCE;
      }

      if(bowtieWidths.size() == 0){
	return errors::NO_FILTER_PROVIDED;
      }

      if(particlesIn.size() == 0){
	//No particles to simulate
	return errors::SUCCESS;
      }

      //Get maximum energy and total number of histories
      double maxE = 0.0;
      unsigned long long nHists = 1;
      for(const detectedPart& part : particlesIn){
	if(maxE < part.state.E)
	  maxE = part.state.E;
	nHists += part.dh;
      }

      //Get maximum bowtie width
      const auto minmaxWidth =
	std::minmax_element(bowtieWidths.cbegin(), bowtieWidths.cend());

      const double minWidth = *minmaxWidth.first;
      const double maxWidth = *minmaxWidth.second;

      //Ensure all widths are greater than 0
      if(minWidth <= 0.0){
	return errors::NEGATIVE_DISTANCE;	
      }

      double filterZtop = filterZbot + maxWidth;

      //Calculate filter size
      double filterSizeX = std::max(coldx1, coldx2);
      double filterSizeY = std::max(coldy1, coldy2);

      if(verbose > 1){
	printf("\n + Simulate bowtie filter:\n"
	       "      filter size   : (%.4f, %.4f) cm\n"
	       "      minimum width : %.4f cm\n"
	       "      maximum width : %.4f cm\n"
	       "      Z up          : %.4f cm\n"
	       "      Z bot         : %.4f cm\n",	       
	       filterSizeX, filterSizeY,
	       minWidth, maxWidth,
	       filterZtop, filterZbot);
      }
      
      //Create a string to store the created geometries
      std::stringstream ss;

      //Create the void world containing the filter and detector geometry
      ss << "# Number of objects:\n " << 3 << std::endl;

      // ** Create the world
      
      //Create a void world to contain all filters and the detector
      double worldSize = maxWidth + filter2col + coldz + 1.0;
      createBaseFilter(2.0*filterSizeX+0.1,  //dx
		       2.0*filterSizeY+0.1,  //dy
		       worldSize,  //dz
		       1,  //Number of subdivisions in X axis
		       ss, //Output stream
		       0,  //Material index
		       "world",  //"Filter" name
		       "void",   //Parent name
		       false,    //Do not print number of objects
		       vector3D<double>(0.0,0.0, filterZtop - maxWidth/2.0 + 0.5)); //Origin

      // ** Create the bowtie filter

      // Material
      std::string errorString;
      int err = penred::penMaterialCreator::createMat(matZ,
						      "bowtie.mat",
						      errorString);
      if(err != 0){
	if(verbose > 0){
	  printf("Unable to create material for bowtie filter.\n"
		 "%s\n"
		 "IRETRN =%d\n",
		 errorString.c_str(),
		 err);
	}
	return errors::UNABLE_TO_CREATE_MATERIAL;
      }

      //Create filter geometry
      createBaseFilter(filterSizeX,  //dx
		       filterSizeY,  //dy
		       minWidth,  //dz
		       bowtieWidths.size(),  //Number of subdivisions in X axis
		       ss, //Output stream
		       1,  //Material index
		       "bowtie",  //Filter name
		       "world",   //Parent name
		       false,    //Do not print number of objects
		       vector3D<double>(0.0,0.0, filterZbot + minWidth/2.0)); //Origin      

      // ** Create the detector

      //Create the geometry
      createBaseFilter(2.0*filterSizeX,  //dx
		       2.0*filterSizeY,  //dy
		       0.05,  //dz
		       1,  //Number of subdivisions in X axis
		       ss, //Output stream
		       1,  //Material index
		       "detector",  //Filter name
		       "world",   //Parent name
		       false,    //Do not print number of objects
		       vector3D<double>(0.0,0.0,filterZbot - filter2col + 0.9*coldz/2.0));//Origin

      
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

      contextConfig.set("materials/bowtie/number", 1);
      contextConfig.set("materials/bowtie/filename", "bowtie.mat");

      //Configure context with no geometry
      pen_parserSection mateInfoSection;
      err = contextSim.configure(maxE, contextConfig, mateInfoSection, verbose >= 3 ? verbose : 1);
      if(err != pen_context::SUCCESS){
	return errors::ERROR_ON_CONTEXT_INITIALIZATION;
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

      //Conform the bowtie filter
      for(unsigned i = 0; i < bowtieWidths.size(); ++i){
	std::string transPath = "transforms/bowtie/conform" + std::to_string(i);
	std::string vergexGroup = "+x-" + std::to_string(i);
	geoConfig.set((transPath + "/index").c_str(), static_cast<int>(i));
	geoConfig.set((transPath + "/vertex-group").c_str(), vergexGroup);
	geoConfig.set((transPath + "/transforms/move/type").c_str(), "TRANSLATION_Z");
	geoConfig.set((transPath + "/transforms/move/index").c_str(), 0);
	geoConfig.set((transPath + "/transforms/move/ds").c_str(), bowtieWidths[i]-minWidth);
      }
      
      //Configure geometry
      if(geometry->configure(geoConfig,verbose >= 3 ? verbose : 1) != PEN_MESHBODY_GEO_SUCCESS){
	printf("Unexpected Error: Unable to construct the geometry. "
	       "Please, report this error\n");
	delete geometry;
	printf("      Geometry string: \n%s\n.", ss.str().c_str());
	return errors::ERROR_ON_GEOMETRY_INITIALIZATION;
      }

      //Check if geometry configuration must be saved
      if(!geoFilename.empty()){
	FILE* fout = fopen(("config_" + geoFilename).c_str(), "w");
	if(fout != nullptr){
	  geoConfig.set("type", "MESH_BODY");
	  geoConfig.set("input-file", geoFilename);	  
	  geoConfig.set("memory-file", false, true);
	  
	  fprintf(fout,"%s", geoConfig.stringify().c_str());
	  fclose(fout);
	}
      }

      //Set geometry in context
      if(contextSim.setGeometry(geometry) != 0){
	printf("Unexpected Error: Unable to set the geometry to context. "
	       "Please, report this error\n");
	delete geometry;
	return errors::ERROR_ON_GEOMETRY_INITIALIZATION;
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
      basicConf.verbose = verbose >= 3 ? verbose : 1;

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

      const double colz1 = filterZbot - filter2col;
      const double colz2 = colz1 - coldz;

      if(verbose > 1){
	printf("\n + Apply collimation after bowtie filter:\n"
	       "      Z up        : %.4f cm\n"
	       "      Z bot       : %.4f cm\n"
	       "      Top aperture: (%.4f x %.4f) cm\n"
	       "      Bot aperture: (%.4f x %.4f) cm\n",
	       colz1, colz2,
	       coldx1, coldy1,
	       coldx2, coldy2);
      }
      
      collimateInVoid(particlesIn,
		      colz1, colz2,
		      coldx1, coldy1,
		      coldx2, coldy2);
      
      return errors::SUCCESS;
    }

    int readerBeamCT::beginSectionFamily(const std::string& pathInSection,
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

    int readerBeamCT::endSectionFamily(const unsigned){

      if(family == 0){
	family = -1;
      }
      else{
	return UNHANDLED;
      }
      return SUCCESS;
    }

    int readerBeamCT::beginSection(const std::string&,
				const unsigned){
      
      if(family == 0){
	filters.emplace_back(13,0.1);
      }
      else{
	return UNHANDLED;
      }
      return SUCCESS;      
    }

    int readerBeamCT::endSection(const unsigned){ return SUCCESS; }
    
    
    int readerBeamCT::storeElement(const std::string& pathInSection,
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
	else if(pathInSection.compare("bowtie/minWidth") == 0){
	  bowtieMinW = element;
	}
	else if(pathInSection.compare("bowtie/maxWidth") == 0){
	  bowtieMaxW = element;
	}
	else if(pathInSection.compare("bowtie/z") == 0){
	  bowtieZ = element;
	}
	else if(pathInSection.compare("bowtie/segments") == 0){
	  bowtieSegments = element;
	}	
	else if(pathInSection.compare("source/zpos") == 0){
	  sourcePositionZ = element;
	}	
	else if(pathInSection.compare("distance/source-isocenter") == 0){
	  source2isocenter = element;
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
	else if(pathInSection.compare("detector/dx") == 0){
	  detectorDx = element;
	}
	else if(pathInSection.compare("detector/dy") == 0){
	  detectorDy = element;
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

    int readerBeamCT::storeString(const std::string& pathInSection,
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

    int generateBeamCT(const pen_parserSection config,
		       const std::vector<detectedPart>& particlesIn,
		       const unsigned verbose){
      
      //This function computes the HVL value from a x-ray configuration.
      //The returned value is in cm
      
      //Read material information from config section
      readerBeamCT reader;
      int err = reader.read(config,verbose);
      if(err != readerBeamCT::SUCCESS){
	return err;
      }

      if(verbose > 1){
	printf("\n + CT characteristics:\n"
	       "     Source to isocenter  : %.4f cm\n"
	       "     Source to detector   : %.4f cm\n"
	       "     Source to filters    : %.4f cm\n"
	       "     Detector X size      : %.4f cm\n"
	       "     Detector Y size      : %.4f cm\n"
	       "     Bowtie minimum width : %.4f cm\n"
	       "     Bowtie maximum width : %.4f cm\n",
	       reader.source2isocenter,
	       reader.source2det,
	       reader.source2filter,
	       reader.detectorDx,
	       reader.detectorDy,
	       reader.bowtieMinW,
	       reader.bowtieMaxW);
      }

      //Check distance between isocenter and detector

      const double phantomDiameter = 32.0;
      const double iso2det = reader.source2det - reader.source2isocenter;
      if(iso2det <= phantomDiameter/2.0){
	if(verbose > 0)
	printf("generateBeamCT:Error:The distance from isocenter to detector"
	       " is lesser than the phantom radius (%.4f cm)\n"
	       " Source to isocenter: %.4f\n"
	       " Source to detector : %.4f\n",
	       phantomDiameter/2.0,
	       reader.source2isocenter, reader.source2det);
	return errors::INVALID_DISTANCE;
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

      //Calculate desired beam aperture in X and Y axis
      double tanBeamApertureX;
      double tanBeamApertureY;

      //X axis
      if(reader.focalSpot < reader.detectorDx)
	tanBeamApertureX = (reader.detectorDx/2.0 - reader.focalSpot/2.0)/reader.source2det;
      else
	tanBeamApertureX = 0.0;

      //Y Axis
      if(reader.focalSpot < reader.detectorDy)
	tanBeamApertureY = (reader.detectorDy/2.0 - reader.focalSpot/2.0)/reader.source2det;
      else
	tanBeamApertureY = 0.0;
	  
      //const double beamSemiApertureAngleX = atan(tanBeamApertureX);
      //const double beamSemiApertureAngleY = atan(tanBeamApertureY);
      

      // ** Anode
      //-----------------

      std::vector<detectedPart> actualParticles;
      unsigned long long nHists = reader.nhists;

      //Check if the anode must be simulated
      if(reader.simAnode){

	//Calculate anode reject angle
	const double detectorRad =
	  sqrt(reader.detectorDx*reader.detectorDx +
	       reader.detectorDy*reader.detectorDy)/2.0;
	double rejectAngle;
	if(reader.focalSpot < 2.0*detectorRad){
	  double rejectAngleTan =
	    (detectorRad + reader.focalSpot/2.0)/reader.source2det;
	  rejectAngle = atan(rejectAngleTan);
	}else{
	  double rejectAngleTan =
	    (reader.focalSpot/2.0)/reader.source2det;
	  rejectAngle = atan(rejectAngleTan);	  
	}
	

	//Simulate anode
	if(verbose > 1){
	  printf("\n + Simulating anode: \n"
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
	  return errors::ERROR_UNABLE_TO_OPEN_FILE;
	}

	printf("   - Printing anode spatial distribution to 'x-ray-anode-spatial.dat'.\n");
	
	if(detectedPart::printSpatialXY(actualParticles,
					"x-ray-anode-spatial.dat",
					200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-anode-spatial.dat'\n");
	  return errors::ERROR_UNABLE_TO_OPEN_FILE;
	}	
	
      }

      //Calculate isocenter
      const vector3D<double> isocenter(0.0,
				       0.0,
				       reader.sourcePositionZ-
				       reader.source2isocenter);

      //Calculate detector position
      const double detectorZ = reader.sourcePositionZ - reader.source2det;

      // ** Collimator
      //-----------------
      
      //Now, collimate the resulting particles from the anode simulation

      //Calculate the first collimator apertures
      const double dcol1 = 0.1;
      const double ds2col1 = reader.source2register + dcol1;
      const double ds2col1Bot = ds2col1 + dcol1;
      
      const double col1z1 = reader.sourcePositionZ - ds2col1;
      const double col1z2 = col1z1 - dcol1;
      
      const double col1dxTop = reader.focalSpot + 2.0*ds2col1*tanBeamApertureX;
      const double col1dyTop = reader.focalSpot + 2.0*ds2col1*tanBeamApertureY;

      const double col1dxBot = reader.focalSpot + 2.0*ds2col1Bot*tanBeamApertureX;
      const double col1dyBot = reader.focalSpot + 2.0*ds2col1Bot*tanBeamApertureY;

      if(verbose > 1){
	printf("\n + Apply anode collimation:\n"
	       "      Z up        : %.4f cm\n"
	       "      Z bot       : %.4f cm\n"
	       "      Top aperture: (%.4f x %.4f) cm\n"
	       "      Bot aperture: (%.4f x %.4f) cm\n"
	       ,col1z1, col1z2,
	       col1dxTop, col1dyTop,
	       col1dxBot, col1dyBot);
      }
      
      collimateInVoid(actualParticles,
		      col1z1, col1z2,
		      col1dxTop, col1dyTop,
		      col1dxBot, col1dyBot);

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
	  return errors::ERROR_UNABLE_TO_OPEN_FILE;
	}

	printf("   - Printing collimated anode spectrum to"
	       " 'x-ray-collimated-anode-spatial.dat'.\n");
	
	if(detectedPart::printSpatialXY(actualParticles,
					"x-ray-collimated-anode-spatial.dat",
					200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-collimated-anode-spatial.dat'\n");
	  printf(" Number of particles: %lu\n",
		 static_cast<unsigned long>(actualParticles.size()));
	  return errors::ERROR_UNABLE_TO_OPEN_FILE;
	}
		
      }

      // ** Create phantom water material
      //--------------------------------
      
      std::string errorString;
      err = penred::penMaterialCreator::createMat(278,
						  "phantom.mat",
						  errorString);
      if(err != 0){
	if(verbose > 0){
	  printf("Unable to create phantom material.\n"
		 "%s\n"
		 "IRETRN =%d\n",
		 errorString.c_str(),
		 err);
	}
	return errors::UNABLE_TO_CREATE_MATERIAL;
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

      const double dcol2 = dcol1;

      const double col2z1 = filtersEnd - dcol1;
      const double col2z2 = col2z1 - dcol2;
      
      const double ds2col2 = reader.sourcePositionZ - col2z1;
      const double ds2col2Bot = ds2col2 + dcol2;

      const double col2dxTop = reader.focalSpot + 2.0*ds2col2*tanBeamApertureX;
      const double col2dyTop = reader.focalSpot + 2.0*ds2col2*tanBeamApertureY;

      const double col2dxBot = reader.focalSpot + 2.0*ds2col2Bot*tanBeamApertureX;
      const double col2dyBot = reader.focalSpot + 2.0*ds2col2Bot*tanBeamApertureY;

      //Avoid simulation if no filters have been provided
      if(reader.filters.size() > 0){
      

	std::string filterGeoFilename;
	if(reader.printGeometries)
	  filterGeoFilename = "filters.geo";
	err = filterWithCollimation(actualParticles,
				    filtersOrigin,
				    dcol2,
				    reader.filters,
				    dcol2,
				    col2dxTop, col2dyTop,
				    col2dxBot, col2dyBot,
				    reader.minE,
				    reader.nThreads,
				    verbose,
				    filterGeoFilename);

	if(err != errors::SUCCESS){
	  if(verbose > 0){
	    printf("generateBeamCT: Error on filters simulations. "
		   "Error code: %d\n", err);
	  }
	  return err;
	}

	if(verbose > 1){

	  //Print filtered spectrum
	  printf("   - Printing filtered spectrum to"
		 " 'x-ray-filtered-spectrum.dat'.\n");
	
	  if(detectedPart::printSpectrum(actualParticles,
					 "x-ray-filtered-spectrum.dat",
					 200, nHists) != 0){
	    printf("Error: Unable to open file 'x-ray-filtered-spectrum.dat'\n");
	    return errors::ERROR_UNABLE_TO_OPEN_FILE;
	  }

	  printf("   - Printing filtered spatial distribution to"
		 " 'x-ray-filtered-spatial.dat'.\n");
	
	  if(detectedPart::printSpatialXY(actualParticles,
					  "x-ray-filtered-spatial.dat",
					  200, nHists) != 0){
	    printf("Error: Unable to open file 'x-ray-filtered-spatial.dat'\n");
	    return errors::ERROR_UNABLE_TO_OPEN_FILE;
	  }	
	}
      }

      // ** Bowtie Initialization
      //-------------------------------------

      //Save particles before bowtie filter
      const std::vector<detectedPart> bowtieIncidentParticles = actualParticles;

      std::vector<double> bowtieWidths(reader.bowtieSegments,reader.bowtieMinW);

      const double bowtieOrigin = col2z2 - dcol2;
      
      const double bowtieEnd =
	bowtieOrigin - reader.bowtieMaxW;
      
      const double dcol3 = dcol1;

      const double col3z1 = bowtieEnd - dcol2;
      //const double col3z2 = col3z1 - dcol3;
      
      const double ds2col3 = reader.sourcePositionZ - col3z1;
      const double ds2col3Bot = ds2col3 + dcol1;

      const double col3dxTop = reader.focalSpot + 2.0*ds2col3*tanBeamApertureX;
      const double col3dyTop = reader.focalSpot + 2.0*ds2col3*tanBeamApertureY;

      const double col3dxBot = reader.focalSpot + 2.0*ds2col3Bot*tanBeamApertureX;
      const double col3dyBot = reader.focalSpot + 2.0*ds2col3Bot*tanBeamApertureY;

      //Calculate bowtie size
      const double bowtieSizeX = std::max(col3dxTop, col3dxBot);
      //const double bowtieSizeY = std::max(col3dyTop, col3dyBot);

      const double bowtieXmin = -bowtieSizeX/2.0;
      const double bowtieXmax =  bowtieSizeX/2.0;
      const double bowtieDx = bowtieSizeX/static_cast<double>(reader.bowtieSegments);

      //Set the number of bins in the estimated profile
      const unsigned profBins = reader.bowtieSegments;

      //Calculate detector limits
      const double detectorXmin = -reader.detectorDx/2.0;
      const double detectorXmax =  reader.detectorDx/2.0;
      const double detectorYmin = -reader.detectorDy/2.0;
      const double detectorYmax =  reader.detectorDy/2.0;      
      
      const double profileDx = reader.detectorDx/static_cast<double>(profBins);
            
      //Generate bowtie segment influence map for each profile bin
      struct segmentInfluence{

	unsigned first, last;
	unsigned nSegments;
	double fracFirst, fracLast;
      };
      std::vector<segmentInfluence> profileInfluenceMap(profBins);

      for(size_t i = 0; i < profBins; ++i){
	//Calculate influence region at bowtie

	//         __  <- source
        //bowtie  /| |
	//  |    / | | 
	//  u__ /__|_| source2det
	//     /|  | |
	// ___/_|__|_|
	// n
	// |
	//Detector

	//Use lines between detector and source
	//
	// x = m*z + b
	//
	const double xLeft  = detectorXmin + static_cast<double>(i)*profileDx;
	const double xRight = detectorXmin + static_cast<double>(i+1)*profileDx;

	const double mLeft =
	  (xLeft - reader.focalSpot/2)/(-reader.source2det);

	const double mRight =
	  (xRight + reader.focalSpot/2)/(-reader.source2det);
	
	const double bLeft  = -reader.focalSpot/2.0 - mLeft*reader.sourcePositionZ;
	const double bRight =  reader.focalSpot/2.0 - mRight*reader.sourcePositionZ;

	//Calculate the X positions, using the lines, where the bowtie is cut
	// x = b + m*z
	
	const double xLeftBowtieBot  = bLeft  + mLeft*bowtieEnd;
	const double xRightBowtieBot = bRight + mRight*bowtieEnd;

	const double xLeftBowtieTop  = bLeft  + mLeft*bowtieOrigin;
	const double xRightBowtieTop = bRight + mRight*bowtieOrigin;

	const double xLeftBowtietMin  = std::min(xLeftBowtieBot,xLeftBowtieTop);
	const double xRightBowtietMax = std::max(xRightBowtieBot,xRightBowtieTop);

	//Set map
	profileInfluenceMap[i].first = 0;
	profileInfluenceMap[i].fracFirst = 1.0;
	profileInfluenceMap[i].last = reader.bowtieSegments;
	profileInfluenceMap[i].fracLast = 1.0;

	if(xLeftBowtietMin >= bowtieXmin){
	  profileInfluenceMap[i].first = (xLeftBowtietMin-bowtieXmin)/bowtieDx;
	  
	  profileInfluenceMap[i].fracFirst =
	    1.0 - (xLeftBowtietMin - (bowtieXmin + profileInfluenceMap[i].first*bowtieDx))/bowtieDx; 
	}
	if(xRightBowtietMax < bowtieXmax){
	  profileInfluenceMap[i].last = (xRightBowtietMax-bowtieXmin)/bowtieDx;
	  profileInfluenceMap[i].fracLast =
	    ((bowtieXmin + profileInfluenceMap[i].last*bowtieDx) - xRightBowtietMax)/bowtieDx;
	}

	//Calculate number of influenced bowtie segments
	profileInfluenceMap[i].nSegments = profileInfluenceMap[i].last - profileInfluenceMap[i].first+1;
      }

      //Create a context to get particle ranges in water
      pen_context context;

      //Set the number of materials to context (1)
      int errmat = context.setMats<pen_material>(1);
      if(errmat != 0){
	printf("Error at context material creation: %d.\n",errmat);
	return -1;
      }
  
      //Get the material
      pen_material& mat = context.getBaseMaterial(0);
	
      //Configure the material
      mat.C1=0.2;
      mat.C2=0.2;
      mat.WCC=1.0e3;
      mat.WCR=1.0e3;

      mat.EABS[PEN_ELECTRON] = 50.0E0;
      mat.EABS[PEN_PHOTON]   = 50.0E0;
      mat.EABS[PEN_POSITRON] = 50.0E0;

      const double EMAX=1.0E9;
      int INFO = 1;
      std::string PMFILEstr[constants::MAXMAT];
      PMFILEstr[0].assign("phantom.mat");
      err = context.init(EMAX,nullptr,INFO,PMFILEstr);
      if(err != 0){
	printf("Error: Unable to configure range context. Please, report this error.\n");
	return -3;
      }

      constexpr const unsigned maxIter = 1000;
      for(unsigned ibowtie = 0; ibowtie < maxIter; ++ibowtie){
      
	// ** Simulate the bowtie
	//-------------------------------------

	//Get particles before the bowtie
	if(ibowtie > 0){
	  actualParticles = bowtieIncidentParticles;
	}

	if(verbose > 1){
	  printf("\n\n ** Bowtie iteration %u\n\n", ibowtie);
	}

	std::string bowtieGeoFilename;
	if(reader.printGeometries)
	  bowtieGeoFilename = "bowtie.geo";
	err = bowtieWithCollimation(actualParticles,
				    bowtieEnd,
				    dcol3,
				    reader.bowtieZ,
				    bowtieWidths,
				    dcol3,
				    col3dxTop, col3dyTop,
				    col3dxBot, col3dyBot,
				    reader.minE,
				    reader.nThreads,
				    verbose,
				    bowtieGeoFilename);

	if(err != errors::SUCCESS){
	  if(verbose > 0){
	    printf("generateBeamCT: Error on bowtie simulation. "
		   "Error code: %d\n", err);
	  }
	  return err;
	}

	//Create a first aproximation of the bowtie filter
	//according to the actual particles


	//First, estimate the profile in the detector
	std::vector<double> profile(profBins,0.0);

	for(const detectedPart& p : actualParticles){

	  //Calculate the X position in the detector
	  const double ds2det = (detectorZ - p.state.Z)/p.state.W; //Notice W sign
	  const double xAtDet = ds2det*p.state.U + p.state.X;
	  const double yAtDet = ds2det*p.state.V + p.state.Y;
	  
	  if(xAtDet <  detectorXmin ||
	     xAtDet >= detectorXmax ||
	     yAtDet <  detectorYmin ||
	     yAtDet >= detectorYmax){
	    //The particle does not reach the detector
	    continue;
	  }
	  
	  //Add the particle contribution to the profile estimation
	  const unsigned ix = (xAtDet - detectorXmin)/profileDx;
	  
	  //Calculate the particle intersection with the phantom
	  //
	  // x**2 + (z-zp)**2 = r**2 = (phantom Radius)**2
	  // x = U*ds + x0
	  // y = V*ds + y0
	  // z = W*ds + z0
	  //
	  //Define z0' = z0-zp
	  //
	  // (U**2 + W**2)*ds**2 + 2*(U*x0 + W*z0)*ds + (x0**2 + z0**2) = (phantom Radius)**2
	  //
	  // A = (U**2 + W**2)
	  // B = 2*(U*x0 + W*z0)
	  // C = (x0**2 + z0**2) - (phantom Radius)**2
	  //
	  // ds = (-B +- sqrt(B**2 - 4*A*C))/(2*A)
	  //
	  // dsInCyl = ds_+ - ds_- = 2*sqrt(B**2 - 4*A*C)/(2A) = sqrt(B**2 - 4*A*C)/A;
	  //

	  const double z0 = p.state.Z - isocenter.z;
	  const double A = p.state.U*p.state.U + p.state.W*p.state.W;
	  const double B = 2.0*(p.state.U*p.state.X + p.state.W*z0);
	  const double C = (p.state.X*p.state.X + z0*z0 - 16*16);

	  const double sqrtArg = B*B - 4.0*A*C;
	  if(sqrtArg > 0.0 && false){
	    const double sqrtVal = sqrt(sqrtArg);
	    const double dsInCyl = sqrtVal/A;

	    //Get gamma range at the particle energy
	    const double mu = 1.0/context.range(p.state.E,PEN_PHOTON,0);

	    profile[ix] += p.state.WGHT*exp(-mu*dsInCyl);
	  }else{
	    profile[ix] += p.state.WGHT;
	  }
	}

	/*
	//Set the origin to the central value
	double maxProfValue = *std::max_element(profile.cbegin(), profile.cend());
	double meanDiff = 0.0;
	for(size_t i = 0; i < profile.size(); ++i){
	  profile[i] = profile[i]/maxProfValue - 1.0;
	  meanDiff += std::fabs(profile[i]);
	}

	meanDiff /= (static_cast<double>(profBins));
	*/
	
	if(verbose > 1){
	  //printf("\n * Mean difference: %.4f%%\n", meanDiff*100.0);

	  //Write the profile to a file
	  std::string profFilename = "profile_" + std::to_string(ibowtie) + ".dat";
	  FILE* fprof = fopen(profFilename.c_str(), "w");
	  if(fprof != nullptr){
	    fprintf(fprof, "# Estimated profile in iteration %u\n", ibowtie);
	    fprintf(fprof, "# Xlow (cm) part/(hist*cm**2) \n");

	    for(size_t i = 0; i < profile.size(); ++i){
	      fprintf(fprof, " %15.5E     %15.5E\n",
		      static_cast<double>(i)*profileDx, profile[i]);
	    }
	    fclose(fprof);
	  }

	  //Write the bowtie to a file
	  std::string bowtieFilename = "bowtie_" + std::to_string(ibowtie) + ".dat";
	  FILE* fbowtie = fopen(bowtieFilename.c_str(), "w");
	  if(fbowtie != nullptr){
	    fprintf(fbowtie, "# Bowtie profile in iteration %u\n", ibowtie);
	    fprintf(fbowtie, "# Xlow (cm) part/(hist*cm**2) \n");

	    for(size_t i = 0; i < profile.size(); ++i){
	      fprintf(fbowtie, " %15.5E     %15.5E\n",
		      static_cast<double>(i)*bowtieDx, bowtieWidths[i]);
	    }
	    fclose(fbowtie);
	  }
	}

	printf("CACA0\n");
	return 0;

	/*
	if(fabs(meanDiff) < 0.02){
	  //Stop the iteration
	  if(verbose > 1){
	    printf("     Agreement achieved\n");
	  }
	  break;
	}
	*/

	//Create bowtie widths displacement vector to be applied
	std::vector<double> bowtieIncrements(reader.bowtieSegments,0.0);

	for(size_t i = 0; i < profBins; ++i){

	  //Get number of influence segments
	  unsigned nSegments = profileInfluenceMap[i].nSegments;

	  //Iterate over bowtie segments with influence in this bin
	  if(nSegments == 1){
	    bowtieIncrements[profileInfluenceMap[i].first] += profile[i]*0.1;
	  }else{

	    double iNSegments = 1.0/static_cast<double>(nSegments);
	    
	    //Apply the modification to the first and last segments
	    bowtieIncrements[profileInfluenceMap[i].first] += profile[i]*profileInfluenceMap[i].fracFirst*iNSegments;
	    bowtieIncrements[profileInfluenceMap[i].last]  += profile[i]*profileInfluenceMap[i].fracLast*iNSegments;
	    
	    for(size_t j = 1; j < profileInfluenceMap[i].nSegments-1; ++j){
	      size_t index = profileInfluenceMap[i].first + j;
	      bowtieIncrements[index] += profile[i];//*iNSegments;
	    }
	  }
	}

	//Apply bowtie increments
	for(size_t i = 0; i < bowtieWidths.size(); ++i){
	  bowtieWidths[i] += bowtieIncrements[i]/static_cast<double>(profile.size())*0.1;
	  if(bowtieWidths[i] < reader.bowtieMinW)
	    bowtieWidths[i] = reader.bowtieMinW;
	  if(bowtieWidths[i] > reader.bowtieMaxW)
	    bowtieWidths[i] = reader.bowtieMaxW;
	}
	
      }
      
      if(verbose > 1){

	//Print filtered spectrum
	printf("   - Printing bowtie filtered spectrum to"
	       " 'x-ray-bowtie-spectrum.dat'.\n");
	
	if(detectedPart::printSpectrum(actualParticles,
				       "x-ray-bowtie-spectrum.dat",
				       200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-bowtie-spectrum.dat'\n");
	  return errors::ERROR_UNABLE_TO_OPEN_FILE;
	}

	printf("   - Printing bowtie filtered spatial distribution to"
	       " 'x-ray-bowtie-spatial.dat'.\n");
	
	if(detectedPart::printSpatialXY(actualParticles,
					"x-ray-bowtie-spatial.dat",
					200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-bowtie-spatial.dat'\n");
	  return errors::ERROR_UNABLE_TO_OPEN_FILE;
	}	
      }
      
      // ** Simulate the phantom
      //-------------------------------------

      err = simCylPhantom("phantom.mat",
			  reader.minE,
			  isocenter,
			  32.0,
			  reader.source2det - reader.source2isocenter,
			  reader.detectorDx,
			  reader.detectorDy,
			  actualParticles,
			  true,
			  verbose,
			  reader.nThreads,
			  "phantom.geo");

      if(err != errors::SUCCESS){
	if(verbose > 0){
	  printf("generateBeamCT: Error on phantom simulation. "
		 "Error code: %d\n", err);
	}
	return err;
      }
      
      if(verbose > 1){

	//Print filtered spectrum
	printf("   - Printing detected spectrum to"
	       " 'x-ray-detected-spectrum.dat'.\n");
	
	if(detectedPart::printSpectrum(actualParticles,
				       "x-ray-detected-spectrum.dat",
				       200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-detected-spectrum.dat'\n");
	  return errors::ERROR_UNABLE_TO_OPEN_FILE;
	}

	printf("   - Printing detected spatial distribution to"
	       " 'x-ray-detected-spatial.dat'.\n");
	
	if(detectedPart::printSpatialXY(actualParticles,
					"x-ray-detected-spatial.dat",
					200, nHists) != 0){
	  printf("Error: Unable to open file 'x-ray-detected-spatial.dat'\n");
	  return errors::ERROR_UNABLE_TO_OPEN_FILE;
	}	
      }

      return errors::SUCCESS;
      
    }
    
  };
};
