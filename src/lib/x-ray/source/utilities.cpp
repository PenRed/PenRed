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
    
  } // namespace xray
} // namespace penred
