//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2025-2026 Vicent Giménez Alventosa
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
                        const double detectorDx,
                        const double detectorDy,
                        const double detectorDz,
                        const double inherentFilterSize,
                        const std::vector<double>& filters,
                        unsigned& initMat,
                        const vector3D<double> detectorPos,
                        const bool constructAnode,
                        const double anodeAngle,
                        const unsigned verbose,
                        FilterGeo* PSFFilterGeo){

      //This function constructs a mesh based geometry of a x-ray device
      //following the specifications provided by the function parameters

      //Calculate source pos
      vector3D<double> sourcePos = detectorPos;
      sourcePos.z += source2det;

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
	if(verbose > 0){
	  penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                   "constructDevice: Error: The minimum distance between "
                                   "inherent filter and first added filter must be %f cm\n"
                                   "    Inherent filter end to source: %f\n"
                                   "    Added filter start to source   : %f\n",
                                   4.0*elementSpacing+2.0*collHeight,
                                   inherentFilter2filters,
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
      const double source2detCollTop = source2filtersEnd + elementSpacing;
      const double source2detCollBot = source2detCollTop + collHeight;

      if(verbose > 1){
        penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                     "\n + X-ray characteristics:\n"
                                     "     Source to detector   : %.4f cm\n"
                                     "     Source to inherent f.: %.4f cm\n"
                                     "     Source to filters    : %.4f cm\n"
                                     "     Detector X size      : %.4f cm\n"
                                     "     Detector Y size      : %.4f cm\n"
                                     "     Detector depth       : %.4f cm\n",
                                     source2det,
                                     source2inherentFilter,
                                     source2filter,
                                     detectorDx,
                                     detectorDy,
                                     detectorDz);
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

      if(PSFFilterGeo != nullptr)
        nBodies += 1;

      out << "# Number of bodies" << std::endl;
      out << " " << nBodies << std::endl;
      out << "#" << std::endl;

      // ** Material index counter
      initMat = std::max(initMat,1u);
      unsigned collMat = initMat; //Collimator material
      unsigned nextMat = initMat+1;
      
      // ** Anode
      //-----------------

      //Check if the anode must be simulated
      if(constructAnode){

	if(verbose > 1){
	  penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                   "\n + Creating anode: \n"
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

      double zorigin = source2filter;
      if(filters.size() > 0){
      
	size_t ifilter = 0;
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
      
      if(PSFFilterGeo){

        //Calculate filter position and size
        vector3D<double> filterCenter = sourcePos;
        std::pair<double, double> PSFFilterFieldSizes;
        
        if(filtersWidth + inherentFilterSize <= 0.0){
          //The PSF/distribution will be recorded just before the inherent filter position

          //Calculate field size and position
          PSFFilterFieldSizes = fieldSize(focalSpot,
                                          detectorDx,
                                          detectorDy,
                                          source2det,
                                          source2inherentFilter);
          filterCenter.z -= source2inherentFilter - 0.05;
        }
        else{
          //The PSF/distribution will be recorded after the last filter
          filterCenter.z -= zorigin + 0.05;
          PSFFilterFieldSizes = filterFieldSizes;
        }
        
        //Create the auxiliary recording filter
        std::string filterName("filter-psf");
        createBaseFilter(1.2*PSFFilterFieldSizes.first,
                         1.2*PSFFilterFieldSizes.second,
                         0.1,
                         1,
                         out,
                         collMat,
                         filterName,
                         "world",
                         false,
                         filterCenter);

        PSFFilterGeo->size.x = 1.2*PSFFilterFieldSizes.first;
        PSFFilterGeo->size.y = 1.2*PSFFilterFieldSizes.first;
        PSFFilterGeo->size.z = 0.1;

        PSFFilterGeo->center = filterCenter;
      }
      
      // ** Detector
      //

      //Create final detector
      vector3D<double> detCenter = sourcePos;
      detCenter.z -= source2det + 0.5;
      createBaseFilter(detectorDx,
		       detectorDy,
		       detectorDz,
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
                       1.2*(source2det + detectorDz),
                       1,
                       out,
                       0,
                       "world",
                       "void",
                       false,
                       worldCenter);

      if(verbose > 1){
        penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                     "\n + X-ray materials:\n"
                                     "     Initial material index   : %u\n"
                                     "     Used                     : %u\n"
                                     "     Next free material index : %u\n",
                                     initMat, nextMat-initMat, nextMat);
      }
      
      initMat = nextMat;
      return errors::SUCCESS;
    }
    
    int constructSimDevice(const pen_parserSection& config,
                           penred::simulation::simulator<pen_context>& simula,
                           const unsigned verbose){

      
      // ** Parse configuration
      
      //Read information from config section
      readerXRayDeviceSimulate reader;
      int err = reader.read(config,verbose);
      if(err != readerXRayDeviceSimulate::SUCCESS){
        return errors::INVALID_CONFIGURATION;
      }

      //Ensure a source is provided
      if(!reader.simAnode && reader.PSFFile.empty() && reader.energyDistribFile.empty()){
        penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                     "simDevice: Error: No source defined. Enable"
                                     " anode simulation, PSF or distribution source.\n");
        return errors::INVALID_CONFIGURATION;
      }

      //If a PSF is created at the detector, force an ideal detector
      if(reader.storeFilteredPSF)
        reader.detectorIdeal = true;

      //Create the simulation configuration
      pen_parserSection simConf;

      // ** Sampling configuration
      
      double maxE = 1.0e6;

      const vector3D<double> sourcePos =
        vector3D<double>(reader.detectorPosition.x,
                         reader.detectorPosition.y,
                         reader.detectorPosition.z + reader.source2det);      

      //Calculate beam radius
      constexpr double pi = 3.141592653589793;
      constexpr double pi05 = pi/2.0;
      constexpr double deg2rad = pi/180.0;
      
      const double anodeAngleRad = deg2rad*reader.anodeAngle;
      const double beamDiameter = reader.focalSpot*tan(pi05-anodeAngleRad);
      const double beamRad = beamDiameter/2.0;


      simConf.set("sources/generic/beam/nhist", static_cast<double>(reader.nHists));
      
      if(reader.simAnode){
	
        const double beamE = reader.kvp*1.0e3;
        maxE = beamE;
        if(maxE > 1.0e6){
          if(verbose > 0){
            penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                         "simDevice: Error: Beam energy must be "
                                         "lesser than 1 MeV for anode simulations.\n");
          }
          return errors::BEAM_ENERGY_TOO_HIGH;
        }

        simConf.set("sources/generic/beam/kpar", "electron");
        
        //Configure a circle spatial source
        simConf.set("sources/generic/beam/spatial/type", "CIRCLE");
        simConf.set("sources/generic/beam/spatial/position/x", sourcePos.x);
        simConf.set("sources/generic/beam/spatial/position/y", sourcePos.y + 10.0);
        simConf.set("sources/generic/beam/spatial/position/z", sourcePos.z);

        simConf.set("sources/generic/beam/spatial/euler/omega",  0.0);
        simConf.set("sources/generic/beam/spatial/euler/theta", 90.0);
        simConf.set("sources/generic/beam/spatial/euler/phi"  ,-90.0);

        simConf.set("sources/generic/beam/spatial/radius", beamRad);

        //Configure a monoenergetic beam
        simConf.set("sources/generic/beam/energy/type", "MONOENERGETIC");
        simConf.set("sources/generic/beam/energy/energy", beamE);

        //Configure direction (-Y)
        simConf.set("sources/generic/beam/direction/type", "SOLID_ANGLE");
        simConf.set("sources/generic/beam/direction/u", 0.0);
        simConf.set("sources/generic/beam/direction/v",-1.0);
        simConf.set("sources/generic/beam/direction/w", 0.0);

        simConf.set("sources/generic/beam/direction/theta0", 0.0);
        simConf.set("sources/generic/beam/direction/theta1", 0.0);
        
        simConf.set("sources/generic/beam/direction/phi0", 0.0);
        simConf.set("sources/generic/beam/direction/dphi", 0.0);
      }
      else{

        if(!reader.energyDistribFile.empty()){
          //Use a distribution

          //Read the distribution
          penred::measurements::results<double, 4> eDistrib;
          std::ifstream fin(reader.energyDistribFile, std::ifstream::in);
          if(!fin){
            if(verbose > 0){
              penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                           "simDevice: Error: Unable to "
                                           "open distribution file '%s'",
                                           reader.energyDistribFile.c_str());
            }
            return errors::ERROR_UNABLE_TO_OPEN_FILE;
          }
  
          int errDis = eDistrib.read(fin);
          if(errDis != 0){
            if(verbose > 0){
              penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                           "simDevice: Error: Unable to "
                                           "read distribution file '%s'.\n  "
                                           "Error code: %d\n",
                                           reader.energyDistribFile.c_str(), errDis);
            }
            return errors::ERROR_INVALID_FILE;
          }

          //Get energy limits
          std::pair<double, double> elimits = eDistrib.readLimits()[0];
          maxE = elimits.second;

          //Energy 4D distribution
          simConf.set("sources/generic/beam/specific/type", "COMBINED_DISTRIBUTIONS");
          simConf.set("sources/generic/beam/specific/distribution/filename", reader.energyDistribFile);
          simConf.set("sources/generic/beam/specific/distribution/origin/x", reader.detectorPosition.x);
          simConf.set("sources/generic/beam/specific/distribution/origin/y", reader.detectorPosition.y);
          simConf.set("sources/generic/beam/specific/distribution/origin/z",
                      reader.detectorPosition.z + reader.source2det);

          //Spatial sampler using the 2D distribution
          simConf.set("sources/generic/beam/spatial/type", "2D_MEASURE");
          simConf.set("sources/generic/beam/spatial/filename", reader.spatialDistribFile);
          simConf.set("sources/generic/beam/spatial/plane", "xy");
          simConf.set("sources/generic/beam/spatial/constant-coordinate",
                      reader.detectorPosition.z + reader.distrib2det);
        }
        else{
          //Use a PSF
          simConf.set("sources/generic/beam/specific/type", "PSF");
          simConf.set("sources/generic/beam/specific/translation/dx", reader.PSFTrans.x);
          simConf.set("sources/generic/beam/specific/translation/dy", reader.PSFTrans.y);
          simConf.set("sources/generic/beam/specific/translation/dz", reader.PSFTrans.z);

          simConf.set("sources/generic/beam/specific/rotation/omega",
                      reader.PSFRotation.x);
          simConf.set("sources/generic/beam/specific/rotation/theta",
                      reader.PSFRotation.y);
          simConf.set("sources/generic/beam/specific/rotation/phi",
                      reader.PSFRotation.z);
        
          pen_parserArray window;
          window.append(0.0);
          window.append(5.0e-2);
          simConf.set("sources/generic/beam/specific/wght-window", window);
          simConf.set("sources/generic/beam/specific/nsplit", 10);

          simConf.set("sources/generic/beam/specific/Emax", 1.0e6);
          simConf.set("sources/generic/beam/specific/filename", reader.PSFFile);
        }
      }

      // ** Geometry
      
      //Create device geometry
      std::stringstream geoStream;
      unsigned deviceNextMat = 1;
      FilterGeo PSFFilterGeo;
      err = constructDevice(geoStream,
                            reader.focalSpot,
                            reader.source2det,
                            reader.source2filter,
                            reader.detectorDx,
                            reader.detectorDy,
                            reader.detectorDz,
                            reader.inherentFilterWidth,
                            reader.filtersWidth,
                            deviceNextMat,
                            reader.detectorPosition,
                            reader.simAnode,
                            reader.anodeAngle,
                            verbose,
                            (reader.storeFilteredPSF || reader.storeFilteredDistrib) ?
                            &PSFFilterGeo : nullptr);

      if(err != 0){
        if(verbose > 0)
          penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                       "simDevice: Error: Unable to create device geometry\n");
        return errors::ERROR_ON_GEOMETRY_INITIALIZATION;
      }

	  //Print geometry file
	  std::ofstream out("device.msh", std::ofstream::out);
	  out << geoStream.str() << std::endl;
	  out.close();
      
      //Save added geometry configuration
      pen_parserSection geoAddConf;
      config.readSubsection("geometry/config", geoAddConf);
      
      //Try to read geometry type
      int nextkdet = 1;
      bool isCombo;
      if(config.read("geometry/config/type", reader.addedGeoType) == INTDATA_SUCCESS){
        //Combined geometry
        isCombo = true;
        simConf.set("geometry/type", "COMBO");
        simConf.set("geometry/geometries/device/priority", 0);
        simConf.set("geometry/geometries/added/priority", 1);

        simConf.set("geometry/geometries/device/config/type", "MESH_BODY");
        simConf.set("geometry/geometries/device/config/input-file", "device.msh");
    
        //Set detectors (kdet)
        simConf.set("geometry/geometries/device/config/kdet/detector", nextkdet++);

        if(reader.inherentFilterWidth > 0.0)
          simConf.set("geometry/geometries/device/config/kdet/inherent-filter", nextkdet++);
        for(unsigned i = 0; i < reader.filtersZ.size(); ++i){
          std::string filterName("filter-");
          filterName += std::to_string(i);
          std::string kdetKey = "geometry/geometries/device/config/kdet/";
          kdetKey += filterName;
          simConf.set(kdetKey, nextkdet++);
        }
        
        //Set PSF filter kdet, if needed
        if(reader.storeFilteredPSF || reader.storeFilteredDistrib){
          simConf.set("geometry/geometries/device/config/kdet/filter-psf", nextkdet++);
        }
        
        simConf.addSubsection("geometry/geometries/added/config", geoAddConf);
        simConf.set("geometry/geometries/added/config/material-shift", int(deviceNextMat)-1);
        
        if(reader.simAnode){
          simConf.set("geometry/geometries/device/config/dsmax/anode", 2.0e-2);
        }
      }else{
        //Device-only geometry
        isCombo = false;
        simConf.set("geometry/type", "MESH_BODY");
        simConf.set("geometry/input-file", "device.msh");
        
        simConf.set("geometry/kdet/detector", nextkdet++);

        if(reader.inherentFilterWidth > 0.0){
          simConf.set("geometry/kdet/inherent-filter", nextkdet++);
        }
        for(unsigned i = 0; i < reader.filtersZ.size(); ++i){
          std::string filterName("filter-");
          filterName += std::to_string(i);
          std::string kdetKey = "geometry/kdet/";
          kdetKey += filterName;
          simConf.set(kdetKey, nextkdet++);
        }

        //Set PSF filter kdet, if needed
        if(reader.storeFilteredPSF || reader.storeFilteredDistrib){
          simConf.set("geometry/kdet/filter-psf", nextkdet++);
        }
        
        if(reader.simAnode){
          simConf.set("geometry/dsmax/anode", 2.0e-2);
        }
      }

      // ** Materials
      int nextMat = 1;

      // Collimators (perfect absorber)

      //Create collimators material file
      std::string errorString;
      err = penred::penMaterialCreator::createMat(82,
                                                  "collimator.mat",
                                                  errorString);

      if(err != 0){
        penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                     "simDevice: Error: Unable to create "
                                     "collimator material: %s\n", errorString.c_str());
        penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                     "IRETRN =%d\n", err);
        return errors::UNABLE_TO_CREATE_MATERIAL;
      }
	
      simConf.set("materials/collimators/number", nextMat++);
      for(unsigned i = 0; i < constants::nParTypes; ++i){
        std::string path("materials/collimators/eabs/");
        path += particleName(i);
        simConf.set(path, 1.0e35);
      }
      simConf.set("materials/collimators/filename", "collimator.mat");

      // Anode
      if(reader.simAnode){

        //Create anode material file
        err = penred::penMaterialCreator::createMat(reader.anodeZ,
                                                    "anode.mat",
                                                    errorString);
        if(err != 0){
          penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                       "simDevice: Error: Unable to create "
                                       "anode material: %s\n", errorString.c_str());
          penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                       "IRETRN =%d\n", err);
          return errors::UNABLE_TO_CREATE_MATERIAL;
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
	
      }

      // Inherent filter
      if(reader.inherentFilterWidth > 0.0){

        //Create inherent filter material
        err = penred::penMaterialCreator::createMat(13,
                                                    "inherent.mat",
                                                    errorString);
        if(err != 0){
          penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                       "simDevice: Error: Unable to create "
                                       "inherent filter material: %s\n",
                                       errorString.c_str());
          penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                       "IRETRN =%d\n", err);
          return errors::UNABLE_TO_CREATE_MATERIAL;
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

      // Filters
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
            penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                         "simDevice: Error: Unable to create "
                                         "filter %lu material: %s\n",
                                         static_cast<unsigned long>(i),
                                         errorString.c_str());
            penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                         "IRETRN =%d\n", err);
            return errors::UNABLE_TO_CREATE_MATERIAL;
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

      // Detector
      simConf.set("materials/detector/number", nextMat++);
      if(reader.detectorIdeal){
        for(unsigned i = 0; i < constants::nParTypes; ++i){
          std::string path("materials/detector/eabs/");
          path += particleName(i);
          simConf.set(path, 1.0e35);
        }
        simConf.set("materials/detector/filename", "collimator.mat");
      }
      else{
        double minSize = std::min(reader.detectorDx, std::min(reader.detectorDx, reader.detectorDz));
        for(unsigned i = 0; i < constants::nParTypes; ++i){
          std::string path("materials/detector/range/");
          path += particleName(i);
          simConf.set(path, minSize/10.0);
        }
        simConf.set("materials/detector/filename", reader.detectorMatFile.c_str());        
      }

      // Added geometry
      if(isCombo){
        //Provided geometry
        for(const readerXRayDeviceSimulate::materialData& mat : reader.addedGeoMats){

          //Add configuration for this material
          std::string prefix = "materials/added_" + mat.name + "/";
          simConf.set((prefix + "number").c_str(), mat.index + int(deviceNextMat) - 1);
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

      // ** VR
      
      if(reader.simAnode){
        simConf.set("VR/IForcing/bremss/particle", "electron");
        simConf.set("VR/IForcing/bremss/interaction", BETAe_HARD_BREMSSTRAHLUNG);
        simConf.set("VR/IForcing/bremss/factor", 400);
        simConf.set("VR/IForcing/bremss/min-weight", 0.1);
        simConf.set("VR/IForcing/bremss/max-weight", 2.0);
        if(isCombo)
          simConf.set("VR/IForcing/bremss/bodies/device_anode", true);
        else
          simConf.set("VR/IForcing/bremss/bodies/anode", true);

        simConf.set("VR/IForcing/innerShell/particle", "electron");
        simConf.set("VR/IForcing/innerShell/interaction", BETAe_HARD_INNER_SHELL);
        simConf.set("VR/IForcing/innerShell/factor", 400);
        simConf.set("VR/IForcing/innerShell/min-weight", 0.1);
        simConf.set("VR/IForcing/innerShell/max-weight", 2.0);
        if(isCombo)
          simConf.set("VR/IForcing/innerShell/bodies/device_anode", true);
        else
          simConf.set("VR/IForcing/innerShell/bodies/anode", true);

        simConf.set("VR/bremss/split4/splitting", 4);
        if(isCombo)
          simConf.set("VR/bremss/split4/bodies/device_anode", true);
        else
          simConf.set("VR/bremss/split4/bodies/anode", true);
      }

      // ** Simulation parameters
      simConf.set("simulation/threads", (int)reader.nThreads);
      simConf.set("simulation/max-time", reader.maxTime);
      if(!reader.dump2Read.empty())
        simConf.set("simulation/dump2read", reader.dump2Read);
      if(!reader.dump2Write.empty()){
        simConf.set("simulation/dump2write", reader.dump2Write);
        simConf.set("simulation/finalDump", true);
        if(reader.dumpTime > 0.0)
          simConf.set("simulation/dump-interval", reader.dumpTime);
        else
          simConf.set("simulation/dump-interval", 3600.0);
      }

      // ** Tallies

      bool simHalted = reader.storeFilteredPSF || reader.storeFilteredDistrib;

      //Check if a PSF is requested
      if(reader.storeFilteredPSF || reader.storeDetectedPSF){
        simConf.set("tallies/psf/type", "PSF");
        simConf.set("tallies/psf/detector", reader.storeFilteredPSF ? nextkdet-1 : 1);
        simConf.set("tallies/psf/emin", reader.minEnergy);
        simConf.set("tallies/psf/emax", 1.0e6);        
        simConf.set("tallies/psf/particles/default", false);        
        simConf.set("tallies/psf/particles/gamma", true);        
        if(!reader.outputPrefix.empty())
          simConf.set("tallies/psf/outputdir", reader.outputPrefix);
      }

      //Check if a distribution is requested
      if(reader.storeFilteredDistrib){
        
        simConf.set("tallies/FilteredDistrib/type", "DETECTION_SPATIAL_DISTRIB");
        
        simConf.set("tallies/FilteredDistrib/spatial/xmin",
                    PSFFilterGeo.center.x - PSFFilterGeo.size.x/2.0);
        simConf.set("tallies/FilteredDistrib/spatial/xmax",
                    PSFFilterGeo.center.x + PSFFilterGeo.size.x/2.0);
        simConf.set("tallies/FilteredDistrib/spatial/nx",
                    (int)reader.detBinsX);

        simConf.set("tallies/FilteredDistrib/spatial/ymin",
                    PSFFilterGeo.center.y - PSFFilterGeo.size.y/2.0);
        simConf.set("tallies/FilteredDistrib/spatial/ymax",
                    PSFFilterGeo.center.y + PSFFilterGeo.size.y/2.0);
        simConf.set("tallies/FilteredDistrib/spatial/ny",
                    (int)reader.detBinsY);

        simConf.set("tallies/FilteredDistrib/spatial/zmin",
                    PSFFilterGeo.center.z - PSFFilterGeo.size.z/2.0);
        simConf.set("tallies/FilteredDistrib/spatial/zmax",
                    PSFFilterGeo.center.z + PSFFilterGeo.size.z/2.0);
        simConf.set("tallies/FilteredDistrib/spatial/nz", 1);

        simConf.set("tallies/FilteredDistrib/energy/nbins", (int)reader.eBins);
        simConf.set("tallies/FilteredDistrib/energy/emin", reader.minEnergy);
        simConf.set("tallies/FilteredDistrib/energy/emax", maxE);
        
        simConf.set("tallies/FilteredDistrib/detector", nextkdet-1);
        simConf.set("tallies/FilteredDistrib/particle", "gamma");
        if(!reader.outputPrefix.empty())
          simConf.set("tallies/FilteredDistrib/outputdir", reader.outputPrefix);
      }

      //Check if the simulation is halted before reaching the detector
      if(!simHalted){
        // Spatial energy detector
        simConf.set("tallies/SpatialDetector/type", "DETECTION_SPATIAL_DISTRIB");

        simConf.set("tallies/SpatialDetector/spatial/nx", (int)reader.detBinsX);
        simConf.set("tallies/SpatialDetector/spatial/xmin",
                    reader.detectorPosition.x-reader.detectorDx/2.0);
        simConf.set("tallies/SpatialDetector/spatial/xmax",
                    reader.detectorPosition.x+reader.detectorDx/2.0);

        simConf.set("tallies/SpatialDetector/spatial/ny", (int)reader.detBinsY);
        simConf.set("tallies/SpatialDetector/spatial/ymin",
                    reader.detectorPosition.y-reader.detectorDy/2.0);
        simConf.set("tallies/SpatialDetector/spatial/ymax",
                    reader.detectorPosition.y+reader.detectorDy/2.0);

        simConf.set("tallies/SpatialDetector/spatial/nz", 1);
        simConf.set("tallies/SpatialDetector/spatial/zmin",
                    reader.detectorPosition.z-reader.detectorDz);
        simConf.set("tallies/SpatialDetector/spatial/zmax",
                    reader.detectorPosition.z+0.1);
      
        simConf.set("tallies/SpatialDetector/detector",  1);
        simConf.set("tallies/SpatialDetector/particle",  "gamma");
        if(!reader.outputPrefix.empty())
          simConf.set("tallies/SpatialDetector/outputdir", reader.outputPrefix);

        // Energy Spectrum
        simConf.set("tallies/SpectrumDetector/type", "DETECTION_SPATIAL_DISTRIB");

        simConf.set("tallies/SpectrumDetector/spatial/nx", 1);
        simConf.set("tallies/SpectrumDetector/spatial/xmin",
                    reader.detectorPosition.x-reader.detectorDx/2.0);
        simConf.set("tallies/SpectrumDetector/spatial/xmax",
                    reader.detectorPosition.x+reader.detectorDx/2.0);

        simConf.set("tallies/SpectrumDetector/spatial/ny", 1);
        simConf.set("tallies/SpectrumDetector/spatial/ymin",
                    reader.detectorPosition.y-reader.detectorDy/2.0);
        simConf.set("tallies/SpectrumDetector/spatial/ymax",
                    reader.detectorPosition.y+reader.detectorDy/2.0);

        simConf.set("tallies/SpectrumDetector/spatial/nz", 1);
        simConf.set("tallies/SpectrumDetector/spatial/zmin",
                    reader.detectorPosition.z-reader.detectorDz);
        simConf.set("tallies/SpectrumDetector/spatial/zmax",
                    reader.detectorPosition.z+0.1);
      
        simConf.set("tallies/SpectrumDetector/energy/nbins", (int)reader.eBins);
        simConf.set("tallies/SpectrumDetector/energy/emin", reader.minEnergy);
        simConf.set("tallies/SpectrumDetector/energy/emax", maxE);
        simConf.set("tallies/SpectrumDetector/detector",  1);
        simConf.set("tallies/SpectrumDetector/particle",  "gamma");
        if(!reader.outputPrefix.empty())
          simConf.set("tallies/SpectrumDetector/outputdir", reader.outputPrefix);        
      }

      //Check user-defined tallies
      std::vector<std::string> userTallies;
      config.ls("tallies", userTallies);
      for(const std::string& tallyName : userTallies){
        pen_parserSection tallyConf;
        std::string key = "tallies/" + tallyName;
        config.readSubsection(key, tallyConf);
        
        std::string safeName = "tallies/added_" + tallyName;
        simConf.addSubsection(safeName.c_str(), tallyConf);

        if(!reader.outputPrefix.empty()){
          key = safeName + "/outputdir";
          simConf.set(key, reader.outputPrefix, true);
        }
      }

      // Print configuration
      out.open("device.in", std::ofstream::out);
	  out << simConf.stringify() << std::endl;
	  out.close();
      
      // Configure simulator
      //***********************
      simula.configure(simConf);
      
      return errors::SUCCESS;
    }

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
	else if(pathInSection.compare("distance/detector") == 0){
	  source2det = element;
	}
	else if(pathInSection.compare("detector/dx") == 0){
	  detectorDx = element;
	}
	else if(pathInSection.compare("detector/dy") == 0){
	  detectorDy = element;
	}
	else if(pathInSection.compare("detector/dz") == 0){
	  detectorDz = element;
	}
	else if(pathInSection.compare("inherent-filter/width") == 0){
	  inherentFilterWidth = element;
	}
	else if(pathInSection.compare("distance/filter") == 0){
	  source2filter = element;
	}
	else if(pathInSection.compare("detector/pos/x") == 0){
	  detectorPosition.x = element;
	}    
	else if(pathInSection.compare("detector/pos/y") == 0){
	  detectorPosition.y = element;
	}    
	else if(pathInSection.compare("detector/pos/z") == 0){
	  detectorPosition.z = element;
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

    int readerXRayDeviceCreate::beginArray(const std::string& /*pathInSection*/,
					   const size_t,
					   const unsigned){
      return errors::UNHANDLED;
    }
      
    int readerXRayDeviceCreate::endArray(const unsigned){

      if(family == -1){
	return errors::SUCCESS;
      }
      return errors::UNHANDLED;
    }

    int readerXRayDeviceCreate::storeArrayElement(const std::string& /*pathInSection*/,
                                                  const pen_parserData& /*element*/,
                                                  const size_t,
                                                  const unsigned){
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
	    penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                     "Error: Unknown atomic number: %s\n", name.c_str());
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
	else if(pathInSection.compare("simulation/dump/time") == 0){
	  dumpTime = element;
	}    
	else if(pathInSection.compare("x-ray/focal-spot") == 0){
	  focalSpot = element;
	}
	else if(pathInSection.compare("x-ray/distance/detector") == 0){
	  source2det = element;
	}
	else if(pathInSection.compare("x-ray/detector/dx") == 0){
	  detectorDx = element;
	}
	else if(pathInSection.compare("x-ray/detector/dy") == 0){
	  detectorDy = element;
	}
	else if(pathInSection.compare("x-ray/detector/dz") == 0){
	  detectorDz = element;
	}
	else if(pathInSection.compare("x-ray/detector/ideal") == 0){
	  detectorIdeal = element;
	}
	else if(pathInSection.compare("x-ray/inherent-filter/width") == 0){
	  inherentFilterWidth = element;
	}
	else if(pathInSection.compare("x-ray/distance/filter") == 0){
	  source2filter = element;
	}
	else if(pathInSection.compare("x-ray/anode/angle") == 0){
	  anodeAngle = element;
	}
	else if(pathInSection.compare("x-ray/anode/z") == 0){
	  anodeZ = element;
	}	
	else if(pathInSection.compare("x-ray/kvp") == 0){
	  kvp = element;
	}
	else if(pathInSection.compare("x-ray/source/distribution/distance") == 0){
	  distrib2det = element;
	}
	else if(pathInSection.compare("psf/detected") == 0){
	  storeDetectedPSF = element;
	}
	else if(pathInSection.compare("psf/filtered") == 0){
	  storeFilteredPSF = element;
	}
	else if(pathInSection.compare("distributions/filtered") == 0){
	  storeFilteredDistrib = element;
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
	if(pathInSection.compare("x-ray/source/psf/path") == 0){	  
	  PSFFile = element;
	}
	else if(pathInSection.compare("x-ray/source/distribution/spatial") == 0){
	  spatialDistribFile = element;
	}
	else if(pathInSection.compare("x-ray/source/distribution/energy") == 0){
	  energyDistribFile = element;
	}
	else if(pathInSection.compare("simulation/output-prefix") == 0){
	  outputPrefix = element;
	}
	else if(pathInSection.compare("simulation/dump/read") == 0){
	  dump2Read = element;
	}
	else if(pathInSection.compare("simulation/dump/write") == 0){
	  dump2Write = element;
	}    
	else if(pathInSection.compare("x-ray/detector/material") == 0){
      detectorMatFile = element;
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
	if(pathInSection.compare("x-ray/detector/position") == 0){
	  if(size != 3){
	    if(verbose > 0){
	      penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                       "Error: Bad position vector (x,y,z).\n"
                                       "      Required coordinates: 3\n"
                                       "      provided coordinates: %lu\n",
                                       static_cast<unsigned long>(size));
	    }
	    return errors::BAD_DIMENSIONS;
	  }
	}
	else if(pathInSection.compare("x-ray/source/psf/translation") == 0){
	  if(size != 3){
	    if(verbose > 0){
	      penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                       "Error: Bad PSF translation vector (dx,dy,dz).\n"
                                       "      Required coordinates: 3\n"
                                       "      provided coordinates: %lu\n",
                                       static_cast<unsigned long>(size));
	    }
	    return errors::BAD_DIMENSIONS;
	  }
	}
	else if(pathInSection.compare("x-ray/source/psf/rotation") == 0){
	  if(size != 3){
	    if(verbose > 0){
	      penred::logs::logger::printf(penred::logs::CONFIGURATION,
                                       "Error: Bad PSF rotation angles (Z,Y,Z).\n"
                                       "      Required angles: 3\n"
                                       "      provided angles: %lu\n",
                                       static_cast<unsigned long>(size));
	    }
	    return errors::BAD_DIMENSIONS;
	  }
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
	if(pathInSection.compare("x-ray/detector/position") == 0){
	  if(pos == 0){
	    detectorPosition.x = element;
	  }
	  else if(pos == 1){
	    detectorPosition.y = element;	    
	  }
	  else if(pos == 2){
	    detectorPosition.z = element;	    
	  }
	  else{
	    return errors::UNHANDLED;
	  }
	}
    else if(pathInSection.compare("x-ray/source/psf/translation") == 0){
	  if(pos == 0){
	    PSFTrans.x = element;
	  }
	  else if(pos == 1){
	    PSFTrans.y = element;	    
	  }
	  else if(pos == 2){
	    PSFTrans.z = element;	    
	  }
	  else{
	    return errors::UNHANDLED;
	  }
	}
    else if(pathInSection.compare("x-ray/source/psf/rotation") == 0){
	  if(pos == 0){
	    PSFRotation.x = element;
	  }
	  else if(pos == 1){
	    PSFRotation.y = element;	    
	  }
	  else if(pos == 2){
	    PSFRotation.z = element;	    
	  }
	  else{
	    return errors::UNHANDLED;
	  }
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
    
    
  } // namespace xray
} // namespace penred
