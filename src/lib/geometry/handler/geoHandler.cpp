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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Olver Gil)
//

#include "geoHandler.hh"

namespace penred{
  namespace geometry{

    penred::errors::Error Handler::setDensity(const unsigned imat, const double dens) noexcept {

      penred::errors::SpecificError<Handler> error;

      if(isConfigured()){
        error.code = GEOMETRY_CONFIGURED;
        error.description = "penred:geometry:Handler:setDensity:Error: Geometry has already been "
          "configured. Changing densities after configuration is not possible.";
        return error;
      }
      
      if(imat < constants::MAXMAT && dens > 0.0)
        densities[imat] = dens;
      else{
        error.code = BAD_VALUE;
        error.description = "penred:geometry:Handler:setDensity:Error: Material index must be lesser "
          "than " + std::to_string(constants::MAXMAT) + " and density greater than zero.";
      }
      return error;
    }
    penred::errors::Error Handler::setDensities(const std::array<double, constants::MAXMAT>& arr) noexcept {

      penred::errors::SpecificError<Handler> error;

      if(isConfigured()){
        error.code = GEOMETRY_CONFIGURED;
        error.description = "penred:geometry:Handler:setDensities:Error: Geometry has already been "
          "configured. Changing densities after configuration is not possible.";
        return error;
      }
      
      for(const double dens : arr){
        if(dens < 0.0){
          error.code = BAD_VALUE;
          error.description = "penred:geometry:Handler:setDensities:Error: Material density must be "
            " greater than zero.";
          return error;
        }
      }
      densities = arr;
      return error;
    }

    penred::errors::Error Handler::configure(const pen_parserSection& config,
                                             const unsigned verbose,
                                             const std::string& prefix){

      //Clear previous geometry
      clear();

      penred::errors::SpecificError<Handler> error;

      //Check if some prefix has been provided
      if(!prefix.empty()){
        //Get geometry section
        int err;
        pen_parserSection geometrySection;
        err = config.readSubsection(prefix.c_str(),geometrySection);
        if(err != INTDATA_SUCCESS){
          error.code = MISSING_CONFIG_SECTION;
          error.description = "penred:geometry:Handler:configure:Error: Unable to read section '" + prefix +
            "' from the provided configuration";
          return error;
        }
        return configure(geometrySection, verbose);
      }

      //Get geometry type
      std::string geoType;
      if(config.read("type",geoType) != INTDATA_SUCCESS){
        error.code = INVALID_GEO_TYPE;
        error.description = "penred:geometry:Handler:configure:Error: No geometry type provided";
        return error;
      }
	
      //Instantiate the specified geometry type
      geometry = std::shared_ptr<wrapper_geometry>(penGeoRegister_create(geoType.c_str()));
      if(!geometry){
        error.code = INVALID_GEO_TYPE;
        error.description = "penred:geometry:Handler:configure:Error: Unknown geometry type '" +
          geoType + "'";
        return error;
      }

      //Create material information section
      pen_parserSection matInfo;
      for(unsigned imat = 0; imat < densities.size(); imat++){
        char key[400];
        sprintf(key,"mat%03d/ID",imat+1);
        matInfo.set(key,(int)imat+1);
        sprintf(key,"mat%03d/density",imat+1);
        matInfo.set(key,densities[imat]);
      }

      //Construct complete configuration section
      pen_parserSection completeConf = config;
      //Ensure the material section is not already set
      completeConf.remove("materials");
      //Add materials info
      completeConf.addSubsection("materials", matInfo);
	
      //Configure the geometry
      penred::errors::Error geoError = geometry->configure(completeConf, verbose);
      if(geoError){
        error.code = GEOMETRY_CONFIG_FAILED;
        error.description = "penred:geometry:Handler:configure: Error configuring geometry.";
        error.setTrace(geoError);
        return error;
      }

      //Set geometry to viewer
      viewer.init(geometry, 0);
	  
      return error;
    }

    penred::errors::Error Handler::configure(const std::string& configString,
                                             const unsigned verbose,
                                             const std::string& prefix){

      penred::errors::SpecificError<Handler> error;

      //Parse configuration file
      pen_parserSection config;
      std::string errorLine;
      unsigned long errorLineNum;
      int err = parseString(configString,config,errorLine,errorLineNum);
  
      if(err != INTDATA_SUCCESS){
        error.code = ERROR_ON_CONFIGURAITON_PARSING;
        error.description = "  Error code: " + std::to_string(err) +  "\n";
        error.description += "  Error message: ";
        error.description += pen_parserError(err);
        error.description += "\n";
        error.description += "  Error located at line " + std::to_string(errorLineNum);
        error.description += ", at text: ";
        error.description += errorLine + "\n";
        return error;
      }
      
      return configure(config, verbose, prefix);
    }
    
    penred::errors::Error Handler::configFromFile(const char* filename,
                                                  const unsigned verbose,
                                                  const std::string& prefix){

      penred::errors::SpecificError<Handler> error;

      //Parse configuration file
      pen_parserSection config;
      std::string errorLine;
      unsigned long errorLineNum;
      int err = parseFile(filename,config,errorLine,errorLineNum);
  
      if(err != INTDATA_SUCCESS){
        error.code = ERROR_ON_CONFIGURAITON_PARSING;
        error.description = "  Error code: " + std::to_string(err) +  "\n";
        error.description += "  Error message: ";
        error.description += pen_parserError(err);
        error.description += "\n";
        error.description += "  Error located at line " + std::to_string(errorLineNum);
        error.description += ", at text: ";
        error.description += errorLine + "\n";
        return error;
      }
      
      return configure(config, verbose, prefix);
    }
    
    penred::errors::Error Handler::createVoxelized(std::shared_ptr<pen_voxelGeo>& voxelGeo,
                                                   const unsigned nx, const unsigned ny, const unsigned nz,
                                                   const double dx, const double dy, const double dz,
                                                   const double ox, const double oy, const double oz,
                                                   const double time, const unsigned granul) const {

      penred::errors::SpecificError<Handler> error;
      if(!isConfigured()){
        error.code = GEOMETRY_NOT_CONFIGURED;
        error.description = "penred:geometry:Handler:createVoxelized:Error: Geometry must be "
          "configured previously.";
        return error;
      }
      
      if(time < 0){
        error.code = BAD_VALUE;
        error.description = "penred:geometry:Handler:createVoxelized:Error: Time cannot be negative.";
        return error;        
      }
      
      if(nx == 0 || ny == 0 || nz == 0){
        error.code = INVALID_DIMENSIONS;
        error.description = "penred:geometry:Handler:createVoxelized:Error: Number of voxels"
          " must be greater than zero in every axis. Provided: (" + std::to_string(nx) +
          "," + std::to_string(ny) + "," + std::to_string(nz) + ").";
        return error;
      }

      if(dx <= 0.0 || dy <= 0.0 || dz <= 0.0){
        error.code = INVALID_DIMENSIONS;
        error.description = "penred:geometry:Handler:createVoxelized:Error: Voxels' dimensions must "
          "be greater than zero in every axis. Provided: (" + std::to_string(dx) + "," +
          std::to_string(dy) + "," + std::to_string(dz) + ").";
        return error;
      }

      if(granul <= 0){
        error.code = BAD_VALUE;
        error.description = "penred:geometry:Handler:createVoxelized:Error: Granularity must be greater"
          " than zero. Provided value: " + std::to_string(granul);
        return error;
      }

      //Create voxelized geometry
      const size_t tvox = nx*ny*nz;
      std::vector<unsigned> voxMats(tvox);
      std::vector<double> voxDensFact(tvox);

      unsigned nvox[3] = {nx,ny,nz};
      double sizes[3] = {dx,dy,dz};

      //Fill voxels
      const double dgranul = static_cast<unsigned>(granul);
      double subSteps[3] = {dx/dgranul,dy/dgranul,dz/dgranul};

      size_t index = 0;
      const size_t nsubPoints = granul*granul*granul;
      pen_particleState state;
      state.PAGE = time;
      for(unsigned k = 0; k < nz; ++k){
        for(unsigned j = 0; j < ny; ++j){
          for(unsigned i = 0; i < nx; ++i){

            unsigned pointsPerMat[constants::MAXMAT];
            for(unsigned imat = 0; imat < constants::MAXMAT; ++imat)
              pointsPerMat[imat] = 0;
	
            double meanDens = 0.0;

            //Set initial state position
            double xInit = static_cast<double>(i)*dx+subSteps[0]*0.5+ox;
            double yInit = static_cast<double>(j)*dy+subSteps[1]*0.5+oy;
            double zInit = static_cast<double>(k)*dz+subSteps[2]*0.5+oz;
            state.X = xInit;
            state.Y = yInit;
            state.Z = zInit;
            for(unsigned k2 = 0; k2 < granul; ++k2){
              state.Z = zInit + static_cast<double>(k2)*subSteps[2];
              for(unsigned j2 = 0; j2 < granul; ++j2){
                state.Y = yInit + static_cast<double>(j2)*subSteps[1];
                for(unsigned i2 = 0; i2 < granul; ++i2){
                  state.X = xInit + static_cast<double>(i2)*subSteps[0];

                  //Locate the state
                  geometry->locate(state);
                  ++pointsPerMat[state.MAT];

                  if(state.MAT > 0){ //Ensure non void material
                    meanDens += densities[state.MAT-1];
                  }
                }
              }
            }

            //Set material
            unsigned maxPoints = 0;
            for(unsigned imat = 1; imat < constants::MAXMAT; ++imat){
              //Skyp void, which is a non valid voxel material
              if(pointsPerMat[imat] > maxPoints){
                maxPoints = pointsPerMat[imat];
                voxMats[index] = imat;
              }
            }

            if(voxMats[index] == 0 || meanDens <= 0.0){
              error.code = VOID_VOXEL;
              error.description = "penred:geometry:Handler:createVoxelized:Error: Voxels can't be set "
                "with neither void material nor null density.";
              return error;
            }
	
            voxDensFact[index] = (meanDens/double(nsubPoints))/densities[voxMats[index]-1];
	
            ++index;
          }
        }
      }

      // Initialize voxel geometry
      //******************************
      voxelGeo = std::make_shared<pen_voxelGeo>();
      penred::errors::Error errorVox = voxelGeo->setVoxels(nvox,sizes,voxMats.data(),voxDensFact.data());
      if(error){
        error.code = SET_VOXELS_ERROR;
        error.description = "penred:geometry:Handler:createVoxelized:Error: Unable to ser voxels.";
        error.setTrace(errorVox);
        return error;
      }
      return error;
    }


    penred::errors::Error Handler::createMesh(const unsigned nx, const unsigned ny, const unsigned nz,
                                              const double dx, const double dy, const double dz,
                                              const double ox, const double oy, const double oz,
                                              const double time, const std::string& output,
                                              const unsigned granul, const unsigned smoothSteps,
                                              const double smoothFactor, const unsigned verbose) const{

      //Define auxiliary structures
      struct vertex{
        vector3D<float> v;
        std::vector<unsigned long> triangles;
        vector3D<float> normal;

        vertex(const vector3D<float>& vIn) : v(vIn){}
      };

      struct meshTriangle{
        vector3D<unsigned long> vi;
        vector3D<unsigned long> ei;
        vector3D<float> c;
        vector3D<float> normal;
      };
      
      if(verbose > 0){
        printf("- Creating 3D mesh from geometry\n");
        fflush(stdout);
      }

      penred::errors::SpecificError<Handler> error;

      //Check values
      if(!isConfigured()){
        error.code = GEOMETRY_NOT_CONFIGURED;
        error.description = "penred:geometry:Handler:createMesh:Error: Geometry must be "
          "configured previously.";
        return error;
      }

      if(time < 0){
        error.code = BAD_VALUE;
        error.description = "penred:geometry:Handler:createMesh:Error: Time cannot be negative.";
        return error;        
      }

      if(nx == 0 || ny == 0 || nz == 0){
        error.code = INVALID_DIMENSIONS;
        error.description = "penred:geometry:Handler:createMesh:Error: Number of voxels"
          " must be greater than zero in every axis. Provided: (" + std::to_string(nx) +
          "," + std::to_string(ny) + "," + std::to_string(nz) + ").";
        return error;
      }

      if(dx <= 0.0 || dy <= 0.0 || dz <= 0.0){
        error.code = INVALID_DIMENSIONS;
        error.description = "penred:geometry:Handler:createMesh:Error: Voxels' dimensions must "
          "be greater than zero in every axis. Provided: (" + std::to_string(dx) + "," +
          std::to_string(dy) + "," + std::to_string(dz) + ").";
        return error;
      }
      const float maxd = std::max(dx, std::max(dy,dz));

      if(granul <= 0){
        error.code = BAD_VALUE;
        error.description = "penred:geometry:Handler:createMesh:Error: Granularity must be greater"
          " than zero. Provided value: " + std::to_string(granul);
        return error;
      }
      
      //Get number of bodies
      unsigned long nBodies = geometry->getBodies();
      //For quadric geometries, remove the last one
      if(std::strcmp(geometry->getType(), "PEN_QUADRIC") == 0){
        nBodies = nBodies > 0 ? nBodies-1 : 0;
      }

      if(nBodies == 0){
        error.code = EMPTY_GEOMETRY;
        error.description = "penred:geometry:Handler:createMesh:Error: "
          "Geometry does not contain any body.";
        return error;
      }
  
      //Create containers for body meshes
      std::vector<std::vector<vertex>> bodyVertex(nBodies);
      std::vector<std::vector<meshTriangle>> bodyTriangles(nBodies);
      std::vector<std::unordered_map<size_t, size_t>> bodyVertexMap(nBodies);
      
      auto toGlobal = [=](const unsigned i,
                          const unsigned j,
                          const unsigned k,
                          const unsigned e) -> vector3D<float>{
        vector3D<float> ret;
        switch(e){
          //Front
        case 0:
          ret.x = (static_cast<float>(i)+0.5)*dx + ox;
          ret.y = static_cast<float>(j)*dy + oy;
          ret.z = static_cast<float>(k)*dz + oz;
          break;
        case 1:
          ret.x = static_cast<float>(i+1)*dx + ox;
          ret.y = (static_cast<float>(j)+0.5)*dy + oy;
          ret.z = static_cast<float>(k)*dz + oz;
          break;
        case 2:
          ret.x = (static_cast<float>(i)+0.5)*dx + ox;
          ret.y = static_cast<float>(j+1)*dy + oy;
          ret.z = static_cast<float>(k)*dz + oz;
          break;
        case 3:
          ret.x = static_cast<float>(i)*dx + ox;
          ret.y = (static_cast<float>(j)+0.5)*dy + oy;
          ret.z = static_cast<float>(k)*dz + oz;
          break;

          //Back
        case 4:
          ret.x = (static_cast<float>(i)+0.5)*dx + ox;
          ret.y = static_cast<float>(j)*dy + oy;
          ret.z = static_cast<float>(k+1)*dz + oz;
          break;
        case 5:
          ret.x = static_cast<float>(i+1)*dx + ox;
          ret.y = (static_cast<float>(j)+0.5)*dy + oy;
          ret.z = static_cast<float>(k+1)*dz + oz;
          break;
        case 6:
          ret.x = (static_cast<float>(i)+0.5)*dx + ox;
          ret.y = static_cast<float>(j+1)*dy + oy;
          ret.z = static_cast<float>(k+1)*dz + oz;
          break;
        case 7:
          ret.x = static_cast<float>(i)*dx + ox;
          ret.y = (static_cast<float>(j)+0.5)*dy + oy;
          ret.z = static_cast<float>(k+1)*dz + oz;
          break;

          //Center
        case 8:
          ret.x = static_cast<float>(i)*dx + ox;
          ret.y = static_cast<float>(j)*dy + oy;
          ret.z = (static_cast<float>(k)+0.5)*dz + oz;
          break;
        case 9:
          ret.x = static_cast<float>(i+1)*dx + ox;
          ret.y = static_cast<float>(j)*dy + oy;
          ret.z = (static_cast<float>(k)+0.5)*dz + oz;
          break;
        case 10:
          ret.x = static_cast<float>(i+1)*dx + ox;
          ret.y = static_cast<float>(j+1)*dy + oy;
          ret.z = (static_cast<float>(k)+0.5)*dz + oz;
          break;
        case 11:
          ret.x = static_cast<float>(i)*dx + ox;
          ret.y = static_cast<float>(j+1)*dy + oy;
          ret.z = (static_cast<float>(k)+0.5)*dz + oz;
          break;

        default:
          ret.x = 0.0;
          ret.y = 0.0;
          ret.z = 0.0;          
        }

        return ret;
      };

      auto checkneighbour = [=](std::unordered_map<unsigned long long, unsigned long>& used,
                                const std::vector<vertex>& bodyV,
                                const int i,
                                const int j,
                                const int k,
                                const unsigned e)
        -> std::unordered_map<unsigned long long, unsigned long>::iterator{

        const vector3D<float> v = toGlobal(i,j,k,e);
          
        for(int kk = k-1; kk <= k+1; ++kk){
          if(kk < 0 || kk >= static_cast<int>(nz)) continue;
          for(int jj = j-1; jj <= j+1; ++jj){
            if(jj < 0 || jj >= static_cast<int>(ny)) continue;
            for(int ii = i-1; ii <= i+1; ++ii){
              if(ii < 0 || ii >= static_cast<int>(nx)) continue;
              if(ii == i && jj == j && kk == k) continue; //Avoid central voxel
              const unsigned long long voxGlob =
                static_cast<unsigned long long>(12) * //12 possible vertex per voxel
                (static_cast<unsigned long long>(ii) +
                 static_cast<unsigned long long>(jj) * static_cast<unsigned long long>(nx) +
                 static_cast<unsigned long long>(kk) * static_cast<unsigned long long>(nx*ny));              
              for(unsigned long long ee = 0; ee < 12; ++ee){
                //Calculate vertex global indexes
                unsigned long long viused = voxGlob + ee;
                auto it = used.find(viused);
                if(it != used.end()){
                  if(v.dist(bodyV[it->second].v)/maxd < 1.0e-3){
                    return it;
                  }
                }
              }
            }
          }
        }
        return used.end();
      };
      
      // Create mesh
      //***************************
  
      //Allocate memory for two slices
      std::vector<unsigned long> voxBodies(nx*ny*2);

      //Calculate sub-step sizes
      double subSteps[3] = {dx/(double)granul,dy/(double)granul,dz/(double)granul};

      // ** Set thread variables
      unsigned int nConcurrency = std::max(std::thread::hardware_concurrency(), 2u);
      std::vector<std::thread> threads;
      threads.reserve(nConcurrency);
      std::atomic<unsigned long> atomicCounter(0);
      //****************************

      if(verbose > 0){
        printf("- Creating auxiliary mesh with %u threads\n", nConcurrency);
        fflush(stdout);
      }
      
      //Fill the first plane with multiple threads
      for(size_t ith = 0; ith < nConcurrency; ++ith){

        threads.push_back(std::thread([=,&atomicCounter,&voxBodies]{

          unsigned j = atomicCounter++;
          while(j < ny){
            for(unsigned i = 0; i < nx; ++i){

              if(granul > 1){
                std::array<unsigned, pen_geoconst::NB> pointsPerBody;
                pointsPerBody.fill(0);

                //Set initial state position
                pen_particleState state;
                state.PAGE = time;
                double xInit = static_cast<double>(i)*dx+subSteps[0]*0.5+ox;
                double yInit = static_cast<double>(j)*dy+subSteps[1]*0.5+oy;
                double zInit = subSteps[2]*0.5+oz;
		
                state.X = xInit;
                state.Y = yInit;
                state.Z = zInit;
                for(unsigned k2 = 0; k2 < granul; ++k2){
                  state.Z = zInit + static_cast<double>(k2)*subSteps[2];
                  for(unsigned j2 = 0; j2 < granul; ++j2){
                    state.Y = yInit + static_cast<double>(j2)*subSteps[1];
                    for(unsigned i2 = 0; i2 < granul; ++i2){
                      state.X = xInit + static_cast<double>(i2)*subSteps[0];

                      //Locate the state
                      geometry->locate(state);
                      ++pointsPerBody[state.IBODY];
                    }
                  }
                }

                //Set material
                unsigned maxPoints = 0;
                unsigned long dominantBody = nBodies;
                for(unsigned ibody = 0; ibody < nBodies; ++ibody){
                  if(pointsPerBody[ibody] > maxPoints){
                    maxPoints = pointsPerBody[ibody];
                    dominantBody = ibody;
                  }
                }
                voxBodies[j*nx + i] = dominantBody;
	
              }else{

                //Set auxiliary state position
                pen_particleState state;
                state.PAGE = time;
                state.X = (static_cast<double>(i)+0.5)*dx+ox;
                state.Y = (static_cast<double>(j)+0.5)*dy+oy;
                state.Z = 0.5*dz+oz;
	  
                //Locate the state
                geometry->locate(state);
                voxBodies[j*nx + i] = state.IBODY;
              }
            }
            j = atomicCounter++;
          }
        }));
      }

      //Join threads
      for(size_t ith = 0; ith < nConcurrency; ++ith){
        threads[ith].join();
      }
      
      //Clear thread vector
      threads.clear();

      //Calculate the mesh
      unsigned long* firstPlane = voxBodies.data();  
      unsigned long* lastPlane = &voxBodies[nx*ny];

      for(unsigned k = 1; k < nz; ++k){

        if(verbose > 0){
          printf("    + Processing Z plane (%u/%u)\n", k, nz);
          fflush(stdout);
        }
    
        //Init atomic counter
        atomicCounter = 0;
    
        //Fill the next plane with multiple threads
        for(size_t ith = 0; ith < nConcurrency; ++ith){
      
          threads.push_back(std::thread([=,&atomicCounter]{

            //Get next plane
            unsigned j = atomicCounter++;
            while(j < ny){
	  
              for(unsigned i = 0; i < nx; ++i){

                if(granul > 1){
                  std::array<unsigned, pen_geoconst::NB> pointsPerBody;
                  pointsPerBody.fill(0);

                  //Set initial state position
                  pen_particleState state;
                  state.PAGE = time;
                  double xInit = static_cast<double>(i)*dx+subSteps[0]*0.5+ox;
                  double yInit = static_cast<double>(j)*dy+subSteps[1]*0.5+oy;
                  double zInit = static_cast<double>(k)*dz+subSteps[2]*0.5+oz;
                  state.X = xInit;
                  state.Y = yInit;
                  state.Z = zInit;
                  for(unsigned k2 = 0; k2 < granul; ++k2){
                    state.Z = zInit + static_cast<double>(k2)*subSteps[2];
                    for(unsigned j2 = 0; j2 < granul; ++j2){
                      state.Y = yInit + static_cast<double>(j2)*subSteps[1];
                      for(unsigned i2 = 0; i2 < granul; ++i2){
                        state.X = xInit + static_cast<double>(i2)*subSteps[0];

                        //Locate the state
                        geometry->locate(state);
                        ++pointsPerBody[state.IBODY];
                      }
                    }
                  }

                  //Set material
                  unsigned maxPoints = 0;
                  unsigned long dominantBody = nBodies;
                  for(unsigned ibody = 0; ibody < nBodies; ++ibody){
                    if(pointsPerBody[ibody] > maxPoints){
                      maxPoints = pointsPerBody[ibody];
                      dominantBody = ibody;
                    }
                  }
                  lastPlane[j*nx + i] = dominantBody;	  
                }else{

                  //Set auxiliary state position
                  pen_particleState state;
                  state.PAGE = time;
                  state.X = (static_cast<double>(i)+0.5)*dx+ox;
                  state.Y = (static_cast<double>(j)+0.5)*dy+oy;
                  state.Z = (static_cast<double>(k)+0.5)*dz+oz;
	  
                  //Locate the state
                  geometry->locate(state);
                  lastPlane[j*nx + i] = state.IBODY;
                }
              }
              j = atomicCounter++;
            }

          }));
        }
        
        //Join threads
        for(size_t ith = 0; ith < nConcurrency; ++ith){
          threads[ith].join();
        }
        	
        //Clear thread vector
        threads.clear();
      
        //Set vertex, avoiding borders

        //Init atomic counter
        atomicCounter = 0;

        for(size_t ith = 0; ith < nConcurrency; ++ith){
          
          threads.push_back(std::thread([&]{
                
            //Iterate over geometry bodies
            unsigned ibody = atomicCounter++;
            while(ibody < nBodies){

              //Iterate over y axis
              for(unsigned j = 0; j < ny-1; ++j){

                //Iterate over x axis
                for(unsigned i = 0; i < nx-1; ++i){

                  //Get vertex body values
                  const unsigned long v1 = firstPlane[j*nx + i];
                  const unsigned long v2 = firstPlane[j*nx + i + 1];
                  const unsigned long v3 = firstPlane[(j+1)*nx + i + 1];
                  const unsigned long v4 = firstPlane[(j+1)*nx + i];
                  const unsigned long v5 = lastPlane[j*nx + i];
                  const unsigned long v6 = lastPlane[j*nx + i + 1];
                  const unsigned long v7 = lastPlane[(j+1)*nx + i + 1];
                  const unsigned long v8 = lastPlane[(j+1)*nx + i];
              
              
                  //Get cube index according the body index match
                  int cubeIndex = 0;
                  if(v1 == ibody) cubeIndex |= (1 << 0);
                  if(v2 == ibody) cubeIndex |= (1 << 1);
                  if(v3 == ibody) cubeIndex |= (1 << 2);
                  if(v4 == ibody) cubeIndex |= (1 << 3);
                  if(v5 == ibody) cubeIndex |= (1 << 4);
                  if(v6 == ibody) cubeIndex |= (1 << 5);
                  if(v7 == ibody) cubeIndex |= (1 << 6);
                  if(v8 == ibody) cubeIndex |= (1 << 7);

                  //Check if it is completely inside or outside
                  if(edgeTable[cubeIndex] == 0) continue;
                
	  
                  //Iterate over cube triangles according to its pattern
                  for(int n = 0; triTable[cubeIndex][n] != -1; n+=3){

                    //Save triangle
                    bodyTriangles[ibody].push_back({
                        //vertex index
                        {i,j,k-1},
                        //e indexes
                        {static_cast<unsigned long>(triTable[cubeIndex][n+2]),
                         static_cast<unsigned long>(triTable[cubeIndex][n+1]),
                         static_cast<unsigned long>(triTable[cubeIndex][n])},
                        //centroid
                        {0.0,0.0,0.0}, //Filled after
                        //normal
                        {0.0,0.0,0.0}  //Filled after
                      });
                  }
                }
              }

              ibody = atomicCounter++;
            }
          }));
        }

        //Join threads
        for(size_t ith = 0; ith < nConcurrency; ++ith){
          threads[ith].join();
        }

        threads.clear();
    
        //Swap planes
        unsigned long* aux = firstPlane;
        firstPlane = lastPlane;
        lastPlane = aux;
    
      }
  
      //Post process vertex
      if(verbose > 0){
        printf("- Post-processing bodies\n");
        fflush(stdout);
      }

      //Init atomic counter
      atomicCounter = 0;

      for(size_t ith = 0; ith < nConcurrency; ++ith){

        threads.push_back(std::thread([&]{
                
          //Iterate over geometry bodies
          unsigned long ibody = atomicCounter++;
          while(ibody < nBodies){

            if(verbose > 0){
              printf("    + Body (%lu/%lu)\n", ibody, nBodies);
              fflush(stdout);
            }

            //Get body vertex vector
            std::vector<vertex>& bodyV = bodyVertex[ibody];        
            std::vector<meshTriangle>& bodyT = bodyTriangles[ibody];        

            //Construct final vertex and complete triangles
            std::unordered_map<unsigned long long, unsigned long> used;
            for(unsigned long itr = 0; itr < bodyT.size(); ++itr){

              meshTriangle& t = bodyT[itr];
          
              //Calculate voxel mesh index
              const unsigned long long voxGlobalIndex =
                static_cast<unsigned long long>(12) * //12 possible vertex per voxel
                (static_cast<unsigned long long>(t.vi.x) +
                 static_cast<unsigned long long>(t.vi.y) * static_cast<unsigned long long>(nx) +
                 static_cast<unsigned long long>(t.vi.z) * static_cast<unsigned long long>(nx*ny));
          
              //Calculate vertex global indexes
              vector3D<unsigned long long> globalIndex(voxGlobalIndex + static_cast<unsigned long long>(t.ei.x),
                                                       voxGlobalIndex + static_cast<unsigned long long>(t.ei.y),
                                                       voxGlobalIndex + static_cast<unsigned long long>(t.ei.z));

              //Check if vertex are already included
              vector3D<unsigned long> vLocalIndex;
          
              //0
              vLocalIndex.x = bodyV.size();
              auto it = used.find(globalIndex.x);
              if(it != used.end()) {
                //Already included
                vLocalIndex.x = it->second;
              }else{
                //Not included, check for equivalent vertex
                it = checkneighbour(used, bodyV, t.vi.x, t.vi.y, t.vi.z, t.ei.x);
                if(it != used.end()) {
                  //Equivalent found
                  vLocalIndex.x = it->second;
                }else{
                  //No equivalent found, add it
                  bodyV.push_back(toGlobal(t.vi.x, t.vi.y, t.vi.z, t.ei.x));
                }
              
                used[globalIndex.x] = vLocalIndex.x;
              }

              //1
              vLocalIndex.y = bodyV.size();
              it = used.find(globalIndex.y);
              if(it != used.end()) {
                //Already included
                vLocalIndex.y = it->second;
              }else{
                //Not included, check for equivalent vertex
                it = checkneighbour(used, bodyV, t.vi.x, t.vi.y, t.vi.z, t.ei.y);
                if(it != used.end()) {
                  //Equivalent found
                  vLocalIndex.y = it->second;
                }else{
                  //No equivalent found, add it
                  bodyV.push_back(toGlobal(t.vi.x, t.vi.y, t.vi.z, t.ei.y));
                }
            
                used[globalIndex.y] = vLocalIndex.y;
              }

              //2
              vLocalIndex.z = bodyV.size();
              it = used.find(globalIndex.z);
              if(it != used.end()) {
                //Already included
                vLocalIndex.z = it->second;
              }else{
                //Not included, check for equivalent vertex
                it = checkneighbour(used, bodyV, t.vi.x, t.vi.y, t.vi.z, t.ei.z);
                if(it != used.end()) {
                  //Equivalent found
                  vLocalIndex.z = it->second;
                }else{
                  //No equivalent found, add it
                  bodyV.push_back(toGlobal(t.vi.x, t.vi.y, t.vi.z, t.ei.z));
                }
            
                used[globalIndex.z] = vLocalIndex.z;
              }

              //Reset indexes according to local vertex vector positions
              t.vi = vLocalIndex;

              //Add this triangle to each used vertex
              bodyV[t.vi.x].triangles.push_back(itr);
              bodyV[t.vi.y].triangles.push_back(itr);
              bodyV[t.vi.z].triangles.push_back(itr);

              //Calculate triangle center
              t.c = (bodyV[t.vi.x].v +
                     bodyV[t.vi.y].v +
                     bodyV[t.vi.z].v)/3.0f;
            }
            
            //Smooth mesh
            for(unsigned is = 0; is < smoothSteps; ++is){
          
              for(vertex& v : bodyV){

                //Smooth this vertex
                vector3D<float> midC{0.0f,0.0f,0.0f};
                for(const unsigned it : v.triangles){
                  midC.add(bodyT[it].c);
                }
	
                midC = midC/static_cast<float>(v.triangles.size());
                const vector3D<float> diff = midC-v.v;
                v.v.add(diff*smoothFactor);
              }

              //Recalculate triangle centers
              for(meshTriangle& t : bodyT){
                t.c = (bodyV[t.vi.x].v +
                       bodyV[t.vi.y].v +
                       bodyV[t.vi.z].v)/3.0f;
              }
            }

            //Compute resulting triangle normals
            for(meshTriangle& t : bodyT){

              t.normal = bodyV[t.vi.y].v - bodyV[t.vi.x].v; //e01
              const vector3D<float> e02 = bodyV[t.vi.z].v   - bodyV[t.vi.x].v;
              t.normal.crossProd(e02);
              t.normal.normalize();
            }    

            //Compute resulting vertex normals
            for(vertex& v : bodyV){
              //Calculate the mean normal for this vertex
              vector3D<float> normal = {0.0f,0.0f,0.0f};
              for(const unsigned it : v.triangles){
                normal.add(bodyT[it].normal);
              }
              normal.normalize();
              v.normal = normal;
            }

            ibody = atomicCounter++;
              
          }
        }));
      }

      //Join threads
      for(size_t ith = 0; ith < nConcurrency; ++ith){
        threads[ith].join();
      }

      threads.clear();      
      
      //Print wavefront file
      FILE* fout = nullptr;
      std::string objFilename = output + ".obj";
      fout = fopen(objFilename.c_str(),"w");
      if(fout == nullptr){
        error.code = INVALID_FILE;
        error.description = "penred:geometry:Handler:createMesh:Error: Unable to create "
          "output object file: " + output + ".obj";
        return error;
      }

      FILE* foutMat = nullptr;
      std::string matFilename = output + ".mtl";
      foutMat = fopen(matFilename.c_str(),"w");
      if(foutMat == nullptr){
        fclose(fout);
        error.code = INVALID_FILE;
        error.description = "penred:geometry:Handler:createMesh:Error: Unable to create "
          "output material file: " + output + ".mtl";
        return error;
      }

      //Specify material file
      fprintf(fout,"mtllib %s\n", matFilename.c_str());
      //Get used materials
      bool usedMats[constants::MAXMAT+1];
      geometry->usedMat(usedMats);

      //Fill materials file
      for(unsigned long b = 1; b < constants::MAXMAT+1; ++b){

        if(usedMats[b]){
          unsigned matIndex = b%195;
    
          //Define material for this object
          fprintf(foutMat,"#Material %lu\n",b);
          fprintf(foutMat,"newmtl %lu\n",b);
          fprintf(foutMat,"Ka 1.0 1.0 1.0\n");
          fprintf(foutMat,"Kd %.6f %.6f %.6f\n\n",palette[matIndex*3],palette[matIndex*3+1],palette[matIndex*3+2]);
        }
      }
      fclose(foutMat);

      //Bodies definition file
      unsigned long nextBodyFirstVertex = 0;
      for(unsigned long b = 0; b < nBodies; ++b){
    
        fprintf(fout,"usemtl %u\n",geometry->getMat(b));
        fprintf(fout,"o %s\n",geometry->getBodyName(b).c_str());
        //Vertex coordinates
        for(const vertex& v : bodyVertex[b]){
          fprintf(fout,"v %e %e %e\n", v.v.x, v.v.y, v.v.z);
        }
        //Vertex normals
        for(const vertex& v : bodyVertex[b]){
          fprintf(fout,"vn %e %e %e\n", v.normal.x, v.normal.y, v.normal.z);
        }

        for(const meshTriangle& f : bodyTriangles[b]){
          fprintf(fout,"f %lu//%lu %lu//%lu %lu//%lu\n",
                  f.vi.x+1+nextBodyFirstVertex, f.vi.x+1+nextBodyFirstVertex,
                  f.vi.y+1+nextBodyFirstVertex, f.vi.y+1+nextBodyFirstVertex,
                  f.vi.z+1+nextBodyFirstVertex, f.vi.z+1+nextBodyFirstVertex);
        }
        nextBodyFirstVertex += bodyVertex[b].size();
      }

      fclose(fout);

      return error;
    }

    //Tables used by Marching Cubes Algorithm
    //these tables came from Paul Baurke's web page at 
    //     http://astronomy.swin.edu.au/~pbourke/modelling/polygonise/

    const int Handler::edgeTable[256]={
      0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
      0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
      0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
      0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
      0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
      0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
      0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
      0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
      0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
      0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
      0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
      0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
      0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
      0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
      0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
      0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
      0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
      0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
      0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
      0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
      0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
      0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
      0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
      0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
      0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
      0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
      0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
      0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
      0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
      0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
      0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
      0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

    const int Handler::triTable[256][16] =
      {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
       {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
       {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
       {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
       {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
       {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
       {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
       {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
       {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
       {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
       {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
       {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
       {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
       {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
       {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
       {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
       {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
       {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
       {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
       {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
       {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
       {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
       {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
       {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
       {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
       {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
       {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
       {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
       {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
       {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
       {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
       {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
       {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
       {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
       {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
       {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
       {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
       {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
       {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
       {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
       {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
       {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
       {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
       {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
       {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
       {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
       {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
       {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
       {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
       {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
       {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
       {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
       {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
       {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
       {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
       {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
       {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
       {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
       {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
       {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
       {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
       {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
       {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
       {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
       {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
       {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
       {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
       {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
       {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
       {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
       {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
       {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
       {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
       {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
       {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
       {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
       {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
       {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
       {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
       {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
       {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
       {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
       {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
       {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
       {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
       {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
       {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
       {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
       {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
       {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
       {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
       {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
       {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
       {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
       {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
       {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
       {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
       {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
       {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
       {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
       {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
       {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
       {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
       {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
       {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
       {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
       {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
       {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
       {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
       {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
       {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
       {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
       {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
       {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
       {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
       {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
       {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
       {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
       {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
       {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
       {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
       {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
       {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
       {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
       {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
       {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
       {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
       {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
       {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
       {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
       {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
       {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
       {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
       {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
       {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
       {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
       {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
       {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
       {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
       {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
       {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
       {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
       {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
       {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
       {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
       {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
       {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
       {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
       {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
       {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
       {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
       {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
       {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
       {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
       {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
       {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
       {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
       {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
       {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
       {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
       {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
       {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
       {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
       {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
       {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
       {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
       {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
       {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
       {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
       {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
       {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
       {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
       {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
       {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
       {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
       {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
       {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
       {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
       {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
       {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
       {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

    //Define palette
    const std::vector<float> Handler::palette = {0.337255,0.752941,0.752941,
                                                 0.92549,0.752941,0.752941,
                                                 0.92549,0.337255,0.752941,
                                                 0.92549,0.92549,0.752941,
                                                 0.337255,0.92549,0.752941,
                                                 0.752941,0.92549,0.752941,
                                                 0.752941,0.92549,0.337255,
                                                 0.752941,0.92549,0.92549,
                                                 0.752941,0.337255,0.92549,
                                                 0.752941,0.752941,0.92549,
                                                 0.337255,0.752941,0.92549,
                                                 0.92549,0.752941,0.92549,
                                                 0.92549,0.752941,0.337255,
                                                 0.580392,0.580392,0.580392,
                                                 0.164706,0.580392,0.580392,
                                                 0.752941,0.580392,0.580392,
                                                 0.752941,0.164706,0.580392,
                                                 0.752941,0.752941,0.580392,
                                                 0.164706,0.752941,0.580392,
                                                 0.580392,0.752941,0.580392,
                                                 0.580392,0.752941,0.164706,
                                                 0.580392,0.752941,0.752941,
                                                 0.580392,0.164706,0.752941,
                                                 0.580392,0.580392,0.752941,
                                                 0.164706,0.580392,0.752941,
                                                 0.752941,0.580392,0.752941,
                                                 0.752941,0.580392,0.164706,
                                                 0.407843,0.407843,0.407843,
                                                 0.996078,0.407843,0.407843,
                                                 0.580392,0.407843,0.407843,
                                                 0.580392,0.996078,0.407843,
                                                 0.580392,0.580392,0.407843,
                                                 0.996078,0.580392,0.407843,
                                                 0.407843,0.580392,0.407843,
                                                 0.407843,0.580392,0.996078,
                                                 0.407843,0.580392,0.580392,
                                                 0.407843,0.996078,0.580392,
                                                 0.407843,0.407843,0.580392,
                                                 0.996078,0.407843,0.580392,
                                                 0.580392,0.407843,0.580392,
                                                 0.580392,0.407843,0.996078,
                                                 0.235294,0.235294,0.235294,
                                                 0.823529,0.235294,0.235294,
                                                 0.407843,0.235294,0.235294,
                                                 0.407843,0.823529,0.235294,
                                                 0.407843,0.407843,0.235294,
                                                 0.823529,0.407843,0.235294,
                                                 0.235294,0.407843,0.235294,
                                                 0.235294,0.407843,0.823529,
                                                 0.235294,0.407843,0.407843,
                                                 0.235294,0.823529,0.407843,
                                                 0.235294,0.235294,0.407843,
                                                 0.823529,0.235294,0.407843,
                                                 0.407843,0.235294,0.407843,
                                                 0.407843,0.235294,0.823529,
                                                 0.0627451,0.0627451,0.0627451,
                                                 0.65098,0.0627451,0.0627451,
                                                 0.235294,0.0627451,0.0627451,
                                                 0.235294,0.65098,0.0627451,
                                                 0.235294,0.235294,0.0627451,
                                                 0.65098,0.235294,0.0627451,
                                                 0.0627451,0.235294,0.0627451,
                                                 0.0627451,0.235294,0.65098,
                                                 0.0627451,0.235294,0.235294,
                                                 0.0627451,0.65098,0.235294,
                                                 0.0627451,0.0627451,0.235294,
                                                 0.65098,0.0627451,0.235294,
                                                 0.235294,0.0627451,0.235294,
                                                 0.235294,0.0627451,0.65098,
                                                 0.894118,0.894118,0.894118,
                                                 0.478431,0.894118,0.894118,
                                                 0.0627451,0.894118,0.894118,
                                                 0.0627451,0.478431,0.894118,
                                                 0.0627451,0.0627451,0.894118,
                                                 0.478431,0.0627451,0.894118,
                                                 0.894118,0.0627451,0.894118,
                                                 0.894118,0.0627451,0.478431,
                                                 0.894118,0.0627451,0.0627451,
                                                 0.894118,0.478431,0.0627451,
                                                 0.894118,0.894118,0.0627451,
                                                 0.478431,0.894118,0.0627451,
                                                 0.0627451,0.894118,0.0627451,
                                                 0.0627451,0.894118,0.478431,
                                                 0.721569,0.721569,0.721569,
                                                 0.305882,0.721569,0.721569,
                                                 0.894118,0.721569,0.721569,
                                                 0.894118,0.305882,0.721569,
                                                 0.894118,0.894118,0.721569,
                                                 0.305882,0.894118,0.721569,
                                                 0.721569,0.894118,0.721569,
                                                 0.721569,0.894118,0.305882,
                                                 0.721569,0.894118,0.894118,
                                                 0.721569,0.305882,0.894118,
                                                 0.721569,0.721569,0.894118,
                                                 0.305882,0.721569,0.894118,
                                                 0.894118,0.721569,0.894118,
                                                 0.894118,0.721569,0.305882,
                                                 0.54902,0.54902,0.54902,
                                                 0.133333,0.54902,0.54902,
                                                 0.721569,0.54902,0.54902,
                                                 0.721569,0.133333,0.54902,
                                                 0.721569,0.721569,0.54902,
                                                 0.133333,0.721569,0.54902,
                                                 0.54902,0.721569,0.54902,
                                                 0.54902,0.721569,0.133333,
                                                 0.54902,0.721569,0.721569,
                                                 0.54902,0.133333,0.721569,
                                                 0.54902,0.54902,0.721569,
                                                 0.133333,0.54902,0.721569,
                                                 0.721569,0.54902,0.721569,
                                                 0.721569,0.54902,0.133333,
                                                 0.376471,0.376471,0.376471,
                                                 0.964706,0.376471,0.376471,
                                                 0.54902,0.376471,0.376471,
                                                 0.54902,0.964706,0.376471,
                                                 0.54902,0.54902,0.376471,
                                                 0.964706,0.54902,0.376471,
                                                 0.376471,0.54902,0.376471,
                                                 0.376471,0.54902,0.964706,
                                                 0.376471,0.54902,0.54902,
                                                 0.376471,0.964706,0.54902,
                                                 0.376471,0.376471,0.54902,
                                                 0.964706,0.376471,0.54902,
                                                 0.54902,0.376471,0.54902,
                                                 0.54902,0.376471,0.964706,
                                                 0.203922,0.203922,0.203922,
                                                 0.792157,0.203922,0.203922,
                                                 0.376471,0.203922,0.203922,
                                                 0.376471,0.792157,0.203922,
                                                 0.376471,0.376471,0.203922,
                                                 0.792157,0.376471,0.203922,
                                                 0.203922,0.376471,0.203922,
                                                 0.203922,0.376471,0.792157,
                                                 0.203922,0.376471,0.376471,
                                                 0.203922,0.792157,0.376471,
                                                 0.203922,0.203922,0.376471,
                                                 0.792157,0.203922,0.376471,
                                                 0.376471,0.203922,0.376471,
                                                 0.376471,0.203922,0.792157,
                                                 0.0313725,0.0313725,0.0313725,
                                                 0.619608,0.0313725,0.0313725,
                                                 0.203922,0.0313725,0.0313725,
                                                 0.203922,0.619608,0.0313725,
                                                 0.203922,0.203922,0.0313725,
                                                 0.619608,0.203922,0.0313725,
                                                 0.0313725,0.203922,0.0313725,
                                                 0.0313725,0.203922,0.619608,
                                                 0.0313725,0.203922,0.203922,
                                                 0.0313725,0.619608,0.203922,
                                                 0.0313725,0.0313725,0.203922,
                                                 0.619608,0.0313725,0.203922,
                                                 0.203922,0.0313725,0.203922,
                                                 0.203922,0.0313725,0.619608,
                                                 0.862745,0.862745,0.862745,
                                                 0.447059,0.862745,0.862745,
                                                 0.0313725,0.862745,0.862745,
                                                 0.0313725,0.447059,0.862745,
                                                 0.0313725,0.0313725,0.862745,
                                                 0.447059,0.0313725,0.862745,
                                                 0.862745,0.0313725,0.862745,
                                                 0.862745,0.0313725,0.447059,
                                                 0.862745,0.0313725,0.0313725,
                                                 0.862745,0.447059,0.0313725,
                                                 0.862745,0.862745,0.0313725,
                                                 0.447059,0.862745,0.0313725,
                                                 0.0313725,0.862745,0.0313725,
                                                 0.0313725,0.862745,0.447059,
                                                 0.690196,0.690196,0.690196,
                                                 0.27451,0.690196,0.690196,
                                                 0.862745,0.690196,0.690196,
                                                 0.862745,0.27451,0.690196,
                                                 0.862745,0.862745,0.690196,
                                                 0.27451,0.862745,0.690196,
                                                 0.690196,0.862745,0.690196,
                                                 0.690196,0.862745,0.27451,
                                                 0.690196,0.862745,0.862745,
                                                 0.690196,0.27451,0.862745,
                                                 0.690196,0.690196,0.862745,
                                                 0.27451,0.690196,0.862745,
                                                 0.862745,0.690196,0.862745,
                                                 0.862745,0.690196,0.27451,
                                                 0.517647,0.517647,0.517647,
                                                 0.101961,0.517647,0.517647,
                                                 0.690196,0.517647,0.517647,
                                                 0.690196,0.101961,0.517647,
                                                 0.690196,0.690196,0.517647,
                                                 0.101961,0.690196,0.517647,
                                                 0.517647,0.690196,0.517647,
                                                 0.517647,0.690196,0.101961,
                                                 0.517647,0.690196,0.690196,
                                                 0.517647,0.101961,0.690196,
                                                 0.517647,0.517647,0.690196,
                                                 0.101961,0.517647,0.690196,
                                                 0.690196,0.517647,0.690196,
                                                 0.690196,0.517647,0.101961,
                                                 0.345098,0.345098,0.345098,
                                                 0.933333,0.345098,0.345098,
                                                 0,0,0};    
    
  } //namespace geometry
} //namespace penred
