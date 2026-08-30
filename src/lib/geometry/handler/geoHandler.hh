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

#ifndef __PEN_GEO_HANDLER__
#define __PEN_GEO_HANDLER__

#include <array>
#include <unordered_map>
#include "pen_geoView.hh"

namespace penred{
  namespace geometry{

    class Handler{

    private:
      //Store geometry
      std::shared_ptr<wrapper_geometry> geometry;

      //Store densities for each material
      std::array<double, constants::MAXMAT> densities;

      //Geometry viewer
      Viewer viewer;
      
    public:

      enum errors{
        SUCCESS = 0,
        INVALID_GEO_TYPE,
        GEOMETRY_CONFIG_FAILED,
        GEOMETRY_NOT_CONFIGURED,
        GEOMETRY_CONFIGURED,
        INVALID_DIMENSIONS,
        BAD_VALUE,
        VOID_VOXEL,
        SET_VOXELS_ERROR,
        MISSING_CONFIG_SECTION,
        ERROR_ON_CONFIGURAITON_PARSING,
        EMPTY_GEOMETRY,
        INVALID_FILE,
      };

      static constexpr const char* errorMessage(const int val) noexcept {
        switch(val){
        case SUCCESS: return "Success";
        case INVALID_GEO_TYPE: return "Invalid geometry type";
        case GEOMETRY_CONFIG_FAILED: return "Geometry configuration failed";
        case GEOMETRY_NOT_CONFIGURED: return "No geometry has been configured";
        case GEOMETRY_CONFIGURED: return "Geometry has already been configured";
        case INVALID_DIMENSIONS: return "Invalid dimensions";
        case BAD_VALUE: return "Bad value";
        case VOID_VOXEL: return "Void voxel";
        case SET_VOXELS_ERROR: return "Unable to set voxels in voxel geometry";
        case MISSING_CONFIG_SECTION: return "Missing configuration section";
        case ERROR_ON_CONFIGURAITON_PARSING: return "Error parsing configuration";
        case INVALID_FILE: return "Invalid file";
        default: return "Unknown error";
        }
      }
      
      //Tables used by Marching Cubes Algorithm
      //these tables came from Paul Baurke's web page at 
      //     http://astronomy.swin.edu.au/~pbourke/modelling/polygonise/

      static const int edgeTable[256];
      static const int triTable[256][16];

      static const std::vector<float> palette;

      Handler(){
        densities.fill(1.0);
      }

      inline penred::errors::Error readGeometry(std::shared_ptr<const wrapper_geometry>& shrGeo) const {
        penred::errors::SpecificError<Handler> error;
        
        //Ensure geometry is configured
        if(!isConfigured()){
          error.code = GEOMETRY_NOT_CONFIGURED;
          error.description = "penred:geometry:Handler:readGeometry:Error: Geometry must be "
            "configured previously.";
          return error;
        }
        shrGeo = geometry;
        return error;
      }

      inline void clear() noexcept {
        geometry.reset();
        viewer.clear();
      }

      penred::errors::Error setDensity(const unsigned imat, const double dens) noexcept;
      penred::errors::Error setDensities(const std::array<double, constants::MAXMAT>& arr) noexcept;

      inline bool isConfigured() const noexcept{
        return geometry ? (geometry->configureStatus() ? false : true) : false;
      }

      penred::errors::Error configure(const pen_parserSection& config,
                                      const unsigned verbose,
                                      const std::string& prefix = "");

      penred::errors::Error configure(const std::string& configString,
                                      const unsigned verbose,
                                      const std::string& prefix = "");
      
      penred::errors::Error configFromFile(const char* filename,
                                           const unsigned verbose,
                                           const std::string& prefix = "");

      inline void set3DResolution(const unsigned nx, const unsigned ny,
                                  const float dx, const float dy,
                                  const float perspective){
        viewer.set3DResolution(nx, ny, dx, dy, perspective);
      }

      penred::errors::Error createVoxelized(std::shared_ptr<pen_voxelGeo>& voxelGeo,
                                            const unsigned nx, const unsigned ny, const unsigned nz,
                                            const double dx, const double dy, const double dz,
                                            const double ox, const double oy, const double oz,
                                            const double t, const unsigned granul) const;

      penred::errors::Error createMesh(const unsigned nx, const unsigned ny, const unsigned nz,
                                       const double dx, const double dy, const double dz,
                                       const double ox, const double oy, const double oz,
                                       const double t, const std::string& output = "geo",
                                       const unsigned granul = 3,
                                       const unsigned smoothSteps = 10,
                                       const double smoothFactor = 0.2,
                                       const unsigned verbose = 0) const;

      //Viewer functions
      inline penred::errors::Error renderX(std::vector<unsigned int>& renderMat,
                                           std::vector<unsigned int>& renderBody,
                                           const float x, const float y, const float z,
                                           const float dy, const float dz,
                                           const unsigned ny, const unsigned nz,
                                           const float t,
                                           const unsigned nthreads = 1) const {
        penred::errors::SpecificError<Handler> error;        
        //Ensure geometry is configured
        if(!isConfigured()){
          error.code = GEOMETRY_NOT_CONFIGURED;
          error.description = "penred:geometry:Handler:renderX:Error: Geometry must be "
            "configured previously.";
          return error;
        }
        renderMat.resize(ny*nz);
        renderBody.resize(ny*nz);
        viewer.renderX(renderMat.data(), renderBody.data(), x, y, z, dy, dz, ny, nz, t, nthreads);
        return error;
      }

      inline penred::errors::Error renderY(std::vector<unsigned int>& renderMat,
                                           std::vector<unsigned int>& renderBody,
                                           const float x, const float y, const float z,
                                           const float dx, const float dz,
                                           const unsigned nx, const unsigned nz,
                                           const float t,
                                           const unsigned nthreads = 1) const {
        penred::errors::SpecificError<Handler> error;        
        //Ensure geometry is configured
        if(!isConfigured()){
          error.code = GEOMETRY_NOT_CONFIGURED;
          error.description = "penred:geometry:Handler:renderX:Error: Geometry must be "
            "configured previously.";
          return error;
        }
        renderMat.resize(nx*nz);
        renderBody.resize(nx*nz);
        viewer.renderY(renderMat.data(), renderBody.data(), x, y, z, dx, dz, nx, nz, t, nthreads);
        return error;
      }

      inline penred::errors::Error renderZ(std::vector<unsigned int>& renderMat,
                                           std::vector<unsigned int>& renderBody,
                                           const float x, const float y, const float z,
                                           const float dx, const float dy,
                                           const unsigned nx, const unsigned ny,
                                           const float t,
                                           const unsigned nthreads = 1) const {
        penred::errors::SpecificError<Handler> error;        
        //Ensure geometry is configured
        if(!isConfigured()){
          error.code = GEOMETRY_NOT_CONFIGURED;
          error.description = "penred:geometry:Handler:renderZ:Error: Geometry must be "
            "configured previously.";
          return error;
        }
        renderMat.resize(nx*ny);
        renderBody.resize(nx*ny);
        viewer.renderZ(renderMat.data(), renderBody.data(), x, y, z, dx, dy, nx, ny, t, nthreads);
        return error;
      }
      
    };
    
  } //namespace geometry

  //Define geometry handler error message function
  template<>
  constexpr const char* errors::errorMessage<geometry::Handler>(const int val) noexcept {
    return geometry::Handler::errorMessage(val);
  }
  
} //namespace penred

#endif
