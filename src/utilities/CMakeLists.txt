
#
#
#    Copyright (C) 2019 Universitat de València - UV
#    Copyright (C) 2019 Universitat Politècnica de València - UPV
#
#    This file is part of PenRed: Parallel Engine for Radiation Energy Deposition.
#
#    PenRed is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PenRed is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with PenRed.  If not, see <https://www.gnu.org/licenses/>. 
#
#
#    contact emails:
#
#        vicent.gimenez.alventosa@gmail.com
#        vicente.gimenez@uv.es
#    
#



#############################
#  UTILITIES
#############################



#       EXECUTABLES         #
#############################


#  Error list generator
#############################
add_executable(createErrors ${PROJECT_SOURCE_DIR}/utilities/errors/createErrors.cpp)


#-----------------------------------------------------------------------------
#  Read material
#############################

# Create read material executable
add_executable(readMat ${PROJECT_SOURCE_DIR}/utilities/penMats/readMat.cpp)

#target libs
target_link_libraries(readMat pen_materials)

# Install readMat executable
install(
  TARGETS readMat
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/penMats
  )


#-----------------------------------------------------------------------------
#  Pen material ranges
#############################

# Create read material executable
add_executable(range ${PROJECT_SOURCE_DIR}/utilities/penMats/range.cpp)

#target libs
target_link_libraries(range pen_materials)

# Install readMat executable
install(
  TARGETS range
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/penMats
  )

#-----------------------------------------------------------------------------
#  registered tallies
#############################

# Create registered tallies executable
add_executable(regTallies ${PROJECT_SOURCE_DIR}/utilities/registers/registeredTallies.cpp)

#target libs
target_link_libraries(regTallies tallies)

# Install registered tallies executable
install(
  TARGETS regTallies
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/registers
  )
#-----------------------------------------------------------------------------
#  registered samplers
#############################

# Create registered tallies executable
add_executable(regSamplers ${PROJECT_SOURCE_DIR}/utilities/registers/registeredSamplers.cpp)

#target libs
target_link_libraries(regSamplers samplers)

# Install registered samplers executable
install(
  TARGETS regSamplers
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/registers
  )
#-----------------------------------------------------------------------------
#  registered geometries
#############################

# Create registered tallies executable
add_executable(regGeometries ${PROJECT_SOURCE_DIR}/utilities/registers/registeredGeometries.cpp)

#target libs
target_link_libraries(regGeometries pen_geos)


# Install registered geometries executable
install(
  TARGETS regGeometries
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/registers
  )
#-----------------------------------------------------------------------------
#  PSF spectrum
#############################

# Create PSF spectrum executable
add_executable(psfSpect ${PROJECT_SOURCE_DIR}/utilities/psf/psf_spectre.cpp)

#target libs
target_link_libraries(psfSpect pen_psf)


# Install registered geometries executable
install(
  TARGETS psfSpect
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/psf
  )
#-----------------------------------------------------------------------------
#  PSF 2 ASCII
#############################

# Create PSF spectrum executable
add_executable(psf2ascii ${PROJECT_SOURCE_DIR}/utilities/psf/psf2ascii.cpp)

#target libs
target_link_libraries(psf2ascii pen_psf)


# Install registered geometries executable
install(
  TARGETS psf2ascii
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/psf
  )
#-----------------------------------------------------------------------------
#  geo to voxel
#############################

# Create registered tallies executable
add_executable(geo2voxel ${PROJECT_SOURCE_DIR}/utilities/geometry/geo2voxel.cpp)

#target libs
target_link_libraries(geo2voxel pen_geos)


# Install registered geometries executable
install(
  TARGETS geo2voxel
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/geometry/voxels
  )
#-----------------------------------------------------------------------------
