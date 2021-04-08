
#
#
#    Copyright (C) 2019-2021 Universitat de València - UV
#    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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
#        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
#        vicente.gimenez@uv.es (Vicente Giménez Gómez)
#    
#



#############################
#  UTILITIES
#############################
include(CMakeDependentOption)

#         OPTIONS           #
#############################

CMAKE_DEPENDENT_OPTION(READ_MAT  "Compile read material utility" ON 
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(MAT_RANGE "Compile material range utility" OFF
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(TALLIES_REGISTER "Compile tallies register" ON
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(SAMPLE_REGISTER "Compile samplers register" ON
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(GEOMETRY_REGISTER "Compile geometries register" ON
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(PHOTON_VR_REGISTER "Compile photon VR register" ON
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(GENERIC_VR_REGISTER "Compile generic VR register" ON
                       "BUILD_UTILITIES" OFF)                       
CMAKE_DEPENDENT_OPTION(PSF_SPECTRE "Compile PSF spectrum generator" ON
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(PSF_2_ASCII "Compile PSF to ASCII conversor" ON
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(GEO_2_VOXEL "Compile conversor to voxelized geometry" ON
                       "BUILD_UTILITIES" OFF)
CMAKE_DEPENDENT_OPTION(LB_SERVER "Compile RUPER_LB basic server program" ON
                       "BUILD_UTILITIES;WITH_LB" OFF)
                       
CMAKE_DEPENDENT_OPTION(DICOM_KERNELS "Compile DICOM kernels utility" ON
                       "BUILD_UTILITIES;WITH_DICOM" OFF)
                       
CMAKE_DEPENDENT_OPTION(MUTREN "Compile mutren" ON
                       "BUILD_UTILITIES" OFF)                       
                       
#       EXECUTABLES         #
#############################

#  Error list generator
#############################
add_executable(createErrors ${PROJECT_SOURCE_DIR}/utilities/errors/createErrors.cpp)  

#-----------------------------------------------------------------------------
#  Read material
#############################

if(READ_MAT)
    # Create read material executable
    add_executable(readMat ${PROJECT_SOURCE_DIR}/utilities/penMats/readMat.cpp)

    #target libs
    target_link_libraries(readMat pen_materials)

    # Install readMat executable
    install(
    TARGETS readMat
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/penMats
    )
endif(READ_MAT)

#-----------------------------------------------------------------------------
#  Pen material ranges
#############################

if(MAT_RANGE)
    # Create read material executable
    add_executable(range ${PROJECT_SOURCE_DIR}/utilities/penMats/range.cpp)

    #target libs
    target_link_libraries(range pen_particles)

    # Install readMat executable
    install(
    TARGETS range
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/penMats
    )
endif(MAT_RANGE)
  
#-----------------------------------------------------------------------------
#  psampler
#############################

# Create psampler executable
add_executable(psampler ${PROJECT_SOURCE_DIR}/utilities/psampler/psampler.cpp)

#target libs
target_link_libraries(psampler pen_particles)
target_link_libraries(psampler pen_timer)

# Install psampler executable
install(
  TARGETS psampler
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/psampler
  )
  
#-----------------------------------------------------------------------------
#  registered tallies
#############################

if(TALLIES_REGISTER)
    # Create registered tallies executable
    add_executable(regTallies ${PROJECT_SOURCE_DIR}/utilities/registers/registeredTallies.cpp)

    #target libs
    target_link_libraries(regTallies tallies)

    # Install registered tallies executable
    install(
    TARGETS regTallies
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/registers
    )
endif(TALLIES_REGISTER)
#-----------------------------------------------------------------------------
#  registered samplers
#############################

if(SAMPLE_REGISTER)
    # Create registered tallies executable
    add_executable(regSamplers ${PROJECT_SOURCE_DIR}/utilities/registers/registeredSamplers.cpp)

    #target libs
    target_link_libraries(regSamplers samplers)

    # Install registered samplers executable
    install(
    TARGETS regSamplers
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/registers
    )
endif(SAMPLE_REGISTER)
#-----------------------------------------------------------------------------
#  registered geometries
#############################

if(GEOMETRY_REGISTER)
    # Create registered tallies executable
    add_executable(regGeometries ${PROJECT_SOURCE_DIR}/utilities/registers/registeredGeometries.cpp)

    #target libs
    target_link_libraries(regGeometries pen_geos)


    # Install registered geometries executable
    install(
    TARGETS regGeometries
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/registers
    )
endif(GEOMETRY_REGISTER)

#-----------------------------------------------------------------------------
#  registered photon VR
#############################

if(PHOTON_VR_REGISTER)
    # Create registered tallies executable
    add_executable(regPhotonVR ${PROJECT_SOURCE_DIR}/utilities/registers/registeredPhotonVR.cpp)

    #target libs
    target_link_libraries(regPhotonVR pen_VR)

    # Install registered tallies executable
    install(
    TARGETS regPhotonVR
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/registers
    )
endif(PHOTON_VR_REGISTER)

#-----------------------------------------------------------------------------
#  registered generic VR
#############################

if(GENERIC_VR_REGISTER)
    # Create registered tallies executable
    add_executable(regGenericVR ${PROJECT_SOURCE_DIR}/utilities/registers/registeredGenericVR.cpp)

    #target libs
    target_link_libraries(regGenericVR pen_VR)

    # Install registered tallies executable
    install(
    TARGETS regGenericVR
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/registers
    )
endif(GENERIC_VR_REGISTER)

#-----------------------------------------------------------------------------
#  PSF spectrum
#############################

if(PSF_SPECTRE)
    # Create PSF spectrum executable
    add_executable(psfSpect ${PROJECT_SOURCE_DIR}/utilities/psf/psf_spectre.cpp)

    #target libs
    target_link_libraries(psfSpect pen_psf)


    # Install registered geometries executable
    install(
    TARGETS psfSpect
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/psf
    )
endif(PSF_SPECTRE)
#-----------------------------------------------------------------------------
#  PSF 2 ASCII
#############################

if(PSF_2_ASCII)
    # Create PSF spectrum executable
    add_executable(psf2ascii ${PROJECT_SOURCE_DIR}/utilities/psf/psf2ascii.cpp)

    #target libs
    target_link_libraries(psf2ascii pen_psf)


    # Install registered geometries executable
    install(
    TARGETS psf2ascii
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/psf
    )
endif(PSF_2_ASCII)
#-----------------------------------------------------------------------------
#  geo to voxel
#############################

if(GEO_2_VOXEL)
    # Create registered tallies executable
    add_executable(geo2voxel ${PROJECT_SOURCE_DIR}/utilities/geometry/geo2voxel.cpp)

    #target libs
    target_link_libraries(geo2voxel pen_geos)


    # Install registered geometries executable
    install(
    TARGETS geo2voxel
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/geometry/voxels
    )
endif(GEO_2_VOXEL)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#  DICOM kernels
#############################

if(DICOM_KERNELS AND WITH_DICOM)

    add_executable(dicomKernels ${PROJECT_SOURCE_DIR}/utilities/DICOM/kernels.cpp)

    #target libs
    target_link_libraries(dicomKernels pen_geos)
    
    # Install DICOM kernels program
    install(
    TARGETS dicomKernels
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/DICOM
    )    
    
endif()

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#  LB server
#############################

if(LB_SERVER AND WITH_LB)
    # Create LB server executable
    add_executable(LBserver ${PROJECT_SOURCE_DIR}/utilities/LB/LBserver.cpp)

    #target libs
    target_link_libraries(LBserver pen_LB)

    # Create stdin LB server executable
    add_executable(stdinServer ${PROJECT_SOURCE_DIR}/utilities/LB/stdinServer.cpp)

    #target libs
    target_link_libraries(stdinServer pen_LB)

    # Install load balance executables
    install(
    TARGETS LBserver
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/LB/server
    )
    install(
    TARGETS stdinServer
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/LB/server
    )    
endif(LB_SERVER AND WITH_LB)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#  MUTREN
#############################

if(MUTREN)
    # Create mutren executable
    add_executable(mutren ${PROJECT_SOURCE_DIR}/utilities/mutren/mutren.cpp)

    target_link_libraries(mutren pen_particles)
    target_link_libraries(mutren pen_geos)
    target_link_libraries(mutren pen_timer)

    # Install mutren executable
    install(
    TARGETS mutren
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/mutren/mutren
    )
endif(MUTREN)
#-----------------------------------------------------------------------------
