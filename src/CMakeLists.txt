
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

 
#----------------------------------------------------------------------------
# Check CMake version
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
option(DEVELOPMENT_WARNINGS "Compile with extra warnings for development " OFF)
option(WITH_DICOM "Build PenRed DICOM " OFF)
option(WITH_MULTI_THREADING "Build PenRed with multi-threading " ON)
option(WITH_MPI "Build PenRed with MPI " OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Project name
project(PenRed LANGUAGES CXX)
#-----------------------------------------------------------------------------

# Compilation options
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # using GCC
  if (DEVELOPMENT_WARNINGS)
  	add_compile_options(-std=c++11 -march=native -Wall -Wpedantic -Wextra -Wshadow -O2 -Wno-unused-result -Wno-deprecated-copy -frecord-gcc-switches)
  else (DEVELOPMENT_WARNINGS)
	add_compile_options(-std=c++11 -march=native -Wall -O2 -frecord-gcc-switches)
  endif (DEVELOPMENT_WARNINGS)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
  add_compile_options(-std=c++11 -xHost -O2 -no-prec-div -fp-model precise -prec-sqrt)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # using Visual Studio C++
endif()

#-----------------------------------------------------------------------------

# Check if DICOM option has been enabled
if(WITH_DICOM)
  add_definitions(-D_PEN_USE_DICOM_=1) #enable DICOM   
endif(WITH_DICOM)

# Check if multi-threading option has been enabled
if(WITH_MULTI_THREADING)
  add_definitions(-D_PEN_USE_THREADS_=1) #enable Threads     
endif(WITH_MULTI_THREADING)

# Check if MPI option has been enabled
if(WITH_MPI)
  add_definitions(-D_PEN_USE_MPI_=1) #enable MPI
endif(WITH_MPI)

# Check platform
if(UNIX)
  add_definitions(-D_PEN_UNIX_=1) #enable DICOM
elseif(WIN32)
  
elseif(APPLE)
  
endif(UNIX)

# Generate error list
add_custom_command(
  OUTPUT
  ${PROJECT_SOURCE_DIR}/config/errors/eenum.h
  ${PROJECT_SOURCE_DIR}/config/errors/emess.h
  DEPENDS ${CMAKE_BINARY_DIR}/utilities/createErrors
  COMMAND ${CMAKE_BINARY_DIR}/utilities/createErrors
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/config/errors
  )
# Create target for error headers
add_custom_target(
  errorHeaders
  DEPENDS ${PROJECT_SOURCE_DIR}/config/errors/eenum.h
  DEPENDS ${PROJECT_SOURCE_DIR}/config/errors/emess.h  
  )

#-----------------------------------------------------------------------------
#############################
#  ADD SUB-FOLDERS
#############################

# Utilities folder
add_subdirectory(utilities)

# PenRed Kernel folder
add_subdirectory(kernel)

# PenRed Libs folder
add_subdirectory(lib)

# PenRed Geometry folder
add_subdirectory(geometry)

# particle generator folder
add_subdirectory(particleGen)

# tallies folder
add_subdirectory(tallies)

# main programs folder
add_subdirectory(mains)

# tests folder
add_subdirectory(tests)

#-----------------------------------------------------------------------------
