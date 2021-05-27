
#
#
#    Copyright (C) 2019-2020 Universitat de València - UV
#    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
#-----------------------------------------------------------------------------
CMAKE_POLICY(SET CMP0054 NEW)
#-----------------------------------------------------------------------------
include(CMakeDependentOption)

option(DEVELOPMENT_WARNINGS "Compile with extra warnings for development " OFF)
option(WITH_DICOM "Build PenRed DICOM " OFF)
option(WITH_MULTI_THREADING "Build PenRed with multi-threading " ON)
option(WITH_MPI "Build PenRed with MPI " OFF)

option(WITH_NATIVE "Compile PenRed with native architecture optimizations" ON)

CMAKE_DEPENDENT_OPTION(WITH_LB "Build PenRed with load balancing support" OFF
                       "WITH_MULTI_THREADING" OFF)

CMAKE_DEPENDENT_OPTION(WITH_SSL "Build TCP with SSL support" OFF
                       "WITH_LB" OFF)
                       
CMAKE_DEPENDENT_OPTION(WITH_HTTP "Build LB with HTTP support" OFF
                       "WITH_LB" OFF)

option(BUILD_TESTS "Compile PenRed tests" OFF)
option(BUILD_UTILITIES "Compile PenRed utilities" ON)
option(BUILD_STATIC "Enable static compilation" OFF)
mark_as_advanced(FORCE BUILD_STATIC)
option(BUILD_C_BINDS "Compile C binds" OFF)
mark_as_advanced(FORCE BUILD_C_BINDS)
option(WITH_DEBUG "Compile with debug flag" OFF)
mark_as_advanced(WITH_DEBUG)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Project name
project(PenRed LANGUAGES CXX)
#-----------------------------------------------------------------------------

if(BUILD_STATIC)
    #SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    #SET(BUILD_SHARED_LIBS OFF)
    SET(CMAKE_EXE_LINKER_FLAGS "-static")
    #SET(CMAKE_MODULE_LINKER_FLAGS "-static")
    #SET(CMAKE_SHARED_LINKER_FLAGS "-static")
    SET(CMAKE_LINK_SEARCH_START_STATIC ON)
    SET(CMAKE_LINK_SEARCH_END_STATIC ON)
    add_definitions(-D_PEN_BUILD_STATIC_=1) #Set Static build
else(BUILD_STATIC)
endif(BUILD_STATIC)

# Compilation options

if (WITH_DEBUG)
	add_compile_options(-g)
endif(WITH_DEBUG)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  # using GCC
  if (WITH_NATIVE)
        add_compile_options(-march=native)
  endif (WITH_NATIVE)
  if (DEVELOPMENT_WARNINGS)
  	add_compile_options(-std=c++11 -Wall -Wpedantic -Wextra -Wshadow -O2 -Wno-unused-result -Wno-deprecated-copy -frecord-gcc-switches)
  else (DEVELOPMENT_WARNINGS)
	add_compile_options(-std=c++11 -Wall -O2 -frecord-gcc-switches)
  endif (DEVELOPMENT_WARNINGS)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  # using Intel C++
  if (WITH_NATIVE)
        add_compile_options(-xHost)
  endif (WITH_NATIVE)
  add_compile_options(-std=c++11 -O2 -no-prec-div -fp-model precise -prec-sqrt)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # using Visual Studio C++
endif()

if (MSVC) 
  # using Visual Studio C++
  add_compile_options(/arch:AVX2)
  # Check if MPI option has been enabled
  if(WITH_MPI)
	add_compile_options(/MP) #enable MPI
  endif(WITH_MPI) 
  # Needed for M_PI to be visible in visual studio.
  ADD_DEFINITIONS(-D_USE_MATH_DEFINES)

  # Removes the warnings about unsafe methods such as strcpy, std::copy,
  # memcmp, etc.
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
  
  # We need a stack size of 8MB (as it is default on linux)
  # See https://github.com/open62541/open62541/issues/1326
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:8388608")
	
endif(MSVC)

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

# Check if load balancing has been enabled
if(WITH_LB)
  add_definitions(-D_PEN_USE_LB_=1) #enable LB
endif(WITH_LB)

# Check if SSL support has been enabled
if(WITH_SSL)
    add_definitions(-D_PEN_USE_SSL_=1) #enable SSL
endif(WITH_SSL)

# Check if HTTP support has been enabled
if(WITH_HTTP)
    add_definitions(-D_PEN_USE_LB_HTTP_=1) #enable HTTP
endif(WITH_HTTP)

# Check platform
if(UNIX)
  add_definitions(-D_PEN_UNIX_=1) #enable Unix
elseif(WIN32)
  
elseif(APPLE)
  
endif(UNIX)
  
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

# VR folder
add_subdirectory(VR)

# main programs folder
add_subdirectory(mains)

# tests folder
if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)

# bindings folder
if(BUILD_C_BINDS)
    add_subdirectory(bindings/C)
endif(BUILD_C_BINDS)
#-----------------------------------------------------------------------------
