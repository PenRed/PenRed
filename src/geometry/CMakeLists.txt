
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
#        vicent.gimenez.alventosa@gmail.com
#        vicente.gimenez@uv.es
#    
#




#############################
#  GEOMETRY 
#############################

if(WITH_MULTI_THREADING)
    #Find threads package
    find_package (Threads REQUIRED)
endif(WITH_MULTI_THREADING)

# Create geo lib
add_library(pen_geos STATIC ${PROJECT_SOURCE_DIR}/geometry/source/pen_geometries.cpp)

#In some compilers the string instantiator method causes the source files to be 
#detected as unused code. Therefore, the code is not included and during the execution 
#the instantiator fails. Linking the library as:
#
# target_link_libraries(executable -Wl,--whole-archive pen_geos -Wl,--no-whole-archive)
#
# instead of:
#
# target_link_libraries(executable pen_geos)
#
# will solve this problem.

target_include_directories(pen_geos PUBLIC ${PROJECT_SOURCE_DIR}/geometry/includes)
target_include_directories(pen_geos PUBLIC ${PROJECT_SOURCE_DIR}/geometry/objects/includes)
target_include_directories(pen_geos PUBLIC ${PROJECT_SOURCE_DIR}/geometry/meshes/includes)
target_include_directories(pen_geos PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)
target_include_directories(pen_geos PUBLIC ${PROJECT_SOURCE_DIR}/kernel/states)
target_include_directories(pen_geos PUBLIC ${PROJECT_SOURCE_DIR}/lib/instantiator)

target_link_libraries (pen_geos pen_errors)
target_link_libraries (pen_geos pen_internalData)
target_link_libraries (pen_geos pen_aux)
target_link_libraries (pen_geos pen_math)
target_link_libraries (pen_geos pen_dump)
target_link_libraries (pen_geos pen_image)

if(WITH_MULTI_THREADING)
    #Link thread libs
    target_link_libraries (pen_geos ${CMAKE_THREAD_LIBS_INIT})
endif(WITH_MULTI_THREADING)

#Check if DICOM option has been enabled
if(WITH_DICOM)
  #dcmtk libs
  target_link_libraries(pen_geos pen_dicom)
endif()
