
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
#  PARTICLE GENERATION
#############################


#        LIBRARIES          #
#############################

if(WITH_MPI)
    #Find MPI package
    find_package(MPI REQUIRED)
endif(WITH_MPI)


# Source Classes lib 
#############################
add_library(samplers STATIC ${PROJECT_SOURCE_DIR}/particleGen/source/pen_samplers.cpp)

target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/includes)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/source)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/kernel/states)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/lib/instantiator)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/specific/includes)

target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/spatial/includes)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/spatial/source)

target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/direction/includes)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/direction/source)

target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/energy/includes)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/energy/source)

target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/time/includes)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/time/source)

target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/specific/includes)
target_include_directories(samplers PUBLIC ${PROJECT_SOURCE_DIR}/particleGen/specific/source)

if(WITH_MPI)
    target_include_directories(samplers PUBLIC ${MPI_INCLUDE_PATH})
endif(WITH_MPI)


target_link_libraries(samplers pen_randoms)
target_link_libraries(samplers pen_math)
target_link_libraries(samplers pen_aux)
target_link_libraries(samplers pen_contexts)
target_link_libraries(samplers pen_internalData)
target_link_libraries(samplers pen_psf)
target_link_libraries(samplers pen_sharedFile)

if(WITH_MPI)
    target_link_libraries(samplers ${MPI_LIBRARIES})
endif(WITH_MPI)
