
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
#         TALLIES
#############################


#        LIBRARIES          #
#############################

if(WITH_MULTI_THREADING)
    #Find threads package
    find_package (Threads REQUIRED)
endif(WITH_MULTI_THREADING)

if(WITH_MPI)
    #Find MPI package
    find_package(MPI REQUIRED)
endif(WITH_MPI)

# Source Classes lib 
#############################
add_library(tallies STATIC ${PROJECT_SOURCE_DIR}/tallies/source/pen_tallies.cpp)

target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/kernel/states)
target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)
target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/tallies/includes)
target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/tallies/source)
target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/lib/instantiator)

target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/tallies/specific/includes)
target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/tallies/specific/source)

target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/tallies/generic/includes)
target_include_directories(tallies PUBLIC ${PROJECT_SOURCE_DIR}/tallies/generic/source)

if(WITH_MPI)
    target_include_directories(tallies PUBLIC ${MPI_INCLUDE_PATH})
endif(WITH_MPI)


target_link_libraries(tallies pen_geos)
target_link_libraries(tallies pen_randoms)
target_link_libraries(tallies pen_math)
target_link_libraries(tallies pen_dump)
target_link_libraries(tallies pen_aux)
target_link_libraries(tallies pen_internalData)
target_link_libraries(tallies pen_splittedFile)
target_link_libraries(tallies pen_psf)
target_link_libraries(tallies pen_grids)

if(WITH_MULTI_THREADING)
    #Link thread libs
    target_link_libraries (tallies ${CMAKE_THREAD_LIBS_INIT})
endif(WITH_MULTI_THREADING)

if(WITH_MPI)
    target_link_libraries(tallies ${MPI_LIBRARIES})
endif(WITH_MPI)
