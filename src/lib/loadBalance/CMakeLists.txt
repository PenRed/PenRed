
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
#  Load Balance 
#############################

if(WITH_MPI)
    #Find MPI package
    find_package(MPI REQUIRED)
    #Find threads package
    find_package (Threads REQUIRED)    
endif(WITH_MPI)


# Create lib
add_library(pen_LB STATIC ${PROJECT_SOURCE_DIR}/lib/loadBalance/loadBalance.cpp)

target_include_directories(pen_LB PUBLIC ${PROJECT_SOURCE_DIR}/lib/loadBalance/)

if(WITH_MPI)
    target_include_directories(pen_LB PUBLIC ${MPI_INCLUDE_PATH})
    target_link_libraries(pen_LB ${MPI_LIBRARIES})
    #Link thread libs
    target_link_libraries (pen_LB ${CMAKE_THREAD_LIBS_INIT})
endif(WITH_MPI)
