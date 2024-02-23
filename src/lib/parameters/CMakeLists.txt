
#
#
#    Copyright (C) 2023 Universitat de València - UV
#    Copyright (C) 2023 Universitat Politècnica de València - UPV
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
#    
#



#############################
#  MUEN
#############################

if(WITH_MULTI_THREADING)
    #Find threads package
    find_package (Threads REQUIRED)
endif(WITH_MULTI_THREADING)

# Create lib
add_library(pen_muen STATIC ${PROJECT_SOURCE_DIR}/lib/parameters/source/pen_muen.cpp)

target_include_directories(pen_muen PUBLIC ${PROJECT_SOURCE_DIR}/lib/parameters/includes)

target_link_libraries(pen_muen pen_particles)
target_link_libraries(pen_muen pen_geos)
target_link_libraries(pen_muen samplers)
target_link_libraries(pen_muen pen_timer)

if(WITH_MULTI_THREADING)
    #Link thread libs
    target_link_libraries (pen_muen ${CMAKE_THREAD_LIBS_INIT})
endif(WITH_MULTI_THREADING)
