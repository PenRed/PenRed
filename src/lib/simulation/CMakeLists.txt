
#
#
#    Copyright (C) 2024 Universitat de València - UV
#    Copyright (C) 2024 Universitat Politècnica de València - UPV
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
#        sanolgi@upvnet.upv.es
#        vicente.gimenez@uv.es (Vicente Giménez Gómez)
#    
#



#############################
#  Simulation lib
#############################

if(WITH_MULTI_THREADING)
    #Find threads package
    find_package (Threads REQUIRED)
endif(WITH_MULTI_THREADING)

# Create lib
add_library(pen_sim STATIC ${PROJECT_SOURCE_DIR}/lib/simulation/source/pen_simulation.cpp)

target_include_directories(pen_sim PUBLIC ${PROJECT_SOURCE_DIR}/lib/simulation/includes)
target_link_libraries(pen_sim pen_geos)
target_link_libraries(pen_sim samplers)
target_link_libraries(pen_sim tallies)
target_link_libraries(pen_sim pen_VR)


#Link thread libs
if(WITH_MULTI_THREADING)
    target_link_libraries(pen_sim ${CMAKE_THREAD_LIBS_INIT})
endif(WITH_MULTI_THREADING)
