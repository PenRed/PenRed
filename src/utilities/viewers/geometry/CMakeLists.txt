
#
#
#    Copyright (C) 2021-2022 Universitat de València - UV
#    Copyright (C) 2021-2022 Universitat Politècnica de València - UPV
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
#        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
#



#############################
#  GEOMETRY VIEWERS
#############################

#Find threads package
find_package (Threads REQUIRED)

# Create geoView lib
add_library(geoViewLib STATIC ${PROJECT_SOURCE_DIR}/utilities/viewers/geometry/source/pen_geoView.cpp)

target_include_directories(geoViewLib PUBLIC ${PROJECT_SOURCE_DIR}/utilities/viewers/geometry/includes)

#target libs
target_link_libraries(geoViewLib pen_geos)
#Link thread libs
target_link_libraries(geoViewLib ${CMAKE_THREAD_LIBS_INIT})

# Create geoView executable
add_executable(geoView ${PROJECT_SOURCE_DIR}/utilities/viewers/geometry/asciiView.cpp)

#target libs
target_link_libraries(geoView geoViewLib)

# Install readMat executable
install(
    TARGETS geoView
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/viewers/geometry/geoView
)

