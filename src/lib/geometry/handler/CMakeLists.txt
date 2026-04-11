
#
#
#    Copyright (C) 2026 Vicent Giménez Alventosa
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
#



#############################
#  GEOMETRY HANDLER
#############################

#Find threads package
find_package (Threads REQUIRED)

# Create geoView lib
add_library(geoHandlerLib STATIC ${PROJECT_SOURCE_DIR}/lib/geometry/handler/geoHandler.cpp)

target_include_directories(geoHandlerLib PUBLIC ${PROJECT_SOURCE_DIR}/lib/geometry/handler)

#target libs
target_link_libraries(geoHandlerLib geoViewLib)
#Link thread libs
target_link_libraries(geoHandlerLib ${CMAKE_THREAD_LIBS_INIT})
