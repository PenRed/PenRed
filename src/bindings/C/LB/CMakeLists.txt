
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
#        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
#        vicente.gimenez@uv.es (Vicente Giménez Gómez)
#    
#

#set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-all-symbols -all_load --whole-file")

# Set TCP lib link flags
##########################
if(NOT BUILD_STATIC)
    target_compile_options(pen_tcp PRIVATE -fPIC)
endif(NOT BUILD_STATIC)

# Set RUPER_LB lib link flags
##############################
if(NOT BUILD_STATIC)
    target_compile_options(pen_LB PRIVATE -fPIC)
endif(NOT BUILD_STATIC)

# Create C warper module
#############################
if(BUILD_STATIC)
    add_library(RUPER_LB_C STATIC server_wraper.cpp)
else(BUILD_STATIC)
    add_library(RUPER_LB_C MODULE server_wraper.cpp)
endif(BUILD_STATIC)
target_link_libraries(RUPER_LB_C pen_LB)

install(
TARGETS RUPER_LB_C
DESTINATION ${PROJECT_SOURCE_DIR}/bindings/C/LB/
)
