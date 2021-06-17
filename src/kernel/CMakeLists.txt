
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



#-----------------------------------------------------------------------------
# Add include directories
include_directories(${PROJECT_SOURCE_DIR}/kernel/includes)
#-----------------------------------------------------------------------------


#############################
#  KERNEL DEPENDENCIES
#############################

# Create pen_errors lib
#############################
add_library(pen_errors STATIC ${PROJECT_SOURCE_DIR}/kernel/errors/pen_errors.cpp)

target_include_directories(pen_errors PUBLIC ${PROJECT_SOURCE_DIR}/kernel/errors)
target_include_directories(pen_errors PUBLIC ${PROJECT_SOURCE_DIR}/config/errors)

# Create pen_randoms lib
#############################
add_library(pen_randoms STATIC ${PROJECT_SOURCE_DIR}/kernel/rands/source/pen_random.cpp)

target_include_directories(pen_randoms PUBLIC ${PROJECT_SOURCE_DIR}/kernel/rands/includes)

target_link_libraries(pen_randoms pen_errors)

# Create pen_grids lib
#############################
add_library(pen_grids STATIC ${PROJECT_SOURCE_DIR}/kernel/grids/source/pen_grids.cpp)

target_include_directories(pen_grids PUBLIC ${PROJECT_SOURCE_DIR}/kernel/grids/includes)
target_include_directories(pen_grids PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

target_link_libraries(pen_grids pen_errors)

# Create pen_contexts lib
#############################
add_library(pen_contexts STATIC ${PROJECT_SOURCE_DIR}/kernel/contexts/source/pen_contexts.cpp)

target_include_directories(pen_contexts PUBLIC ${PROJECT_SOURCE_DIR}/kernel/contexts/includes)

target_link_libraries(pen_contexts pen_materials)

# Create pen_particles lib
#############################
add_library(pen_particles STATIC ${PROJECT_SOURCE_DIR}/kernel/particles/source/pen_particles.cpp)

target_include_directories(pen_particles PUBLIC ${PROJECT_SOURCE_DIR}/kernel/particles/includes)

target_link_libraries(pen_particles pen_contexts)

# Create pen_timer lib
#############################
add_library(pen_timer STATIC ${PROJECT_SOURCE_DIR}/kernel/timer/pen_timer.cpp)

target_include_directories(pen_timer PUBLIC ${PROJECT_SOURCE_DIR}/kernel/timer)

#############################
#  ADD SUB-FOLDERS
#############################

# data folder
add_subdirectory(parsers)

# materials folder
add_subdirectory(materials)
