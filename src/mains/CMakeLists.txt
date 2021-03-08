
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
#  MAIN PROGRAMS
#############################

# Create penmain
add_executable(pen_main ${PROJECT_SOURCE_DIR}/mains/source/pen_main.cpp)

target_include_directories(pen_main PUBLIC ${PROJECT_SOURCE_DIR}/mains/includes)
target_include_directories(pen_main PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

target_link_libraries(pen_main pen_particles)
target_link_libraries(pen_main pen_geos)
target_link_libraries(pen_main samplers)
target_link_libraries(pen_main tallies)
target_link_libraries(pen_main pen_VR)
target_link_libraries(pen_main pen_timer)


#-----------------------------------------------------------------------------
# Install penmain executable
install(
  TARGETS pen_main
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/mains/penmain
  )
#-----------------------------------------------------------------------------
