
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



#############################
#  MAIN PROGRAMS
#############################

# Create penmain
add_executable(pen_main ${PROJECT_SOURCE_DIR}/mains/source/pen_main.cpp)

target_include_directories(pen_main PUBLIC ${PROJECT_SOURCE_DIR}/mains/includes)
target_include_directories(pen_main PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

target_link_libraries(pen_main pen_particles)
if (MSVC) 
	target_link_libraries(pen_main pen_geos)
	target_link_libraries(pen_main samplers)
	target_link_libraries(pen_main tallies)
	target_link_libraries(pen_main pen_VR)
	set_target_properties(pen_main PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /WHOLEARCHIVE:pen_geos /WHOLEARCHIVE:samplers /WHOLEARCHIVE:tallies /WHOLEARCHIVE:pen_VR")
else(MSVC)
	target_link_libraries(pen_main -Wl,--whole-archive pen_geos -Wl,--no-whole-archive )
	target_link_libraries(pen_main -Wl,--whole-archive samplers -Wl,--no-whole-archive )
	target_link_libraries(pen_main -Wl,--whole-archive tallies -Wl,--no-whole-archive )
	target_link_libraries(pen_main -Wl,--whole-archive pen_VR -Wl,--no-whole-archive )
endif(MSVC)
target_link_libraries(pen_main pen_timer)


#-----------------------------------------------------------------------------
# Install penmain executable
install(
  TARGETS pen_main
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/mains/penmain
  )
#-----------------------------------------------------------------------------
