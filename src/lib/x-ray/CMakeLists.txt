
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
#        vicent.gimenez.alventosa@gmail.com
#    
#




#############################
#  X-RAY LIB
#############################

set(SOURCE
${CMAKE_CURRENT_SOURCE_DIR}/source/collimator.cpp
${CMAKE_CURRENT_SOURCE_DIR}/source/filter.cpp
${CMAKE_CURRENT_SOURCE_DIR}/source/anode.cpp
${CMAKE_CURRENT_SOURCE_DIR}/source/phantom.cpp
${CMAKE_CURRENT_SOURCE_DIR}/source/utilities.cpp
${CMAKE_CURRENT_SOURCE_DIR}/source/device.cpp)

# Create lib
add_library(pen_xray STATIC ${SOURCE})

target_include_directories(pen_xray PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/includes)
target_link_libraries(pen_xray pen_sim)

if (MSVC)
	target_link_libraries(pen_xray pen_VR)
	set_target_properties(pen_xray PROPERTIES LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /WHOLEARCHIVE:pen_VR")
else(MSVC)

	if(APPLE) #Check if is an Apple OS
		target_link_libraries(pen_xray -Wl,-force_load pen_VR)
	else(APPLE) #By default assume a Linux OS
		target_link_libraries(pen_xray -Wl,--whole-archive pen_VR -Wl,--no-whole-archive )
	endif(APPLE)

endif(MSVC)
