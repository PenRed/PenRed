
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

# Create lib
add_library(pen_xray STATIC ${PROJECT_SOURCE_DIR}/lib/x-ray/source/x-ray.cpp ${PROJECT_SOURCE_DIR}/lib/x-ray/source/anode.cpp)

target_include_directories(pen_xray PUBLIC ${PROJECT_SOURCE_DIR}/lib/x-ray/includes)
target_link_libraries(pen_xray pen_sim)
target_link_libraries(pen_xray -Wl,--whole-archive pen_VR -Wl,--no-whole-archive )
