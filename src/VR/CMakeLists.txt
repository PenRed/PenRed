
#
#
#    Copyright (C) 2020-2021 Universitat de València - UV
#    Copyright (C) 2020-2021 Universitat Politècnica de València - UPV
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
#         VR
#############################


# Source Classes lib 
#############################
add_library(pen_VR STATIC ${PROJECT_SOURCE_DIR}/VR/source/pen_vr.cpp)

#In some compilers the string instantiator method causes the source files to be 
#detected as unused code. Therefore, the code is not included and during the execution 
#the instantiator fails. Linking the library as:
#
# target_link_libraries(executable -Wl,--whole-archive pen_VR -Wl,--no-whole-archive)
#
# instead of:
#
# target_link_libraries(executable pen_VR)
#
# will solve this problem.

target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/kernel/states)
target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/VR/includes)
target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/VR/source)
target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/lib/instantiator)

target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/VR/specific/includes)
target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/VR/specific/source)

target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/VR/generic/includes)
target_include_directories(pen_VR PUBLIC ${PROJECT_SOURCE_DIR}/VR/generic/source)

target_link_libraries(pen_VR pen_randoms)
target_link_libraries(pen_VR pen_math)
target_link_libraries(pen_VR pen_aux)
target_link_libraries(pen_VR pen_internalData)
target_link_libraries(pen_VR pen_grids)

