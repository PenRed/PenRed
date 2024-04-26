
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
#  SAMPLER TESTS
#############################


#       EXECUTABLES         #
#############################

add_executable(sample_t0 ${CMAKE_CURRENT_SOURCE_DIR}/test0.cpp)

target_link_libraries(sample_t0 PUBLIC pen_randoms)

#-----------------------------------------------------------------------------
install(
  TARGETS sample_t0
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/samplers
  )
#-----------------------------------------------------------------------------

add_executable(sample_t1 ${CMAKE_CURRENT_SOURCE_DIR}/test1.cpp)

target_link_libraries(sample_t1 PUBLIC pen_randoms)

#-----------------------------------------------------------------------------
install(
  TARGETS sample_t1
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/samplers
  )  
#-----------------------------------------------------------------------------
