
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
#  MEASUREMENT TESTS
#############################


#       EXECUTABLES         #
#############################

add_executable(measure_t1 ${CMAKE_CURRENT_SOURCE_DIR}/test1.cpp)

target_include_directories(measure_t1 PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

#-----------------------------------------------------------------------------
install(
  TARGETS measure_t1
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/measurement
  )  
#-----------------------------------------------------------------------------

add_executable(measure_t2 ${CMAKE_CURRENT_SOURCE_DIR}/test2.cpp)

target_include_directories(measure_t2 PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

#-----------------------------------------------------------------------------
install(
  TARGETS measure_t2
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/measurement
  )
#-----------------------------------------------------------------------------

add_executable(measure_t3 ${CMAKE_CURRENT_SOURCE_DIR}/test3.cpp)

target_include_directories(measure_t3 PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

#-----------------------------------------------------------------------------
install(
  TARGETS measure_t3
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/measurement
  )
#-----------------------------------------------------------------------------

add_executable(measure_t4 ${CMAKE_CURRENT_SOURCE_DIR}/test4.cpp)

target_include_directories(measure_t4 PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

#-----------------------------------------------------------------------------
install(
  TARGETS measure_t4
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/measurement
  )
#-----------------------------------------------------------------------------

add_executable(measure_t5 ${CMAKE_CURRENT_SOURCE_DIR}/test5.cpp)

target_include_directories(measure_t5 PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

#-----------------------------------------------------------------------------
install(
  TARGETS measure_t5
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/measurement
  )
#-----------------------------------------------------------------------------
