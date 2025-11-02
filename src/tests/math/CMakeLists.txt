
#
#
#    Copyright (C) 2025 Vicent Giménez Alventosa
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
#  MATH LIBRARY TESTS
#############################


#       EXECUTABLES         #
#############################

add_executable(math_t1 ${CMAKE_CURRENT_SOURCE_DIR}/math1.cpp)

target_include_directories(math_t1 PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

#-----------------------------------------------------------------------------
install(
  TARGETS math_t1
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/math
  )  
#-----------------------------------------------------------------------------
