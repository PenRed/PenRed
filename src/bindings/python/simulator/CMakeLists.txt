
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

find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

# Get python version
execute_process(
    COMMAND ${Python3_EXECUTABLE}-config --extension-suffix
    OUTPUT_VARIABLE PYTHON_EXTENSION_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Create an extension module
add_library(simulator MODULE ${CMAKE_CURRENT_SOURCE_DIR}/simulator.cpp)
target_link_libraries(simulator PUBLIC pybind11::module)
target_link_libraries(simulator PUBLIC pen_sim)
set_target_properties(simulator PROPERTIES
      OUTPUT_NAME simulator${PYTHON_EXTENSION_SUFFIX}
      SUFFIX ""
      PREFIX "")

#-----------------------------------------------------------------------------
# Install penred executable
install(
  TARGETS simulator
  DESTINATION ${PROJECT_SOURCE_DIR}/bindings/python/pyPenred/pyPenred
  )
#-----------------------------------------------------------------------------
