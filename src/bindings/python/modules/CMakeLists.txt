
#
#
#    Copyright (C) 2024 Universitat de València - UV
#    Copyright (C) 2024 Universitat Politècnica de València - UPV
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

find_package(Python3 COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

# If pybind11 is not installed, download it from github
if(NOT pybind11_FOUND)

        # Ensure the correct Python version is used with pybind11
        set(PYTHON_INCLUDE_DIR ${Python3_INCLUDE_DIRS})
        set(PYTHON_LIBRARIES ${Python3_LIBRARIES})

        set(PYBIND11_FINDPYTHON ON CACHE BOOL "" FORCE)
        set(PYBIND11_PYTHON_VERSION ${Python3_VERSION} CACHE STRING "" FORCE)

        include(FetchContent)

        # Download pybind11
        FetchContent_Declare(
                pybind11
                GIT_REPOSITORY https://github.com/pybind/pybind11.git
        )

        # Add pybind11 to the project
        FetchContent_MakeAvailable(pybind11)

endif(NOT pybind11_FOUND)

# Get python version
execute_process(
    COMMAND ${Python3_EXECUTABLE}-config --extension-suffix
    OUTPUT_VARIABLE PYTHON_EXTENSION_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(DEFINED PYTHON_MODULE_EXTENSION)
    set(MODULE_SUFFIX ${PYTHON_MODULE_EXTENSION})
else()
    set(MODULE_SUFFIX ${PYTHON_EXTENSION_SUFFIX})
endif()

# Create extension modules
add_library(simulator MODULE ${CMAKE_CURRENT_SOURCE_DIR}/simulator.cpp)
target_link_libraries(simulator PUBLIC pybind11::module)
target_link_libraries(simulator PUBLIC pen_sim)
if(XRAY)
    target_link_libraries(simulator PUBLIC pen_xray)
endif(XRAY)
set_target_properties(simulator PROPERTIES
      OUTPUT_NAME simulation${MODULE_SUFFIX}
      SUFFIX ""
      PREFIX "")
      
add_library(pyPSF MODULE ${CMAKE_CURRENT_SOURCE_DIR}/psf.cpp)
target_link_libraries(pyPSF PUBLIC pybind11::module)
target_link_libraries(pyPSF PUBLIC pen_psf)
set_target_properties(pyPSF PROPERTIES
      OUTPUT_NAME psf${MODULE_SUFFIX}
      SUFFIX ""
      PREFIX "")

#-----------------------------------------------------------------------------
# Install python modules
install(
  TARGETS simulator
  DESTINATION ${PROJECT_SOURCE_DIR}/bindings/python/pyPenred/pyPenred
  )
  
install(
  TARGETS pyPSF
  DESTINATION ${PROJECT_SOURCE_DIR}/bindings/python/pyPenred/pyPenred
  )
#-----------------------------------------------------------------------------
