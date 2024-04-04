
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
#        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
#    
#

                       
# PENELOPE MATERIAL UTILITIES #
###############################

#-----------------------------------------------------------------------------
#  Read material
#############################

    # Create read material executable
    add_executable(readMat ${CMAKE_CURRENT_SOURCE_DIR}/readMat.cpp)

    #target libs
    target_link_libraries(readMat pen_materials)

    # Install readMat executable
    install(
    TARGETS readMat
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/penMats
    )

#-----------------------------------------------------------------------------
#  Pen material ranges
#############################

    # Create read material executable
    add_executable(range ${CMAKE_CURRENT_SOURCE_DIR}/range.cpp)

    #target libs
    target_link_libraries(range pen_particles)

    # Install readMat executable
    install(
    TARGETS range
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/penMats
    )
  
#-----------------------------------------------------------------------------
#  Pen material creator
#############################
if(NOT MSVC)
    # Create read material executable
    add_executable(createMat ${CMAKE_CURRENT_SOURCE_DIR}/createMat.cpp)

    #target libs
    target_link_libraries(createMat penelope_material_db)

    # Install readMat executable
    install(
    TARGETS createMat
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/penMats
    )
endif(NOT MSVC)

#-----------------------------------------------------------------------------
#  Pen material creator from list
###################################

    # Create read material executable
    add_executable(createListedMat ${CMAKE_CURRENT_SOURCE_DIR}/createListedMat.cpp)

    #target libs
    target_link_libraries(createListedMat penelope_material_db)

    # Install readMat executable
    install(
    TARGETS createListedMat
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/penMats
    )

