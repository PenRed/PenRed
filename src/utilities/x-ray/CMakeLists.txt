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

    # Create x-ray anode simulation executable
    add_executable(anodeSim ${PROJECT_SOURCE_DIR}/utilities/x-ray/simAnode.cpp)

    target_link_libraries(anodeSim pen_xray)

    # Install anode simulation executable
    install(
    TARGETS anodeSim
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/x-ray/anodeSim
    )

    # Create x-ray basic filter generator executable
    add_executable(createBaseFilter ${PROJECT_SOURCE_DIR}/utilities/x-ray/createBaseFilter.cpp)

    target_link_libraries(createBaseFilter pen_xray)

    # Install anode simulation executable
    install(
    TARGETS createBaseFilter
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/x-ray/createBaseFilter
    )

    # Create x-ray basic collimator generator executable
    add_executable(createBaseCollimator ${PROJECT_SOURCE_DIR}/utilities/x-ray/createBaseCollimator.cpp)

    target_link_libraries(createBaseCollimator pen_xray)

    # Install anode simulation executable
    install(
    TARGETS createBaseCollimator
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/x-ray/createBaseCollimator
    )

    # Create x-ray utility executable
    add_executable(x-ray ${CMAKE_CURRENT_SOURCE_DIR}/x-ray.cpp)

    target_link_libraries(x-ray pen_xray)

    # Install anode simulation executable
    install(
    TARGETS x-ray
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/x-ray
    )
