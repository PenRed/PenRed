
#
#
#    Copyright (C) 2021 Universitat de València - UV
#    Copyright (C) 2021 Universitat Politècnica de València - UPV
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
#        vicente.gimenez@uv.es (Vicente Giménez Gómez)
#        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
#

# Create C warper shared library
#############################

if(NOT BUILD_STATIC)
    target_compile_options(geoViewLib PRIVATE -fPIC)
    target_compile_options(pen_geos PRIVATE -fPIC)
    target_compile_options(pen_internalData PRIVATE -fPIC)
    target_compile_options(pen_dump PRIVATE -fPIC)
    target_compile_options(pen_aux PRIVATE -fPIC)
    target_compile_options(pen_math PRIVATE -fPIC)
    if(WITH_DICOM)
        target_compile_options(pen_image PRIVATE -fPIC)
        target_compile_options(pen_dicom PRIVATE -fPIC)
    endif(WITH_DICOM)
endif(NOT BUILD_STATIC)


if(BUILD_STATIC)
    add_library(geoView_C STATIC geoView_C.cpp)
else(BUILD_STATIC)
    add_library(geoView_C SHARED geoView_C.cpp)
endif(BUILD_STATIC)

target_link_libraries(geoView_C geoViewLib)

install(
TARGETS geoView_C
DESTINATION ${PROJECT_SOURCE_DIR}/bindings/C/viewers/geometry
)
