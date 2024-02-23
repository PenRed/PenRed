
#
#
#    Copyright (C) 2019 Universitat de València - UV
#    Copyright (C) 2019 Universitat Politècnica de València - UPV
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



#####################################
#  DICOM LOADER 
#####################################

add_definitions(-D_PEN_USE_DICOM_=1) #enable DICOM   
find_package(DCMTK REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${DCMTK_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${DCMTK_LINK_FLAGS})
include_directories(SYSTEM ${DCMTK_INCLUDE_DIRS})      

# Create DICOM lib
add_library(pen_dicom STATIC ${PROJECT_SOURCE_DIR}/kernel/parsers/DICOM/source/pen_dicom.cpp)

target_include_directories(pen_dicom PUBLIC ${PROJECT_SOURCE_DIR}/kernel/parsers/DICOM/includes)
target_include_directories(pen_dicom PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

# Link DCMTK lib
target_link_libraries(pen_dicom ${DCMTK_LIBRARIES})
