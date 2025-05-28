
#
#
#    Copyright (C) 2019-2024 Universitat de València - UV
#    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
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
#  INTERNAL DATA FORMAT AND PARSERS 
#####################################

# Create geo lib
add_library(pen_internalData STATIC ${CMAKE_CURRENT_SOURCE_DIR}/source/pen_parser.cpp ${CMAKE_CURRENT_SOURCE_DIR}/source/pen_reader.cpp ${CMAKE_CURRENT_SOURCE_DIR}/source/pen_data.cpp)

target_include_directories(pen_internalData PUBLIC ${PROJECT_SOURCE_DIR}/kernel/parsers/internalData/includes)

target_link_libraries(pen_internalData PUBLIC pen_logger)
