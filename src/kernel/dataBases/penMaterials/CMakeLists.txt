
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




##############################
#  PENELOPE MATERIAL DATABASE
##############################

# Create the program to create the data base
add_executable(penelope_material_db_creator ${CMAKE_CURRENT_SOURCE_DIR}/source/createDB.cpp)

# Create a directory to store the generated database
add_custom_target(makeDirDB ALL COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/database)

# Run the DB creator if needed
add_custom_command(
	OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/database/database.hh
	COMMAND penelope_material_db_creator ${PROJECT_SOURCE_DIR}/../dataBases/penmaterials/filesInDB.txt ${PROJECT_SOURCE_DIR}/../dataBases/penmaterials/pdfiles ${CMAKE_CURRENT_BINARY_DIR}/database
	DEPENDS penelope_material_db_creator
	DEPENDS makeDirDB
)
#add_custom_target(genDB DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/databse/database.hh)

# Create lib
add_library(penelope_material_db STATIC ${CMAKE_CURRENT_SOURCE_DIR}/source/materialCreator.cpp ${CMAKE_CURRENT_BINARY_DIR}/database/database.hh)
#add_dependencies(penelope_material_db genDB)

target_include_directories(penelope_material_db PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/includes)
target_include_directories(penelope_material_db PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/database)
