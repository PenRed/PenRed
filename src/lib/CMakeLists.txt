
#
#
#    Copyright (C) 2019-2022 Universitat de València - UV
#    Copyright (C) 2019-2022 Universitat Politècnica de València - UPV
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
#    
#



#############################
#  LIBS
#############################


#############################
#  ADD SUB-FOLDERS
#############################

# math folder
add_subdirectory(math)

# dump folder
add_subdirectory(dump)

# phase space file folder
add_subdirectory(psf)

# splitted file folder
add_subdirectory(splittedFile)

# shared file folder
add_subdirectory(sharedFile)

# auxiliar folder
add_subdirectory(auxiliar)

# TCP folder
add_subdirectory(tcp)

# image folder
add_subdirectory(image)

# load balance folder
if(WITH_LB)
    add_subdirectory(loadBalance)
endif(WITH_LB)

# muen folder
add_subdirectory(parameters)

# simulation folder
add_subdirectory(simulation)

# x-ray folder
add_subdirectory(x-ray)
