
#
#
#    Copyright (C) 2020 Universitat de València - UV
#    Copyright (C) 2020 Universitat Politècnica de València - UPV
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
#  TCP lib
#############################
                       
if(WITH_SSL)
    find_package(OpenSSL REQUIRED)
endif(WITH_SSL)

#Find threads package
find_package (Threads REQUIRED)    

# Create lib
add_library(pen_tcp STATIC ${PROJECT_SOURCE_DIR}/lib/tcp/tcp_cs.cpp)

target_include_directories(pen_tcp PUBLIC ${PROJECT_SOURCE_DIR}/lib/tcp/)
target_include_directories(pen_tcp PUBLIC ${PROJECT_SOURCE_DIR}/lib/tcp/ASIO/include/)

#Link thread libs
target_link_libraries(pen_tcp ${CMAKE_THREAD_LIBS_INIT})

if(WITH_SSL)
    target_include_directories(pen_tcp PUBLIC ${OPENSSL_INCLUDE_DIR})
    target_link_libraries(pen_tcp ${OPENSSL_LIBRARIES})
endif(WITH_SSL)
