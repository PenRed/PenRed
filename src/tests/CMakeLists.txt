
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



#############################
#  LIBS
#############################

if(WITH_MULTI_THREADING)
    #Find threads package
    find_package (Threads REQUIRED)
endif(WITH_MULTI_THREADING)


#############################
#  TESTS
#############################


#       EXECUTABLES         #
#############################

#   Internal data tests
#############################

### Test 1

add_executable(internalData_t1 ${PROJECT_SOURCE_DIR}/tests/internalData/test1.cpp)

target_include_directories(internalData_t1 PUBLIC ${PROJECT_SOURCE_DIR}/utilities/parsers/includes)

#target libs
target_link_libraries(internalData_t1 pen_internalData)

# Install executables
install(
  TARGETS internalData_t1
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/internalData
  )
  
  
### Test 2

add_executable(internalData_t2 ${PROJECT_SOURCE_DIR}/tests/internalData/test2.cpp)

target_include_directories(internalData_t2 PUBLIC ${PROJECT_SOURCE_DIR}/utilities/parsers/includes)

#target libs
target_link_libraries(internalData_t2 pen_internalData)

# Install executables
install(
  TARGETS internalData_t2
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/internalData
  )


##########################################################

#   Sampling tests
#############################

add_executable(generalSampling_t1 ${PROJECT_SOURCE_DIR}/tests/stateSampler/generalSampling.cpp)

target_include_directories(generalSampling_t1 PUBLIC ${PROJECT_SOURCE_DIR}/utilities/parsers/includes)
target_include_directories(generalSampling_t1 PUBLIC ${PROJECT_SOURCE_DIR}/geometry/includes)

#target libs
target_link_libraries(generalSampling_t1 samplers)
target_link_libraries(generalSampling_t1 pen_geos)


#-----------------------------------------------------------------------------
install(
  TARGETS generalSampling_t1
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/particleGenerators
  )  
#-----------------------------------------------------------------------------

#   Dump tests
#############################

add_executable(dump_t1 ${PROJECT_SOURCE_DIR}/tests/dump/dumpTest.cpp)

#target libs
target_link_libraries(dump_t1 pen_dump)

#-----------------------------------------------------------------------------
install(
  TARGETS dump_t1
  DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/dump
  )  
  
#-----------------------------------------------------------------------------

#   Splitted files tests
#############################

# Example 1
#-----------

if(WITH_MULTI_THREADING)

    add_executable(split_t1 ${PROJECT_SOURCE_DIR}/tests/splittedFiles/split_t1.cpp)

    #target libs
    target_link_libraries(split_t1 pen_splittedFile)
    #Link thread libs
    target_link_libraries (split_t1 ${CMAKE_THREAD_LIBS_INIT})


#-----------------------------------------------------------------------------
    install(
        TARGETS split_t1
        DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/splittedFiles
    )  
#-----------------------------------------------------------------------------

endif(WITH_MULTI_THREADING)

#   shared file tests
#############################

# Example 1
#-----------

if(WITH_MULTI_THREADING)

    add_executable(shared_t1 ${PROJECT_SOURCE_DIR}/tests/sharedFile/test1.cpp)

    #target libs
    target_link_libraries(shared_t1 pen_sharedFile)
    #Link thread libs
    target_link_libraries (shared_t1 ${CMAKE_THREAD_LIBS_INIT})


#-----------------------------------------------------------------------------
    install(
        TARGETS shared_t1 
        DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/sharedFile
    )  
#-----------------------------------------------------------------------------

endif(WITH_MULTI_THREADING)

#-----------------------------------------------------------------------------

#   PSF tests
#############################

# Example 1
#-----------

add_executable(psf_t1 ${PROJECT_SOURCE_DIR}/tests/psf/writeRead.cpp)

#target libs
target_link_libraries(psf_t1 samplers)

#-----------------------------------------------------------------------------
install(
    TARGETS psf_t1 
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/psf
)  
#-----------------------------------------------------------------------------

# Example 2
#-----------

if(WITH_MULTI_THREADING)

    add_executable(psf_t2 ${PROJECT_SOURCE_DIR}/tests/psf/multithreadRead.cpp)

    #target libs
    target_link_libraries(psf_t2 samplers)
    target_link_libraries(psf_t2 pen_splittedFile)
    #Link thread libs
    target_link_libraries(psf_t2 ${CMAKE_THREAD_LIBS_INIT})


#-----------------------------------------------------------------------------
    install(
        TARGETS psf_t2
        DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/psf
    )  
#-----------------------------------------------------------------------------

endif(WITH_MULTI_THREADING)

#-----------------------------------------------------------------------------

#   DICOM tests
#############################

if(WITH_DICOM)

    # Example 1
    #-----------

    add_executable(dicom_t1 ${PROJECT_SOURCE_DIR}/tests/DICOM/test_dicom.cpp)

    #target libs
    target_link_libraries(dicom_t1 pen_dicom)

    #-----------------------------------------------------------------------------
    install(
    TARGETS dicom_t1
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/DICOM
    )  
    #-----------------------------------------------------------------------------

    # Example 2
    #-----------

    add_executable(dicom_t2 ${PROJECT_SOURCE_DIR}/tests/DICOM/test_dicom2voxel.cpp)

    #target libs
    target_link_libraries(dicom_t2 pen_geos)

    #-----------------------------------------------------------------------------
    install(
    TARGETS dicom_t2
    DESTINATION ${PROJECT_SOURCE_DIR}/compiled/tests/DICOM
    )  
    #-----------------------------------------------------------------------------

endif()

add_subdirectory(geometry)
