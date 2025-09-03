
#
#
#    Copyright (C) 2019 Universitat de València - UV
#    Copyright (C) 2019 Universitat Politècnica de València - UPV
#    Copyright (C) 2024-2025 Vicent Giménez Alventosa
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

# Add penRed DICOM lib
add_library(pen_dicom STATIC ${PROJECT_SOURCE_DIR}/kernel/parsers/DICOM/source/pen_dicom.cpp)

if(NOT FORCE_DICOM_DOWNLOAD)
    # Find the zlib package
    find_package(ZLIB QUIET)
endif()

if(FORCE_DICOM_DOWNLOAD OR NOT ZLIB_FOUND)

    message(STATUS "ZLIB not found, fetching it...")

    include(ExternalProject)

    # Define installation directories for ZLib
    set(ZLIB_INSTALL_DIR "${CMAKE_BINARY_DIR}/zlib_install")
    
    # Set CMake comon arguments
    set(ZLIB_CMAKE_COMMON_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=${ZLIB_INSTALL_DIR}   # Install location
        -DBUILD_SHARED_LIBS:BOOL=OFF                        # Build static library
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}       # Match the project's build type
        -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=${CMAKE_POSITION_INDEPENDENT_CODE} # To handle python modules and normal compilations
    )
    
    # macOS-specific configuration
    if (APPLE)
        if (DEFINED CMAKE_OSX_ARCHITECTURES AND NOT CMAKE_OSX_ARCHITECTURES STREQUAL "")
            list(APPEND ZLIB_CMAKE_COMMON_ARGS
                -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
            )
        endif()
    endif()
    

    # External project for ZLib
    ExternalProject_Add(
        zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG v1.3.1  # You can specify a version tag
        INSTALL_DIR "${ZLIB_INSTALL_DIR}"                # Installation directory
        CMAKE_CACHE_ARGS ${ZLIB_CMAKE_COMMON_ARGS}       # Arguments passed to the CMake configuration step
    )

    # Get the ZLib include and library paths for later use
    ExternalProject_Get_Property(zlib install_dir)
    set(ZLIB_INCLUDE_DIR "${install_dir}/include")
    if(WIN32)
        if (CMAKE_BUILD_TYPE STREQUAL "Release")
            set(ZLIB_LIBRARY "${install_dir}/lib/zlibstatic.lib")
        elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
            set(ZLIB_LIBRARY "${install_dir}/lib/zlibstaticd.lib")
        else()
            set(ZLIB_LIBRARY "${install_dir}/lib/zlibstatic.lib")
        endif()
    else()
        set(ZLIB_LIBRARY "${install_dir}/lib/libz.a")
    endif()

    # Print ZLib paths to verify setup
    message(STATUS "Set ZLib include directory: ${ZLIB_INCLUDE_DIR}")
    message(STATUS "Set ZLib static library: ${ZLIB_LIBRARY}")

    set(ZLIB_LIBRARIES ${ZLIB_LIBRARY})

else()
    message(STATUS "ZLIB found: ${ZLIB_INCLUDE_DIRS}; ${ZLIB_LIBRARIES}")
endif()

if(NOT FORCE_DICOM_DOWNLOAD)
    # Find the dcmtk package
    find_package(DCMTK QUIET)
endif()

# Check if it is found
if(FORCE_DICOM_DOWNLOAD OR NOT DCMTK_FOUND)
    message(STATUS "DCMTK library not found, trying to download it from github")

    # Set parallel level to build DCMTK
    if(NOT DEFINED CMAKE_BUILD_PARALLEL_LEVEL)
        set(CMAKE_BUILD_PARALLEL_LEVEL 4)
    endif()

    include(ExternalProject)

    # Define installation and build directories for dcmtk
    set(DCMTK_INSTALL_DIR "${CMAKE_BINARY_DIR}/dcmtk_install")

    # Select the desired DCMTK modules to compile
    set(DCMTK_MODULES
        ofstd
        oflog
        dcmdata
        dcmimgle            # Enable support for basic DICOM image processing
        dcmimage            # Enable support for advanced DICOM image processing
        dcmjpeg             # Enable support for JPEG DICOM images
        dcmjpls             # Enable support for JPEG-LS DICOM images
        # dcmtls            # Disable support for TLS/SSL
        # dcmnet            # Disable network capabilities
        dcmsr               # Enable support for Structured Reports
        # dcmwlm            # Disable worklist management (Modality worklist)
        # dcmqrdb           # Disable support for DB
        # dcmpstat          # Disable support for print management
        dcmrt              # Enable support for DICOM RT-Radiation Therapy
        dcmiod             # Enable support for IOD-Information Object Definitions
        dcmfg
        dcmseg             # Enable support for DICOM segmentation storage
        # dcmtract         # Disable support for tractography storage
        dcmpmap            # Enable support for parametric maps
        dcmect
    )

    message(STATUS "Selected DCMTK modules: ${DCMTK_MODULES}")

    # Define a command to set the parallel level
    if(MSVC)
        set(PARALLEL_FLAG /maxcpucount:${CMAKE_BUILD_PARALLEL_LEVEL})
    else()
        set(PARALLEL_FLAG -j${CMAKE_BUILD_PARALLEL_LEVEL})
    endif()

    set(DCMTK_CMAKE_COMMON_ARGS
        -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=${CMAKE_POSITION_INDEPENDENT_CODE} # To handle python modules and normal compilations
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}       # Match the project's build type
        -DCMAKE_INSTALL_PREFIX:PATH=${DCMTK_INSTALL_DIR}
        -DDCMTK_MODULES:STRING=${DCMTK_MODULES}
        -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
        -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
        -DBUILD_SHARED_LIBS:BOOL=OFF                             # Build static libraries
        -DBUILD_APPS:BOOL=OFF                                    # Disable apps build
        -DDCMTK_WITH_ZLIB:BOOL=ON
        -DDCMTK_WITH_ICU:BOOL=OFF
        -DDCMTK_WITH_OPENJPEG:BOOL=OFF
        -DDCMTK_WITH_OPENSSL:BOOL=OFF
        -DDCMTK_ENABLE_STL:BOOL=ON
        -DDCMTK_DEFAULT_DICT:STRING=builtin                      # Embedded dictionary        
    )

    # macOS-specific configuration
    if (APPLE)
        list(APPEND DCMTK_CMAKE_COMMON_ARGS
            -DDCMTK_WITH_ICONV:BOOL=OFF # Use built-in iconv
        )
        if (DEFINED CMAKE_OSX_ARCHITECTURES AND NOT CMAKE_OSX_ARCHITECTURES STREQUAL "")
            list(APPEND DCMTK_CMAKE_COMMON_ARGS
                -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
            )
        endif()
    endif()

    # Add DCMTK as an external project
    ExternalProject_Add(
        dcmtk
        GIT_REPOSITORY https://github.com/DCMTK/dcmtk.git
        GIT_TAG        DCMTK-3.6.9  # DCMTK version
        INSTALL_DIR "${DCMTK_INSTALL_DIR}"                 # Installation directory
        CMAKE_CACHE_ARGS ${DCMTK_CMAKE_COMMON_ARGS}
        BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> -- ${PARALLEL_FLAG}
    )

    # Get DCMTK install directory
    ExternalProject_Get_Property(dcmtk install_dir)

    # Add DCMTK include and library paths based on OS
    set(DCMTK_LOCAL_LIBRARY_PATHS ${install_dir}/lib64 ${install_dir}/lib)
    target_link_directories(pen_dicom PUBLIC ${DCMTK_LOCAL_LIBRARY_PATHS})
    set(DCMTK_INCLUDE_DIRS "${install_dir}/include")
    if(WIN32)
            #"Ws2_32" #Required for "gethostname" definition
            #"Iphlpapi.lib" #Required for "GetAdaptersInfo" definition
        set(DCMTK_LIBRARIES
            "dcmdata.lib"
            "oflog.lib"
            "ofstd.lib"
            "dcmdata.lib"
            "i2d.lib"
            "dcmimgle.lib"
            "dcmimage.lib"
            "dcmjpeg.lib"
            "ijg8.lib"
            "ijg12.lib"
            "ijg16.lib"
            "dcmjpls.lib"
            "dcmtkcharls.lib"
            # "dcmtls.lib"
            # "dcmnet.lib"
            "dcmsr.lib"
            "cmr.lib"
            # "dcmdsig.lib"
            # "dcmwlm.lib"
            # "dcmqrdb.lib"
            # "dcmpstat.lib"
            "dcmrt.lib"
            "dcmiod.lib"
            "dcmfg.lib"
            "dcmseg.lib"
            # "dcmtract.lib"
            "dcmpmap.lib"
            "dcmect.lib"
            "Ws2_32"         #Required windows specific libraries
            "Iphlpapi.lib"
            "Netapi32.lib"
        )
    elseif(UNIX)
        set(DCMTK_LIBRARIES
            "libdcmdata.a"
            "liboflog.a"
            "libofstd.a"
            "libi2d.a"
            "libdcmimgle.a"
            "libdcmimage.a"
            "libdcmjpeg.a"
            "libijg8.a"
            "libijg12.a"
            "libijg16.a"
            "libdcmjpls.a"
            "libdcmtkcharls.a"
            # "libdcmtls.a"
            # "libdcmnet.a"
            "libdcmsr.a"
            "libcmr.a"
            # "libdcmdsig.a"
            # "libdcmwlm.a"
            # "libdcmqrdb.a"
            # "libdcmpstat.a"
            "libdcmrt.a"
            "libdcmiod.a"
            "libdcmfg.a"
            "libdcmseg.a"
            # "libdcmtract.a"
            "libdcmpmap.a"
            "libdcmect.a"
        )
    endif()

    # Create a target to link DCMTK after it's built
    add_dependencies(pen_dicom dcmtk)

    message(STATUS "Set DCMTK include directory: ${DCMTK_INCLUDE_DIRS}")
    message(STATUS "Set DCMTK libraries: ${DCMTK_LIBRARIES}")

    # If zlib has been fetched, add it as dcmtk dependencies
    if(TARGET zlib)
        add_dependencies(dcmtk zlib)
    endif()
else()
    message(STATUS "DCMTK library found: ${DCMTK_INCLUDE_DIRS}")
endif()

target_include_directories(pen_dicom PUBLIC ${PROJECT_SOURCE_DIR}/kernel/parsers/DICOM/includes)
target_include_directories(pen_dicom PUBLIC ${PROJECT_SOURCE_DIR}/kernel/includes)

# Handle dirent header depending on OS
if (WIN32)
    message(STATUS "Set windows dirent")
    target_include_directories(pen_dicom PUBLIC ${PROJECT_SOURCE_DIR}/kernel/parsers/DICOM/includes/winDirent)
else ()
    message(STATUS "Use standard Unix dirent")
endif ()

# Link DCMTK lib
target_include_directories(pen_dicom PUBLIC ${DCMTK_INCLUDE_DIRS})
target_link_libraries(pen_dicom ${DCMTK_LIBRARIES} ${ZLIB_LIBRARIES})

# Link logger
target_link_libraries(pen_dicom pen_logger)
