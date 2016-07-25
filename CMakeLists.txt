project(GenomicsDB)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")

#Build parameters
set(DISABLE_MPI False CACHE BOOL "Disable use of any MPI compiler/libraries")
set(DISABLE_OPENMP False CACHE BOOL "Disable OpenMP")
set(BUILD_DISTRIBUTABLE_LIBRARY False CACHE BOOL "Build the TileDB/GenomicsDB library with minimal runtime dependencies")
set(BUILD_JAVA False CACHE BOOL "Build Java/JNI interface for combined VCF records")
set(HTSLIB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/htslib" CACHE PATH "Path to htslib source directory")
set(HTSLIB_INSTALL_DIR "" CACHE PATH "Path to htslib install directory")
set(OPENSSL_PREFIX_DIR "" CACHE PATH "OpenSSL installation prefix")
set(TILEDB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/TileDB" CACHE PATH "Path to TileDB source directory")
set(TILEDB_INSTALL_DIR "" CACHE PATH "Path to TileDB install directory")
set(USE_LIBCSV False CACHE BOOL "Disable library components that import data from csv files")
set(LIBCSV_DIR "" CACHE PATH "Path to libcsv header and library")
set(HTSJDK_JAR_PATH "" CACHE PATH "Path to pre-built htsjdk JAR")
set(DO_PROFILING False CACHE BOOL "Collect some stats during execution - doesn't add much overhead")
set(DO_MEMORY_PROFILING False CACHE BOOL "Collect memory consumption in parts of the combine gVCF program - high overhead")

#Platform check modules
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)

if(BUILD_DISTRIBUTABLE_LIBRARY)
    cmake_minimum_required(VERSION 3.4)
    set(DISABLE_MPI True)
    set(DISABLE_OPENMP True)
    #For the GNU compiler, link in static gcc libraries
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    endif()
endif()

#Warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-result")

#Check C++ 2011 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++11 test_cpp_2011)
if(NOT test_cpp_2011)
    message(FATAL_ERROR "Your compiler does not support C++ 2011, exiting")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

#MPI library/compiler etc
if(NOT DISABLE_MPI)
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

if(NOT DISABLE_OPENMP)
    find_package(OpenMPv4)
    if(NOT OPENMPV4_FOUND)
        set(DISABLE_OPENMP False)
    else()
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

add_definitions(-D_FILE_OFFSET_BITS=64)  #large file support
add_definitions(-DHTSDIR=1) #htslib is a mandatory requirement
add_definitions(-DDUPLICATE_CELL_AT_END=1) #mandatory

#RapidJSON
find_package(RapidJSON REQUIRED)
include_directories(${RAPIDJSON_INCLUDE_DIR})

#htslib
find_package(htslib REQUIRED)
include_directories(${HTSLIB_INCLUDE_DIR})

#TileDB
find_package(TileDB REQUIRED)
include_directories(${TILEDB_INCLUDE_DIR})

#librt
find_library(LIBRT_LIBRARY rt)

#JNI
if(BUILD_JAVA)
    find_package(JNI REQUIRED)
    include_directories(${JNI_INCLUDE_DIRS})
endif()

#libcsv
if(NOT BUILD_DISTRIBUTABLE_LIBRARY)
    if(USE_LIBCSV)
        find_package(libcsv REQUIRED)
    else()
        find_package(libcsv)
    endif()
endif()
if(LIBCSV_FOUND)
    include_directories(${LIBCSV_INCLUDE_DIR})
    set(USE_LIBCSV True)
    add_definitions(-DUSE_LIBCSV)
endif()

#Collect stats
if(DO_PROFILING)
    message(STATUS "Enabling stats collection")
    add_definitions(-DDO_PROFILING=1)
endif()

#Collect memory consumption stats while producing combined VCF records
if(DO_MEMORY_PROFILING)
    message(STATUS "Enabling memory consumption profiling while producing combined VCF records")
    add_definitions(-DDO_MEMORY_PROFILING=1)
endif()

#Produce verbose output
if(VERBOSE)
    add_definitions(-DVERBOSE=${VERBOSE})
    message(STATUS "Programs will produce output with verbosity=${VERBOSE}")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g3 -gdwarf-3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -fvisibility=hidden")
set(CMAKE_CXX_FLAGS_COVERAGE "-DDEBUG -g3 -gdwarf-3 --coverage")

#Clean all
add_custom_target(clean-all ${CMAKE_MAKE_PROGRAM} clean)
if(TILEDB_SOURCE_DIR)
    add_custom_target(clean-TileDB ${CMAKE_MAKE_PROGRAM} clean -C ${TILEDB_SOURCE_DIR})
    add_dependencies(clean-all clean-TileDB)
endif()
if(HTSLIB_SOURCE_DIR)
    add_custom_target(clean-htslib ${CMAKE_MAKE_PROGRAM} clean -C ${HTSLIB_SOURCE_DIR})
    add_dependencies(clean-all clean-htslib)
endif()

include_directories (
    include/genomicsdb
    include/loader
    include/query_operations
    include/utils
    include/vcf
)
if(BUILD_JAVA)
    include_directories(src/java/JNI/include)
endif()

add_subdirectory(src)
add_subdirectory(example)
