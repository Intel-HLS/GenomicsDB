# Determine compiler flags for TileDB
# Once done this will define
# TILEDB_FOUND - TileDB found

#Build if TileDB source directory specified
set(TileDB_Debug_BUILD_FLAGS "DEBUG=1")
#Build type
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(TileDB_BUILD_TYPE "debug")
else()
    set(TileDB_BUILD_TYPE "release")
endif()
#OpenMP
if(DISABLE_OPENMP)
    set(GNU_PARALLEL 0)
else()
    set(GNU_PARALLEL 1)
endif()
#OpenSSL
if(OPENSSL_PREFIX_DIR AND NOT OPENSSL_ROOT_DIR)
    set(OPENSSL_ROOT_DIR "${OPENSSL_PREFIX_DIR}")
endif()

if(BUILD_DISTRIBUTABLE_LIBRARY)
   set(OPENSSL_USE_STATIC_LIBS True)
endif()

find_package(OpenSSL REQUIRED)

list(APPEND CMAKE_ARGS "-DGNU_PARALLEL=${GNU_PARALLEL} -DCMAKE_BUILD_TYPE=${TileDB_BUILD_TYPE}")

#Build as external project
include(ExternalProject)
ExternalProject_Add(
    TileDB
    DOWNLOAD_COMMAND ""
    SOURCE_DIR "${TILEDB_SOURCE_DIR}"
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CMAKE_ARGS "${CMAKE_ARGS}"
    INSTALL_COMMAND ""
)
ExternalProject_Add_Step(
    TileDB Make
    COMMAND ${CMAKE_MAKE_COMMAND} OPENSSL_PREFIX_DIR=${OPENSSL_ROOT_DIR} GNU_PARALLEL=${GNU_PARALLEL} BUILD=${TileDB_BUILD_TYPE}
)

set(TileDB_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/dependencies/TileDB/core/include/c_api")
set(TileDB_LIBRARIES "${CMAKE_SHARED_LIBRARY_PREFIX}tiledb${CMAKE_SHARED_LIBRARY_PREFIX}")
include_directories(${TileDB_INCLUDE_DIRS})
#    find_path(TILEDB_INCLUDE_DIR NAMES tiledb.h HINTS "${TILEDB_SOURCE_DIR}/core/include/c_api")
#    set(TILEDB_LIBRARY "${TILEDB_SOURCE_DIR}/core/lib/${TileDB_BUILD_TYPE}/libtiledb.a")
#    add_custom_target(${TILEDB_LIBRARY} ${CMAKE_MAKE_COMMAND} TileDB)
#    find_package_handle_standard_args(TileDB "Could not find TileDB headers ${DEFAULT_MSG}" TILEDB_INCLUDE_DIR)
