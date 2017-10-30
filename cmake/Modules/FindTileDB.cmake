# Determine compiler flags for TileDB
# Once done this will define
# TILEDB_FOUND - TileDB found

#Disable Master Catalog in TileDB
list(APPEND TILEDB_CMAKE_ARGS "-DENABLE_MASTER_CATALOG:BOOL=False")

#Zlib
find_package(ZLIB REQUIRED)
if(ZLIB_ROOT)
    list(APPEND TILEDB_CMAKE_ARGS "-DZLIB_ROOT:PATH=${ZLIB_ROOT}")
endif()

#OpenSSL
if(OPENSSL_PREFIX_DIR AND NOT OPENSSL_ROOT_DIR)
    set(OPENSSL_ROOT_DIR "${OPENSSL_PREFIX_DIR}")
endif()
if(OPENSSL_ROOT_DIR)
    list(APPEND TILEDB_CMAKE_ARGS "-DOPENSSL_ROOT_DIR:PATH=${OPENSSL_ROOT_DIR}")
endif()
if(BUILD_DISTRIBUTABLE_LIBRARY)
    set(OPENSSL_USE_STATIC_LIBS True)
    list(APPEND TILEDB_CMAKE_ARGS "-DOPENSSL_USE_STATIC_LIBS:BOOL=True")
endif()
find_package(OpenSSL REQUIRED) #now performed inside TileDB

#libuuid
if(LIBUUID_DIR)
    list(APPEND TILEDB_CMAKE_ARGS "-DLIBUUID_DIR:PATH=${LIBUUID_DIR}")
endif()
find_package(libuuid REQUIRED)

include(FindPackageHandleStandardArgs)

#Build if TileDB source directory specified
if(TILEDB_SOURCE_DIR)
    #OpenMP
    if(DISABLE_OPENMP)
        set(TILEDB_USE_OPENMP False)
    else()
        set(TILEDB_USE_OPENMP True)
    endif()
    list(APPEND TILEDB_CMAKE_ARGS "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}")
    list(APPEND TILEDB_CMAKE_ARGS "-DUSE_OPENMP:BOOL=${TILEDB_USE_OPENMP}")
    list(APPEND TILEDB_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}")

    #Build as external project
    include(ExternalProject)
    ExternalProject_Add(
        TileDB
        DOWNLOAD_COMMAND ""
        SOURCE_DIR "${TILEDB_SOURCE_DIR}"
        CMAKE_ARGS "${TILEDB_CMAKE_ARGS}"
        BINARY_DIRECTORY "${CMAKE_BINARY_DIR}"
        )
    #since the headers and libraries are available only after make completes in the TileDB dir
    set(TILEDB_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
    set(TILEDB_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/libtiledb.a")
    #find_path(TILEDB_INCLUDE_DIR NAMES tiledb.h HINTS "${CMAKE_INSTALL_PREFIX}/include")
    #add_custom_target(${TILEDB_LIBRARY} ${CMAKE_MAKE_COMMAND} TileDB)
    #find_package_handle_standard_args(TileDB "Could not find TileDB headers ${DEFAULT_MSG}" TILEDB_INCLUDE_DIR)
else()
    find_path(TILEDB_INCLUDE_DIR NAMES "tiledb.h" HINTS "${TILEDB_INSTALL_DIR}/include")
    find_library(TILEDB_LIBRARY NAMES libtiledb.a tiledb HINTS "${TILEDB_INSTALL_DIR}" "${TILEDB_INSTALL_DIR}/lib64" "${TILEDB_INSTALL_DIR}/lib")
    find_package_handle_standard_args(TileDB "Could not find TileDB headers and/or libraries ${DEFAULT_MSG}" TILEDB_INCLUDE_DIR TILEDB_LIBRARY)
endif()
