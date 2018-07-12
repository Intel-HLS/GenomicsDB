# Determine compiler flags for TileDB
# Once done this will define
# TILEDB_FOUND - TileDB found

#Disable Master Catalog in TileDB
set(ENABLE_MASTER_CATALOG False)
set(TILEDB_VERBOSE True)

#Zlib
find_package(ZLIB REQUIRED)

#OpenSSL
if(OPENSSL_PREFIX_DIR AND NOT OPENSSL_ROOT_DIR)
    set(OPENSSL_ROOT_DIR "${OPENSSL_PREFIX_DIR}")
endif()
if(BUILD_DISTRIBUTABLE_LIBRARY)
    set(OPENSSL_USE_STATIC_LIBS True)
endif()
find_package(OpenSSL REQUIRED) #now performed inside TileDB

#libuuid
find_package(libuuid REQUIRED)

include(FindPackageHandleStandardArgs)

#Build if TileDB source directory specified
if(TILEDB_SOURCE_DIR)
    #OpenMP
    if(DISABLE_OPENMP)
        set(USE_OPENMP False CACHE BOOL "Disable OpenMP" FORCE)
    else()
        set(USE_OPENMP True CACHE BOOL "Enable OpenMP" FORCE)
    endif()

    find_path(TILEDB_INCLUDE_DIR NAMES tiledb.h HINTS "${TILEDB_SOURCE_DIR}/core/include/c_api")
    find_package_handle_standard_args(TileDB "Could not find TileDB headers ${DEFAULT_MSG}" TILEDB_INCLUDE_DIR)
    add_subdirectory(${TILEDB_SOURCE_DIR} EXCLUDE_FROM_ALL)
else()
    find_path(TILEDB_INCLUDE_DIR NAMES "tiledb.h" HINTS "${TILEDB_INSTALL_DIR}/include")
    find_library(TILEDB_LIBRARY NAMES libtiledb.a tiledb HINTS "${TILEDB_INSTALL_DIR}" "${TILEDB_INSTALL_DIR}/lib64" "${TILEDB_INSTALL_DIR}/lib")
    find_package_handle_standard_args(TileDB "Could not find TileDB headers and/or libraries ${DEFAULT_MSG}" TILEDB_INCLUDE_DIR TILEDB_LIBRARY)
endif()
