# Determine compiler flags for libdbi
# Once done this will define
# LIBDBI_FOUND - libdbi found

find_path(LIBDBI_INCLUDE_DIR NAMES dbi.h HINTS "${LIBDBI_DIR}/include/dbi")

find_library(LIBPGSQL_DRIVER_LIBRARY NAMES dbdpgsql HINTS "${LIBDBI_DIR}/lib/${CMAKE_C_LIBRARY_ARCHITECTURE}/dbd")

find_library(LIBDBI_DEV_LIBRARY NAMES dbi)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(libdbi "Could not find libdbi headers and/or libraries ${DEFAULT_MSG}" LIBDBI_INCLUDE_DIR LIBDBI_DEV_LIBRARY LIBPGSQL_DRIVER_LIBRARY)
