# Determine compiler flags for libdbi
# Once done this will define
# LIBDBI_FOUND - libdbi found

find_path(LIBDBI_INCLUDE_DIR NAMES dbi.h HINTS "/usr/include/dbi")
find_library(LIBDBI_DEV_LIBRARY NAMES dbi HINTS "/usr/lib/x86_64-linux-gnu")
find_library(LIBDBI1_LIBRARY NAMES dbi HINTS "/usr/lib/x86_64-linux-gnu")
find_library(LIBPGSQL_DRIVER_LIBRARY NAMES dbdpgsql HINTS "/usr/lib/x86_64-linux-gnu/dbd")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libcsv "Could not find libcsv headers and/or libraries ${DEFAULT_MSG}" LIBDBI_INCLUDE_DIR LIBDBI_DEV_LIBRARY LIBDBI1_LIBRARY LIBPGSQL_DRIVER_LIBRARY)
