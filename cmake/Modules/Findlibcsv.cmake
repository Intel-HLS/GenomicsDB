# Determine compiler flags for libcsv
# Once done this will define
# LIBCSV_FOUND - libcsv found

find_path(LIBCSV_INCLUDE_DIR NAMES csv.h HINTS "${LIBCSV_DIR}")
find_library(LIBCSV_LIBRARY NAMES csv HINTS "${LIBCSV_DIR}" "${LIBCSV_DIR}/.libs")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libcsv "Could not find libcsv headers and/or libraries ${DEFAULT_MSG}" LIBCSV_INCLUDE_DIR LIBCSV_LIBRARY)
