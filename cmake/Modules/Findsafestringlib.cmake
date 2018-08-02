# Determine compiler flags for libuuid
# Once done this will define
# SAFESTRINGLIB_FOUND - libuuid found

find_path(SAFESTRINGLIB_INCLUDE_DIR NAMES safe_mem_lib.h HINTS "${SAFESTRINGLIB_DIR}/include")
find_library(SAFESTRINGLIB_LIBRARY NAMES libsafestring.a safestring HINTS "${SAFESTRINGLIB_DIR}/lib64"
  "${SAFESTRINGLIB_DIR}/lib" "${SAFESTRINGLIB_DIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(safestringlib "Could not find safestring headers and/or libraries ${DEFAULT_MSG}" SAFESTRINGLIB_INCLUDE_DIR
  SAFESTRINGLIB_LIBRARY)
