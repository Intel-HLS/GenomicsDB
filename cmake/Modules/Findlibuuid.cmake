# Determine compiler flags for libuuid
# Once done this will define
# LIBUUID_FOUND - libuuid found

find_path(LIBUUID_INCLUDE_DIR NAMES uuid/uuid.h HINTS "${LIBUUID_DIR}/include" "${LIBUUID_DIR}")

find_library(LIBUUID_LIBRARY NAMES uuid HINTS "${LIBUUID_DIR}/lib64" "${LIBUUID_DIR}/lib" "${LIBUUID_DIR}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(libuuid "Could not find libuuid headers and/or libraries ${DEFAULT_MSG}" LIBUUID_INCLUDE_DIR LIBUUID_LIBRARY)
