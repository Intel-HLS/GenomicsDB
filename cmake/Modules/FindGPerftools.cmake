# Determine compiler flags for gperftools
# Once done this will define
# GPERFTOOLS_FOUND - gperftools found

find_path(GPERFTOOLS_INCLUDE_DIR NAMES gperftools/profiler.h HINTS "${GPERFTOOLS_DIR}" "${GPERFTOOLS_DIR}/include")
find_library(GPERFTOOLS_PROFILER_LIBRARY NAMES profiler HINTS "${GPERFTOOLS_DIR}" "${GPERFTOOLS_DIR}/lib64" "${GPERFTOOLS_DIR}/lib")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(gperftools "Could not find gperftools headers and/or libraries ${DEFAULT_MSG}" GPERFTOOLS_INCLUDE_DIR GPERFTOOLS_PROFILER_LIBRARY)
