# Determine path to IPP libraries for Zlib
# Once done this will define
# LIBIPP_FOUND - IPP found

find_library(LIBIPPDC_LIBRARY NAMES libippdc.a ippdc HINTS "${IPPROOT}/lib/intel64/")
find_library(LIBIPPS_LIBRARY NAMES libipps.a ipps HINTS "${IPPROOT}/lib/intel64/")
find_library(LIBIPPCORE_LIBRARY NAMES libippcore.a ippcore HINTS "${IPPROOT}/lib/intel64/")
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPP "Could not find IPP libraries ${DEFAULT_MSG}" LIBIPPDC_LIBRARY LIBIPPS_LIBRARY LIBIPPCORE_LIBRARY)
