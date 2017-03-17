# Determine compiler flags for htslib
# Once done this will define
# HTSLIB_FOUND - htslib found

include(FindPackageHandleStandardArgs)
#Build if htslib source directory specified
if(HTSLIB_SOURCE_DIR)
    set(HTSLIB_Debug_BUILD_FLAGS "DEBUG=1")  #will be picked if compiling in debug mode
    include(ExternalProject)
    if(APPLE AND BUILD_DISTRIBUTABLE_LIBRARY)
        set(HTSLIB_EXTRA_CFLAGS -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET})
    endif()
    ExternalProject_Add(
        htslib
        DOWNLOAD_COMMAND ""
        SOURCE_DIR "${HTSLIB_SOURCE_DIR}"
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND $(MAKE) ${HTSLIB_${CMAKE_BUILD_TYPE}_BUILD_FLAGS} CPPFLAGS=${HTSLIB_EXTRA_CFLAGS}
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
        )
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/vcf.h HINTS "${HTSLIB_SOURCE_DIR}")
    set(HTSLIB_LIBRARY "${HTSLIB_SOURCE_DIR}/libhts.a")
    find_package_handle_standard_args(htslib "Could not find htslib headers ${DEFAULT_MSG}" HTSLIB_INCLUDE_DIR)
else()
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/vcf.h HINTS "${HTSLIB_INSTALL_DIR}")
    find_library(HTSLIB_LIBRARY NAMES libhts.a hts HINTS "${HTSLIB_INSTALL_DIR}")
    find_package_handle_standard_args(htslib "Could not find htslib headers and/or libraries ${DEFAULT_MSG}" HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY)
endif()
