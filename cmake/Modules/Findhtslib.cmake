# Determine compiler flags for htslib
# Once done this will define
# HTSLIB_FOUND - htslib found

include(FindPackageHandleStandardArgs)
#Build if htslib source directory specified
if(HTSLIB_SOURCE_DIR)
    set(HTSLIB_Debug_CFLAGS "-Wall -fPIC -DDEBUG  -g3 -gdwarf-3")  #will be picked if compiling in debug mode
    set(HTSLIB_Coverage_CFLAGS "${HTSLIB_Debug_CFLAGS}")
    set(HTSLIB_Release_CFLAGS " -Wall -fPIC -O3")
    set(HTSLIB_Debug_LDFLAGS "-g3 -gdwarf-3")
    set(HTSLIB_Coverage_LDFLAGS "${HTSLIB_Debug_LDFLAGS}")
    set(HTSLIB_Release_LDFLAGS "")
    include(ExternalProject)
    if(APPLE AND BUILD_DISTRIBUTABLE_LIBRARY)
        set(HTSLIB_EXTRA_CFLAGS -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET})
    endif()
    #Cross compiling for MacOSX
    if((NOT (CMAKE_SYSTEM_NAME STREQUAL CMAKE_HOST_SYSTEM_NAME)) AND APPLE)
        set(HTSLIB_OSXCROSS_COMPILE_FLAGS LIBS=${OSXCROSS_LIBS} CPPFLAGS=${OSXCROSS_CPPFLAGS} --host=${CMAKE_SYSTEM_PROCESSOR}-${CMAKE_SYSTEM})
    endif()
    ExternalProject_Add(
        htslib
        DOWNLOAD_COMMAND ""
        SOURCE_DIR "${HTSLIB_SOURCE_DIR}"
        UPDATE_COMMAND "autoreconf"
        PATCH_COMMAND ""
        CONFIGURE_COMMAND ${HTSLIB_SOURCE_DIR}/configure CFLAGS=${HTSLIB_${CMAKE_BUILD_TYPE}_CFLAGS} LDFLAGS=${HTSLIB_${CMAKE_BUILD_TYPE}_LDFLAGS}
            CC=${CMAKE_C_COMPILER} AR=${CMAKE_AR} RANLIB=${CMAKE_RANLIB}
            ${HTSLIB_OSXCROSS_COMPILE_FLAGS}
            --disable-lzma --disable-bz2 --disable-libcurl
        BUILD_COMMAND ${CMAKE_COMMAND} -E make_directory cram
            COMMAND ${CMAKE_COMMAND} -E make_directory test
            COMMAND $(MAKE) -f ${HTSLIB_SOURCE_DIR}/Makefile VPATH=${HTSLIB_SOURCE_DIR} SOURCE_DIR=${HTSLIB_SOURCE_DIR}
            AR=${CMAKE_AR}
        #BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
        )
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/vcf.h HINTS "${HTSLIB_SOURCE_DIR}" CMAKE_FIND_ROOT_PATH_BOTH)
    ExternalProject_Get_Property(htslib BINARY_DIR)
    set(HTSLIB_DIR_IN_BUILD_DIR "${BINARY_DIR}")
    set(HTSLIB_LIBRARY "${BINARY_DIR}/libhts.a")
    find_package_handle_standard_args(htslib "Could not find htslib headers ${DEFAULT_MSG}" HTSLIB_INCLUDE_DIR)
else()
    find_path(HTSLIB_INCLUDE_DIR NAMES htslib/vcf.h HINTS "${HTSLIB_INSTALL_DIR}")
    find_library(HTSLIB_LIBRARY NAMES libhts.a hts HINTS "${HTSLIB_INSTALL_DIR}")
    find_package_handle_standard_args(htslib "Could not find htslib headers and/or libraries ${DEFAULT_MSG}" HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY)
endif()
