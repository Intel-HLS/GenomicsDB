# Determine compiler flags for libdbi
# Once done this will define
# LIBDBI_FOUND - libdbi found

find_path(LIBDBI_INCLUDE_DIR NAMES dbi/dbi.h HINTS "${LIBDBI_DIR}/include")

find_library(LIBPGSQL_DRIVER_LIBRARY NAMES dbdpgsql HINTS "${LIBDBI_DIR}/lib/dbd")

find_library(LIBDBI_DEV_LIBRARY NAMES dbi HINTS "${LIBDBI_DIR}/lib")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(libdbi "Could not find libdbi headers and/or libraries ${DEFAULT_MSG}" LIBDBI_INCLUDE_DIR LIBDBI_DEV_LIBRARY LIBPGSQL_DRIVER_LIBRARY)

if(LIBDBI_FOUND)
    include(CheckCSourceCompiles)
    file(READ ${CMAKE_SOURCE_DIR}/cmake/Modules/libdbi_test_program.c LIBDBI_C_TEST_SOURCE)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "-I${LIBDBI_INCLUDE_DIR}")
    set(SAFE_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES ${LIBPGSQL_DRIVER_LIBRARY} ${LIBDBI_DEV_LIBRARY})
    check_c_source_compiles("${LIBDBI_C_TEST_SOURCE}" LIBDBI_TEST_PROGRAM_COMPILES)
    if(NOT LIBDBI_TEST_PROGRAM_COMPILES)
        message(STATUS "libdbi headers and libraries found; however, test program fails to compile. GenomicsDB requires libdbi >= 0.9.0.")
        message(STATUS "Mapping DB support disabled")
        unset(LIBDBI_FOUND)
    endif()
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_LIBRARIES "${SAFE_CMAKE_REQUIRED_LIBRARIES}")
endif()
