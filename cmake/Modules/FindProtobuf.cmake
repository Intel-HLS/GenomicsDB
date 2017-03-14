# Following the standard CMake FindProtobuf module
# Determine compiler flags for protobuf
# Once done this will define
# PROTOBUF_LIBRARY_FOUND - protobuf found

find_path(PROTOBUF_INCLUDE_DIRS NAMES google/protobuf/service.h PATHS "${PROTOBUF_INCLUDE_DIR}" "${PROTOBUF_LIBRARY}/include")
if(PROTOBUF_STATIC_LINKING OR BUILD_DISTRIBUTABLE_LIBRARY)
    find_library(PROTOBUF_LIBRARIES NAMES libprotobuf.a protobuf PATHS "${PROTOBUF_LIBRARY}" "${PROTOBUF_LIBRARY}/lib64" "${PROTOBUF_LIBRARY}/lib")
else()
    find_library(PROTOBUF_LIBRARIES NAMES protobuf PATHS "${PROTOBUF_LIBRARY}" "${PROTOBUF_LIBRARY}/lib64" "${PROTOBUF_LIBRARY}/lib")
endif()
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Protobuf "Could not find Protobuf headers and/or libraries\n${DEFAULT_MSG}" PROTOBUF_INCLUDE_DIRS PROTOBUF_LIBRARIES)
