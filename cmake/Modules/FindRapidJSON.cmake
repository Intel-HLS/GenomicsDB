# Determine compiler flags for RapidJSON
# Once done this will define
# RAPIDJSON_FOUND - RapidJSON found
# RapidJSON_C_FLAGS
# RapidJSON_CXX_FLAGS

#If specified by user
if(RAPIDJSON_INCLUDE_DIR)
    set(RAPIDJSON_INCLUDE_DIR_HINT "${RAPIDJSON_INCLUDE_DIR}")
    unset(RAPIDJSON_INCLUDE_DIR CACHE)
endif()
find_path(RAPIDJSON_INCLUDE_DIR NAMES rapidjson/rapidjson.h HINTS "${RAPIDJSON_INCLUDE_DIR_HINT}"
    PATHS "${CMAKE_SOURCE_DIR}/dependencies/RapidJSON/include" CMAKE_FIND_ROOT_PATH_BOTH)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RapidJSON "Could not find RapidJSON header: rapidjson/rapidjson.h - specify the variable RAPIDJSON_INCLUDE_DIR to point to the directory <RapidJSON>/include" RAPIDJSON_INCLUDE_DIR)
