set(PROTOBUF_PROTO_FILES
    genomicsdb_callsets_mapping.proto
    genomicsdb_export_config.proto
    genomicsdb_import_config.proto
    genomicsdb_vid_mapping.proto
    )
if(PROTOBUF_REGENERATE)
    PROTOBUF_GENERATE_CPP(PROTOBUF_GENERATED_CXX_SRCS PROTOBUF_GENERATED_CXX_HDRS ${PROTOBUF_PROTO_FILES})
    set(PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
        CACHE INTERNAL "Path to protocol buffers generated C++ headers")
    add_custom_target(PROTOBUF_GENERATED_CXX_TARGET DEPENDS ${PROTOBUF_GENERATED_CXX_SRCS} ${PROTOBUF_GENERATED_CXX_HDRS})
    set(PROTOBUF_GENERATED_JAVA_SRCS
        ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/com/intel/genomicsdb/GenomicsDBCallsetsMapProto.java
        ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/com/intel/genomicsdb/GenomicsDBExportConfiguration.java
        ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/com/intel/genomicsdb/GenomicsDBImportConfiguration.java
        ${GENOMICSDB_MAVEN_BUILD_DIR}/protobuf/com/intel/genomicsdb/GenomicsDBVidMapProto.java
        )
else()
    set(PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/src/main/protobuf-generated/cpp/include")
    set(PROTOBUF_GENERATED_CXX_SRCS "")
    set(PROTOBUF_GENERATED_CXX_HDRS "")
    foreach(PROTO_FILE ${PROTOBUF_PROTO_FILES})
        get_filename_component(CURR_FILENAME ${PROTO_FILE} NAME_WE)
        list(APPEND PROTOBUF_GENERATED_CXX_SRCS "protobuf-generated/cpp/src/${CURR_FILENAME}.pb.cc")
        list(APPEND PROTOBUF_GENERATED_CXX_HDRS "${PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS}/${CURR_FILENAME}.pb.h")
    endforeach()
    set(PROTOBUF_GENERATED_JAVA_SRCS
        src/main/protobuf-generated/java/com/intel/genomicsdb/GenomicsDBCallsetsMapProto.java
        src/main/protobuf-generated/java/com/intel/genomicsdb/GenomicsDBExportConfiguration.java
        src/main/protobuf-generated/java/com/intel/genomicsdb/GenomicsDBImportConfiguration.java
        src/main/protobuf-generated/java/com/intel/genomicsdb/GenomicsDBVidMapProto.java
        )
endif()
set(PROTOBUF_GENERATED_CXX_SRCS ${PROTOBUF_GENERATED_CXX_SRCS} CACHE INTERNAL "Protocol buffers generated C++ sources")
set(PROTOBUF_GENERATED_CXX_HDRS ${PROTOBUF_GENERATED_CXX_HDRS} CACHE INTERNAL "Protocol buffers generated C++ headers")
set(PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS ${PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS} CACHE INTERNAL "Directory containing Protocol buffers generated C++ headers")
set(PROTOBUF_GENERATED_JAVA_SRCS ${PROTOBUF_GENERATED_JAVA_SRCS} CACHE INTERNAL "Protocol buffers generated Java sources")
