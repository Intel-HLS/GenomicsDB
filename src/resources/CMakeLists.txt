PROTOBUF_GENERATE_CPP(PROTOBUF_GENERATED_CXX_SRCS PROTOBUF_GENERATED_CXX_HDRS
    genomicsdb_callsets_mapping.proto
    genomicsdb_export_config.proto
    genomicsdb_import_config.proto
    genomicsdb_vid_mapping.proto)
set(PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}
    CACHE INTERNAL "Path to protocol buffers generated C++ headers")
set(PROTOBUF_GENERATED_CXX_SRCS ${PROTOBUF_GENERATED_CXX_SRCS} CACHE INTERNAL "Protocol buffers generated C++ sources")
set(PROTOBUF_GENERATED_CXX_HDRS ${PROTOBUF_GENERATED_CXX_HDRS} CACHE INTERNAL "Protocol buffers generated C++ headers")
add_custom_target(PROTOBUF_GENERATED_CXX_TARGET DEPENDS ${PROTOBUF_GENERATED_CXX_SRCS} ${PROTOBUF_GENERATED_CXX_HDRS})
