set(GenomicsDB_library_sources 
    cpp/src/query_operations/variant_operations.cc
    cpp/src/query_operations/broad_combined_gvcf.cc
    cpp/src/genomicsdb/variant_cell.cc
    cpp/src/genomicsdb/variant_storage_manager.cc
    cpp/src/genomicsdb/variant_field_data.cc
    cpp/src/genomicsdb/variant_array_schema.cc
    cpp/src/genomicsdb/variant_field_handler.cc
    cpp/src/genomicsdb/variant.cc
    cpp/src/genomicsdb/variant_query_config.cc
    cpp/src/genomicsdb/query_variants.cc
    cpp/src/genomicsdb/genomicsdb_columnar_field.cc
    cpp/src/genomicsdb/genomicsdb_iterators.cc
    cpp/src/genomicsdb/genomicsdb_multid_vector_field.cc
    cpp/src/loader/tiledb_loader_text_file.cc
    cpp/src/loader/load_operators.cc
    cpp/src/loader/genomicsdb_importer.cc
    cpp/src/loader/tiledb_loader_file_base.cc
    cpp/src/loader/tiledb_loader.cc
    cpp/src/utils/command_line.cc
    cpp/src/utils/memory_measure.cc
    cpp/src/utils/histogram.cc
    cpp/src/utils/json_config.cc
    cpp/src/utils/vid_mapper_pb.cc
    cpp/src/utils/lut.cc
    cpp/src/utils/known_field_info.cc
    cpp/src/utils/vid_mapper.cc
    cpp/src/utils/vid_mapper_sql.cc
    cpp/src/utils/timer.cc
    cpp/src/vcf/vcf_adapter.cc
    cpp/src/vcf/genomicsdb_bcf_generator.cc
    cpp/src/vcf/vcf2binary.cc
    )

include_directories(${PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS})
set(GenomicsDB_library_sources
    ${GenomicsDB_library_sources}
    ${PROTOBUF_GENERATED_CXX_SRCS}
    )

if(PROTOBUF_REGENERATE)
    #Must be set here - see https://cmake.org/cmake/help/v3.3/command/set_source_files_properties.html
    set_source_files_properties(${PROTOBUF_GENERATED_CXX_SRCS} ${PROTOBUF_GENERATED_CXX_HDRS} PROPERTIES GENERATED True)
endif()

if(BUILD_JAVA)
    set(GenomicsDB_library_sources ${GenomicsDB_library_sources}
        jni/src/genomicsdb_GenomicsDBImporter.cc
        jni/src/genomicsdb_GenomicsDBQueryStream.cc
        jni/src/genomicsdb_jni_init.cc
        )
endif()

#Compile sources with PIC
add_library(GenomicsDB_library_object_files OBJECT ${GenomicsDB_library_sources})
add_dependencies(GenomicsDB_library_object_files TileDB)
if(PROTOBUF_REGENERATE)
    add_dependencies(GenomicsDB_library_object_files PROTOBUF_GENERATED_CXX_TARGET)
endif()
set_property(TARGET GenomicsDB_library_object_files PROPERTY POSITION_INDEPENDENT_CODE ON)

#Create the shared and static libraries
add_library(genomicsdb STATIC $<TARGET_OBJECTS:GenomicsDB_library_object_files>)
#Link in other libraries
add_library(tiledbgenomicsdb SHARED $<TARGET_OBJECTS:GenomicsDB_library_object_files>)
if(TILEDB_SOURCE_DIR)
    add_dependencies(tiledbgenomicsdb TileDB)
endif()
if(HTSLIB_SOURCE_DIR)
    add_dependencies(tiledbgenomicsdb htslib)
endif()
target_link_libraries(tiledbgenomicsdb ${HTSLIB_LIBRARY} ${TILEDB_LIBRARY} ${GENOMICSDB_EXTERNAL_DEPENDENCIES_LIBRARIES})
if(LIBRT_LIBRARY)
    target_link_libraries(tiledbgenomicsdb ${LIBRT_LIBRARY})
endif()

if(LIBCSV_FOUND)
    target_link_libraries(tiledbgenomicsdb ${LIBCSV_LIBRARY})
endif()
#if(BUILD_DISTRIBUTABLE_LIBRARY)
    #target_link_libraries(tiledbgenomicsdb ${OPENSSL_LIBRARIES})
#endif()

install(TARGETS genomicsdb tiledbgenomicsdb LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
