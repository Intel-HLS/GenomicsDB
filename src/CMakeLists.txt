set(GenomicsDB_library_sources 
    ./genomicsdb/variant.cc
    ./genomicsdb/variant_query_config.cc
    ./genomicsdb/variant_field_handler.cc
    ./genomicsdb/variant_field_data.cc
    ./genomicsdb/variant_array_schema.cc
    ./genomicsdb/variant_storage_manager.cc
    ./genomicsdb/query_variants.cc
    ./genomicsdb/variant_cell.cc
    ./loader/tiledb_loader.cc
    ./loader/tiledb_loader_text_file.cc
    ./loader/tiledb_loader_file_base.cc
    ./loader/load_operators.cc
    ./query_operations/variant_operations.cc
    ./query_operations/broad_combined_gvcf.cc
    ./utils/memory_measure.cc
    ./utils/command_line.cc
    ./utils/histogram.cc
    ./utils/libtiledb_variant.cc
    ./utils/lut.cc
    ./utils/timer.cc
    ./utils/known_field_info.cc
    ./utils/json_config.cc
    ./utils/vid_mapper.cc
    ./vcf/vcf_adapter.cc
    ./vcf/vcf2binary.cc
    ./vcf/jni_bcf_reader.cc
    )

if(BUILD_JAVA)
    set(GenomicsDB_library_sources ${GenomicsDB_library_sources} ./java/JNI/src/genomicsdb_GenomicsDBQueryStream.cc)
endif()

#Compile sources with PIC
add_library(GenomicsDB_library_object_files OBJECT ${GenomicsDB_library_sources})
set_property(TARGET GenomicsDB_library_object_files PROPERTY POSITION_INDEPENDENT_CODE ON)

#Create the shared and static libraries
add_library(genomicsdb STATIC $<TARGET_OBJECTS:GenomicsDB_library_object_files>)
add_library(tiledbgenomicsdb SHARED $<TARGET_OBJECTS:GenomicsDB_library_object_files>)
target_link_libraries(tiledbgenomicsdb ${HTSLIB_LIBRARY} ${TILEDB_LIBRARY})
if(LIBCSV_FOUND)
    target_link_libraries(tiledbgenomicsdb ${LIBCSV_LIBRARY})
endif()
if(BUILD_DISTRIBUTABLE_LIBRARY)
    target_link_libraries(tiledbgenomicsdb ${OPENSSL_LIBRARIES})
endif()

install(TARGETS genomicsdb tiledbgenomicsdb LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

add_subdirectory(java)
