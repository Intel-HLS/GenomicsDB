set(GenomicsDB_java_sources
    genomicsdb/GenomicsDBException.java
    genomicsdb/GenomicsDBFeatureReader.java
    genomicsdb/GenomicsDBQueryStream.java
    genomicsdb/GenomicsDBUtils.java
    )

if(BUILD_JAVA) 
    #Copy libtiledbgenomicsdb into genomicsdb directory
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/genomicsdb)
    set(copy_libtiledbgenomicsdb_path genomicsdb/$<TARGET_FILE_NAME:tiledbgenomicsdb>)
    add_custom_target(copy_libtiledbgenomicsdb
        COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:tiledbgenomicsdb> ${copy_libtiledbgenomicsdb_path}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    add_dependencies(copy_libtiledbgenomicsdb tiledbgenomicsdb)
    #Create JAR with Java files
    find_package(Java 1.8 REQUIRED)
    include(UseJava)
    #htsjdk dependency
    find_package(htsjdk) 
    set(CMAKE_JAVA_INCLUDE_PATH "${CMAKE_JAVA_INCLUDE_PATH}:${HTSJDK_JAR_PATH}")
    #JAR name
    set(CMAKE_JAVA_TARGET_OUTPUT_NAME "genomicsdb")
    add_jar(genomicsdb_jar SOURCES ${GenomicsDB_java_sources})
    get_property(genomicsdb_jar_path TARGET genomicsdb_jar PROPERTY JAR_FILE)
    #Add the native library
    add_custom_command(TARGET genomicsdb_jar POST_BUILD
        COMMAND jar uf ${genomicsdb_jar_path} ${copy_libtiledbgenomicsdb_path}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
       )
   add_dependencies(genomicsdb_jar copy_libtiledbgenomicsdb)
   install(FILES ${genomicsdb_jar_path} DESTINATION lib)
endif()
