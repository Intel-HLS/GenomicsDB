project(GenomicsDB)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

#git tag
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

#Build parameters
set(GENOMICSDB_RELEASE_VERSION "0.6.4" CACHE STRING "GenomicsDB release version") #used in Maven builds
set(GENOMICSDB_VERSION "${GENOMICSDB_RELEASE_VERSION}-${GIT_COMMIT_HASH}" CACHE STRING "GenomicsDB full version string")
set(DISABLE_MPI False CACHE BOOL "Disable use of any MPI compiler/libraries")
set(DISABLE_OPENMP False CACHE BOOL "Disable OpenMP")
set(BUILD_DISTRIBUTABLE_LIBRARY False CACHE BOOL "Build the TileDB/GenomicsDB library with minimal runtime dependencies")
set(BUILD_JAVA False CACHE BOOL "Build Java/JNI interface for combined VCF records")
set(HTSLIB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/htslib" CACHE PATH "Path to htslib source directory")
set(HTSLIB_INSTALL_DIR "" CACHE PATH "Path to htslib install directory")
set(OPENSSL_PREFIX_DIR "" CACHE PATH "OpenSSL installation prefix")
set(GENOMICSDB_PROTOBUF_VERSION "3.0.2" CACHE STRING "Version of Google Protocol Buffer library")
set(PROTOBUF_LIBRARY "" CACHE PATH "Path to protobuf headers and library")
set(PROTOBUF_STATIC_LINKING True CACHE BOOL "Statically link Protobuf libraries")
set(PROTOBUF_REGENERATE False CACHE BOOL "Regenerate protocol buffers C++ files")
set(TILEDB_SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/TileDB" CACHE PATH "Path to TileDB source directory")
set(TILEDB_INSTALL_DIR "" CACHE PATH "Path to TileDB install directory")
set(USE_LIBCSV False CACHE BOOL "Disable library components that import data from csv files")
set(LIBCSV_DIR "" CACHE PATH "Path to libcsv header and library")
set(DO_PROFILING False CACHE BOOL "Collect some stats during execution - doesn't add much overhead")
set(DO_MEMORY_PROFILING False CACHE BOOL "Collect memory consumption in parts of the combine gVCF program - high overhead")
set(GENOMICSDB_MAVEN_BUILD_DIR ${CMAKE_BINARY_DIR}/target CACHE PATH "Path to maven build directory")
set(LIBDBI_DIR "" CACHE PATH "Path to libdbi install directory")
set(MAVEN_QUIET False CACHE BOOL "Do not print mvn messages")
set(GENOMICSDB_SPARK_PROFILE "" CACHE STRING "Profile to choose Apache Spark version") # used in Maven builds

#See https://cmake.org/Wiki/CMake_RPATH_handling#Common_questions
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH True)
enable_testing()

#Platform check modules
include(CheckIncludeFileCXX)
include(CheckCXXSymbolExists)

if(BUILD_DISTRIBUTABLE_LIBRARY)
    cmake_minimum_required(VERSION 3.4)
    set(DISABLE_MPI True)
    set(DISABLE_OPENMP True)
    #For the GNU compiler, link in static gcc libraries
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
    endif()
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
endif()

#Check C++ 2011 support
include(CheckAndSetStdCXX2011Flag)
CHECK_AND_SET_STD_CXX_2011_FLAG(test_cpp_2011)
if(NOT test_cpp_2011)
    message(FATAL_ERROR "Your compiler does not support C++ 2011, exiting")
endif()

#Warnings
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-result")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-unused-variable")
endif()

#MPI library/compiler etc
if(NOT DISABLE_MPI)
    find_package(MPI REQUIRED)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
else()
    add_definitions(-DDISABLE_MPI=1)
endif()

if(NOT DISABLE_OPENMP)
    find_package(OpenMPv4)
    if(NOT OPENMPV4_FOUND)
        set(DISABLE_OPENMP True)
    else()
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()
if(DISABLE_OPENMP)
    add_definitions(-DDISABLE_OPENMP=1)
endif()

add_definitions(-D_FILE_OFFSET_BITS=64)  #large file support
add_definitions(-DHTSDIR=1) #htslib is a mandatory requirement
add_definitions(-DDUPLICATE_CELL_AT_END=1) #mandatory
add_definitions(-DGENOMICSDB_VERSION=\"${GENOMICSDB_VERSION}\")

#RapidJSON
find_package(RapidJSON REQUIRED)
include_directories(${RAPIDJSON_INCLUDE_DIR})

#htslib
find_package(htslib REQUIRED)
include_directories(${HTSLIB_INCLUDE_DIR})

#TileDB
find_package(TileDB REQUIRED)
include_directories(${TILEDB_INCLUDE_DIR})

#librt
find_library(LIBRT_LIBRARY rt)

#pgsql driver and dbi libs
find_package(libdbi)
if(LIBDBI_FOUND)
    include_directories(${LIBDBI_INCLUDE_DIR})
    add_definitions(-DLIBDBI)
endif()

#Protobuf library
find_package(ProtobufWrapper REQUIRED)
include_directories(${PROTOBUF_INCLUDE_DIRS})

#JNI
if(BUILD_JAVA)
    find_package(Java 1.8 REQUIRED)
    find_package(JNI REQUIRED)
    include_directories(${JNI_INCLUDE_DIRS})
    include_directories(src/main/jni/include)
    include(UseJava)
endif()

#libcsv
if(NOT BUILD_DISTRIBUTABLE_LIBRARY)
    if(USE_LIBCSV)
        find_package(libcsv REQUIRED)
    else()
        find_package(libcsv)
    endif()
endif()
if(LIBCSV_FOUND)
    include_directories(${LIBCSV_INCLUDE_DIR})
    set(USE_LIBCSV True)
    add_definitions(-DUSE_LIBCSV)
endif()

#Collect stats
if(DO_PROFILING)
    message(STATUS "Enabling stats collection")
    add_definitions(-DDO_PROFILING=1)
endif()

#Collect memory consumption stats while producing combined VCF records
if(DO_MEMORY_PROFILING)
    message(STATUS "Enabling memory consumption profiling while producing combined VCF records")
    add_definitions(-DDO_MEMORY_PROFILING=1)
endif()

#Produce verbose output
if(VERBOSE)
    add_definitions(-DVERBOSE=${VERBOSE})
    message(STATUS "Programs will produce output with verbosity=${VERBOSE}")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -g3 -gdwarf-3")
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -O3 -fvisibility=hidden")
set(CMAKE_CXX_FLAGS_COVERAGE "-DDEBUG -g3 -gdwarf-3 --coverage")

#Clean all
add_custom_target(clean-all ${CMAKE_MAKE_PROGRAM} clean)
if(TILEDB_SOURCE_DIR)
    add_custom_target(clean-TileDB ${CMAKE_MAKE_PROGRAM} clean -C ${TILEDB_SOURCE_DIR})
    add_dependencies(clean-all clean-TileDB)
endif()
if(HTSLIB_SOURCE_DIR)
    add_custom_target(clean-htslib ${CMAKE_MAKE_PROGRAM} clean -C ${HTSLIB_SOURCE_DIR})
    add_dependencies(clean-all clean-htslib)
endif()

include_directories (
    src/main/cpp/include/genomicsdb
    src/main/cpp/include/loader
    src/main/cpp/include/query_operations
    src/main/cpp/include/utils
    src/main/cpp/include/vcf
    src/test/cpp/include
)

function(build_GenomicsDB_executable_common target)
    #add_executable(${target} src/${target}.cc)
    if(TILEDB_SOURCE_DIR)
        add_dependencies(${target} TileDB)
    endif()
    if(HTSLIB_SOURCE_DIR)
        add_dependencies(${target} htslib)
    endif()
    target_link_libraries(${target} genomicsdb ${TILEDB_LIBRARY} ${HTSLIB_LIBRARY})
    if(LIBCSV_FOUND)
        target_link_libraries(${target} ${LIBCSV_LIBRARY})
    endif()
    if(LIBDBI_FOUND)
        target_link_libraries(${target} ${LIBPGSQL_DRIVER_LIBRARY} ${LIBDBI_DEV_LIBRARY})
    endif()
    if(LIBRT_LIBRARY)
        target_link_libraries(${target} ${LIBRT_LIBRARY})
    endif()
    target_link_libraries(${target} ${PROTOBUF_LIBRARIES} ${OPENSSL_LIBRARIES} ${ZLIB_LIBRARIES})
    install(TARGETS ${target} RUNTIME DESTINATION bin)
endfunction()

function(build_GenomicsDB_executable target)
    add_executable(${target} src/${target}.cc)
    build_GenomicsDB_executable_common(${target})
endfunction()

add_subdirectory(src)
include_directories(${PROTOBUF_GENERATED_CXX_HDRS_INCLUDE_DIRS})
add_subdirectory(tools)
add_subdirectory(example)

if(BUILD_JAVA)
    set(JAVA_SCALA_SOURCES
        ./src/main/java/com/intel/genomicsdb/GenomicsDBConfiguration.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBQueryStream.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBJavaSparkFactory.java
        ./src/main/java/com/intel/genomicsdb/SilentByteBufferStream.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBImporter.java
        ./src/main/java/com/intel/genomicsdb/SampleInfo.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBTimer.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBPartitionInfo.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBFeatureReader.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBInputSplit.java
        ./src/main/java/com/intel/genomicsdb/MultiChromosomeIterator.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBException.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBImporterStreamWrapper.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBInputFormat.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBRecordReader.java
        ./src/main/java/com/intel/genomicsdb/ChromosomeInterval.java
        ./src/main/java/com/intel/genomicsdb/GenomicsDBUtils.java
        ./src/test/java/com/intel/genomicsdb/GenomicsDBVidMappingProtoSpec.java
        ./src/test/java/com/intel/genomicsdb/GenomicsDBImportConfigurationSpec.java
        ./src/test/java/com/intel/genomicsdb/GenomicsDBImporterSpec.java
        ./src/main/scala/com/intel/genomicsdb/GenomicsDBPartition.scala
        ./src/main/scala/com/intel/genomicsdb/GenomicsDBScalaSparkFactory.scala
        ./src/main/scala/com/intel/genomicsdb/GenomicsDBContext.scala
        ./src/main/scala/com/intel/genomicsdb/GenomicsDBRDD.scala
        )
    if(PROTOBUF_REGENERATE)
        #Must be set here - see https://cmake.org/cmake/help/v3.3/command/set_source_files_properties.html
        set_source_files_properties(${PROTOBUF_GENERATED_JAVA_SRCS} PROPERTIES GENERATED True)
        #Protobuf files are generated by Maven
        set(MAVEN_PROTOBUF_REGENERATE_ARGS "-Dprotoc.filepath=${PROTOBUF_PROTOC_EXECUTABLE}" "-Dprotobuf.regenerate=true")
    else()
        set(JAVA_SCALA_SOURCES ${JAVA_SCALA_SOURCES} ${PROTOBUF_GENERATED_JAVA_SRCS})
        set(MAVEN_PROTOBUF_REGENERATE_ARGS "")
    endif()
    set(MAVEN_QUIET_ARGS "")
    if(MAVEN_QUIET)
        set(MAVEN_QUIET_ARGS "--quiet")
    endif()
    set(MAVEN_ARGS
        ${MAVEN_QUIET_ARGS}
        -Dgenomicsdb.version=${GENOMICSDB_RELEASE_VERSION}
        -Dgenomicsdb_source_directory=${CMAKE_SOURCE_DIR}
        -Dgenomicsdb_build_directory=${GENOMICSDB_MAVEN_BUILD_DIR}
        -Dgenomicsdb_lib_directory=$<TARGET_FILE_DIR:tiledbgenomicsdb>
        -Dprotobuf.version=${GENOMICSDB_PROTOBUF_VERSION}
        ${MAVEN_PROTOBUF_REGENERATE_ARGS})
    #Maven build - depends on dynamic library
    add_custom_command(
        OUTPUT ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}.jar ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}-jar-with-dependencies.jar
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/pom.xml ${CMAKE_BINARY_DIR}/pom.xml
        COMMAND mvn versions:set ${MAVEN_QUIET_ARGS} -DnewVersion=${GENOMICSDB_RELEASE_VERSION} -P${GENOMICSDB_SPARK_PROFILE}
        COMMAND mvn package -DskipTests
          ${MAVEN_ARGS}
        DEPENDS tiledbgenomicsdb ${JAVA_SCALA_SOURCES} pom.xml
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    install(FILES ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}.jar ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}-jar-with-dependencies.jar DESTINATION bin)

    add_test(NAME javatests
        COMMAND mvn test
          ${MAVEN_ARGS}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    set(GENOMICSDB_EXAMPLE_SOURCES
        ${CMAKE_SOURCE_DIR}/example/java/TestGenomicsDB.java
        ${CMAKE_SOURCE_DIR}/example/java/TestBufferStreamGenomicsDBImporter.java
        ${CMAKE_SOURCE_DIR}/example/java/TestGenomicsDBImporterWithMergedVCFHeader.java)

    add_jar(genomicsdb-${GENOMICSDB_RELEASE_VERSION}-examples
        SOURCES ${GENOMICSDB_EXAMPLE_SOURCES}
        INCLUDE_JARS ${GENOMICSDB_MAVEN_BUILD_DIR}/genomicsdb-${GENOMICSDB_RELEASE_VERSION}-jar-with-dependencies.jar
        OUTPUT_DIR ${GENOMICSDB_MAVEN_BUILD_DIR})

    add_test(NAME CI_tests
        COMMAND /usr/bin/python ${CMAKE_SOURCE_DIR}/tests/run.py
        ${CMAKE_BINARY_DIR}
        ${CMAKE_INSTALL_PREFIX})

endif()
