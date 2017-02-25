
##########
# Macros #
##########

# --- Release Version --- #
RELEASE_VERSION = 0.4.0

# --- Build Flags --- #
# Large file support
LFS_CFLAGS = -D_FILE_OFFSET_BITS=64

CFLAGS=-Wall -Wno-reorder -Wno-unknown-pragmas -Wno-unused-variable \
           -Wno-unused-but-set-variable -Wno-unused-result
JAVA_BUILD_FLAGS=
# LINKFLAGS appear before the object file list
# in the link command (e.g. -fopenmp, -O3)
LINKFLAGS:=

# LDFLAGS appear after the list of object files (-lz etc)
LDFLAGS:=
ifdef MAXIMIZE_STATIC_LINKING
  LINKFLAGS+=-static-libgcc -static-libstdc++
  LDFLAGS+=-Wl,-Bstatic -lcrypto -Wl,-Bdynamic
else
  LDFLAGS+=-lcrypto
endif
LDFLAGS+= -lz -lrt
SHARED_LIBRARY_EXTENSION:=so
SHARED_LIBRARY_FLAGS:=-shared

# --- OS-specific Build Flags ---#
OS := $(shell uname)
#Only build shared library on MacOS
ifeq ($(OS), Darwin)
  OPENSSL_PREFIX_DIR?=/usr/local/opt/openssl/
  CFLAGS=-mmacosx-version-min=10.9
  LINKFLAGS:=
  ifdef MAXIMIZE_STATIC_LINKING
    LDFLAGS:=$(OPENSSL_PREFIX_DIR)/lib/libcrypto.a -lz
  else
    LDFLAGS:=-L$(OPENSSL_PREFIX_DIR)/lib -lcrypto -lz
  endif
  SHARED_LIBRARY_EXTENSION=dylib
  SHARED_LIBRARY_FLAGS:=-dynamiclib -mmacosx-version-min=10.9
  DISABLE_OPENMP:=1
endif

GNU_PARALLEL=0
ifdef OPENMP
  ifndef DISABLE_OPENMP
    CFLAGS+=-fopenmp
    GNU_PARALLEL=1
    LINKFLAGS+=-fopenmp
  endif
endif

# --- Debug/Release/Verbose mode handler --- #
BUILD ?= release
VERBOSE ?= 0
HTSLIB_BUILD=

ifeq ($(BUILD),debug)
  CFLAGS+= -gdwarf-3 -g3 -DDEBUG
  LINKFLAGS+=-gdwarf-3 -g3
  TILEDB_BUILD:=debug
  HTSLIB_BUILD="DEBUG=1"
  JAVA_BUILD_FLAGS+=-g
endif

ifeq ($(BUILD),debug-coverage)
  CFLAGS+= -gdwarf-3 -g3 -DDEBUG --coverage
  LINKFLAGS+=-gdwarf-3 -g3 --coverage
  TILEDB_BUILD:=debug
  HTSLIB_BUILD="DEBUG=1"
endif

ifeq ($(BUILD),release)
  CFLAGS += -DNDEBUG -O3 -fvisibility=hidden
  LINKFLAGS+=-O3
  TILEDB_BUILD:=release
endif

ifdef DISABLE_MPI
  CC = gcc
  CXX = g++
  CFLAGS += -DDISABLE_MPI
else
  # MPI compiler for C++
  ifdef MPIPATH
    CC  = $(MPIPATH)/mpicc
    CXX = $(MPIPATH)/mpicxx
  else
    MPIPATH=
    CC  = mpicc
    CXX = mpicxx
  endif
endif

ifdef DISABLE_OPENMP
  CFLAGS+=-DDISABLE_OPENMP
endif

CPPFLAGS=-std=c++11 -fPIC $(LFS_CFLAGS) $(CFLAGS)

# In the current version, this is mandatory
CPPFLAGS += -DDUPLICATE_CELL_AT_END

# --- TileDB Source --- #
ifndef TILEDB_DIR
  TILEDB_DIR=dependencies/TileDB
endif
CPPFLAGS+=-I$(TILEDB_DIR)/core/include/c_api
ifeq ($(OS), Darwin)
  LDFLAGS:= $(TILEDB_DIR)/core/lib/$(TILEDB_BUILD)/libtiledb.a $(LDFLAGS)
else
  LDFLAGS:= -Wl,-Bstatic -L$(TILEDB_DIR)/core/lib/$(TILEDB_BUILD) -ltiledb -Wl,-Bdynamic $(LDFLAGS)
endif

# --- Htslib Source --- #
HTSLIB_BUILD_NUM_THREADS ?= 1
HTSLIB_EXTRA_CFLAGS=
ifndef HTSDIR
  HTSDIR=dependencies/htslib
endif
ifdef HTSDIR
  CPPFLAGS+=-I$(HTSDIR) -DHTSDIR
  ifeq ($(OS), Darwin)
    LDFLAGS:=$(HTSDIR)/libhts.a $(LDFLAGS)
    HTSLIB_EXTRA_CFLAGS=-mmacosx-version-min=10.9
  else
    LDFLAGS:=-Wl,-Bstatic -L$(HTSDIR) -lhts -Wl,-Bdynamic $(LDFLAGS)
  endif
endif

# --- RapidJSON Source (header only library) --- #
ifndef RAPIDJSON_INCLUDE_DIR
  RAPIDJSON_INCLUDE_DIR=dependencies/RapidJSON/include
endif
CPPFLAGS+=-I$(RAPIDJSON_INCLUDE_DIR)

#libcsv - optional, but required if csvs need to be imported
ifdef LIBCSV_DIR
  CPPFLAGS+=-DUSE_LIBCSV -I$(LIBCSV_DIR)
  LDFLAGS+=-L$(LIBCSV_DIR)/.libs -L$(LIBCSV_DIR)/lib -lcsv
else
  ifdef USE_LIBCSV
    CPPFLAGS+=-DUSE_LIBCSV
    LDFLAGS+=-lcsv
  endif
endif

# --- JNI flag - optional, but required if the JNI library is needed --- #
ifdef JNI_FLAGS
  CPPFLAGS+=$(JNI_FLAGS)
endif

# --- BigMPI Flags (optional) --- #
ifdef USE_BIGMPI
  CPPFLAGS+=-I$(USE_BIGMPI)/src -DUSE_BIGMPI
  LDFLAGS+=-L$(USE_BIGMPI)/src -lbigmpi
endif

ifdef DO_PROFILING
  CPPFLAGS+=-DDO_PROFILING
endif

ifdef DO_MEMORY_PROFILING
  CPPFLAGS+=-DDO_MEMORY_PROFILING
endif

# --- Google Performance Tools Library (optional) --- #
ifdef USE_GPERFTOOLS
  ifdef GPERFTOOLSDIR
    CPPFLAGS+=-DUSE_GPERFTOOLS -I$(GPERFTOOLSDIR)/include
    LDFLAGS += -Wl,-Bstatic -L$(GPERFTOOLSDIR)/lib -lprofiler -Wl,-Bdynamic  -lunwind
  endif
endif

ifdef VERBOSE
  CPPFLAGS+= -DVERBOSE=$(VERBOSE)
endif

LDFLAGS+=-lprotobuf

# --- Directories --- #

GENOMICSDB_OBJ_DIR=./obj
GENOMICSDB_BIN_DIR=./bin

# --- Header directories --- #
GENOMICSDB_LIBRARY_INCLUDE_DIRS:=\
  src/main/cpp/include/genomicsdb \
  src/main/cpp/include/loader \
  src/main/cpp/include/query_operations \
  src/main/cpp/include/utils \
  src/main/cpp/include/vcf \
  src/main/jni/include \
  example/include \
  tools/include

CPPFLAGS+=$(GENOMICSDB_LIBRARY_INCLUDE_DIRS:%=-I%)

# 'vpath' to know which directories to search for sources
vpath %.cc src/main/cpp/src/genomicsdb:src/main/cpp/src/loader:src/main/cpp/src/query_operations:src/main/cpp/src/utils:src/main/cpp/src/vcf:src/main/jni/src:example/src:tools/src

EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
vpath %.h %(subst $(SPACE),:,$(GENOMICSDB_LIBRARY_INCLUDE_DIRS)) :

# --- Files --- #

GENOMICSDB_LIBRARY_SOURCES:= \
			    vcf_adapter.cc \
			    json_config.cc \
			    vid_mapper.cc \
			    vid_mapper_pb.cc \
			    libtiledb_variant.cc \
			    variant_cell.cc \
			    variant_query_config.cc \
			    variant_field_handler.cc \
			    variant_field_data.cc \
			    variant.cc \
			    histogram.cc \
			    lut.cc \
			    known_field_info.cc \
			    vcf2binary.cc \
			    command_line.cc \
			    variant_array_schema.cc \
			    tiledb_loader.cc \
			    broad_combined_gvcf.cc \
			    variant_operations.cc \
			    load_operators.cc \
			    variant_storage_manager.cc \
			    query_variants.cc \
			    tiledb_loader_file_base.cc \
			    tiledb_loader_text_file.cc \
			    genomicsdb_bcf_generator.cc \
			    timer.cc \
			    memory_measure.cc \
			    genomicsdb_importer.cc \
			    genomicsdb_import_config.pb.cc \
			    genomicsdb_export_config.pb.cc \
			    genomicsdb_vid_mapping.pb.cc \
			    genomicsdb_callsets_mapping.pb.cc

ifdef BUILD_JAVA

  # --- Jars --- #
  GENOMICSDB_JAR = target/genomicsdb-$(RELEASE_VERSION).jar
  GENOMICSDB_JAR_WITH_DEPS = target/genomicsdb-$(RELEASE_VERSION)-jar-with-dependencies.jar

  # --- Java/Scala Source Directories --- #
  SCALA_SRCDIR = src/main/scala/com/intel/genomicsdb
  SCALA_TESTDIR = src/test/scala/com/intel/genomicsdb
  SCALA_SRC_SUBDIRS = $(wildcard $(SCALA_SRCDIR)/*)
  SCALA_TEST_SUBDIRS = $(wildcard $(SCALA_TESTDIR)/*)
  JAVA_SRCDIR = src/main/java/com/intel/genomicsdb
  JAVA_TESTDIR = src/test/java/com/intel/genomicsdb
  JAVA_SRC_SUBDIRS = $(wildcard $(JAVA_SRCDIR)/*)
  JAVA_TEST_SUBDIRS = $(wildcard $(JAVA_TESTDIR)/*)

  # --- C++ Sources --- #
  GENOMICSDB_LIBRARY_SOURCES:= $(GENOMICSDB_LIBRARY_SOURCES) \
  	genomicsdb_GenomicsDBQueryStream.cc \
  	genomicsdb_GenomicsDBImporter.cc \
    genomicsdb_jni_init.cc

  # --- Java/Scala Sources --- #
  SCALA_SRC = $(wildcard $(foreach D,$(SCALA_SRC_SUBDIRS),$D/*.scala))
  SCALA_SRC += $(wildcard $(SCALA_SRCDIR)/*.scala)
  SCALA_TEST_SRC = $(wildcard $(foreach D,$(SCALA_TEST_SUBDIRS),$D/*.scala))
  SCALA_TEST_SRC += $(wildcard $(SCALA_TESTDIR)/*.scala)

  JAVA_SRC = $(wildcard $(foreach D,$(JAVA_SRC_SUBDIRS),$D/*.java))
  JAVA_SRC += $(wildcard $(JAVA_SRCDIR)/*.java)
  JAVA_TEST_SRC = $(wildcard $(foreach D,$(JAVA_TEST_SUBDIRS),$D/*.java))
  JAVA_TEST_SRC += $(wildcard $(JAVA_TESTDIR)/*.java)

  # --- Spark Submit shell script --- #
  GENOMICSDB_SPARK_SUBMIT_SCRIPT=src/resources/genomicsdb-spark-submit.sh
endif

GENOMICSDB_EXAMPLE_SOURCES:= \
    			create_tiledb_workspace.cc \
			    vcf2tiledb.cc \
			    vcfdiff.cc \
			    example_libtiledb_variant_driver.cc \
			    vcf_histogram.cc \
			    gt_mpi_gather.cc \
			    test_genomicsdb_bcf_generator.cc \
			    test_genomicsdb_importer.cc

ALL_GENOMICSDB_SOURCES := $(GENOMICSDB_LIBRARY_SOURCES) $(GENOMICSDB_EXAMPLE_SOURCES)

GENOMICSDB_LIBRARY_OBJ_FILES := $(patsubst %.cc, $(GENOMICSDB_OBJ_DIR)/%.o, $(GENOMICSDB_LIBRARY_SOURCES))
GENOMICSDB_EXAMPLE_OBJ_FILES := $(patsubst %.cc, $(GENOMICSDB_OBJ_DIR)/%.o, $(GENOMICSDB_EXAMPLE_SOURCES))
GENOMICSDB_EXAMPLE_BIN_FILES := $(patsubst %.cc, $(GENOMICSDB_BIN_DIR)/%, $(GENOMICSDB_EXAMPLE_SOURCES))

ALL_GENOMICSDB_OBJ_FILES:=$(GENOMICSDB_LIBRARY_OBJ_FILES) $(GENOMICSDB_EXAMPLE_OBJ_FILES)
ALL_GENOMICSDB_HEADER_DEPENDENCIES = $(ALL_GENOMICSDB_OBJ_FILES:%.o=%.d)

GENOMICSDB_STATIC_LIBRARY:=$(GENOMICSDB_BIN_DIR)/libgenomicsdb.a
GENOMICSDB_SHARED_LIBRARY_BASENAME:=libtiledbgenomicsdb.$(SHARED_LIBRARY_EXTENSION)
GENOMICSDB_SHARED_LIBRARY:=$(GENOMICSDB_BIN_DIR)/$(GENOMICSDB_SHARED_LIBRARY_BASENAME)

#Put GENOMICSDB_STATIC_LIBRARY as first component of LDFLAGS
ifeq ($(OS), Darwin)
  LDFLAGS:=$(GENOMICSDB_BIN_DIR)/libgenomicsdb.a $(LDFLAGS)
else
  LDFLAGS:=-Wl,-Bstatic -L$(GENOMICSDB_BIN_DIR) -lgenomicsdb -Wl,-Bdynamic $(LDFLAGS)
endif

###################
# General Targets #
###################

.PHONY: all genomicsdb_library clean clean-dependencies clean-all \
        TileDB_library TileDB_clean htslib_library htslib_clean

ALL_BUILD_TARGETS:= genomicsdb_library
ifndef DISABLE_MPI
  ALL_BUILD_TARGETS += $(GENOMICSDB_EXAMPLE_BIN_FILES)
endif
ifdef BUILD_JAVA
  ALL_BUILD_TARGETS += $(GENOMICSDB_JAR)
endif

all: $(ALL_BUILD_TARGETS)

genomicsdb_library: $(GENOMICSDB_STATIC_LIBRARY) $(GENOMICSDB_SHARED_LIBRARY)

clean:
	rm -rf $(GENOMICSDB_BIN_DIR)/* $(GENOMICSDB_OBJ_DIR)/*
	mvn clean -Dgenomicsdb.version=$(RELEASE_VERSION)

clean-dependencies: TileDB_clean htslib_clean

clean-all: clean clean-dependencies

#TileDB library
TileDB_library:
	$(MAKE) -C $(TILEDB_DIR) MPIPATH=$(MPIPATH) BUILD=$(TILEDB_BUILD) GNU_PARALLEL=$(GNU_PARALLEL) \
	  OPENSSL_PREFIX_DIR=$(OPENSSL_PREFIX_DIR)

TileDB_clean:
	$(MAKE) -C $(TILEDB_DIR) clean

$(TILEDB_DIR)/core/lib/$(TILEDB_BUILD)/libtiledb.a:
	$(MAKE) -C $(TILEDB_DIR) MPIPATH=$(MPIPATH) BUILD=$(TILEDB_BUILD) GNU_PARALLEL=$(GNU_PARALLEL)

#htslib library
htslib_library:
	$(MAKE) -C $(HTSDIR) $(HTSLIB_BUILD) CPPFLAGS=$(HTSLIB_EXTRA_CFLAGS)

htslib_clean:
	$(MAKE) -C $(HTSDIR) clean

$(HTSDIR)/libhts.a:
	$(MAKE) -C $(HTSDIR) $(HTSLIB_BUILD) CPPFLAGS=$(HTSLIB_EXTRA_CFLAGS)

# --- Compilation and dependency generation --- #

-include $(ALL_GENOMICSDB_HEADER_DEPENDENCIES)

.PRECIOUS: $(GENOMICSDB_OBJ_DIR)/%.o

#All object files
$(GENOMICSDB_OBJ_DIR)/%.o: %.cc
	@mkdir -p $(GENOMICSDB_OBJ_DIR) 
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) -MMD -c $< -o $@

#Library
$(GENOMICSDB_STATIC_LIBRARY): $(GENOMICSDB_LIBRARY_OBJ_FILES)
	@mkdir -p $(GENOMICSDB_BIN_DIR)
	@echo "Creating static library $@"
	@ar rcs $@ $^

$(GENOMICSDB_SHARED_LIBRARY): $(GENOMICSDB_LIBRARY_OBJ_FILES) \
    $(TILEDB_DIR)/core/lib/$(TILEDB_BUILD)/libtiledb.a $(HTSDIR)/libhts.a
	@mkdir -p $(GENOMICSDB_BIN_DIR)
	@echo "Creating dynamic library $@"
	@$(CXX) $(LINKFLAGS) $(SHARED_LIBRARY_FLAGS) -o $@ $^ $(LDFLAGS)

$(GENOMICSDB_JAR): $(GENOMICSDB_SHARED_LIBRARY) $(JAVA_SRC) $(JAVA_TEST_SRC) $(SCALA_SRC) \
                   $(SCALA_TEST_SRC)
	@echo "Creating GenomicsDB jar file"
	@mvn package -DskipTests -Dgenomicsdb.version=$(RELEASE_VERSION)
	@echo "Copying GenomicsDB jars to bin/"
	@cp $(GENOMICSDB_JAR) $(GENOMICSDB_JAR_WITH_DEPS) $(GENOMICSDB_BIN_DIR)
	@cp $(GENOMICSDB_SPARK_SUBMIT_SCRIPT) bin/


# --- GenomicsDB examples --- #

#Linking
$(GENOMICSDB_BIN_DIR)/%: $(GENOMICSDB_OBJ_DIR)/%.o $(GENOMICSDB_STATIC_LIBRARY) \
    $(TILEDB_DIR)/core/lib/$(TILEDB_BUILD)/libtiledb.a $(HTSDIR)/libhts.a
	@mkdir -p $(GENOMICSDB_BIN_DIR)
	@echo "Linking example $@"
	@$(CXX) $(LINKFLAGS) -o $@ $< $(LDFLAGS)
