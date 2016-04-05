##########
# Macros #
##########

# Large file support
LFS_CFLAGS = -D_FILE_OFFSET_BITS=64

CFLAGS=-Wall -Wno-reorder -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable
#LINKFLAGS appear before the object file list in the link command (e.g. -fopenmp, -O3)
LINKFLAGS=
#LDFLAGS appear after the list of object files (-lz etc)
LDFLAGS:=-lz -lrt -lcrypto

ifdef OPENMP
  CFLAGS+=-fopenmp
endif
LINKFLAGS+=-fopenmp

# --- Debug/Release/Verbose mode handler --- #
BUILD ?= debug
VERBOSE ?= 0
HTSLIB_BUILD=

ifeq ($(BUILD),debug)
  CFLAGS+= -gdwarf-3 -g3 -DDEBUG
  LINKFLAGS+=-gdwarf-3 -g3
  TILEDB_BUILD:=debug
  HTSLIB_BUILD="DEBUG=1"
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

# MPI compiler for C++
ifdef MPIPATH
    CC  = $(MPIPATH)/mpicc
    CXX = $(MPIPATH)/mpicxx
else
    MPIPATH=
    CC  = mpicc
    CXX = mpicxx
endif
CPPFLAGS=-std=c++11 -fPIC $(LFS_CFLAGS) $(CFLAGS)

#In the current version, this is mandatory
CPPFLAGS += -DDUPLICATE_CELL_AT_END

#TileDB source
TILEDB_BUILD_NUM_THREADS ?= 16
ifndef TILEDB_DIR
    TILEDB_DIR=dependencies/TileDB
endif
CPPFLAGS+=-I$(TILEDB_DIR)/core/include/c_api
LDFLAGS:= -Wl,-Bstatic -L$(TILEDB_DIR)/core/lib/$(TILEDB_BUILD) -ltiledb -Wl,-Bdynamic $(LDFLAGS)

#htslib
HTSLIB_BUILD_NUM_THREADS ?= 8
ifndef HTSDIR
    HTSDIR=dependencies/htslib
endif
ifdef HTSDIR
    CPPFLAGS+=-I$(HTSDIR) -DHTSDIR
    LDFLAGS:=-Wl,-Bstatic -L$(HTSDIR) -lhts -Wl,-Bdynamic $(LDFLAGS)
endif

#RapidJSON - header only library
ifndef RAPIDJSON_INCLUDE_DIR
    RAPIDJSON_INCLUDE_DIR=dependencies/RapidJSON/include
endif
CPPFLAGS+=-I$(RAPIDJSON_INCLUDE_DIR)

#libcsv - optional, but required if csvs need to be imported
ifdef LIBCSV_DIR
    CPPFLAGS+=-DUSE_LIBCSV -I$(LIBCSV_DIR)
    LDFLAGS+=-L$(LIBCSV_DIR)/.libs -lcsv
else
    ifdef USE_LIBCSV
	CPPFLAGS+=-DUSE_LIBCSV
	LDFLAGS+=-lcsv
    endif
endif

#BigMPI - optional
ifdef USE_BIGMPI
    CPPFLAGS+=-I$(USE_BIGMPI)/src -DUSE_BIGMPI
    LDFLAGS+=-L$(USE_BIGMPI)/src -lbigmpi
endif

ifdef DO_PROFILING
    CPPFLAGS+=-DDO_PROFILING
endif

#Google performance tools library - optional
ifdef USE_GPERFTOOLS
    ifdef GPERFTOOLSDIR
	CPPFLAGS+=-DUSE_GPERFTOOLS -I$(GPERFTOOLSDIR)/include
	LDFLAGS += -Wl,-Bstatic -L$(GPERFTOOLSDIR)/lib -lprofiler -Wl,-Bdynamic  -lunwind
    endif
endif

ifdef VERBOSE
    CPPFLAGS+= -DVERBOSE=$(VERBOSE)
endif

# --- Directories --- #

GENOMICSDB_OBJ_DIR=./obj
GENOMICSDB_BIN_DIR=./bin

#Header directories
GENOMICSDB_LIBRARY_INCLUDE_DIRS:=include/genomicsdb include/loader include/query_operations include/utils include/vcf \
    example/include
CPPFLAGS+=$(GENOMICSDB_LIBRARY_INCLUDE_DIRS:%=-I%)

#Using vpath to let Makefile know which directories to search for sources
#For sources
vpath %.cc src/genomicsdb:src/loader:src/query_operations:src/utils:src/vcf:example/src

EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
vpath %.h %(subst $(SPACE),:,$(GENOMICSDB_LIBRARY_INCLUDE_DIRS)) :

# --- Files --- #

GENOMICSDB_LIBRARY_SOURCES:= \
			    vcf_adapter.cc \
			    json_config.cc \
			    vid_mapper.cc \
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
			    tiledb_loader_text_file.cc

GENOMICSDB_EXAMPLE_SOURCES:= \
    			    create_tiledb_workspace.cc \
			    gt_verifier.cc \
			    example_ga4gh_scan_operator.cc \
			    vcf2tiledb.cc \
			    vcfdiff.cc \
			    example_libtiledb_variant_driver.cc \
			    vcf_histogram.cc \
			    gt_profile_query_processor.cc \
			    gt_example_query_processor.cc \
			    gt_mpi_gather.cc

ALL_GENOMICSDB_SOURCES := $(GENOMICSDB_LIBRARY_SOURCES) $(GENOMICSDB_EXAMPLE_SOURCES)

GENOMICSDB_LIBRARY_OBJ_FILES := $(patsubst %.cc, $(GENOMICSDB_OBJ_DIR)/%.o, $(GENOMICSDB_LIBRARY_SOURCES))

GENOMICSDB_EXAMPLE_OBJ_FILES := $(patsubst %.cc, $(GENOMICSDB_OBJ_DIR)/%.o, $(GENOMICSDB_EXAMPLE_SOURCES))
GENOMICSDB_EXAMPLE_BIN_FILES := $(patsubst %.cc, $(GENOMICSDB_BIN_DIR)/%, $(GENOMICSDB_EXAMPLE_SOURCES))

ALL_GENOMICSDB_OBJ_FILES:=$(GENOMICSDB_LIBRARY_OBJ_FILES) $(GENOMICSDB_EXAMPLE_OBJ_FILES)
ALL_GENOMICSDB_HEADER_DEPENDENCIES = $(ALL_GENOMICSDB_OBJ_FILES:%.o=%.d)

GENOMICSDB_STATIC_LIBRARY:=$(GENOMICSDB_BIN_DIR)/libgenomicsdb.a
GENOMICSDB_SHARED_LIBRARY:=$(GENOMICSDB_BIN_DIR)/libgenomicsdb.so

#Put GENOMICSDB_STATIC_LIBRARY as first component of LDFLAGS
LDFLAGS:=-Wl,-Bstatic -L$(GENOMICSDB_BIN_DIR) -lgenomicsdb -Wl,-Bdynamic $(LDFLAGS)

###################
# General Targets #
###################

.PHONY: all genomicsdb_library clean clean-dependencies clean-all TileDB_library TileDB_clean htslib_library htslib_clean

all: genomicsdb_library $(GENOMICSDB_EXAMPLE_BIN_FILES)

genomicsdb_library: $(GENOMICSDB_STATIC_LIBRARY) $(GENOMICSDB_SHARED_LIBRARY)

clean:
	rm -rf $(GENOMICSDB_BIN_DIR)/* $(GENOMICSDB_OBJ_DIR)/*

clean-dependencies: TileDB_clean htslib_clean

clean-all: clean clean-dependencies

#TileDB library
TileDB_library:
	make -C $(TILEDB_DIR) MPIPATH=$(MPIPATH) BUILD=$(TILEDB_BUILD) -j $(TILEDB_BUILD_NUM_THREADS)

TileDB_clean:
	make -C $(TILEDB_DIR) clean

$(TILEDB_DIR)/core/lib/$(TILEDB_BUILD)/libtiledb.a:
	make -C $(TILEDB_DIR) MPIPATH=$(MPIPATH) BUILD=$(TILEDB_BUILD) -j $(TILEDB_BUILD_NUM_THREADS)

#htslib library
htslib_library:
	make -C $(HTSDIR) $(HTSLIB_BUILD) -j $(HTSLIB_BUILD_NUM_THREADS)

htslib_clean:
	make -C $(HTSDIR) clean

$(HTSDIR)/libhts.a:
	make -C $(HTSDIR) $(HTSLIB_BUILD) -j $(HTSLIB_BUILD_NUM_THREADS)

# --- Compilation and dependency genration --- #

-include $(ALL_GENOMICSDB_HEADER_DEPENDENCIES)

#All object files
$(GENOMICSDB_OBJ_DIR)/%.o: %.cc
	@mkdir -p $(GENOMICSDB_OBJ_DIR) 
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) -MMD -c $< -o $@

#Library
$(GENOMICSDB_BIN_DIR)/libgenomicsdb.a: $(GENOMICSDB_LIBRARY_OBJ_FILES)
	@mkdir -p $(GENOMICSDB_BIN_DIR)
	@echo "Creating static library $@"
	@ar rcs $@ $^

$(GENOMICSDB_BIN_DIR)/libgenomicsdb.so: $(GENOMICSDB_LIBRARY_OBJ_FILES)
	@mkdir -p $(GENOMICSDB_BIN_DIR)
	@echo "Creating dynamic library $@"
	@$(CXX) -shared -o $@ $^


#GenomicsDB examples

#Linking
$(GENOMICSDB_BIN_DIR)/%: $(GENOMICSDB_OBJ_DIR)/%.o $(GENOMICSDB_STATIC_LIBRARY) \
    $(TILEDB_DIR)/core/lib/$(TILEDB_BUILD)/libtiledb.a $(HTSDIR)/libhts.a
	@mkdir -p $(GENOMICSDB_BIN_DIR)
	@echo "Linking example $@"
	@$(CXX) $(LINKFLAGS) -o $@ $< $(LDFLAGS)
