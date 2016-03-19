##########
# Macros #
##########

OS := $(shell uname)

# Large file support
LFS_CFLAGS = -D_FILE_OFFSET_BITS=64

CFLAGS=-Wall -Wno-reorder -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable
#LINKFLAGS appear before the object file list in the link command (e.g. -fopenmp, -O3)
LINKFLAGS=
#LDFLAGS appear after the list of object files (-lz etc)
LDFLAGS=-lz -lrt -lcrypto

ifdef OPENMP
  CFLAGS+=-fopenmp
endif
LINKFLAGS+=-fopenmp

# --- Debug/Release/Verbose mode handler --- #
BUILD ?= debug
VERBOSE ?= 0

ifeq ($(BUILD),debug)
  CFLAGS+= -g -gdwarf-2 -g3 -DDEBUG
  LINKFLAGS+=-g -gdwarf-2 -g3
endif

ifeq ($(BUILD),release)
  CFLAGS += -DNDEBUG -O3 -fvisibility=hidden
  LINKFLAGS+=-O3
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

ifndef TILEDB_DIR
    TILEDB_DIR=dependencies/TileDB
endif
CPPFLAGS+=-I$(TILEDB_DIR)/core/include/c_api

ifndef HTSDIR
    HTSDIR=dependencies/htslib
endif

ifdef HTSDIR
    CPPFLAGS+=-I$(HTSDIR) -DHTSDIR
    LDFLAGS+=-Wl,-Bstatic -L$(HTSDIR) -lhts -Wl,-Bdynamic
endif

#In the current version, this is mandatory
CPPFLAGS += -DDUPLICATE_CELL_AT_END

ifndef RAPIDJSON_INCLUDE_DIR
    RAPIDJSON_INCLUDE_DIR=dependencies/RapidJSON/include
endif
CPPFLAGS+=-I$(RAPIDJSON_INCLUDE_DIR)

ifdef USE_BIGMPI
    CPPFLAGS+=-I$(USE_BIGMPI)/src -DUSE_BIGMPI
    LDFLAGS+=-L$(USE_BIGMPI)/src -lbigmpi
endif

ifdef DO_PROFILING
    CPPFLAGS+=-DDO_PROFILING
endif

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

#Variant dir
VARIANT_INCLUDE_DIR = variant/include
VARIANT_SRC_DIR = variant/src
VARIANT_OBJ_DIR = variant/obj
VARIANT_BIN_DIR = variant/bin

#Variant examples
VARIANT_EXAMPLE_INCLUDE_DIR = variant/example/include
VARIANT_EXAMPLE_SRC_DIR = variant/example/src
VARIANT_EXAMPLE_OBJ_DIR = variant/example/obj
VARIANT_EXAMPLE_BIN_DIR = variant/example/bin

TILEDB_LDFLAGS= -Wl,-Bstatic -L$(TILEDB_DIR)/core/lib/$(BUILD) -ltiledb -Wl,-Bdynamic 
GENOMICSDB_LDFLAGS= -Wl,-Bstatic -L$(VARIANT_BIN_DIR)/ -ltiledb_variant -Wl,-Bdynamic

# --- Paths --- #

VARIANT_INCLUDE_PATHS = -I$(VARIANT_INCLUDE_DIR)
VARIANT_EXAMPLE_INCLUDE_PATHS = -I$(VARIANT_EXAMPLE_INCLUDE_DIR)

# --- Files --- #

VARIANT_INCLUDE := $(wildcard $(VARIANT_INCLUDE_DIR)/*.h)
VARIANT_SRC := $(wildcard $(VARIANT_SRC_DIR)/*.cc)
VARIANT_OBJ := $(patsubst $(VARIANT_SRC_DIR)/%.cc, $(VARIANT_OBJ_DIR)/%.o, $(VARIANT_SRC))
VARIANT_EXAMPLE_INCLUDE := $(wildcard $(VARIANT_EXAMPLE_INCLUDE_DIR)/*.h)
VARIANT_EXAMPLE_SRC := $(wildcard $(VARIANT_EXAMPLE_SRC_DIR)/*.cc)
VARIANT_EXAMPLE_OBJ := $(patsubst $(VARIANT_EXAMPLE_SRC_DIR)/%.cc, $(VARIANT_EXAMPLE_OBJ_DIR)/%.o, $(VARIANT_EXAMPLE_SRC))
VARIANT_EXAMPLE_BIN := $(patsubst $(VARIANT_EXAMPLE_SRC_DIR)/%.cc, $(VARIANT_EXAMPLE_BIN_DIR)/%, $(VARIANT_EXAMPLE_SRC))

###################
# General Targets #
###################

.PHONY: clean variant variant_example clean_variant clean_variant_example

all: variant variant_example

variant: $(VARIANT_OBJ) $(VARIANT_BIN_DIR)/libtiledb_variant.a

variant_example: $(VARIANT_EXAMPLE_BIN) $(VARIANT_EXAMPLE_OBJ)

clean: clean_variant clean_variant_example TileDB_clean

$(TILEDB_DIR)/core/lib/$(BUILD)/libtiledb.a:
	make -C $(TILEDB_DIR) MPIPATH=$(MPIPATH) BUILD=$(BUILD) -j 16

TileDB_clean:
	make -C $(TILEDB_DIR) clean

$(HTSDIR)/libhts.a:
	make -C $(HTSDIR) -j 16

###############
# Variant specific part of TileDB #
###############

# --- Compilation and dependency genration --- #

-include $(VARIANT_OBJ:.o=.d)

$(VARIANT_OBJ_DIR)/%.o: $(VARIANT_SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(VARIANT_INCLUDE_PATHS) -c $< -o $@
	@$(CXX) $(CPPFLAGS) -MM $(VARIANT_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

$(VARIANT_BIN_DIR)/libtiledb_variant.a: $(VARIANT_OBJ)
	@test -d $(VARIANT_BIN_DIR) || mkdir -p $(VARIANT_BIN_DIR)
	ar rcs $@ $^

clean_variant:
	rm -f $(VARIANT_OBJ_DIR)/* $(VARIANT_BIN_DIR)/*

####################
# Variant examples #
####################

# --- Compilation and dependency genration --- #

-include $(VARIANT_EXAMPLE_OBJ:.o=.d)

$(VARIANT_EXAMPLE_OBJ_DIR)/%.o: $(VARIANT_EXAMPLE_SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(VARIANT_INCLUDE_PATHS) $(VARIANT_EXAMPLE_INCLUDE_PATHS) -c $< -o $@
	@$(CXX) $(CPPFLAGS) -MM $(VARIANT_INCLUDE_PATHS) $(VARIANT_EXAMPLE_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

#Linking
$(VARIANT_EXAMPLE_BIN_DIR)/%: $(VARIANT_EXAMPLE_OBJ_DIR)/%.o $(VARIANT_BIN_DIR)/libtiledb_variant.a \
    $(TILEDB_DIR)/core/lib/$(BUILD)/libtiledb.a $(HTSDIR)/libhts.a
	@test -d $(VARIANT_EXAMPLE_BIN_DIR) || mkdir -p $(VARIANT_EXAMPLE_BIN_DIR)
	$(CXX) $(LINKFLAGS) -o $@ $< $(GENOMICSDB_LDFLAGS) $(TILEDB_LDFLAGS) $(LDFLAGS)

clean_variant_example:
	rm -f $(VARIANT_EXAMPLE_OBJ_DIR)/* $(VARIANT_EXAMPLE_BIN_DIR)/*
