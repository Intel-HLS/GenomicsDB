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
    CC  = mpicc
    CXX = mpicxx
endif
CPPFLAGS=-std=c++11 -fPIC $(LFS_CFLAGS) $(CFLAGS)

ifdef HTSDIR
    CPPFLAGS+=-I$(HTSDIR) -DHTSDIR
    LDFLAGS+=-Wl,-Bstatic -L$(HTSDIR) -lhts -Wl,-Bdynamic
    include $(HTSDIR)/htslib.mk
endif

#In the current version, this is mandatory
CPPFLAGS += -DDUPLICATE_CELL_AT_END

ifdef RAPIDJSON_INCLUDE_DIR
    CPPFLAGS+=-I$(RAPIDJSON_INCLUDE_DIR)
endif

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

STATIC_LINK_TILEDB_LIBRARY=-Wl,-Bstatic -L$(CORE_BIN_DIR)/ -lcore -Wl,-Bdynamic
STATIC_LINK_VARIANT_LIBRARY=-Wl,-Bstatic -L$(VARIANT_BIN_DIR)/ -ltiledb_variant -Wl,-Bdynamic

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

clean: clean_variant clean_variant_example

###############
# Variant specific part of TileDB #
###############

# --- Compilation and dependency genration --- #

-include $(VARIANT_OBJ:.o=.d)

$(VARIANT_OBJ_DIR)/%.o: $(VARIANT_SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(CORE_INCLUDE_PATHS) \
                $(MPI_INCLUDE_PATHS) $(VARIANT_INCLUDE_PATHS) -c $< -o $@
	@$(CXX) $(CPPFLAGS) -MM $(CORE_INCLUDE_PATHS) $(VARIANT_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

$(VARIANT_BIN_DIR)/libtiledb_variant.a: $(CORE_OBJ) $(VARIANT_OBJ)
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
	@$(CXX) $(CPPFLAGS) $(CORE_INCLUDE_PATHS) \
                $(MPI_INCLUDE_PATHS) $(VARIANT_INCLUDE_PATHS) $(VARIANT_EXAMPLE_INCLUDE_PATHS) -c $< -o $@
	@$(CXX) $(CPPFLAGS) -MM $(CORE_INCLUDE_PATHS) $(VARIANT_INCLUDE_PATHS) $(VARIANT_EXAMPLE_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

#Linking
$(VARIANT_EXAMPLE_BIN_DIR)/%: $(VARIANT_EXAMPLE_OBJ_DIR)/%.o $(VARIANT_BIN_DIR)/libtiledb_variant.a
	@test -d $(VARIANT_EXAMPLE_BIN_DIR) || mkdir -p $(VARIANT_EXAMPLE_BIN_DIR)
	$(CXX) $(LINKFLAGS) -o $@ $< $(STATIC_LINK_VARIANT_LIBRARY) $(LDFLAGS)

clean_variant_example:
	rm -f $(VARIANT_EXAMPLE_OBJ_DIR)/* $(VARIANT_EXAMPLE_BIN_DIR)/*
