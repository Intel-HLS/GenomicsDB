##########
# Macros #
##########

OS := $(shell uname)

# Configuration flags
CONFIG_FLAGS =

# Use of mmap function for reading
USE_MMAP =

ifeq ($(USE_MMAP),)
  USE_MMAP = 1
endif

ifeq ($(USE_MMAP),1)
  CONFIG_FLAGS += -D_TILEDB_USE_MMAP
endif

# Large file support
LFS_CFLAGS = -D_FILE_OFFSET_BITS=64

#FIXME: should remove the no-sign-compare flag
CFLAGS=-Wall -Wno-reorder -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare
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

# Parallel sort
GNU_PARALLEL =

ifeq ($(GNU_PARALLEL),)
  GNU_PARALLEL = 1
endif

ifeq ($(GNU_PARALLEL),1)
  CFLAGS += -DGNU_PARALLEL
endif

# MPI compiler for C++
ifdef MPIPATH
    CC  = $(MPIPATH)/mpicc
    CXX = $(MPIPATH)/mpicxx
else
    CC  = mpicc
    CXX = mpicxx
endif
CPPFLAGS=-std=c++11 -fPIC \
      $(LFS_CFLAGS) $(CFLAGS) $(CONFIG_FLAGS)

ifdef HTSDIR
    CPPFLAGS+=-I$(HTSDIR) -DHTSDIR
    LDFLAGS+=-Wl,-Bstatic -L$(HTSDIR) -lhts -Wl,-Bdynamic
endif

CPPFLAGS += -DDUPLICATE_CELL_AT_END

ifdef RAPIDJSON_INCLUDE_DIR
    CPPFLAGS+=-I$(RAPIDJSON_INCLUDE_DIR)
endif

ifdef USE_BIGMPI
    CPPFLAGS+=-I$(USE_BIGMPI)/src -DUSE_BIGMPI
    LDFLAGS+=-L$(USE_BIGMPI)/src -lbigmpi
endif

ifeq ($(VERBOSE),0)
  CFLAGS += -DNVERBOSE
endif

ifneq ($(VERBOSE),0)
  CFLAGS += -DVERBOSE=$(VERBOSE)
endif

# --- Set library path to Google Test shared objects --- #
LDFLAGS += -L$(PWD)/3rdparty/gtest/lib
LDFLAGS += -Wl,-R$(PWD)/3rdparty/gtest/lib `$$ORIGIN`

ifdef DO_PROFILING
    CPPFLAGS+=-DDO_PROFILING
endif

ifdef USE_GPERFTOOLS
    GPERFTOOLSDIR ?= /home/karthikg/softwares/gperftools-2.2/install/
    CPPFLAGS+=-DUSE_GPERFTOOLS -I$(GPERFTOOLSDIR)/include
    LDFLAGS += -Wl,-Bstatic -L$(GPERFTOOLSDIR)/lib -lprofiler -Wl,-Bdynamic  -lunwind
endif

ifdef VERBOSE
    CPPFLAGS+= -DVERBOSE=$(VERBOSE)
endif

# --- Directories --- #
# Directories for the core code of TileDB
CORE_INCLUDE_DIR = core/include
CORE_INCLUDE_SUBDIRS = $(wildcard core/include/*)
CORE_SRC_DIR = core/src
CORE_SRC_SUBDIRS = $(wildcard core/src/*)
CORE_OBJ_DEB_DIR = core/obj/debug
CORE_BIN_DEB_DIR = core/bin/debug
ifeq ($(BUILD),debug)
  CORE_OBJ_DIR = $(CORE_OBJ_DEB_DIR)
  CORE_BIN_DIR = $(CORE_BIN_DEB_DIR)
endif
CORE_OBJ_REL_DIR = core/obj/release
CORE_BIN_REL_DIR = core/bin/release
ifeq ($(BUILD),release)
  CORE_OBJ_DIR = $(CORE_OBJ_REL_DIR)
  CORE_BIN_DIR = $(CORE_BIN_REL_DIR)
endif
CORE_LIB_DEB_DIR = core/lib/debug
ifeq ($(BUILD),debug)
  CORE_LIB_DIR = $(CORE_LIB_DEB_DIR)
endif
CORE_LIB_REL_DIR = core/lib/release
ifeq ($(BUILD),release)
  CORE_LIB_DIR = $(CORE_LIB_REL_DIR)
endif

#Variant dir
VARIANT_INCLUDE_DIR = variant/include
VARIANT_SRC_DIR = variant/src
VARIANT_OBJ_DIR = variant/obj
VARIANT_BIN_DIR = variant/bin
VARIANT_EXAMPLE_INCLUDE_DIR = variant/example/include
VARIANT_EXAMPLE_SRC_DIR = variant/example/src
VARIANT_EXAMPLE_OBJ_DIR = variant/example/obj
VARIANT_EXAMPLE_BIN_DIR = variant/example/bin

# Directories for the examples
EXAMPLES_INCLUDE_DIR = examples/include
EXAMPLES_SRC_DIR = examples/src
EXAMPLES_OBJ_DEB_DIR = examples/obj/debug
EXAMPLES_BIN_DEB_DIR = examples/bin/debug
ifeq ($(BUILD),debug)
  EXAMPLES_OBJ_DIR = $(EXAMPLES_OBJ_DEB_DIR)
  EXAMPLES_BIN_DIR = $(EXAMPLES_BIN_DEB_DIR)
endif
EXAMPLES_OBJ_REL_DIR = examples/obj/release
EXAMPLES_BIN_REL_DIR = examples/bin/release
ifeq ($(BUILD),release)
  EXAMPLES_OBJ_DIR = $(EXAMPLES_OBJ_REL_DIR)
  EXAMPLES_BIN_DIR = $(EXAMPLES_BIN_REL_DIR)
endif

# Directories for TileDB tests
TEST_SRC_SUBDIRS = $(wildcard test/src/*)
TEST_SRC_DIR = test/src
TEST_OBJ_DIR = test/obj
TEST_BIN_DIR = test/bin

# Directory for Doxygen documentation
DOXYGEN_DIR = doxygen
DOXYGEN_MAINPAGE = $(DOXYGEN_DIR)/mainpage.dox

# Directories for the MPI files - not necessary if mpicxx used.
MPI_INCLUDE_DIR := .
MPI_LIB_DIR := .

STATIC_LINK_CORE_LIBRARY=-Wl,-Bstatic -L$(CORE_BIN_DIR)/ -lcore -Wl,-Bdynamic
STATIC_LINK_VARIANT_LIBRARY=-Wl,-Bstatic -L$(VARIANT_BIN_DIR)/ -ltiledb_variant -Wl,-Bdynamic

# --- Paths --- #
CORE_INCLUDE_PATHS = $(addprefix -I, $(CORE_INCLUDE_SUBDIRS))
TEST_INCLUDE_PATHS = $(addprefix -I, $(CORE_INCLUDE_SUBDIRS))

EXAMPLES_INCLUDE_PATHS = -I$(EXAMPLES_INCLUDE_DIR)
VARIANT_INCLUDE_PATHS = -I$(VARIANT_INCLUDE_DIR)
VARIANT_EXAMPLE_INCLUDE_PATHS = -I$(VARIANT_EXAMPLE_INCLUDE_DIR)
MPI_INCLUDE_PATHS = -I$(MPI_INCLUDE_DIR)
MPI_LIB_PATHS = -L$(MPI_LIB_DIR)

# --- File Extensions --- #
ifeq ($(OS), Darwin)
  SHLIB_EXT = dylib
else
  SHLIB_EXT = so
endif

# --- Files --- #

# Files of the TileDB core
CORE_INCLUDE := $(foreach D,$(CORE_INCLUDE_SUBDIRS),$D/*.h)
CORE_SRC := $(wildcard $(foreach D,$(CORE_SRC_SUBDIRS),$D/*.cc))
CORE_OBJ := $(patsubst $(CORE_SRC_DIR)/%.cc, $(CORE_OBJ_DIR)/%.o, $(CORE_SRC))

VARIANT_INCLUDE := $(wildcard $(VARIANT_INCLUDE_DIR)/*.h)
VARIANT_SRC := $(wildcard $(VARIANT_SRC_DIR)/*.cc)
VARIANT_OBJ := $(patsubst $(VARIANT_SRC_DIR)/%.cc, $(VARIANT_OBJ_DIR)/%.o, $(VARIANT_SRC))
VARIANT_EXAMPLE_INCLUDE := $(wildcard $(VARIANT_EXAMPLE_INCLUDE_DIR)/*.h)
VARIANT_EXAMPLE_SRC := $(wildcard $(VARIANT_EXAMPLE_SRC_DIR)/*.cc)
VARIANT_EXAMPLE_OBJ := $(patsubst $(VARIANT_EXAMPLE_SRC_DIR)/%.cc, $(VARIANT_EXAMPLE_OBJ_DIR)/%.o, $(VARIANT_EXAMPLE_SRC))
VARIANT_EXAMPLE_BIN := $(patsubst $(VARIANT_EXAMPLE_SRC_DIR)/%.cc, $(VARIANT_EXAMPLE_BIN_DIR)/%, $(VARIANT_EXAMPLE_SRC))

# Files of the examples
EXAMPLES_INCLUDE := $(wildcard $(EXAMPLES_INCLUDE_DIR)/*.h)
EXAMPLES_SRC := $(wildcard $(EXAMPLES_SRC_DIR)/*.cc)
EXAMPLES_OBJ := $(patsubst $(EXAMPLES_SRC_DIR)/%.cc,\
                             $(EXAMPLES_OBJ_DIR)/%.o, $(EXAMPLES_SRC))
EXAMPLES_BIN := $(patsubst $(EXAMPLES_SRC_DIR)/%.cc,\
                             $(EXAMPLES_BIN_DIR)/%, $(EXAMPLES_SRC))

# Files of the TileDB tests
TEST_SRC := $(wildcard $(foreach D,$(TEST_SRC_SUBDIRS),$D/*.cc))
TEST_OBJ := $(patsubst $(TEST_SRC_DIR)/%.cc, $(TEST_OBJ_DIR)/%.o, $(TEST_SRC))

# Files for the HTML version of the Manpages
MANPAGES_MAN := $(wildcard $(MANPAGES_MAN_DIR)/*)
MANPAGES_HTML := $(patsubst $(MANPAGES_MAN_DIR)/%,\
                            $(MANPAGES_HTML_DIR)/%.html, $(MANPAGES_MAN))

###################
# General Targets #
###################

.PHONY: core examples check doc clean_core \
        clean_check clean_examples \
        clean variant example clean_variant clean_variant_example

all: core libtiledb examples variant variant_example

core: $(CORE_OBJ)

libtiledb: core $(CORE_LIB_DIR)/libtiledb.$(SHLIB_EXT) $(CORE_LIB_DIR)/libtiledb.a

examples: core $(EXAMPLES_OBJ) $(EXAMPLES_BIN)

html: $(MANPAGES_HTML)

doc: doxyfile.inc html

ifdef HTSDIR
include $(HTSDIR)/htslib.mk
endif

variant: $(VARIANT_OBJ) $(VARIANT_BIN_DIR)/libtiledb_variant.a

variant_example: $(VARIANT_EXAMPLE_BIN) $(VARIANT_EXAMPLE_OBJ)

check: libtiledb $(TEST_BIN_DIR)/tiledb_test
	@echo "Running TileDB tests"
	@$(TEST_BIN_DIR)/tiledb_test

clean: clean_core clean_libtiledb \
       clean_check clean_doc clean_examples clean_variant clean_variant_example

########
# Core #
########

# --- Compilation and dependency genration --- #

-include $(CORE_OBJ:%.o=%.d)

$(CORE_OBJ_DIR)/%.o: $(CORE_SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(CORE_INCLUDE_PATHS) \
                $(MPI_INCLUDE_PATHS) -c $< -o $@
	@$(CXX) $(CPPFLAGS) -MM $(CORE_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

# --- Cleaning --- #

clean_core:
	@echo 'Cleaning core'
	@rm -rf $(CORE_OBJ_DEB_DIR)/* $(CORE_OBJ_REL_DIR)/* \
                $(CORE_BIN_DEB_DIR)/* $(CORE_BIN_REL_DIR)/*

#############
# libtiledb #
#############

-include $(CORE_OBJ:%.o=%.d)

# --- Linking --- #

ifeq ($(0S), Darwin)
  SHLIB_FLAGS = -dynamiclib
else
  SHLIB_FLAGS = -shared
endif

ifeq ($(SHLIB_EXT), so)
  SONAME = -Wl,-soname=libtiledb.so
else
  SONAME =
endif

$(CORE_LIB_DIR)/libtiledb.$(SHLIB_EXT): $(CORE_OBJ)
	@mkdir -p $(CORE_LIB_DIR)
	@echo "Creating dynamic library libtiledb.$(SHLIB_EXT)"
	@$(CXX) $(SHLIB_FLAGS) $(SONAME) -o $@ $^ $(LDFLAGS)

$(CORE_LIB_DIR)/libtiledb.a: $(CORE_OBJ)
	@mkdir -p $(CORE_LIB_DIR)
	@echo "Creating static library libtiledb.a"
	@ar rcs $(CORE_LIB_DIR)/libtiledb.a $^

# --- Cleaning --- #

clean_libtiledb:
	@echo "Cleaning libtiledb.$(SHLIB_EXT)"
	@rm -rf $(CORE_LIB_DEB_DIR)/* $(CORE_LIB_REL_DIR)/*

##############
#  Examples  #
##############

# --- Compilation and dependency genration --- #

-include $(EXAMPLES_OBJ:.o=.d)

$(EXAMPLES_OBJ_DIR)/%.o: $(EXAMPLES_SRC_DIR)/%.cc
	@mkdir -p $(EXAMPLES_OBJ_DIR)
	@echo "Compiling $<"
	@$(CXX) $(EXAMPLES_INCLUDE_PATHS) $(CORE_INCLUDE_PATHS) -c $< \
         $(ZLIB) $(OPENSSLLIB) -o $@
	@$(CXX) -MM $(EXAMPLES_INCLUDE_PATHS) \
                    $(CORE_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

# --- Linking --- #

$(EXAMPLES_BIN_DIR)/%: $(EXAMPLES_OBJ_DIR)/%.o $(CORE_OBJ)
	@mkdir -p $(EXAMPLES_BIN_DIR)
	@echo "Creating $@"
	@$(CXX) $(LINKFLAGS) -o $@ $^ $(LDFLAGS) 

# --- Cleaning --- #

clean_examples:
	@echo 'Cleaning examples'
	@rm -f $(EXAMPLES_OBJ_DEB_DIR)/* $(EXAMPLES_OBJ_REL_DIR)/* \
               $(EXAMPLES_BIN_DEB_DIR)/* $(EXAMPLES_BIN_REL_DIR)/*

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

################
# TileDB Tests #
################

# --- Compilation and dependency genration --- #

-include $(TEST_OBJ:.o=.d)

$(TEST_OBJ_DIR)/%.o: $(TEST_SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	@echo "Compiling $<"
	@$(CXX) $(TEST_INCLUDE_PATHS) -c $< -o $@
	@$(CXX) -MM $(TEST_INCLUDE_PATHS) \
                    $(CORE_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

# --- Linking --- #

$(TEST_BIN_DIR)/tiledb_test: $(TEST_OBJ) $(CORE_OBJ)
	@mkdir -p $(TEST_BIN_DIR)
	@echo "Creating test_cmd"
	@$(CXX) $(LDFLAGS) $(OPENMP_LIB_PATHS) $(OPENMP_LIB) \
			$(MPI_LIB_PATHS) $(MPI_LIB) \
      -o $@ $^ $(ZLIB) $(OPENSSLLIB) -lgtest -lgtest_main

# --- Cleaning --- #

clean_check:
	@echo "Cleaning check"
	@rm -rf $(TEST_OBJ_DIR) $(TEST_BIN_DIR)

################################
# Documentation (with Doxygen) #
################################

doxyfile.inc: $(CORE_INCLUDE)
	@echo 'Creating Doxygen documentation'
	@echo INPUT = $(DOXYGEN_DIR)/mainpage.dox $(CORE_INCLUDE) \
	    > doxyfile.inc
	@echo FILE_PATTERNS = *.h >> doxyfile.inc
	@doxygen Doxyfile.mk > Doxyfile.log 2>&1

# --- Cleaning --- #

clean_doc:
	@echo "Cleaning documentation"
	@rm -f doxyfile.inc

