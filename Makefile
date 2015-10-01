##########
# Macros #
##########

OS := $(shell uname)

# Large file support
LFS_CFLAGS = -D_FILE_OFFSET_BITS=64

CFLAGS = -fopenmp
#LINKFLAGS appear before the object file list in the link command (e.g. -fopenmp, -O3)
LINKFLAGS=-fopenmp
#LDFLAGS appear after the list of object files (-lz etc)
LDFLAGS=-lz

# --- Debug/Release mode handler --- #
BUILD ?= debug 

ifeq ($(BUILD),debug)
# CFLAGS += -DDEBUG -Wall -O0 -g
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

# --- Compilers --- #

# C++ compiler
# CXX = g++ 

# MPI compiler for C++
#MPIPATH = #/opt/mpich/dev/intel/default/bin/
CC  = $(MPIPATH)mpicc
CXX = $(MPIPATH)mpicxx
CPPFLAGS=-std=c++11 -fPIC -fvisibility=hidden \
      $(LFS_CFLAGS) $(CFLAGS)

#HTSDIR=../../htslib

ifdef HTSDIR
  CPPFLAGS+=-I$(HTSDIR) -DHTSDIR
  LDFLAGS+=-Wl,-Bstatic -L$(HTSDIR) -lhts -Wl,-Bdynamic
endif

#RAPIDJSON_INCLUDE_DIR ?= /home/karthikg/softwares/setup_files/rapidjson/include/
ifdef RAPIDJSON_INCLUDE_DIR
    CPPFLAGS+=-I$(RAPIDJSON_INCLUDE_DIR)
else
   $(error Variable RAPIDJSON_INCLUDE_DIR is not defined) 
endif

SOFLAGS=-shared -Wl,-soname=

ifdef DO_PROFILING
    CPPFLAGS+=-DDO_PROFILING
endif

ifdef USE_GPERFTOOLS
    GPERFTOOLSDIR ?= /home/karthikg/softwares/gperftools-2.2/install/
    CPPFLAGS+=-DUSE_GPERFTOOLS --no-inline -I$(GPERFTOOLSDIR)/include
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

# Directories for the command-line-based frontend of TileDB
TILEDB_CMD_INCLUDE_DIR = tiledb_cmd/include
TILEDB_CMD_SRC_DIR = tiledb_cmd/src
TILEDB_CMD_OBJ_DEB_DIR = tiledb_cmd/obj/debug
TILEDB_CMD_BIN_DEB_DIR = tiledb_cmd/bin/debug
ifeq ($(BUILD),debug)
  TILEDB_CMD_OBJ_DIR = $(TILEDB_CMD_OBJ_DEB_DIR)
  TILEDB_CMD_BIN_DIR = $(TILEDB_CMD_BIN_DEB_DIR)
endif
TILEDB_CMD_OBJ_REL_DIR = tiledb_cmd/obj/release
TILEDB_CMD_BIN_REL_DIR = tiledb_cmd/bin/release
ifeq ($(BUILD),release)
  TILEDB_CMD_OBJ_DIR = $(TILEDB_CMD_OBJ_REL_DIR)
  TILEDB_CMD_BIN_DIR = $(TILEDB_CMD_BIN_REL_DIR)
endif

#Variant dir
VARIANT_INCLUDE_DIR = variant/include
VARIANT_SRC_DIR = variant/src
VARIANT_OBJ_DIR = variant/obj
VARIANT_BIN_DIR = variant/bin
EXAMPLE_INCLUDE_DIR = variant/example/include
EXAMPLE_SRC_DIR = variant/example/src
EXAMPLE_OBJ_DIR = variant/example/obj
EXAMPLE_BIN_DIR = variant/example/bin

# Directories of Google Test
GTEST_DIR = gtest
GTEST_INCLUDE_DIR = gtest/include
GTEST_SRC_DIR = gtest/src
GTEST_OBJ_DIR = gtest/obj
GTEST_BIN_DIR = gtest/bin

# Directories for TileDB tests
TEST_SRC_DIR = test/src
TEST_OBJ_DIR = test/obj
TEST_BIN_DIR = test/bin

# Directories for Linear Algebra applications
LA_INCLUDE_DIR = la/include
LA_SRC_DIR = la/src
LA_OBJ_DIR = la/obj
LA_BIN_DIR = la/bin

# Directories for distributed applications
RVMA_INCLUDE_DIR = rvma/include
RVMA_SRC_DIR     = rvma/src
RVMA_OBJ_DIR     = rvma/obj
RVMA_BIN_DIR     = rvma/bin

# Directory for Doxygen documentation
DOXYGEN_DIR = doxygen

# Manpages directories
MANPAGES_MAN_DIR = manpages/man
MANPAGES_HTML_DIR = manpages/html

# Directories for the MPI files - not necessary if mpicxx used.
MPI_INCLUDE_DIR := .
MPI_LIB_DIR := .

# Directories for the OpenMP files
OPENMP_INCLUDE_DIR = .
OPENMP_LIB_DIR = .

STATIC_LINK_CORE_LIBRARY=-Wl,-Bstatic -L$(CORE_BIN_DIR)/ -lcore -Wl,-Bdynamic
STATIC_LINK_VARIANT_LIBRARY=-Wl,-Bstatic -L$(VARIANT_BIN_DIR)/ -ltiledb_variant -Wl,-Bdynamic

# --- Paths --- #
CORE_INCLUDE_PATHS = $(addprefix -I, $(CORE_INCLUDE_SUBDIRS))
TILEDB_CMD_INCLUDE_PATHS = -I$(TILEDB_CMD_INCLUDE_DIR)
LA_INCLUDE_PATHS = -I$(LA_INCLUDE_DIR)
VARIANT_INCLUDE_PATHS = -I$(VARIANT_INCLUDE_DIR)
EXAMPLE_INCLUDE_PATHS = -I$(EXAMPLE_INCLUDE_DIR)
MPI_INCLUDE_PATHS = -I$(MPI_INCLUDE_DIR)
MPI_LIB_PATHS = -L$(MPI_LIB_DIR)
OPENMP_INCLUDE_PATHS = -L$(OPENMP_INCLUDE_DIR)
OPENMP_LIB_PATHS = -L$(OPENMP_LIB_DIR)

# --- Libs --- #
MPI_LIB = -lmpi
OPENMP_LIB = -fopenmp 

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

# Files of the TileDB command-line-based frontend
TILEDB_CMD_INCLUDE := $(wildcard $(TILEDB_CMD_INCLUDE_DIR)/*.h)
TILEDB_CMD_SRC := $(wildcard $(TILEDB_CMD_SRC_DIR)/*.cc)
TILEDB_CMD_OBJ := $(patsubst $(TILEDB_CMD_SRC_DIR)/%.cc,\
                             $(TILEDB_CMD_OBJ_DIR)/%.o, $(TILEDB_CMD_SRC))
TILEDB_CMD_BIN := $(patsubst $(TILEDB_CMD_SRC_DIR)/%.cc,\
                             $(TILEDB_CMD_BIN_DIR)/%, $(TILEDB_CMD_SRC)) 

VARIANT_INCLUDE := $(wildcard $(VARIANT_INCLUDE_DIR)/*.h)
VARIANT_SRC := $(wildcard $(VARIANT_SRC_DIR)/*.cc)
VARIANT_OBJ := $(patsubst $(VARIANT_SRC_DIR)/%.cc, $(VARIANT_OBJ_DIR)/%.o, $(VARIANT_SRC))
EXAMPLE_INCLUDE := $(wildcard $(EXAMPLE_INCLUDE_DIR)/*.h)
EXAMPLE_SRC := $(wildcard $(EXAMPLE_SRC_DIR)/*.cc)
EXAMPLE_OBJ := $(patsubst $(EXAMPLE_SRC_DIR)/%.cc, $(EXAMPLE_OBJ_DIR)/%.o, $(EXAMPLE_SRC))
EXAMPLE_BIN := $(patsubst $(EXAMPLE_SRC_DIR)/%.cc, $(EXAMPLE_BIN_DIR)/%, $(EXAMPLE_SRC))

# Files of the Google Test
GTEST_INCLUDE := $(wildcard $(GTEST_INCLUDE_DIR)/*.h)
GTEST_OBJ := $(patsubst $(GTEST_SRC_DIR)/%.cc, $(GTEST_OBJ_DIR)/%.o,\
                        $(GTEST_SRC))

# Files of the TileDB tests
TEST_SRC := $(wildcard $(TEST_SRC_DIR)/*.cc)
TEST_OBJ := $(patsubst $(TEST_SRC_DIR)/%.cc, $(TEST_OBJ_DIR)/%.o, $(TEST_SRC))

# Files of the Linear Algebra applications
LA_SRC := $(wildcard $(LA_SRC_DIR)/*.cc)
LA_OBJ := $(patsubst $(LA_SRC_DIR)/%.cc, $(LA_OBJ_DIR)/%.o, $(LA_SRC))
LA_BIN := $(patsubst $(LA_SRC_DIR)/%.cc, $(LA_BIN_DIR)/%, $(LA_SRC))

# Files of the distributed applications
RVMA_SRC := $(wildcard $(RVMA_SRC_DIR)/*.c)
RVMA_OBJ := $(patsubst $(RVMA_SRC_DIR)/%.c, $(RVMA_OBJ_DIR)/%.o, $(RVMA_SRC))
RVMA_BIN := $(patsubst $(RVMA_SRC_DIR)/%.c, $(RVMA_BIN_DIR)/%, $(RVMA_SRC))

# Files for the HTML version of the Manpages
MANPAGES_MAN := $(wildcard $(MANPAGES_MAN_DIR)/*)
MANPAGES_HTML := $(patsubst $(MANPAGES_MAN_DIR)/%,\
                            $(MANPAGES_HTML_DIR)/%.html, $(MANPAGES_MAN)) 

###################
# General Targets #
###################

.PHONY: core gtest test doc clean_core clean_gtest variant example \
        clean_test clean_tiledb_cmd clean_la clean_variant clean_example clean

all: core libtiledb tiledb_cmd gtest test variant example 

core: $(CORE_OBJ) 

libtiledb: core $(CORE_LIB_DIR)/libtiledb.$(SHLIB_EXT)

tiledb_cmd: core $(TILEDB_CMD_OBJ) $(TILEDB_CMD_BIN)

la: core $(LA_OBJ) $(LA_BIN_DIR)/example_transpose

rvma: core $(RVMA_OBJ) #$(RVMA_BIN_DIR)/simple_test

html: $(MANPAGES_HTML)

doc: doxyfile.inc html

ifdef HTSDIR
include $(HTSDIR)/htslib.mk
endif

variant: $(VARIANT_OBJ) $(VARIANT_BIN_DIR)/libtiledb_variant.a

example: $(EXAMPLE_BIN) $(EXAMPLE_OBJ)

gtest: $(GTEST_OBJ_DIR)/gtest-all.o

test: $(TEST_OBJ)

clean: clean_core clean_libtiledb clean_tiledb_cmd clean_gtest \
       clean_test clean_doc clean_variant clean_example 

########
# Core #
########

# --- Compilation and dependency genration --- #

-include $(CORE_OBJ:%.o=%.d)

$(CORE_OBJ_DIR)/%.o: $(CORE_SRC_DIR)/%.cc
	@mkdir -p $(dir $@) 
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(CORE_INCLUDE_PATHS) $(OPENMP_INCLUDE_PATHS) \
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
	@echo "Creating libtiledb.$(SHLIB_EXT)"
	@$(CXX) $(SHLIB_FLAGS) $(SONAME) -o $@ $^ $(LDFLAGS)
	ar rcs $(CORE_LIB_DIR)/libtiledb.a $^

# --- Cleaning --- #

clean_libtiledb:
	@echo "Cleaning libtiledb.$(SHLIB_EXT)"
	@rm -rf $(CORE_LIB_DEB_DIR)/* $(CORE_LIB_REL_DIR)/*

##############
# TileDB_cmd #
##############

# --- Compilation and dependency genration --- #

-include $(TILEDB_CMD_OBJ:.o=.d)

$(TILEDB_CMD_OBJ_DIR)/%.o: $(TILEDB_CMD_SRC_DIR)/%.cc
	@mkdir -p $(TILEDB_CMD_OBJ_DIR)
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(TILEDB_CMD_INCLUDE_PATHS) $(CORE_INCLUDE_PATHS) -c $< -o $@
	@$(CXX) $(CPPFLAGS) -MM $(TILEDB_CMD_INCLUDE_PATHS) \
                    $(CORE_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

# --- Linking --- #

$(TILEDB_CMD_BIN_DIR)/%: $(TILEDB_CMD_OBJ_DIR)/%.o $(CORE_OBJ)
	@mkdir -p $(TILEDB_CMD_BIN_DIR)
	@echo "Creating $@"
	@$(CXX) $(OPENMP_LIB_PATHS) $(OPENMP_LIB) $(MPI_LIB_PATHS) $(MPI_LIB) \
                -o $@ $^ $(LDFLAGS)

# --- Cleaning --- #

clean_tiledb_cmd:
	@echo 'Cleaning tiledb_cmd'
	@rm -f $(TILEDB_CMD_OBJ_DEB_DIR)/* $(TILEDB_CMD_OBJ_REL_DIR)/* \
               $(TILEDB_CMD_BIN_DEB_DIR)/* $(TILEDB_CMD_BIN_REL_DIR)/*
 
######
# LA #
######

# --- Compilation and dependency genration --- #

# -include $(LA_OBJ:.o=.d)

# $(LA_OBJ_DIR)/%.o: $(LA_SRC_DIR)/%.cc
# 	@test -d $(LA_OBJ_DIR) || mkdir -p $(LA_OBJ_DIR)
# 	@echo "Compiling $<"
# 	@$(CXX) $(CPPFLAGS) $(LA_INCLUDE_PATHS) $(CORE_INCLUDE_PATHS) -c $< -o $@
# 	@$(CXX) $(CPPFLAGS) -MM $(CORE_INCLUDE_PATHS) $(LA_INCLUDE_PATHS) $< > $(@:.o=.d)
# 	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
# 	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
# 	@rm -f $(@:.o=.d.tmp)

# --- Linking --- #

# $(LA_BIN_DIR)/example_transpose: $(LA_OBJ) $(CORE_OBJ)
#	@mkdir -p $(LA_BIN_DIR)
#	@echo "Creating example_transpose"
#	@$(CXX) $(OPENMP_LIB_PATHS) $(OPENMP_LIB) $(MPI_LIB_PATHS) $(MPI_LIB) \
#               -o $@ $^

# --- Cleaning --- #

# clean_la:
#	@echo 'Cleaning la'
#	@rm -f $(LA_OBJ_DIR)/* $(LA_BIN_DIR)/* 

###############
# Variant specific part of TileDB #
###############

# --- Compilation and dependency genration --- #

-include $(VARIANT_OBJ:.o=.d)

$(VARIANT_OBJ_DIR)/%.o: $(VARIANT_SRC_DIR)/%.cc
	@mkdir -p $(dir $@) 
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(CORE_INCLUDE_PATHS) $(OPENMP_INCLUDE_PATHS) \
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

########
# RVMA #
########

# --- Compilation and dependency genration --- #

# -include $(RVMA_OBJ:.o=.d)

# $(RVMA_OBJ_DIR)/%.o: $(RVMA_SRC_DIR)/%.c
# 	@test -d $(RVMA_OBJ_DIR) || mkdir -p $(RVMA_OBJ_DIR)
# 	@echo "Compiling $<"
# 	@$(CXX) $(CPPFLAGS) $(RVMA_INCLUDE_PATHS) $(CORE_INCLUDE_PATHS) -c $< -o $@
# 	@$(CXX) $(CPPFLAGS) -MM $(CORE_INCLUDE_PATHS) $(RVMA_INCLUDE_PATHS) $< > $(@:.o=.d)
# 	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
# 	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
# 	@rm -f $(@:.o=.d.tmp)

# $(RVMA_BIN_DIR)/simple_test: $(RVMA_OBJ) $(CORE_OBJ)
# 	@mkdir -p $(RVMA_BIN_DIR)
# 	@echo "Creating simple_test"
# 	@$(CXX) $(OPENMP_LIB_PATHS) $(OPENMP_LIB) $(MPI_LIB_PATHS) $(MPI_LIB) \
#                -o $@ $^

# --- Cleaning --- #

# clean_rvma:
#	@echo 'Cleaning RVMA'
#	@rm -f $(RVMA_OBJ_DIR)/* $(RVMA_BIN_DIR)/*

############
# Examples #
############

# --- Compilation and dependency genration --- #

-include $(EXAMPLE_OBJ:.o=.d)

$(EXAMPLE_OBJ_DIR)/%.o: $(EXAMPLE_SRC_DIR)/%.cc
	@mkdir -p $(dir $@) 
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(CORE_INCLUDE_PATHS) $(OPENMP_INCLUDE_PATHS) \
                $(MPI_INCLUDE_PATHS) $(VARIANT_INCLUDE_PATHS) $(EXAMPLE_INCLUDE_PATHS) -c $< -o $@
	@$(CXX) $(CPPFLAGS) -MM $(CORE_INCLUDE_PATHS) $(VARIANT_INCLUDE_PATHS) $(EXAMPLE_INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

#Linking
$(EXAMPLE_BIN_DIR)/%: $(EXAMPLE_OBJ_DIR)/%.o $(VARIANT_BIN_DIR)/libtiledb_variant.a
	@test -d $(EXAMPLE_BIN_DIR) || mkdir -p $(EXAMPLE_BIN_DIR)
	$(CXX) $(LINKFLAGS) -o $@ $< $(STATIC_LINK_VARIANT_LIBRARY) $(LDFLAGS)

clean_example:
	rm -f $(EXAMPLE_OBJ_DIR)/* $(EXAMPLE_BIN_DIR)/* 

###############
# Google Test #
###############

# --- Compilation --- #

$(GTEST_OBJ_DIR)/gtest-all.o: gtest/src/gtest-all.cc \
                              $(wildcard gtest/include/gtest/*.h)
	@test -d $(GTEST_OBJ_DIR) || mkdir -p $(GTEST_OBJ_DIR)
	@echo "Compiling $<"
	@$(CXX) -isystem $(GTEST_INCLUDE_DIR) -I$(GTEST_DIR) \
                -pthread -c $< -o $@

# --- Cleaning --- #

clean_gtest:
	@echo "Cleaning gtest"
	@rm -f $(GTEST_OBJ_DIR)/* $(GTEST_BIN_DIR)/* 

################
# TileDB Tests #
################

# Coming up soon...

# --- Cleaning --- #

clean_test:
	@echo "Cleaning test"
	@rm -f $(TEST_OBJ_DIR)/* $(TEST_BIN_DIR)/* 

################################
# Documentation (with Doxygen) #
################################

doxyfile.inc: $(CORE_INCLUDE) $(TILEDB_CMD_INCLUDE) $(LA_INCLUDE)
	@echo 'Creating Doxygen documentation'
	@echo INPUT = $(DOXYGEN_DIR)/mainpage.dox $(CORE_INCLUDE) \
                      $(TILEDB_CMD_INCLUDE) $(LA_INCLUDE) > doxyfile.inc
	@echo FILE_PATTERNS = *.h >> doxyfile.inc
	@doxygen Doxyfile.mk > Doxyfile.log 2>&1

$(MANPAGES_HTML_DIR)/%.html: $(MANPAGES_MAN_DIR)/%
	@echo 'Converting $< to HTML'
	@mkdir -p $(MANPAGES_HTML_DIR)
	@man2html $< > $@

# --- Cleaning --- #

clean_doc:
	@echo "Cleaning documentation"
	@rm -f doxyfile.inc
	@rm -f $(MANPAGES_HTML)

