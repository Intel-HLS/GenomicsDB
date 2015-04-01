##########
# Macros #
##########

# --- Compiler --- #
CXX = g++
CPPFLAGS=-std=c++11
LDFLAGS=

#HTSDIR=../../htslib

ifdef HTSDIR
  CPPFLAGS+=-I$(HTSDIR) -DHTSDIR
  LDFLAGS+=-Wl,-Bstatic -L$(HTSDIR) -lhts -Wl,-Bdynamic
endif

CPPFLAGS += -fPIC
SOFLAGS=-shared -Wl,-soname,

LINKFLAGS=
ifdef DEBUG
  CPPFLAGS+= -g -gdwarf-2 -g3
  LINKFLAGS+=-g -gdwarf-2 -g3
endif
ifdef SPEED
  CPPFLAGS+=-O3 -DNDEBUG
  LINKFLAGS+=-O3
endif

ifdef DO_PROFILING
  GPERFTOOLSDIR=/home/karthikg/softwares/gperftools-2.2/install/
  CPPFLAGS+=-DDO_PROFILING --no-inline -I$(GPERFTOOLSDIR)/include
  LDFLAGS += -Wl,-Bstatic -L$(GPERFTOOLSDIR)/lib -lprofiler -Wl,-Bdynamic  -lunwind 
endif


# --- Directories --- #
CORE_INCLUDE_DIR = core/include
CORE_SRC_DIR = core/src
CORE_OBJ_DIR = core/obj
CORE_BIN_DIR = core/bin
VARIANT_INCLUDE_DIR = variant/include
VARIANT_SRC_DIR = variant/src
VARIANT_OBJ_DIR = variant/obj
VARIANT_BIN_DIR = variant/bin
EXAMPLE_INCLUDE_DIR = example/include
EXAMPLE_SRC_DIR = example/src
EXAMPLE_OBJ_DIR = example/obj
EXAMPLE_BIN_DIR = example/bin
GTEST_DIR = gtest
GTEST_INCLUDE_DIR = gtest/include
GTEST_SRC_DIR = gtest/src
GTEST_OBJ_DIR = gtest/obj
GTEST_BIN_DIR = gtest/bin
TEST_SRC_DIR = test/src
TEST_OBJ_DIR = test/obj
TEST_BIN_DIR = test/bin
DOC_DIR = doc

STATIC_LINK_CORE_LIBRARY=-Wl,-Bstatic -L$(CORE_BIN_DIR)/ -lcore -Wl,-Bdynamic
STATIC_LINK_VARIANT_LIBRARY=-Wl,-Bstatic -L$(VARIANT_BIN_DIR)/ -ltiledb_variant -Wl,-Bdynamic

# --- Paths --- #
INCLUDE_PATHS = -I$(CORE_INCLUDE_DIR) -I$(VARIANT_INCLUDE_DIR)

# --- Files --- #
CORE_INCLUDE := $(wildcard $(CORE_INCLUDE_DIR)/*.h)
CORE_SRC := $(wildcard $(CORE_SRC_DIR)/*.cc)
CORE_OBJ := $(patsubst $(CORE_SRC_DIR)/%.cc, $(CORE_OBJ_DIR)/%.o, $(CORE_SRC))
VARIANT_INCLUDE := $(wildcard $(VARIANT_INCLUDE_DIR)/*.h)
VARIANT_SRC := $(wildcard $(VARIANT_SRC_DIR)/*.cc)
VARIANT_OBJ := $(patsubst $(VARIANT_SRC_DIR)/%.cc, $(VARIANT_OBJ_DIR)/%.o, $(VARIANT_SRC))
EXAMPLE_INCLUDE := $(wildcard $(EXAMPLE_INCLUDE_DIR)/*.h)
EXAMPLE_SRC := $(wildcard $(EXAMPLE_SRC_DIR)/*.cc)
EXAMPLE_OBJ := $(patsubst $(EXAMPLE_SRC_DIR)/%.cc, $(EXAMPLE_OBJ_DIR)/%.o, $(EXAMPLE_SRC))
EXAMPLE_BIN := $(patsubst $(EXAMPLE_SRC_DIR)/%.cc, $(EXAMPLE_BIN_DIR)/%, $(EXAMPLE_SRC))
GTEST_INCLUDE := $(wildcard $(GTEST_INCLUDE_DIR)/*.h)
GTEST_OBJ := $(patsubst $(GTEST_SRC_DIR)/%.cc, $(GTEST_OBJ_DIR)/%.o, $(GTEST_SRC))
TEST_SRC := $(wildcard $(TEST_SRC_DIR)/*.cc)
TEST_OBJ := $(patsubst $(TEST_SRC_DIR)/%.cc, $(TEST_OBJ_DIR)/%.o, $(TEST_SRC))

###################
# General Targets #
###################

.PHONY: core example gtest test doc doc_doxygen clean_core clean_example clean_gtest clean_test clean

all: core variant example gtest test 

ifdef HTSDIR
include $(HTSDIR)/htslib.mk
endif

core: $(CORE_OBJ) $(CORE_BIN_DIR)/libcore.a

variant: $(VARIANT_OBJ) $(VARIANT_BIN_DIR)/libtiledb_variant.so $(VARIANT_BIN_DIR)/libtiledb_variant.a

example: $(EXAMPLE_BIN)

gtest: $(GTEST_OBJ_DIR)/gtest-all.o

test: $(TEST_OBJ)

doc: doxyfile.inc

clean: clean_core clean_variant clean_example clean_gtest clean_test

clean_lib:
	rm -f $(CORE_BIN_DIR)/* $(VARIANT_BIN_DIR)/*

###############
# Core TileDB #
###############

# --- Compilation and dependency genration --- #

-include $(CORE_OBJ:.o=.d)

$(CORE_OBJ_DIR)/%.o: $(CORE_SRC_DIR)/%.cc
	@test -d $(CORE_OBJ_DIR) || mkdir -p $(CORE_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(INCLUDE_PATHS) -c $< -o $@
	@$(CXX) -MM $(CPPFLAGS) $(INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

$(CORE_BIN_DIR)/libcore.a: $(CORE_OBJ)
	@test -d $(CORE_BIN_DIR) || mkdir -p $(CORE_BIN_DIR)
	ar rcs $@ $^

clean_core:
	rm -f $(CORE_OBJ_DIR)/* $(CORE_BIN_DIR)/* 

###############
# Variant specific part of TileDB #
###############

# --- Compilation and dependency genration --- #

-include $(VARIANT_OBJ:.o=.d)

$(VARIANT_OBJ_DIR)/%.o: $(VARIANT_SRC_DIR)/%.cc
	@test -d $(VARIANT_OBJ_DIR) || mkdir -p $(VARIANT_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(INCLUDE_PATHS) -c $< -o $@
	@$(CXX) -MM $(CPPFLAGS) $(INCLUDE_PATHS) $< > $(@:.o=.d)
	@mv -f $(@:.o=.d) $(@:.o=.d.tmp)
	@sed 's|.*:|$@:|' < $(@:.o=.d.tmp) > $(@:.o=.d)
	@rm -f $(@:.o=.d.tmp)

$(VARIANT_BIN_DIR)/libtiledb_variant.a: $(CORE_OBJ) $(VARIANT_OBJ)
	@test -d $(VARIANT_BIN_DIR) || mkdir -p $(VARIANT_BIN_DIR)
	ar rcs $@ $^

$(VARIANT_BIN_DIR)/libtiledb_variant.so: $(CORE_OBJ) $(VARIANT_OBJ)
	@test -d $(VARIANT_BIN_DIR) || mkdir -p $(VARIANT_BIN_DIR)
	$(CXX) $(SOFLAGS)libtiledb_variant.so -o $@ $^

clean_variant:
	rm -f $(VARIANT_OBJ_DIR)/* $(VARIANT_BIN_DIR)/* 

############
# Examples #
############

# --- Compilation and dependency genration --- #

-include $(EXAMPLE_OBJ:.o=.d)

$(EXAMPLE_OBJ_DIR)/%.o: $(EXAMPLE_SRC_DIR)/%.cc
	@test -d $(EXAMPLE_OBJ_DIR) || mkdir -p $(EXAMPLE_OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(INCLUDE_PATHS) -I $(EXAMPLE_INCLUDE_DIR) -c $< -o $@
	@$(CXX) -MM $(CPPFLAGS) $(INCLUDE_PATHS) -I $(EXAMPLE_INCLUDE_DIR) $< > $(@:.o=.d)
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
# Google test #
###############

$(GTEST_OBJ_DIR)/gtest-all.o: gtest/src/gtest-all.cc $(wildcard gtest/include/gtest/*.h)
	@test -d $(GTEST_OBJ_DIR) || mkdir -p $(GTEST_OBJ_DIR)
	$(CXX) -isystem $(GTEST_INCLUDE_DIR) -I$(GTEST_DIR) -pthread -c $< -o $@

clean_gtest:
	rm -f $(GTEST_OBJ_DIR)/* $(GTEST_BIN_DIR)/* 

#########
# Tests #
#########

# Coming up soon...

#########################
# Documentation doxygen #
#########################

doxyfile.inc: $(CORE_INCLUDE)
	@echo INPUT         =  $(DOC_DIR)/mainpage.dox $(CORE_INCLUDE) > doxyfile.inc
	@echo FILE_PATTERNS =  *.h >> doxyfile.inc
	doxygen Doxyfile.mk

# LIB_PATHS = /usr/local/lib/libspatialindex.so
# LIBS = -lpqxx -lpthread
