file(GLOB_RECURSE CPP_TEST_SOURCES ./cpp/*.cc)
build_GenomicsDB_test_executable(runAllGTests ${CPP_TEST_SOURCES})
