find_package(GTest)
if(GTEST_FOUND)
    enable_testing()
    set(CPP_TEST_SOURCES
        main_testall.cc
        test_mapping_data_loader.cc)
    message("IN TEST: ${CPP_TEST_SOURCES}") 
    add_executable(runAllGTests ${CPP_TEST_SOURCES})
    target_link_libraries(runAllGTests ${GTEST_BOTH_LIBRARIES})
    build_GenomicsDB_executable_common(runAllGTests)
endif()
