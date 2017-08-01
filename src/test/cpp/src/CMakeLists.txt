find_package(GTest)
if(GTEST_FOUND)
    set(CPP_TEST_SOURCES
        main_testall.cc
        )
    if(LIBDBI_FOUND)
        set(CPP_TEST_SOURCES
            ${CPP_TEST_SOURCES}
            test_mapping_data_loader.cc)
    endif()
    add_executable(runAllGTests ${CPP_TEST_SOURCES})
    target_link_libraries(runAllGTests ${GTEST_BOTH_LIBRARIES})
    build_GenomicsDB_executable_common(runAllGTests)
    add_test(NAME All_GTests COMMAND runAllGTests)
endif()
