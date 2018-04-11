if(NOT DISABLE_MPI)
    build_GenomicsDB_executable(test_genomicsdb_bcf_generator)
    build_GenomicsDB_executable(test_genomicsdb_importer)
endif()
