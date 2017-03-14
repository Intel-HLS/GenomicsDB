include_directories("include")

if(NOT DISABLE_MPI)
    build_GenomicsDB_executable(create_tiledb_workspace)
    build_GenomicsDB_executable(gt_mpi_gather)
    build_GenomicsDB_executable(vcf2tiledb)
    build_GenomicsDB_executable(vcfdiff)
    build_GenomicsDB_executable(vcf_histogram)
endif()
