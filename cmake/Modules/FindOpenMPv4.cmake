# Determine if OpenMP specification v4 is supported by the C/C++ compiler
# Once done this will define
# OPENMPV4_FOUND - OpenMP v4 found
# OpenMP_C_FLAGS
# OpenMP_CXX_FLAGS

include(CheckCSourceCompiles)
find_package(OpenMP QUIET)
set(OpenMPv4_C_TEST_SOURCE
"
#include <stdio.h>
#include <stdlib.h>

/*Define custom reduction operation*/
int main()
{
    int i = 0;
    int A[10];
    int sum = 0;
#pragma omp declare reduction ( sum_up : int : omp_out += omp_in ) initializer(omp_priv = 0)
#pragma omp parallel for default(shared) reduction(sum_up : sum)
    for(i=0;i<10;++i)
        sum += A[i];
    return 0;
}
")

if(NOT OPENMP_FOUND)
    set(OPENMPV4_FOUND False)
else()
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${OpenMP_C_FLAGS}")
    check_c_source_compiles("${OpenMPv4_C_TEST_SOURCE}" OPENMPV4_FOUND)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
endif()
