cmake_minimum_required(VERSION 3.0.0)
project(LADELtest VERSION 0.1.0)

SET(TEST_SRC
    src/run_all_tests.c
    src/test_scale.c
    src/test_matvec.c
    src/test_matmat.c
    src/test_upper_diag.c
    src/test_etree.c
    src/test_postorder.c
    src/test_transpose.c
    src/test_col_counts.c
    src/test_permutation.c
    src/test_ldl.c
    src/test_rank1_mod.c
    src/test_row_mod.c
    src/test_add.c
    src/test_submatrix.c)

add_executable(ladel_run_all_tests
    ${TEST_SRC}
)

target_include_directories(ladel_run_all_tests
    PRIVATE include
    include/Minunit
)

target_link_libraries(ladel_run_all_tests ladel)

if(UNIX)
target_link_libraries(ladel_run_all_tests m)
endif(UNIX)