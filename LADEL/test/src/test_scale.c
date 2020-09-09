#include "minunit.h"
#include "ladel_types.h"
#include "ladel_scale.h"
#include "ladel_global.h"

#define NROW 4
#define NCOL 5
#define NZMAX 11
#define TOL 1e-8

static ladel_sparse_matrix *M;
/*
M = [1.2  0  3.4  0    0.5;
     0   -2  1.1  0    0;
     3.6  0  1.5  1.6  0;
     0   -3  0    1.7 -0.5;]
*/

void scale_test_setup(void)
{
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UNSYMMETRIC, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 2; M->p[2] = 4; M->p[3] = 7; M->p[4] = 9; M->p[5] = 11;
    M->i[0] = 0; M->i[1] = 2; M->i[2] = 1; M->i[3] = 3; M->i[4] = 0; M->i[5] = 1; M->i[6] = 2; M->i[7] = 2; M->i[8] = 3; M->i[9] = 0; M->i[10] = 3;
    M->x[0] = 1.2; M->x[1] = 3.6; M->x[2] = -2; M->x[3] = -3; M->x[4] = 3.4; M->x[5] = 1.1; M->x[6] = 1.5; M->x[7] = 1.6; M->x[8] = 1.7; M->x[9] = 0.5; M->x[10] = -0.5;
}

void scale_test_teardown(void)
{
    ladel_sparse_free(M);
}

MU_TEST(test_scale_rows)
{
    ladel_double row_scale[NROW] = {2, 0.5, -1, -3};
    ladel_scale_rows(M, row_scale);
    mu_assert_double_eq(M->x[0], 2.4, TOL);
    mu_assert_double_eq(M->x[1], -3.6, TOL);
    mu_assert_double_eq(M->x[2], -1, TOL);
    mu_assert_double_eq(M->x[3], 9, TOL);
    mu_assert_double_eq(M->x[4], 6.8, TOL);
    mu_assert_double_eq(M->x[5], 0.55, TOL);
    mu_assert_double_eq(M->x[6], -1.5, TOL);
    mu_assert_double_eq(M->x[7], -1.6, TOL);
    mu_assert_double_eq(M->x[8], -5.1, TOL);
    mu_assert_double_eq(M->x[9], 1, TOL);
    mu_assert_double_eq(M->x[10], 1.5, TOL);
}


MU_TEST(test_scale_columns)
{
    ladel_double col_scale[NCOL] = {-2, 0.5, 1.5, 5, 10};
    ladel_scale_columns(M, col_scale);
    mu_assert_double_eq(M->x[0], -2.4, TOL);
    mu_assert_double_eq(M->x[1], -7.2, TOL);
    mu_assert_double_eq(M->x[2], -1, TOL);
    mu_assert_double_eq(M->x[3], -1.5, TOL);
    mu_assert_double_eq(M->x[4], 5.1, TOL);
    mu_assert_double_eq(M->x[5], 1.65, TOL);
    mu_assert_double_eq(M->x[6], 2.25, TOL);
    mu_assert_double_eq(M->x[7], 8, TOL);
    mu_assert_double_eq(M->x[8], 8.5, TOL);
    mu_assert_double_eq(M->x[9], 5, TOL);
    mu_assert_double_eq(M->x[10], -5, TOL);
}

MU_TEST(test_scale_scalar)
{
    ladel_double s = 100;
    ladel_scale_scalar(M, s);
    mu_assert_double_eq(M->x[0], 120, TOL);
    mu_assert_double_eq(M->x[1], 360, TOL);
    mu_assert_double_eq(M->x[2], -200, TOL);
    mu_assert_double_eq(M->x[3], -300, TOL);
    mu_assert_double_eq(M->x[4], 340, TOL);
    mu_assert_double_eq(M->x[5], 110, TOL);
    mu_assert_double_eq(M->x[6], 150, TOL);
    mu_assert_double_eq(M->x[7], 160, TOL);
    mu_assert_double_eq(M->x[8], 170, TOL);
    mu_assert_double_eq(M->x[9], 50, TOL);
    mu_assert_double_eq(M->x[10], -50, TOL);
}

MU_TEST(test_infinity_norm_rows)
{
    ladel_double norms[NROW];
    ladel_infinity_norm_rows(M, norms);
    mu_assert_double_eq(norms[0], 3.4, TOL);
    mu_assert_double_eq(norms[1], 2, TOL);
    mu_assert_double_eq(norms[2], 3.6, TOL);
    mu_assert_double_eq(norms[3], 3.0, TOL);
}

MU_TEST(test_infinity_norm_columns)
{
    ladel_double norms[NCOL];
    ladel_infinity_norm_columns(M, norms);
    mu_assert_double_eq(norms[0], 3.6, TOL);
    mu_assert_double_eq(norms[1], 3, TOL);
    mu_assert_double_eq(norms[2], 3.4, TOL);
    mu_assert_double_eq(norms[3], 1.7, TOL);
    mu_assert_double_eq(norms[4], 0.5, TOL);
}

MU_TEST_SUITE(suite_scale) 
{
    MU_SUITE_CONFIGURE(NULL, NULL, scale_test_setup, scale_test_teardown);

    MU_RUN_TEST(test_scale_rows);
    MU_RUN_TEST(test_scale_columns);
    MU_RUN_TEST(test_scale_scalar);
    MU_RUN_TEST(test_infinity_norm_rows);
    MU_RUN_TEST(test_infinity_norm_columns);
}