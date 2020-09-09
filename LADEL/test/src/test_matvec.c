#include "minunit.h"
#include "ladel_global.h"
#include "ladel_types.h"
#include "ladel_matvec.h"

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

void matvec_suite_setup(void)
{
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UNSYMMETRIC, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 2; M->p[2] = 4; M->p[3] = 7; M->p[4] = 9; M->p[5] = 11;
    M->i[0] = 0; M->i[1] = 2; M->i[2] = 1; M->i[3] = 3; M->i[4] = 0; M->i[5] = 1; M->i[6] = 2; M->i[7] = 2; M->i[8] = 3; M->i[9] = 0; M->i[10] = 3;
    M->x[0] = 1.2; M->x[1] = 3.6; M->x[2] = -2; M->x[3] = -3; M->x[4] = 3.4; M->x[5] = 1.1; M->x[6] = 1.5; M->x[7] = 1.6; M->x[8] = 1.7; M->x[9] = 0.5; M->x[10] = -0.5;
}

void matvec_suite_teardown(void)
{
    ladel_sparse_free(M);
}

MU_TEST(test_matvec)
{
    ladel_double x[NCOL] = {1.5, -2, -1, 5, -10};
    ladel_double y[NROW];
    ladel_matvec(M, x, y, TRUE);
    mu_assert_double_eq(y[0], -6.6, TOL);
    mu_assert_double_eq(y[1], 2.9, TOL);
    mu_assert_double_eq(y[2], 11.9, TOL);
    mu_assert_double_eq(y[3], 19.5, TOL);

    ladel_double temp[NROW];
    ladel_int index;
    for (index = 0; index < NROW; index++) temp[index] = y[index];
    ladel_matvec(M, x, y, FALSE);
    mu_assert_double_eq(y[0], 2*temp[0], TOL);
    mu_assert_double_eq(y[1], 2*temp[1], TOL);
    mu_assert_double_eq(y[2], 2*temp[2], TOL);
    mu_assert_double_eq(y[3], 2*temp[3], TOL);
}

MU_TEST(test_tpose_matvec)
{
    ladel_double x[NROW] = {1, -1.5, -2, 100};
    ladel_double y[NCOL];
    ladel_tpose_matvec(M, x, y, TRUE);
    mu_assert_double_eq(y[0], -6, TOL);
    mu_assert_double_eq(y[1], -297, TOL);
    mu_assert_double_eq(y[2], -1.25, TOL);
    mu_assert_double_eq(y[3], 166.8, TOL);
    mu_assert_double_eq(y[4], -49.5, TOL);

    ladel_double temp[NCOL];
    ladel_int index;
    for (index = 0; index < NCOL; index++) temp[index] = y[index];
    ladel_tpose_matvec(M, x, y, FALSE);
    mu_assert_double_eq(y[0], 2*temp[0], TOL);
    mu_assert_double_eq(y[1], 2*temp[1], TOL);
    mu_assert_double_eq(y[2], 2*temp[2], TOL);
    mu_assert_double_eq(y[3], 2*temp[3], TOL);
    mu_assert_double_eq(y[4], 2*temp[4], TOL);
}

MU_TEST_SUITE(suite_matvec)
{
    MU_SUITE_CONFIGURE(matvec_suite_setup, matvec_suite_teardown, NULL, NULL);
    MU_RUN_TEST(test_matvec);
    MU_RUN_TEST(test_tpose_matvec);
}