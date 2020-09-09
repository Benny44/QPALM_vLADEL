#include "minunit.h"
#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_transpose.h"

#define NROW 4
#define NCOL 5
#define NZMAX 7
#define TOL 1e-8

static ladel_work *work;
static ladel_sparse_matrix *M;

void transpose_suite_setup(void)
{
    work = ladel_workspace_allocate(NROW);
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UNSYMMETRIC, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 2; M->p[2] = 3; M->p[3] = 5; M->p[4] = 6; M->p[5] = 7; 
    M->i[0] = 0; M->i[1] = 2; M->i[2] = 1; M->i[3] = 2; M->i[4] = 3; M->i[5] = 0; M->i[6] = 3; 
    M->x[0] = 1; M->x[1] = 3; M->x[2] = 2; M->x[3] = 5; M->x[4] = 6; M->x[5] = 4; M->x[6] = 7; 
}

void transpose_suite_teardown(void)
{
    ladel_workspace_free(work);
    ladel_sparse_free(M);
}

MU_TEST(test_transpose_with_values)
{
    ladel_sparse_matrix *M_transpose = ladel_transpose(M, TRUE, work);
    ladel_int p_sol[NROW+1] = {0, 2, 3, 5, 7};
    ladel_int i_sol[NZMAX] = {0, 3, 1, 0, 2, 2, 4};
    ladel_double x_sol[NZMAX] = {1, 4, 2, 3, 5, 6, 7};

    ladel_int index;
    for (index = 0; index < NROW+1; index++)
    {
        mu_assert_long_eq(M_transpose->p[index], p_sol[index]);
    }
    for (index = 0; index < NZMAX; index++)
    {
        mu_assert_long_eq(M_transpose->i[index], i_sol[index]);
        mu_assert_double_eq(M_transpose->x[index], x_sol[index], TOL);
    }

    ladel_sparse_free(M_transpose);
}

MU_TEST(test_transpose_no_values)
{
    ladel_sparse_matrix *M_transpose = ladel_transpose(M, FALSE, work);
    ladel_int p_sol[NROW+1] = {0, 2, 3, 5, 7};
    ladel_int i_sol[NZMAX] = {0, 3, 1, 0, 2, 2, 4};

    ladel_int index;
    for (index = 0; index < NROW+1; index++)
    {
        mu_assert_long_eq(M_transpose->p[index], p_sol[index]);
    }
    for (index = 0; index < NZMAX; index++)
    {
        mu_assert_long_eq(M_transpose->i[index], i_sol[index]);
    }

    mu_assert_long_eq(M_transpose->values, FALSE);
    mu_assert_false(M_transpose->x);

    ladel_sparse_free(M_transpose);
}

MU_TEST_SUITE(suite_transpose)
{
    MU_SUITE_CONFIGURE(transpose_suite_setup, transpose_suite_teardown, NULL, NULL);
    MU_RUN_TEST(test_transpose_with_values);
    MU_RUN_TEST(test_transpose_no_values);
}