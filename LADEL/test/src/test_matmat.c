#include "minunit.h"
#include "ladel.h"

#define NROW 4
#define NCOL 5
#define NZMAX 9
#define TOL 1e-8

static ladel_work *work;
static ladel_sparse_matrix *M;
/*
M = [1.2  0  0    0    0.5;
     0   -2  1.1  0    0;
     3.6  0  1.5  0  0;
     0   -3  0    1.7 -0.5;]
*/

void matmat_suite_setup(void)
{
    work = ladel_workspace_allocate(NROW);
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UNSYMMETRIC, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 2; M->p[2] = 4; M->p[3] = 6; M->p[4] = 7; M->p[5] = 9;
    M->i[0] = 0; M->i[1] = 2; M->i[2] = 1; M->i[3] = 3; M->i[4] = 1; M->i[5] = 2; M->i[6] = 3; M->i[7] = 0; M->i[8] = 3;
    M->x[0] = 1.2; M->x[1] = 3.6; M->x[2] = -2; M->x[3] = -3; M->x[4] = 1.1; M->x[5] = 1.5; M->x[6] = 1.7; M->x[7] = 0.5; M->x[8] = -0.5;
}

void matmat_suite_teardown(void)
{
    ladel_sparse_free(M);
    ladel_workspace_free(work);
}

MU_TEST(test_mat_diag_mat_transpose)
{
    ladel_sparse_matrix *M_transpose = ladel_transpose(M, TRUE, work);
    ladel_double diag[NCOL] = {1, 2, 3, 4, 5};
    ladel_sparse_matrix *MMt = ladel_mat_diag_mat_transpose(M, M_transpose, diag, work);
    mu_assert_true(MMt != NULL);

    /* The resulting matrix is unsorted, so we check correctness through some common operations */
    
    /* Test the resulting mat_vec */
    ladel_double x[NROW] = {1.5, -2, -1, 5};
    ladel_double y[NROW], y_sol[NROW] = {-6.535, 31.79, -23.13, 128.175};
    ladel_symmetric_matvec(MMt, x, y, TRUE);
    
    ladel_int index;
    for (index = 0; index < NROW; index++)
        mu_assert_double_eq(y[index], y_sol[index], TOL);

    /* Test the resulting solve */
    y_sol[0] = 9.418703474431208e-01;
    y_sol[1] = -5.483886607270303e-01;
    y_sol[2] = -1.194498239652705e-01;
    y_sol[3] =  4.140863960736211e-01;

    ladel_symbolics *sym = ladel_symbolics_alloc(NROW);
    ladel_factor *LD;
    ladel_factorize(MMt, sym, NO_ORDERING, &LD, work);
    ladel_dense_solve(LD, x, y, work);
    
    for (index = 0; index < NROW; index++)
        mu_assert_double_eq(y[index], y_sol[index], TOL);

    ladel_symbolics_free(sym);
    ladel_factor_free(LD);

    ladel_sparse_free(M_transpose);
    ladel_sparse_free(MMt);
}

MU_TEST_SUITE(suite_matmat)
{
    MU_SUITE_CONFIGURE(matmat_suite_setup, matmat_suite_teardown, NULL, NULL);
    MU_RUN_TEST(test_mat_diag_mat_transpose);
}