#include "minunit.h"
#include "ladel.h"

#define NROW 3
#define NCOL 3
#define TOL 1e-8

static ladel_work *work;
static ladel_sparse_matrix *A;
static ladel_sparse_matrix *B;
/*
A = [1 2 0;
     2 3 4;
     0 4 5];
B = [10 30 20;
     30 10  0;
     20  0 10];
*/

void add_suite_setup(void)
{
    work = ladel_workspace_allocate(NCOL);
}

void add_test_setup(void)
{
    A = ladel_sparse_alloc(NROW, NCOL, 7, UNSYMMETRIC, TRUE, FALSE);
    B = ladel_sparse_alloc(NROW, NCOL, 7, UNSYMMETRIC, TRUE, FALSE);
    A->p[0] = 0; A->p[1] = 2; A->p[2] = 5; A->p[3] = 7;  
    A->i[0] = 0; A->i[1] = 1; A->i[2] = 0; A->i[3] = 1; A->i[4] = 2; A->i[5] = 1; A->i[6] = 2; 
    A->x[0] = 1; A->x[1] = 2; A->x[2] = 2; A->x[3] = 3; A->x[4] = 4; A->x[5] = 4; A->x[6] = 5;   
    B->p[0] = 0; B->p[1] = 3; B->p[2] = 5; B->p[3] = 7;  
    B->i[0] = 0; B->i[1] = 1; B->i[2] = 2; B->i[3] = 0; B->i[4] = 1; B->i[5] = 0; B->i[6] = 2; 
    B->x[0] = 10; B->x[1] = 30; B->x[2] = 20; B->x[3] = 30; B->x[4] = 10; B->x[5] = 20; B->x[6] = 10;     
}

void add_suite_teardown(void)
{
    ladel_workspace_free(work);
}

void add_test_teardown(void)
{
    A = ladel_sparse_free(A);
    B = ladel_sparse_free(B);
}

MU_TEST(test_add_unsymmetric_matrices)
{
    ladel_double alpha = 4.0, beta = 0.2;
    ladel_sparse_matrix *C = ladel_add_matrices(alpha, A, beta, B, work);
    /* The resulting matrix is unsorted, so we check correctness through some common operations */
    
    /* Test the resulting mat_vec */
    ladel_double x[NROW] = {1.5, -2, 1};
    ladel_double y[NROW], y_sol[NROW] = {-15, 9, -4};
    ladel_matvec(C, x, y, TRUE);
    
    ladel_int index;
    for (index = 0; index < NROW; index++)
        mu_assert_double_eq(y[index], y_sol[index], TOL);

    ladel_sparse_free(C);
}

MU_TEST(test_add_symmetric_matrices)
{
    ladel_double alpha = 4.0, beta = 0.2;
    ladel_to_upper_diag(A);
    ladel_to_upper_diag(B);
    ladel_sparse_matrix *C = ladel_add_matrices(alpha, A, beta, B, work);
    /* The resulting matrix is unsorted, so we check correctness through some common operations */
    
    /* Test the resulting mat_vec */
    ladel_double x[NROW] = {1.5, -2, 1};
    ladel_double y[NROW], y_sol[NROW] = {-15, 9, -4};
    ladel_symmetric_matvec(C, x, y, TRUE);
    
    ladel_int index;
    for (index = 0; index < NROW; index++)
        mu_assert_double_eq(y[index], y_sol[index], TOL);

    /* Test the resulting solve */
    y_sol[0] = -3.018092105263158e-01;
    y_sol[1] = 2.623355263157895e-01;
    y_sol[2] = -9.046052631578948e-02;

    ladel_symbolics *sym = ladel_symbolics_alloc(NCOL);
    ladel_factor *LD;
    ladel_factorize(C, sym, NO_ORDERING, &LD, work);
    ladel_dense_solve(LD, x, y, work);
    
    for (index = 0; index < NCOL; index++)
        mu_assert_double_eq(y[index], y_sol[index], TOL);

    ladel_symbolics_free(sym);
    ladel_factor_free(LD);

    ladel_sparse_free(C);
}

MU_TEST_SUITE(suite_add)
{
    MU_SUITE_CONFIGURE(add_suite_setup, add_suite_teardown, add_test_setup, add_test_teardown);
    MU_RUN_TEST(test_add_unsymmetric_matrices);
    MU_RUN_TEST(test_add_symmetric_matrices);
}