#include "minunit.h"
#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_permutation.h"
#include "ladel_transpose.h"
#include "ladel_copy.h"

#define NCOL 6
#define NROW NCOL
#define NZMAX 10
#define TOL 1e-8

static ladel_work *work;
static ladel_double x[NCOL], y[NCOL];
static ladel_sparse_matrix *M, *Mpp, *Mpp_ref;

void permutation_suite_setup(void)
{
    work = ladel_workspace_allocate(NCOL);
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 1; M->p[2] = 2; M->p[3] = 4; M->p[4] = 5; M->p[5] = 7; M->p[6] = 10;
    M->i[0] = 0; M->i[1] = 1; M->i[2] = 0; M->i[3] = 2; M->i[4] = 3; M->i[5] = 1; M->i[6] = 4; 
    M->i[7] = 0; M->i[8] = 4; M->i[9] = 5; 
    M->x[0] = 1; M->x[1] = 2; M->x[2] = 10; M->x[3] = 3; M->x[4] = 4; M->x[5] = -5; M->x[6] = 5; 
    M->x[7] = 13; M->x[8] = 44; M->x[9] = 6; 

    Mpp = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE, FALSE);
    Mpp_ref = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE, FALSE);
    Mpp_ref->p[0] = 0; Mpp_ref->p[1] = 1; Mpp_ref->p[2] = 3; Mpp_ref->p[3] = 5; Mpp_ref->p[4] = 6; Mpp_ref->p[5] = 8; Mpp_ref->p[6] = 10;
    Mpp_ref->i[0] = 0; Mpp_ref->i[1] = 0; Mpp_ref->i[2] = 1; Mpp_ref->i[3] = 1; Mpp_ref->i[4] = 2; Mpp_ref->i[5] = 3; Mpp_ref->i[6] = 2; 
    Mpp_ref->i[7] = 4; Mpp_ref->i[8] = 4; Mpp_ref->i[9] = 5; 
    Mpp_ref->x[0] = 3; Mpp_ref->x[1] = 10; Mpp_ref->x[2] = 1; Mpp_ref->x[3] = 13; Mpp_ref->x[4] = 6; Mpp_ref->x[5] = 4; Mpp_ref->x[6] = 44; 
    Mpp_ref->x[7] = 5; Mpp_ref->x[8] = -5; Mpp_ref->x[9] = 2; 
}

void permutation_suite_teardown(void)
{
    ladel_workspace_free(work);
    ladel_sparse_free(M);
    ladel_sparse_free(Mpp);
    ladel_sparse_free(Mpp_ref);
}

void permutation_test_setup(void)
{
    x[0] = 0; x[1] = 1; x[2] = 2; x[3] = 3; x[4] = 4; x[5] = 5; 
}

MU_TEST(test_permute_vector)
{
    ladel_int permutation_vector[NCOL] = {2, 0, 5, 3, 4, 1};
    ladel_permute_vector(x, permutation_vector, NCOL, y);
    mu_assert_double_eq(y[0], x[2], TOL);
    mu_assert_double_eq(y[1], x[0], TOL);
    mu_assert_double_eq(y[2], x[5], TOL);
    mu_assert_double_eq(y[3], x[3], TOL);
    mu_assert_double_eq(y[4], x[4], TOL);
    mu_assert_double_eq(y[5], x[1], TOL);
}

MU_TEST(test_permute_inverse_vector)
{
    ladel_int inverse_permutation_vector[NCOL] = {2, 0, 5, 3, 4, 1};
    ladel_inverse_permute_vector(x, inverse_permutation_vector, NCOL, y);
    mu_assert_double_eq(y[2], x[0], TOL);
    mu_assert_double_eq(y[0], x[1], TOL);
    mu_assert_double_eq(y[5], x[2], TOL);
    mu_assert_double_eq(y[3], x[3], TOL);
    mu_assert_double_eq(y[4], x[4], TOL);
    mu_assert_double_eq(y[1], x[5], TOL);
}


MU_TEST(test_permute_symmetric_matrix)
{
    ladel_int permutation_vector[NCOL] = {2, 0, 5, 3, 4, 1};
    ladel_permute_symmetric_matrix(M, permutation_vector, Mpp, work);
    
    /* Sort entries after permuting by doing a double transpose */
    ladel_sparse_matrix *Mpp_transpose = ladel_transpose(Mpp, TRUE, work);
    ladel_sparse_free(Mpp);
    Mpp = ladel_transpose(Mpp_transpose, TRUE, work);
    ladel_sparse_free(Mpp_transpose);

    ladel_int index;
    for (index = 0; index < NCOL+1; index++)
    {
        mu_assert_long_eq(Mpp->p[index], Mpp_ref->p[index]);
    }
    for (index = 0; index < NZMAX; index++)
    {
        mu_assert_long_eq(Mpp->i[index], Mpp_ref->i[index]);
        mu_assert_double_eq(Mpp->x[index], Mpp_ref->x[index], TOL);
    }
}

MU_TEST(test_permute_symmetric_matrix_without_permutation_vector)
{
    ladel_permute_symmetric_matrix(M, NULL, Mpp, work);
    ladel_int index;
    for (index = 0; index < NCOL+1; index++)
    {
        mu_assert_long_eq(Mpp->p[index], M->p[index]);
    }
    for (index = 0; index < NZMAX; index++)
    {
        mu_assert_long_eq(Mpp->i[index], M->i[index]);
        mu_assert_double_eq(Mpp->x[index], M->x[index], TOL);
    }
}

MU_TEST_SUITE(suite_permutation)
{
    MU_SUITE_CONFIGURE(permutation_suite_setup, permutation_suite_teardown, permutation_test_setup, NULL);
    MU_RUN_TEST(test_permute_vector);
    MU_RUN_TEST(test_permute_inverse_vector);
    MU_RUN_TEST(test_permute_symmetric_matrix);
    MU_RUN_TEST(test_permute_symmetric_matrix_without_permutation_vector);
}