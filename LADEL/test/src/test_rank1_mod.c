#include "minunit.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel_types.h"
#include "ladel.h"
#include "ladel_rank1_mod.h"
#include "ladel_debug_print.h"

#define NCOL 8
#define NROW 8
#define NZMAX 24
#define MAX_SIZE_SET 8
#define TOL 1e-8

static ladel_work *work;
static ladel_sparse_matrix *M, *Mbasis;
static ladel_sparse_matrix *W;
static ladel_factor *LD = NULL;
static ladel_symbolics *sym = NULL;

void rank1_mod_suite_setup(void)
{
    work = ladel_workspace_allocate(NCOL);
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 1; M->p[2] = 2; M->p[3] = 4; M->p[4] = 7; M->p[5] = 8; M->p[6] = 10; M->p[7] = 14; M->p[8] = 16;
    M->i[0] = 0; M->i[1] = 1; M->i[2] = 0; M->i[3] = 2; M->i[4] = 1; M->i[5] = 2; M->i[6] = 3; M->i[7] = 4;  
    M->i[8] = 4; M->i[9] = 5; M->i[10] = 2; M->i[11] = 3; M->i[12] = 5; M->i[13] = 6; M->i[14] = 5; M->i[15] = 7; 
    M->x[0] = 1.484233718441293e+00; M->x[1] = 8.503811060900474e-01; M->x[2] = 7.946847833340259e-02; 
    M->x[3] = 4.093923468942169e-01; M->x[4] = 6.944674301359678e-02; M->x[5] = 1.761109237578858e-01; 
    M->x[6] = 3.766406539600695e-01; M->x[7] = 9.168193399034030e-01; M->x[8] = 9.238873678854942e-01; 
    M->x[9] = 1.897902116880709e+00; M->x[10] = 6.052733699027846e-01; M->x[11] = 2.665692902440686e-01; 
    M->x[12] = 7.650155176644616e-02; M->x[13] = 1.151758087181832e+00; M->x[14] = 7.767464464874704e-01; 
    M->x[15] = 6.404488288848778e-01; 

    W = ladel_sparse_alloc(NROW, 1, 3, UNSYMMETRIC, TRUE, TRUE);
    W->p[0] = 0; W->p[1] = 3; 
    W->nz[0] = 3;
    W->i[0] = 3; W->i[1] = 5; W->i[2] = 7;
    W->x[0] = 1.418863386272153e-01; W->x[1] = 4.217612826262750e-01; W->x[2] = 9.157355251890671e-01;  


    sym = ladel_symbolics_alloc(NCOL); 

    Mbasis = ladel_sparse_alloc(NROW, NCOL, 36, UPPER, FALSE, FALSE);
    Mbasis->p[0] = 0; Mbasis->p[1] = 1; Mbasis->p[2] = 3; Mbasis->p[3] = 6; Mbasis->p[4] = 10; 
    Mbasis->p[5] = 15; Mbasis->p[6] = 21; Mbasis->p[7] = 28; Mbasis->p[8] = 36; 
    Mbasis->i[0] = 0; Mbasis->i[1] = 0; Mbasis->i[2] = 1; Mbasis->i[3] = 0; Mbasis->i[4] = 1; 
    Mbasis->i[5] = 2; Mbasis->i[6] = 0; Mbasis->i[7] = 1; Mbasis->i[8] = 2; Mbasis->i[9] = 3; 
    Mbasis->i[10] = 0; Mbasis->i[11] = 1; Mbasis->i[12] = 2; Mbasis->i[13] = 3; Mbasis->i[14] = 4; 
    Mbasis->i[15] = 0; Mbasis->i[16] = 1; Mbasis->i[17] = 2; Mbasis->i[18] = 3; Mbasis->i[19] = 4; 
    Mbasis->i[20] = 5; Mbasis->i[21] = 0; Mbasis->i[22] = 1; Mbasis->i[23] = 2; Mbasis->i[24] = 3; 
    Mbasis->i[25] = 4; Mbasis->i[26] = 5; Mbasis->i[27] = 6; Mbasis->i[28] = 0; Mbasis->i[29] = 1; 
    Mbasis->i[30] = 2; Mbasis->i[31] = 3; Mbasis->i[32] = 4; Mbasis->i[33] = 5; Mbasis->i[34] = 6; 
    Mbasis->i[35] = 7;  
}

void rank1_mod_suite_teardown(void)
{
    ladel_workspace_free(work);
    ladel_sparse_free(M);
    ladel_sparse_free(W);
    ladel_sparse_free(Mbasis);
    ladel_symbolics_free(sym);
}

void rank1_mod_test_setup(void)
{

}

void rank1_mod_test_teardown(void)
{
    LD = ladel_factor_free(LD);
}

MU_TEST(test_set_union)
{
    ladel_int set1_vals[MAX_SIZE_SET] = {1, 2, 5, 7};
    ladel_set *set1 = ladel_malloc(1, sizeof(ladel_set));
    ladel_set_set(set1, set1_vals, 4, MAX_SIZE_SET);

    ladel_int set2_vals[MAX_SIZE_SET] = {1, 2, 7};
    ladel_set *set2 = ladel_malloc(1, sizeof(ladel_set));
    ladel_set_set(set2, set2_vals, 3, MAX_SIZE_SET);

    ladel_set *dif = work->set_preallocated1;
    dif->size_set = 0;

    ladel_int offset[8], insertions[8];

    ladel_int status; 
    status = ladel_set_union(set1, set2, dif, offset, insertions, 0);
    mu_assert_long_eq(status, SET_HAS_NOT_CHANGED);
    mu_assert_long_eq(set1->size_set, 4);
    mu_assert_long_eq(dif->size_set, 0);
    
    ladel_int set3_vals[MAX_SIZE_SET] = {1, 4, 6, 10, 11};
    ladel_set *set3 = ladel_malloc(1, sizeof(ladel_set));
    ladel_set_set(set3, set3_vals, 5, MAX_SIZE_SET);

    status = ladel_set_union(set1, set3, dif, offset, insertions, 0);
    mu_assert_long_eq(status, SET_HAS_CHANGED);
    mu_assert_long_eq(set1->size_set, 8);
    
    mu_assert_long_eq(set1_vals[0], 1);
    mu_assert_long_eq(set1_vals[1], 2);
    mu_assert_long_eq(set1_vals[2], 4);
    mu_assert_long_eq(set1_vals[3], 5);
    mu_assert_long_eq(set1_vals[4], 6);
    mu_assert_long_eq(set1_vals[5], 7);
    mu_assert_long_eq(set1_vals[6], 10);
    mu_assert_long_eq(set1_vals[7], 11);

    mu_assert_long_eq(dif->size_set, 4);
    mu_assert_long_eq(dif->set[0], 4);
    mu_assert_long_eq(dif->set[1], 6);
    mu_assert_long_eq(dif->set[2], 10);
    mu_assert_long_eq(dif->set[3], 11);

    ladel_int set4_vals[MAX_SIZE_SET] = {14};
    ladel_set *set4 = ladel_malloc(1, sizeof(ladel_set));
    ladel_set_set(set4, set4_vals, 1, MAX_SIZE_SET);

    status = ladel_set_union(set1, set4, dif, offset, insertions, 0);
    mu_assert_long_eq(status, MAX_SET_SIZE_EXCEEDED);   

    ladel_free(set1);
    ladel_free(set2);
    ladel_free(set3);
    ladel_free(set4);
}

MU_TEST(test_rank1_mod)
{
    ladel_int status, index;
    status = ladel_factorize_advanced(M, sym, NO_ORDERING, &LD, Mbasis, work);
    mu_assert_long_eq(status, SUCCESS);

    ladel_double rhs[8] = {7.922073295595544e-01, 9.594924263929030e-01, 6.557406991565868e-01,
                            3.571167857418955e-02, 8.491293058687771e-01, 9.339932477575505e-01,
                            6.787351548577735e-01, 7.577401305783334e-01};
    ladel_double sol[8] = {1.626739458964396e+01, 1.404806681196643e+00, -2.938574982347923e+02,
                            -3.385740249326764e+00, 6.516216869838241e+02, -6.457174936667400e+02,
                            1.986908367068955e+02, 7.843195055030803e+02};

    ladel_double x[8];
    ladel_dense_solve(LD, rhs, x, work);
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol[index], TOL);

    ladel_double sol_mod[8] = {4.021591066933866e-01, 1.257689300408343e+00, 2.457694262215502e+00,
                                -1.584275766315217e+00, 3.426556805737886e+00, -2.481259429010251e+00, 
                                -1.707843620820627e-01, 2.602541252030759e+00};
    
    status = ladel_rank1_update(LD, sym, W, 0, 1.0, UPDATE, work);
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, rhs, x, work);
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol_mod[index], TOL);

    status = ladel_rank1_update(LD, sym, W, 0, 1.0, DOWNDATE, work);
    mu_assert_long_eq(status, SUCCESS);
    ladel_dense_solve(LD, rhs, x, work);
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol[index], TOL);
}

MU_TEST_SUITE(suite_rank1_mod)
{
    MU_SUITE_CONFIGURE(rank1_mod_suite_setup, rank1_mod_suite_teardown, rank1_mod_test_setup, rank1_mod_test_teardown);
    MU_RUN_TEST(test_set_union);
    MU_RUN_TEST(test_rank1_mod);
}