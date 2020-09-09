#include "minunit.h"
#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_etree.h"
#include "ladel_postorder.h"

static ladel_work *work;
static ladel_sparse_matrix *M;
static ladel_symbolics *sym;
#define NROW 11
#define NCOL 11
#define NZMAX 43

void postorder_test_setup(void) 
{
    work = ladel_workspace_allocate(NCOL);
    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 3; M->p[2] = 6; M->p[3] = 10; M->p[4] = 13; M->p[5] = 16; 
    M->p[6] = 21; M->p[7] = 24; M->p[8] = 29; M->p[9] = 31; M->p[10] = 37; M->p[11] = 43;
    
    M->i[0] = 0; M->i[1] = 5; M->i[2] = 6; M->i[3] = 1; M->i[4] = 2; M->i[5] = 7; M->i[6] = 1; 
    M->i[7] = 2; M->i[8] = 9; M->i[9] = 10; M->i[10] = 3; M->i[11] = 5; M->i[12] = 9; M->i[13] = 4; 
    M->i[14] = 7; M->i[15] = 10; M->i[16] = 0; M->i[17] = 3; M->i[18] = 5; M->i[19] = 8; M->i[20] = 9;
    M->i[21] = 0; M->i[22] = 6; M->i[23] = 10; M->i[24] = 1; M->i[25] = 4; M->i[26] = 7; M->i[27] = 9; 
    M->i[28] = 10; M->i[29] = 5; M->i[30] = 8; M->i[31] = 2; M->i[32] = 3; M->i[33] = 5; M->i[34] = 7; 
    M->i[35] = 9; M->i[36] = 10; M->i[37] = 2; M->i[38] = 4; M->i[39] = 6; M->i[40] = 7; M->i[41] = 9; M->i[42] = 10;

    M->x[0] = 12; M->x[1] = -13; M->x[2] = 0; M->x[3] = -2; M->x[4] = 6; M->x[5] = 9; M->x[6] = 10; 
    M->x[7] = -9; M->x[8] = 7; M->x[9] = 6; M->x[10] = -14; M->x[11] = -16; M->x[12] = 0; M->x[13] = 19; 
    M->x[14] = -7; M->x[15] = 3; M->x[16] = -11; M->x[17] = 10; M->x[18] = -10; M->x[19] = 0; M->x[20] = 8;
    M->x[21] = 16; M->x[22] = 19; M->x[23] = 2; M->x[24] = -15; M->x[25] = -14; M->x[26] = -10; M->x[27] = 14; 
    M->x[28] = -10; M->x[29] = 13; M->x[30] = -11; M->x[31] = 18; M->x[32] = -6; M->x[33] = -12; M->x[34] = -10; 
    M->x[35] = 5; M->x[36] = -1; M->x[37] = -6; M->x[38] = 14; M->x[39] = 3; M->x[40] = 2; M->x[41] = 17; M->x[42] = -9;

    sym = ladel_symbolics_alloc(NCOL);
}

void postorder_test_teardown(void)
{
    ladel_workspace_free(work);
    ladel_sparse_free(M);
    ladel_symbolics_free(sym);
}

MU_TEST(test_postorder)
{
    ladel_int postorder_ref[NCOL] = {1, 2, 4, 7, 0, 3, 5, 6, 8, 9, 10};
    ladel_etree(M, sym, work);
    ladel_postorder(M, sym, work);
    
    ladel_int i;
    for (i = 0; i < NCOL; i++)
    {
        mu_assert_long_eq(sym->postorder[i], postorder_ref[i]);
    }
}

MU_TEST_SUITE(suite_postorder)
{
    MU_SUITE_CONFIGURE(NULL, NULL, postorder_test_setup, postorder_test_teardown);
    MU_RUN_TEST(test_postorder);
}