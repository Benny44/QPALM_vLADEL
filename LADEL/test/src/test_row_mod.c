#include "minunit.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel_types.h"
#include "ladel.h"
#include "ladel_row_mod.h"
#include "ladel_debug_print.h"

#define NCOL 6
#define NROW 6
#define NZMAX 21
#define TOL 1e-8

static ladel_work *work;
static ladel_sparse_matrix *M, *Mbasis, *new_row;
static ladel_double new_diag;
static ladel_factor *LD;
static ladel_symbolics *sym;

void row_mod_suite_setup(void)
{
    work = ladel_workspace_allocate(NCOL);

    M = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE, FALSE);
    M->p[0] = 0; M->p[1] = 1; M->p[2] = 3; M->p[3] = 6; M->p[4] = 10; M->p[5] = 15; M->p[6] = 16; 
    M->i[0] = 0; M->i[1] = 0; M->i[2] = 1; M->i[3] = 0; M->i[4] = 1; M->i[5] = 2; M->i[6] = 0; 
    M->i[7] = 1; M->i[8] = 2; M->i[9] = 3; M->i[10] = 0; M->i[11] = 1; M->i[12] = 2; M->i[13] = 3; 
    M->i[14] = 4; M->i[15] = 5;
    M->x[0] = 2.055098610936521e-01; M->x[1] = 2.572819890455917e-01; M->x[2] = 2.417367221483280e-02; 
    M->x[3] = 6.893733655105202e-01; M->x[4] = 4.565212273616481e-01; M->x[5] = 8.270829803188939e-02; 
    M->x[6] = 4.517532293767544e-01; M->x[7] = 3.272091623419543e-01; M->x[8] = 9.357072685394627e-01; 
    M->x[9] = 3.748943405611161e-01; M->x[10] = 6.398162776709431e-01; M->x[11] = 3.281149448958269e-01; 
    M->x[12] = 5.614789725017673e-01; M->x[13] = 4.335189702823217e-01; M->x[14] = 9.223654191512636e-01; 
    M->x[15] = 1.0;

    Mbasis = ladel_sparse_alloc(NROW, NCOL, NZMAX, UPPER, TRUE, FALSE);
    Mbasis->p[0] = 0; Mbasis->p[1] = 1; Mbasis->p[2] = 3; Mbasis->p[3] = 6; Mbasis->p[4] = 10; Mbasis->p[5] = 15; Mbasis->p[6] = 21; 
    Mbasis->i[0] = 0; Mbasis->i[1] = 0; Mbasis->i[2] = 1; Mbasis->i[3] = 0; Mbasis->i[4] = 1; Mbasis->i[5] = 2; Mbasis->i[6] = 0; 
    Mbasis->i[7] = 1; Mbasis->i[8] = 2; Mbasis->i[9] = 3; Mbasis->i[10] = 0; Mbasis->i[11] = 1; Mbasis->i[12] = 2; Mbasis->i[13] = 3; 
    Mbasis->i[14] = 4; Mbasis->i[15] = 0; Mbasis->i[16] = 1; Mbasis->i[17] = 2; Mbasis->i[18] = 3; Mbasis->i[19] = 4; Mbasis->i[20] = 5; 
    Mbasis->x[0] = 2.055098610936521e-01; Mbasis->x[1] = 2.572819890455917e-01; Mbasis->x[2] = 2.417367221483280e-02; 
    Mbasis->x[3] = 6.893733655105202e-01; Mbasis->x[4] = 4.565212273616481e-01; Mbasis->x[5] = 8.270829803188939e-02; 
    Mbasis->x[6] = 4.517532293767544e-01; Mbasis->x[7] = 3.272091623419543e-01; Mbasis->x[8] = 9.357072685394627e-01; 
    Mbasis->x[9] = 3.748943405611161e-01; Mbasis->x[10] = 6.398162776709431e-01; Mbasis->x[11] = 3.281149448958269e-01; 
    Mbasis->x[12] = 5.614789725017673e-01; Mbasis->x[13] = 4.335189702823217e-01; Mbasis->x[14] = 9.223654191512636e-01; 
    Mbasis->x[15] = 1.0; Mbasis->x[16] = 1.0; Mbasis->x[17] = 1.0; Mbasis->x[18] = 1.0; Mbasis->x[19] = 1.0; Mbasis->x[20] = 1.0; 

    new_row = ladel_sparse_alloc(NROW, 1, NROW, UNSYMMETRIC, TRUE, TRUE);
    new_row->p[0] = 0; new_row->p[1] = 5;
    new_row->nz[0] = 5;
    new_row->i[0] = 0; new_row->i[1] = 1; new_row->i[2] = 2; new_row->i[3] = 3; new_row->i[4] = 4;   
    new_row->x[0] = 8.746015987152932e-01; new_row->x[1] = 7.532514851657419e-01; new_row->x[2] = 2.226659543716385e-01; 
    new_row->x[3] = 1.380149620432186e-01; new_row->x[4] = 2.208883363414058e-01;   


    new_diag = -1.3;  

    sym = ladel_symbolics_alloc(NCOL);
}

void row_mod_suite_teardown(void)
{
    ladel_sparse_free(M);
    ladel_sparse_free(Mbasis);
    ladel_symbolics_free(sym);
    ladel_sparse_free(new_row);
    ladel_workspace_free(work);
    LD = ladel_factor_free(LD);
}

void row_mod_test_teardown(void)
{
    
}

MU_TEST(test_row_add_at_the_end)
{
    /*Add a row to the end*/
    ladel_int status, index;
    status = ladel_factorize_advanced(M, sym, NO_ORDERING, &LD, Mbasis, work);
    mu_assert_long_eq(status, SUCCESS);

    ladel_double rhs[6] = {5.679645458145579e-01, 6.630438231677197e-01, 5.274175204016522e-01, 
                            5.634939173447451e-02, 3.631625949531299e-01, 0};
    ladel_double sol[6] = {-4.192429648798019e-01, -2.456611750539314e+00, 1.398638846414751e-01,
                            1.636082172326339e+00, 7.043291031616438e-01, 0};

    ladel_double x[6];
    ladel_dense_solve(LD, rhs, x, work);
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol[index], TOL);

    /*Add row to the end*/
    rhs[5] = 9.677274091080609e-01;
    status = ladel_row_add(LD, sym, 5, new_row, 0, new_diag, work);
    mu_assert_long_eq(status, SUCCESS);

    ladel_double sol_mod[6] = {2.609664443874901e+00, 3.027486501422937e+00, -9.998957722639763e-01, 
                                -2.530237213301512e+00, -1.203546137577670e+00, 2.121111804411014e+00};
    ladel_dense_solve(LD, rhs, x, work);
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol_mod[index], TOL);
}

MU_TEST(test_row_del)
{
    /*Delete the third row */ 
    ladel_int status, index;

    status = ladel_row_del(LD, sym, 2, work);
    mu_assert_long_eq(status, SUCCESS);

    ladel_double rhs[6] = {5.679645458145579e-01, 6.630438231677197e-01, 5.274175204016522e-01, 
                            5.634939173447451e-02, 3.631625949531299e-01, 9.677274091080609e-01};

    ladel_double sol[6] = {4.509218534985005e+00, 4.780597529712503e+00, 5.274175204016522e-01,
                            -1.025501171241090e+01, -5.435956090455947e-01, 3.878165798320269e+00};
    ladel_double x[6];
    ladel_dense_solve(LD, rhs, x, work);
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol[index], TOL);
}

MU_TEST(test_row_add_in_middle)
{
    /* Add the third row back */
    new_row->x[0] = 6.893733655105202e-01; new_row->x[1] = 4.565212273616481e-01; 
    new_row->x[2] = 9.357072685394627e-01; new_row->x[3] = 5.614789725017673e-01;
    new_row->x[4] = 2.226659543716385e-01;
    new_diag = 8.270829803188939e-02;
    new_row->i[0] = 0; new_row->i[1] = 1; new_row->i[2] = 3; new_row->i[3] = 4; new_row->i[4] = 5; 
    new_row->p[1] = 5;
    new_row->nz[0] = 5;

    ladel_int status, index;
    status = ladel_row_add(LD, sym, 2, new_row, 0, new_diag, work);
    mu_assert_long_eq(status, SUCCESS);

    ladel_double rhs[6] = {5.679645458145579e-01, 6.630438231677197e-01, 5.274175204016522e-01, 
                            5.634939173447451e-02, 3.631625949531299e-01, 9.677274091080609e-01};

    ladel_double sol[6] = {2.609664443874901e+00, 3.027486501422937e+00, -9.998957722639763e-01, 
                                -2.530237213301512e+00, -1.203546137577670e+00, 2.121111804411014e+00};
    
    ladel_double x[6];
    ladel_dense_solve(LD, rhs, x, work);
    for (index = 0; index < NCOL; index++) mu_assert_double_eq(x[index], sol[index], TOL);
}

MU_TEST_SUITE(suite_row_mod)
{
    MU_SUITE_CONFIGURE(row_mod_suite_setup, row_mod_suite_teardown, NULL, row_mod_test_teardown);
    MU_RUN_TEST(test_row_add_at_the_end);
    MU_RUN_TEST(test_row_del);
    MU_RUN_TEST(test_row_add_in_middle); /*Note, this test depends on test_row_del*/
}