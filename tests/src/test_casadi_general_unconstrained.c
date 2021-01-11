#include "minunit.h"
#include "qpalm.h"

#define N 2
#define M 0
#define ANZMAX 0
#define QNZMAX 2
#define TOL 1e-12

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
solver_common *c;
solver_common common;


void casadi_general_unconstrained_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-12;
    settings->eps_rel = 1e-12;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;
    data->q = c_calloc(N,sizeof(c_float));
    data->q[0] = -0.7; data->q[1] = -2.3;   
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmax = c_calloc(M,sizeof(c_float));

    // solver_common common;
    c = &common;
    solver_sparse *A = ladel_sparse_alloc(M, N, ANZMAX, UNSYMMETRIC, TRUE, FALSE);
    solver_sparse *Q = ladel_sparse_alloc(N, N, QNZMAX, UPPER, TRUE, FALSE);

    c_float *Qx;
    c_int *Qi, *Qp;
    Qx = Q->x;
    Qp = Q->p;
    Qi = Q->i;
    Qx[0] = 1.0; Qx[1] = 1.0;   
    Qp[0] = 0; Qp[1] = 1; Qp[2] = 2; 
    Qi[0] = 0; Qi[1] = 1;  

    ladel_to_upper_diag(Q);

    data->A = A;
    data->Q = Q;
}

void casadi_general_unconstrained_suite_teardown(void) {
    c_free(settings);
    // Clean setup
    data->Q = ladel_sparse_free(data->Q);
    data->A = ladel_sparse_free(data->A);
    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
}

void casadi_general_unconstrained_test_teardown(void) {
    qpalm_cleanup(work);
}


MU_TEST(test_casadi_general_unconstrained) {
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], 0.7, TOL);
    mu_assert_double_eq(work->solution->x[1], 2.3, TOL);
}

MU_TEST_SUITE(suite_casadi_general_unconstrained) {
    MU_SUITE_CONFIGURE(casadi_general_unconstrained_suite_setup, casadi_general_unconstrained_suite_teardown, NULL, casadi_general_unconstrained_test_teardown);
    MU_RUN_TEST(test_casadi_general_unconstrained);
}       
