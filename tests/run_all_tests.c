#include "minunit.h"
#include "test_lin_alg.h"
#include "test_basic_qp.h"
#include "test_prim_inf_qp.h"
#include "test_dua_inf_qp.h"
#include "test_degen_hess.h"
#include "test_solver_interface.h"
#include "test_nonconvex_qp.h"
#include "test_update.h"
#include "test_validate.h"
#include "test_error_handling.h"
#include "test_ls_qp.h"
#include "test_medium_qp.h"
#include "test_casadi_general_unconstrained.h"
#include "test_casadi_general_convex_sparse.h"
#include "test_update_after_setup.h"


int main(){
    MU_INITIALIZE();
    MU_RUN_SUITE(suite_validation);
    MU_RUN_SUITE(suite_error_handling);
    MU_RUN_SUITE(suite_lin_alg);
    MU_RUN_SUITE(suite_solver);
    MU_RUN_SUITE(suite_basic_qp);
    MU_RUN_SUITE(suite_medium_qp);
    MU_RUN_SUITE(suite_prim_inf_qp);
    MU_RUN_SUITE(suite_dua_inf_qp);
    MU_RUN_SUITE(suite_degen_hess);
    #ifdef COMPILE_NONCONVEX
    MU_RUN_SUITE(suite_nonconvex);
    #endif
    MU_RUN_SUITE(suite_update);
    MU_RUN_SUITE(suite_ls_qp);
    MU_RUN_SUITE(suite_casadi_general_unconstrained);
    MU_RUN_SUITE(suite_casadi_general_convex_sparse);
    MU_RUN_SUITE(suite_update_after_setup);
    MU_REPORT();
    
    return minunit_fail; /* =0 if all tests passed, >0 otherwise */
	    
}