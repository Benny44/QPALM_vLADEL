#include "minunit.h"
#include "test_scale.h"
#include "test_matvec.h"
#include "test_matmat.h"
#include "test_upper_diag.h"
#include "test_etree.h"
#include "test_postorder.h"
#include "test_transpose.h"
#include "test_col_counts.h"
#include "test_permutation.h"
#include "test_ldl.h"
#include "test_rank1_mod.h"
#include "test_row_mod.h"
#include "test_add.h"
#include "test_submatrix.h"

int main(){
    MU_INITIALIZE();
    MU_RUN_SUITE(suite_scale);
    MU_RUN_SUITE(suite_matvec);
    MU_RUN_SUITE(suite_matmat);
    MU_RUN_SUITE(suite_add);
    MU_RUN_SUITE(suite_submatrix);
    MU_RUN_SUITE(suite_upper_diag);
    MU_RUN_SUITE(suite_etree);
    MU_RUN_SUITE(suite_postorder);
    MU_RUN_SUITE(suite_transpose);
    MU_RUN_SUITE(suite_col_counts);
    MU_RUN_SUITE(suite_permutation);
    MU_RUN_SUITE(suite_ldl);
    MU_RUN_SUITE(suite_rank1_mod);
    MU_RUN_SUITE(suite_row_mod);
    MU_REPORT();
    return minunit_fail;
}