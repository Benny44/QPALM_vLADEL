#include "ladel_types.h"
#include "ladel_constants.h"
#include "ladel_global.h"
#include "ladel_permutation.h"
#include "ladel_etree.h"
#include "ladel_postorder.h"
#include "ladel_col_counts.h"
#include "ladel_debug_print.h"

#ifdef DAMD
#include "amd.h"
#endif /*DAMD*/

ladel_int ladel_ldl_symbolic(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_sparse_matrix *Mpp, ladel_work* work)
{
    if (!M || !sym || !Mpp || !work) return FAIL;

    ladel_sparse_matrix *Mwork = M;
    if (ordering_method == AMD)
    {
        #ifdef DAMD
        ladel_int status;
        double Info [AMD_INFO];

        #ifdef DLONG
        status = amd_l_order(M->ncol, M->p, M->i, sym->p, NULL, Info);
        #else /*DLONG*/
        status = amd_order(M->ncol, M->p, M->i, sym->p, NULL, Info);
        #endif
        if (status != AMD_OK) return FAIL;
        
        #else /*DAMD*/
        sym->p = ladel_free(sym->p);
        #endif
    } else if (ordering_method == GIVEN_ORDERING)
    {
        /*do nothing, sym->p already contains the permutation*/
    } else if (ordering_method == NO_ORDERING)
    {
        sym->p = ladel_free(sym->p);
    }
    
    if (sym->p)
    {
        ladel_permute_symmetric_matrix(M, sym->p, Mpp, work);
        Mwork = Mpp;
        ladel_invert_permutation_vector(sym->p, sym->pinv, M->ncol);
    }

    #ifdef SIMPLE_COL_COUNTS
    ladel_etree_and_col_counts(Mwork, sym, work);
    #else
    ladel_etree(Mwork, sym, work);
    ladel_postorder(Mwork, sym, work);
    ladel_col_counts(Mwork, sym, work);
    #endif /* SIMPLE_COL_COUNTS */

    return SUCCESS;
}