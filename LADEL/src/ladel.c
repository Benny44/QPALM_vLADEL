#include "ladel_types.h"
#include "ladel_constants.h"
#include "ladel_global.h"
#include "ladel_ldl_symbolic.h"
#include "ladel_ldl_numeric.h"
#include "ladel_permutation.h"
#include "ladel_etree.h"
#include "ladel_debug_print.h"
#include "ladel.h"

ladel_int ladel_factorize(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_work* work)
{
    ladel_diag d;
    d.diag_size = 0;
    return ladel_factorize_with_diag(M, d, sym, ordering_method, LD, work);
}

ladel_int ladel_factorize_with_diag(ladel_sparse_matrix *M, ladel_diag d, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_work* work)
{
    if (!M || !sym || !work) return FAIL;

    ladel_int ok_symbolic, ok_numeric;
    ladel_sparse_matrix *Mpp;
    
    if (ordering_method != NO_ORDERING) Mpp = ladel_sparse_alloc(M->nrow, M->ncol, M->nzmax, M->symmetry, M->values, FALSE);
    else Mpp = M;

    if (!Mpp) return FAIL;
    ok_symbolic = ladel_ldl_symbolic(M, sym, ordering_method, Mpp, work);
    if (ok_symbolic == FAIL) return FAIL;

    *LD = ladel_factor_allocate(sym);
    if (!*LD)
    {
        if (ordering_method != NO_ORDERING) ladel_sparse_free(Mpp);
        return FAIL;
    }
    if (!Mpp) return FAIL;
    ok_numeric = ladel_ldl_numeric_with_diag(Mpp, d, sym, *LD, work);

    if (ordering_method != NO_ORDERING) ladel_sparse_free(Mpp);
    if (ok_symbolic && ok_numeric) return SUCCESS;
    else return FAIL;
}

ladel_int ladel_factorize_advanced(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_sparse_matrix *Mbasis, ladel_work* work)
{
    ladel_diag d;
    d.diag_size = 0;
    return ladel_factorize_advanced_with_diag(M, d, sym, ordering_method, LD, Mbasis, work);
}

ladel_int ladel_factorize_advanced_with_diag(ladel_sparse_matrix *M, ladel_diag d, ladel_symbolics *sym, ladel_int ordering_method, ladel_factor **LD, ladel_sparse_matrix *Mbasis, ladel_work* work)
{
    if (!M || !sym || !Mbasis || !work) return FAIL;

    ladel_int ok_symbolic, ok_numeric;
    ladel_sparse_matrix *Mpp;
    
    if (ordering_method != NO_ORDERING) Mpp = ladel_sparse_alloc(Mbasis->nrow, Mbasis->ncol, Mbasis->nzmax, Mbasis->symmetry, Mbasis->values, FALSE);
    else Mpp = Mbasis;

    if (!Mpp) return FAIL;
    ok_symbolic = ladel_ldl_symbolic(Mbasis, sym, ordering_method, Mpp, work);
    *LD = ladel_factor_allocate(sym);
    if (!*LD)
    {
        if (ordering_method != NO_ORDERING) ladel_sparse_free(Mpp);
        return FAIL;
    }
    if (sym->p)
    {
        ladel_sparse_free(Mpp);
        Mpp = ladel_sparse_alloc(M->nrow, M->ncol, M->nzmax, M->symmetry, M->values, FALSE);
        ladel_permute_symmetric_matrix(M, sym->p, Mpp, work);
    } else
    {
        Mpp = M;
    }
    
    ladel_etree(Mpp, sym, work);

    ok_numeric = ladel_ldl_numeric_with_diag(Mpp, d, sym, *LD, work);

    if (ordering_method != NO_ORDERING) ladel_sparse_free(Mpp);
    if (ok_symbolic && ok_numeric) return SUCCESS;
    else return FAIL;
}

ladel_int ladel_factorize_with_prior_basis(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work)
{
    ladel_diag d;
    d.diag_size = 0;
    return ladel_factorize_with_prior_basis_with_diag(M, d, sym, LD, work);
}

ladel_int ladel_factorize_with_prior_basis_with_diag(ladel_sparse_matrix *M, ladel_diag d, ladel_symbolics *sym, ladel_factor *LD, ladel_work* work)
{
    if (!M || !sym || !LD || !work) return FAIL;

    ladel_int ok_numeric;
    ladel_sparse_matrix *Mpp;

    if (sym->p)
    {
        Mpp = ladel_sparse_alloc(M->nrow, M->ncol, M->nzmax, M->symmetry, M->values, FALSE);
        ladel_permute_symmetric_matrix(M, sym->p, Mpp, work);
    } else
    {
        Mpp = M;
    }

    ladel_etree(Mpp, sym, work);
    ok_numeric = ladel_ldl_numeric_with_diag(Mpp, d, sym, LD, work);

    if (sym->p) ladel_sparse_free(Mpp);
    return ok_numeric;
}

ladel_int ladel_dense_solve(const ladel_factor *LD, const ladel_double *rhs, ladel_double *y, ladel_work* work)
{
    if (!LD || !rhs || !y || !work) return FAIL;
    
    ladel_sparse_matrix *L = LD->L;
    ladel_double *Dinv = LD->Dinv;
    ladel_int index, row, ncol = LD->L->ncol;
    if (LD->p) for (row = 0; row < ncol; row++) y[row] = rhs[LD->p[row]];
    else       for (row = 0; row < ncol; row++) y[row] = rhs[row];

    for (row = 0; row < ncol; row++)
    {
        for (index = L->p[row]; index < L->p[row]+L->nz[row]; index++)
        {
            y[L->i[index]] -= L->x[index]*y[row];
        }
    }
    for (row = 0; row < ncol; row++) y[row] *= Dinv[row];
    for (row = ncol-1; row >= 0; row--)
    {
        for (index = L->p[row]; index < L->p[row]+L->nz[row]; index++)
        {
            y[row] -= L->x[index]*y[L->i[index]];
        }
    } 

    if (LD->p)
    {
        ladel_double *temp = work->array_double_all_zeros_ncol1;
        for (row = 0; row < ncol; row++) temp[row] = y[row];
        for (row = 0; row < ncol; row++) 
        {
            y[LD->p[row]] = temp[row];
            temp[row] = 0.0; /*reset to keep the workspace consistent*/
        }
    }
    return SUCCESS;
}