#include "ladel_global.h"
#include "ladel_types.h"
#include "ladel_constants.h"
#include "ladel_copy.h"
#include "ladel_permutation.h"
#include <stdlib.h>
#include <stdio.h>

#ifdef MATLAB
#include "mex.h"
void* ladel_calloc(ladel_int n, size_t size)
{
    void *m = mxCalloc(LADEL_MAX(n, 1), size);
    mexMakeMemoryPersistent(m);
    return m;
}

void *ladel_malloc(ladel_int n, size_t size) 
{
    void *m = mxMalloc(LADEL_MAX(n, 1) * size);
    mexMakeMemoryPersistent(m);
    return m;
}

void* ladel_realloc(void *p, ladel_int n, size_t size, ladel_int *status) 
{
    void *p_new = mxRealloc(p, LADEL_MAX(n, 1) * size);
    *status = (p_new != NULL);
    mexMakeMemoryPersistent(p_new);
    return ((*status) ? p_new : p);
}

void *ladel_free(void* p) 
{
    if (p) mxFree(p);
    return NULL;
}

#elif defined PYTHON
#include <Python.h>
void *ladel_malloc(ladel_int n, size_t size) 
{
    void *m = PyMem_Malloc(LADEL_MAX(n, 1) * size);
    return m;
} 
void *ladel_free(void* p) 
{
    if (p) PyMem_Free(p);
    return NULL;
}
void* ladel_realloc(void *p, ladel_int n, size_t size, ladel_int *status) 
{
    void *p_new = PyMem_Realloc(p, LADEL_MAX(n, 1) * size);
    *status = (p_new != NULL);
    return ((*status) ? p_new : p);
}
void* ladel_calloc(ladel_int n, size_t size)
{
    #if PY_MAJOR_VERSION >= 3
    void *m = PyMem_Calloc(LADEL_MAX(n, 1), size);
    #else
    void *m = PyMem_Malloc(num * size);
    memset(m, 0, num * size);
    #endif
    return m;
}

#else
void *ladel_malloc(ladel_int n, size_t size) 
{
    return (malloc(LADEL_MAX(n, 1) * size));
}

void *ladel_calloc(ladel_int n, size_t size)
{
    return (calloc(LADEL_MAX(n, 1), size));
}

void *ladel_free(void* p) 
{
    if (p) free(p);
    return NULL;
}

void *ladel_realloc(void *p, ladel_int n, size_t size, ladel_int *status)
{
    void *p_new;
    p_new = realloc(p, LADEL_MAX(n, 1) * size);
    *status = (p_new != NULL);
    return ((*status) ? p_new : p);
}

#endif /*MATLAB*/

ladel_sparse_matrix *ladel_sparse_free(ladel_sparse_matrix *M)
{
    if (!M) return NULL;
    ladel_free(M->p);
    ladel_free(M->i);
    ladel_free(M->x);
    ladel_free(M->nz);
    return ((ladel_sparse_matrix *) ladel_free(M));
}

ladel_sparse_matrix *ladel_sparse_alloc(ladel_int nrow, ladel_int ncol, 
                                                ladel_int nzmax, ladel_int symmetry,
                                                ladel_int values, ladel_int nz)
{
    ladel_sparse_matrix *M = (ladel_sparse_matrix *) ladel_calloc(1, sizeof(ladel_sparse_matrix));
    if (!M) return NULL;
    M->nrow = nrow;
    M->ncol = ncol;
    M->nzmax = nzmax;
    M->values = values;
    M->symmetry = symmetry;
    M->p = (ladel_int *) ladel_malloc(ncol+1, sizeof(ladel_int));
    M->i = (ladel_int *) ladel_malloc(nzmax, sizeof(ladel_int));
    M->x = values ? (ladel_double *) ladel_malloc(nzmax, sizeof(ladel_double)) : NULL;
    M->nz = nz ? (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int)) : NULL;
    if (!M->p || !M->i || (values && !M->x) || (nz && !M->nz)) M = ladel_sparse_free(M);
    return M;
}

ladel_sparse_matrix *ladel_sparse_alloc_empty(ladel_int nrow, ladel_int ncol, 
                                                ladel_int symmetry, ladel_int values, 
                                                ladel_int nz)
{
    ladel_sparse_matrix *M = (ladel_sparse_matrix *) ladel_calloc(1, sizeof(ladel_sparse_matrix));
    if (!M) return NULL;
    M->nrow = nrow;
    M->ncol = ncol;
    M->nzmax = 0;
    M->values = values;
    M->symmetry = symmetry;
    M->p = (ladel_int *) ladel_calloc(ncol+1, sizeof(ladel_int));
    M->i = (ladel_int *) ladel_malloc(1, sizeof(ladel_int));
    M->x = values ? (ladel_double *) ladel_malloc(1, sizeof(ladel_double)) : NULL;
    M->nz = nz ? (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int)) : NULL;
    if (!M->p || !M->i || (values && !M->x) || (nz && !M->nz)) M = ladel_sparse_free(M);
    return M;
}

ladel_int ladel_sparse_realloc(ladel_sparse_matrix* M, ladel_int nzmax)
{
    ladel_int status_i, status, status_x = SUCCESS;
    if (!M) return FAIL;
    if (nzmax <= 0) nzmax = M->p[M->ncol];
    M->i = (ladel_int *) ladel_realloc(M->i, nzmax, sizeof(ladel_int), &status_i);
    if (M->values) M->x = (ladel_double *) ladel_realloc(M->x, nzmax, sizeof(ladel_double), &status_x);
    status = status_i && status_x;
    if (status == SUCCESS) M->nzmax = nzmax;
    return status;
}

ladel_symbolics *ladel_symbolics_free(ladel_symbolics *sym)
{
    if (!sym) return NULL;
    ladel_free(sym->etree);
    ladel_free(sym->postorder);
    ladel_free(sym->col_counts);
    ladel_free(sym->p);
    ladel_free(sym->pinv);
    ladel_free(sym->pattern);
    ladel_free(sym->nodes);
    return (ladel_symbolics *) ladel_free(sym);
}

ladel_symbolics *ladel_symbolics_alloc(ladel_int ncol)
{
    ladel_symbolics *sym = (ladel_symbolics *) ladel_calloc(1, sizeof(ladel_symbolics));
    if (!sym) return NULL;
    sym->ncol = ncol;
    sym->etree = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->postorder = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->col_counts = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->p = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->pinv = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->pattern = (ladel_int *) ladel_malloc(ncol, sizeof(ladel_int));
    sym->nodes = (ladel_int *) ladel_calloc(ncol, sizeof(ladel_int));
    if (!sym->etree || !sym->postorder || !sym->col_counts || !sym->pattern || !sym->nodes) 
        sym = ladel_symbolics_free(sym);
    return sym;
}

ladel_factor *ladel_factor_free(ladel_factor *LD)
{
    if (!LD) return NULL;
    ladel_sparse_free(LD->L);
    ladel_free(LD->D);
    ladel_free(LD->Dinv);
    ladel_free(LD->p);
    ladel_free(LD->pinv);
    return (ladel_factor *) ladel_free(LD);
}

ladel_factor *ladel_factor_allocate(ladel_symbolics *sym)
{
    ladel_factor *LD = (ladel_factor *) ladel_calloc(1, sizeof(ladel_factor));
    if (!LD || !sym) return NULL;
    ladel_int ncol = LD->ncol = sym->ncol;
    LD->L = ladel_sparse_alloc(ncol, ncol, sym->col_counts[ncol-1], UNSYMMETRIC, TRUE, TRUE);
    LD->D = ladel_malloc(ncol, sizeof(ladel_double));
    LD->Dinv = ladel_malloc(ncol, sizeof(ladel_double));
    if (!LD->L || !LD->D || !LD->Dinv)
    {
        ladel_factor_free(LD);
        return NULL;
    }
    if (sym->p)
    {
        LD->p = ladel_malloc(ncol, sizeof(ladel_int));
        LD->pinv = ladel_malloc(ncol, sizeof(ladel_int));
        if (!LD->p || !LD->pinv)
        {
            ladel_factor_free(LD);
            return NULL;
        }
        ladel_int_vector_copy(sym->p, ncol, LD->p);
        ladel_int_vector_copy(sym->pinv, ncol, LD->pinv);
    } else 
    {
        LD->p = NULL;
        LD->pinv = NULL;
    }
    return LD; 
}

ladel_set *ladel_set_free(ladel_set *set)
{
    if (!set) return NULL;
    ladel_free(set->set);
    return (ladel_set *) ladel_free(set);
}

ladel_set *ladel_set_allocate(ladel_int max_size)
{
    ladel_set* set = ladel_malloc(1, sizeof(ladel_set));
    if (!set) return NULL;
    set->set = ladel_malloc(max_size, sizeof(ladel_int));
    if (!set->set) 
    {
        ladel_set_free(set);
        return NULL;
    }
    set->max_size_set = max_size;
    return set;
}

void ladel_set_set(ladel_set *set, ladel_int *set_vals, ladel_int size_set, ladel_int max_size_set)
{
    set->set = set_vals;
    set->size_set = size_set;
    set->max_size_set = max_size_set;
}

ladel_work *ladel_workspace_free(ladel_work* work)
{
    if (!work) return NULL;
    ladel_set_free(work->set_preallocated1);
    ladel_set_free(work->set_preallocated2);
    ladel_set_free(work->set_preallocated3);
    ladel_free(work->set_unallocated_values1);
    ladel_free(work->set_unallocated_values2);
    ladel_free(work->set_unallocated_values3);
    ladel_free(work->array_int_ncol1);
    ladel_free(work->array_int_ncol2);
    ladel_free(work->array_int_ncol3);
    ladel_free(work->array_int_ncol4);
    ladel_free(work->array_int_ncol_flag);
    ladel_free(work->array_double_all_zeros_ncol1);
    ladel_free(work->array_double_ncol1);
    return (ladel_work *) ladel_free(work);
}

ladel_work *ladel_workspace_allocate(ladel_int ncol)
{
    ladel_work *work = ladel_malloc(1, sizeof(ladel_work));
    if (!work) return NULL;
    work->set_preallocated1 = ladel_set_allocate(ncol);
    work->set_preallocated2 = ladel_set_allocate(ncol);
    work->set_preallocated3 = ladel_set_allocate(ncol);
    work->set_unallocated_values1 = ladel_malloc(1, sizeof(ladel_set));
    work->set_unallocated_values2 = ladel_malloc(1, sizeof(ladel_set));
    work->set_unallocated_values3 = ladel_malloc(1, sizeof(ladel_set));
    work->array_int_ncol1 = ladel_malloc(ncol, sizeof(ladel_int));
    work->array_int_ncol2 = ladel_malloc(ncol, sizeof(ladel_int));
    work->array_int_ncol3 = ladel_malloc(ncol, sizeof(ladel_int));
    work->array_int_ncol4 = ladel_malloc(ncol, sizeof(ladel_int));
    work->array_int_ncol_flag = ladel_calloc(ncol, sizeof(ladel_int));
    work->flag = 1;
    work->array_double_ncol1 = ladel_malloc(ncol, sizeof(ladel_double));
    work->array_double_all_zeros_ncol1 = ladel_calloc(ncol, sizeof(ladel_double));

    if (!work->set_preallocated1 || !work->set_preallocated2 || !work->set_preallocated3 
        || !work->set_unallocated_values1 || !work->set_unallocated_values2 || !work->set_unallocated_values3
        || !work->array_int_ncol1 || !work->array_int_ncol2 || !work->array_int_ncol3 
        || !work->array_int_ncol4 || !work->array_double_all_zeros_ncol1) 
    {
        ladel_workspace_free(work);
        return NULL;   
    }
    return work;
}