#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_copy.h"

ladel_sparse_matrix *ladel_sparse_allocate_and_copy(ladel_sparse_matrix *M)
{
    ladel_sparse_matrix *M_copy = ladel_sparse_alloc(M->nrow, M->ncol, M->nzmax, M->symmetry, M->values, (ladel_int) M->nz);
    ladel_sparse_copy(M, M_copy);
    return M_copy;
}

void ladel_sparse_copy(ladel_sparse_matrix *M, ladel_sparse_matrix *M_copy)
{
    if (!M || !M_copy)
    {
        M_copy = NULL;
    } else
    {
        M_copy->ncol = M->ncol;
        M_copy->nrow = M->nrow;
        M_copy->nzmax = M->nzmax;
        M_copy->symmetry = M->symmetry;
        M_copy->values = M->values;
        ladel_int index;
        for (index = 0; index < M->ncol+1; index++) M_copy->p[index] = M->p[index];
        
        if (M->nz) for (index = 0; index < M->ncol; index++) M_copy->nz[index] = M->nz[index];
        else M_copy->nz = ladel_free(M_copy->nz);
            
        for (index = 0; index < M->nzmax; index++)
        {
            M_copy->i[index] = M->i[index];
            if (M->values) M_copy->x[index] = M->x[index];
        }
    }
}

void ladel_int_vector_copy(ladel_int *x, ladel_int size, ladel_int *y)
{
    ladel_int index;
    for (index = 0; index < size; index++)
    {
        y[index] = x[index];
    }
}

void ladel_double_vector_copy(ladel_double *x, ladel_int size, ladel_double *y)
{
    ladel_int index;
    for (index = 0; index < size; index++)
    {
        y[index] = x[index];
    }
}