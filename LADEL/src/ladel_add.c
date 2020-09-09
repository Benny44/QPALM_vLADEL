#include "ladel_constants.h"
#include "ladel_global.h"
#include "ladel_types.h"
#include "ladel_add.h"

ladel_sparse_matrix *ladel_add_matrices(ladel_double alpha, ladel_sparse_matrix* A, ladel_double beta, ladel_sparse_matrix *B, ladel_work *work)
{
    return ladel_add_matrices_advanced(alpha, A, beta, B, TRUE, work);
}

ladel_sparse_matrix *ladel_add_matrices_pattern(ladel_sparse_matrix* A, ladel_sparse_matrix *B, ladel_work *work)
{
    return ladel_add_matrices_advanced(0, A, 0, B, FALSE, work);
}

ladel_sparse_matrix *ladel_add_matrices_advanced(ladel_double alpha, ladel_sparse_matrix* A, ladel_double beta, ladel_sparse_matrix *B, ladel_int values, ladel_work *work)
{
    /* TODO: for different symmetries this will fail. */
    if (!A || !B) return NULL;
    ladel_double *temp = work->array_double_all_zeros_ncol1;
    ladel_int *touched = work->array_int_ncol_flag;
    ladel_int index, col, row, C_nnz = 0;
    ladel_int C_nrow = LADEL_MAX(A->nrow, B->nrow);
    ladel_int C_ncol = LADEL_MAX(A->ncol, B->ncol);
    ladel_int C_symmetry = (A->symmetry == B->symmetry) ? A->symmetry : UNSYMMETRIC;
    ladel_int C_values = values && (A->values || B->values);

    for (col = 0; col < C_ncol; col++)
    {
        work->flag++;
        LADEL_FOR(index, A, col)
        {
            row = A->i[index];
            if (touched[row] != work->flag) 
            {
                touched[row] = work->flag;
                C_nnz++;
            }

        }
        LADEL_FOR(index, B, col)
        {
            row = B->i[index];
            if (touched[row] != work->flag)
            {
                touched[row] = work->flag;
                C_nnz++;
            } 
        }
    }

    ladel_sparse_matrix *C = ladel_sparse_alloc(C_nrow, C_ncol, C_nnz, C_symmetry, C_values, FALSE);
    if (!C) return NULL;

    C_nnz = 0;
    C->p[0] = 0;
    for (col = 0; col < C_ncol; col++)
    {
        work->flag++;
        LADEL_FOR(index, A, col)
        {
            row = A->i[index];
            if (touched[row] != work->flag) 
            {
                touched[row] = work->flag;
                C->i[C_nnz] = row;
                C_nnz++;
            }
            if (C_values) temp[row] += A->values ? alpha*A->x[index] : 0;
        }
        LADEL_FOR(index, B, col)
        {
            row = B->i[index];
            if (touched[row] != work->flag)
            {
                touched[row] = work->flag;
                C->i[C_nnz] = row;
                C_nnz++;
            } 
            if (C_values) temp[row] += B->values ? beta*B->x[index] : 0;
        }
        C->p[col+1] = C_nnz;
        LADEL_FOR(index, C, col)
        {
            if (C_values)
            {
                row = C->i[index];
                C->x[index] = temp[row];
                temp[row] = 0;
            } 
        }   
    }  
    return C;
}