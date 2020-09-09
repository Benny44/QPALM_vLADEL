#include "ladel_constants.h"
#include "ladel_global.h"
#include "ladel_types.h"
#include "ladel_copy.h"
#include "ladel_matmat.h"

ladel_sparse_matrix *ladel_mat_mat_transpose(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_work *work)
{
    return ladel_mat_mat_transpose_advanced(M, M_transpose, NULL, TRUE, work);
}

ladel_sparse_matrix *ladel_mat_mat_transpose_pattern(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_work *work)
{
    return ladel_mat_mat_transpose_advanced(M, M_transpose, NULL, FALSE, work);
}

ladel_sparse_matrix *ladel_mat_diag_mat_transpose(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_double *diag, ladel_work *work)
{
    return ladel_mat_mat_transpose_advanced(M, M_transpose, diag, TRUE, work);
}

ladel_sparse_matrix *ladel_mat_mat_transpose_advanced(ladel_sparse_matrix *M, ladel_sparse_matrix *M_transpose, ladel_double *diag, ladel_int values, ladel_work *work)
{
    if (!M || !M_transpose || !work) return NULL;
    
    ladel_int col, row, row2, index, index2, MMt_nnz = 0;
    ladel_int *touched = work->array_int_ncol_flag;
    ladel_double *MMt_col = work->array_double_all_zeros_ncol1;

    /* Count nnz in M*M_tranpose */
    for (col = 0; col < M_transpose->ncol; col++)
    {
        work->flag++; /*this makes it so that touched==work->flag is false everywhere*/ 
        LADEL_FOR(index, M_transpose, col)
        {
            row = M_transpose->i[index];
            LADEL_FOR(index2, M, row)
            {
                row2 = M->i[index2];
                if (row2 > col) break;
                if (touched[row2] != work->flag)
                {
                    touched[row2] = work->flag;
                    MMt_nnz++; 
                } 
            }
        }
    }

    ladel_sparse_matrix *MMt = ladel_sparse_alloc(M->nrow, M->nrow, MMt_nnz, UPPER, values && M->values, FALSE);
    if (!MMt) return NULL;

    /* Compute M*diag*M_transpose */
    if (MMt->values) for (index = 0; index < MMt_nnz; index++) MMt->x[index] = 0;
    MMt->p[0] = 0;
    MMt_nnz = -1;

    for (col = 0; col < M_transpose->ncol; col++)
    {
        work->flag++;
        LADEL_FOR(index, M_transpose, col)
        {
            row = M_transpose->i[index];
            LADEL_FOR(index2, M, row)
            {
                row2 = M->i[index2];
                
                if (row2 > col) break;
                if (touched[row2] != work->flag)
                {
                    MMt_nnz++;
                    touched[row2] = work->flag;
                    MMt->i[MMt_nnz] = row2;
                } 
                if (MMt->values) 
                    MMt_col[row2] += (diag) ? M->x[index2]*diag[row]*M_transpose->x[index] : M->x[index2]*M_transpose->x[index];
            }
        }
        MMt->p[col+1] = MMt_nnz+1;
        
        if (MMt->values)
        {
            LADEL_FOR(index, MMt, col)
            {
                MMt->x[index] = MMt_col[MMt->i[index]];
                MMt_col[MMt->i[index]] = 0;
            }
        }
    }

    return MMt;
}
