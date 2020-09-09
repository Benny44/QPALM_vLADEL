#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_constants.h"

ladel_sparse_matrix *ladel_transpose(ladel_sparse_matrix *M, ladel_int values, ladel_work* work)
{
    if (!M) return NULL;
    ladel_int *col_pointers, index;
    if (work) col_pointers = work->array_int_ncol1;
    else col_pointers = ladel_malloc(M->nrow, sizeof(ladel_int));
    
    for (index = 0; index < M->nrow; index++) col_pointers[index] = 0;

    ladel_sparse_matrix *M_transpose = ladel_sparse_alloc(M->ncol, M->nrow, M->nzmax, -M->symmetry, values && M->values, FALSE);
    
    if (!M_transpose) 
    {
        if (!work) ladel_free(col_pointers);
        return NULL; 
    }

    ladel_int col, new_index, prev_col_count;
    for (col = 0; col < M->ncol; col++)
        LADEL_FOR(index, M, col)
            col_pointers[M->i[index]]++;
    
    M_transpose->p[0] = 0;
    for (col = 1; col < M_transpose->ncol; col++)
    {
        prev_col_count = col_pointers[col-1];
        col_pointers[col] += prev_col_count;
        M_transpose->p[col] = prev_col_count;
        col_pointers[col-1] = M_transpose->p[col-1];
    } 
    M_transpose->p[M_transpose->ncol] = col_pointers[M_transpose->ncol-1];
    col_pointers[M_transpose->ncol-1] = M_transpose->p[M_transpose->ncol-1];

    for (col = 0; col < M->ncol; col++)
    {
        LADEL_FOR(index, M, col)
        {
            new_index = col_pointers[M->i[index]]++;
            M_transpose->i[new_index] = col;
            if (M_transpose->values) M_transpose->x[new_index] = M->x[index];
        }
    }

    if (!work) ladel_free(col_pointers);
    return M_transpose;
}
