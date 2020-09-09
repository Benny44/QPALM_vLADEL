#include "ladel_types.h"
#include "ladel_constants.h"
#include "ladel_global.h"
#include <stdlib.h>

void ladel_permute_vector(ladel_double *x, ladel_int *p, ladel_int size, ladel_double *y)
{
    ladel_int index;
    for (index = 0; index < size; index++) y[index] = x[p[index]];
}

void ladel_inverse_permute_vector(ladel_double *x, ladel_int *pinv, ladel_int size, ladel_double *y)
{
    ladel_int index;
    for (index = 0; index < size; index++) y[pinv[index]] = x[index];
}

void ladel_invert_permutation_vector(ladel_int *p, ladel_int *pinv, ladel_int size)
{
    ladel_int index;
    for (index = 0; index < size; index++) pinv[p[index]] = index;
}

static int ladel_int_compare(const ladel_int *a, const ladel_int *b) 
{ 
    if (*a > *b) return 1;
    else return -1;
} 

void ladel_permute_sparse_vector(ladel_sparse_matrix *x, ladel_int col, ladel_int *p, ladel_work *work)
{
    ladel_int row, index, index_temp, xnz = x->p[col+1] - x->p[col];
    ladel_double *temp = work->array_double_all_zeros_ncol1;
    
    /* In a relatively dense case, don't sort but just do a dense permutation.
    The dense permutation takes O(nrow) computations.
    The sparse permutation takes O(xnz*log(xnz)) computations, because of the sort. */
    if (xnz > x->nrow / 5)
    {
        LADEL_FOR(index, x, col)
        {
            row = p[x->i[index]];
            temp[row] = x->x[index]; 
        }
        index = x->p[col]; 
        for (index_temp = 0; index_temp < x->nrow; index_temp++)
        {
            if (temp[index_temp] != 0.0)
            {
                x->i[index] = index_temp;
                x->x[index] = temp[index_temp];
                temp[index_temp] = 0.0; /*reset to keep the workspace consistent*/
                index++; 
            }
        }
    }
    else
    {
        LADEL_FOR(index, x, col)
        {
            row = p[x->i[index]];
            x->i[index] = row;
            temp[row] = x->x[index]; 
        }
        qsort(x->i + x->p[col], xnz, sizeof(ladel_int), (int (*) (const void *, const void *)) ladel_int_compare);
        LADEL_FOR(index, x, col)
        {
            row = x->i[index];
            x->x[index] = temp[row];
            temp[row] = 0.0;
        }
    }   
}

void ladel_permute_symmetric_matrix(ladel_sparse_matrix *M, ladel_int *p, ladel_sparse_matrix *Mpp, ladel_work* work)
{
    if (!M || !Mpp) return;
    if (!p) 
    {
        ladel_sparse_copy(M, Mpp);
    } else
    {
        ladel_int col, pcol, prow, index, pindex, prev_col_count, ncol = M->ncol;
        ladel_int *col_counts = work->array_int_ncol1, *pinv = work->array_int_ncol2;
        for (index = 0; index < ncol; index++) col_counts[index] = 0;
        for (col = 0; col < ncol; col++) pinv[p[col]] = col;
        for (col = 0; col < ncol; col++)
        {
            pcol = pinv[col];
            LADEL_FOR(index, M, col)
            {
                prow = pinv[M->i[index]];
                col_counts[LADEL_MAX(pcol, prow)]++;
            }
        }
        Mpp->p[0] = 0;
        for (col = 1; col < ncol; col++)
        {
            prev_col_count = col_counts[col-1];
            Mpp->p[col] = prev_col_count;
            col_counts[col] += prev_col_count; 
            col_counts[col-1] = Mpp->p[col-1]; 
        }
        Mpp->p[ncol] = col_counts[ncol-1];
        col_counts[ncol-1] = Mpp->p[ncol-1];

        for (col = 0; col < ncol; col++)
        {
            pcol = pinv[col];
            LADEL_FOR(index, M, col)
            {
                prow = pinv[M->i[index]];
                if (pcol < prow)
                {
                    pindex = col_counts[prow]++;
                    Mpp->i[pindex] = pcol;
                } else 
                {
                    pindex = col_counts[pcol]++;
                    Mpp->i[pindex] = prow;
                }
                if (M->values) Mpp->x[pindex] = M->x[index]; 
            }
        }
    } 
}