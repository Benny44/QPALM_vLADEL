#include "ladel.h"
#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel_pattern.h"
#include "ladel_rank1_mod.h"
#include "ladel_row_mod.h"
#include "ladel_permutation.h"
#include "ladel_debug_print.h"
#include <math.h>

ladel_int ladel_row_add(ladel_factor *LD, ladel_symbolics *sym, ladel_int row_in_L, ladel_sparse_matrix *W, ladel_int col_in_W, ladel_double diag, ladel_work *work)
{
    if (!LD || !sym || !W || !work) return FAIL;
    
    ladel_int start, index_in_pattern, ncol = sym->ncol, row, index, index2, Wnz, status;
    ladel_sparse_matrix* L = LD->L;
    ladel_double *Dinv = LD->Dinv;
    ladel_int *etree = sym->etree;
    ladel_double temp, d22 = diag;
    ladel_double *l12 = work->array_double_all_zeros_ncol1;

    ladel_set *l31_pattern = work->set_preallocated1;
    l31_pattern->size_set = 0;
    ladel_set *set_L31 = work->set_unallocated_values2;
    ladel_set *difference = work->set_preallocated2;
    ladel_int *offset = work->array_int_ncol1;
    ladel_int *insertions = work->array_int_ncol2;
    
    if (W->nz == NULL) Wnz = W->p[col_in_W+1] - W->p[col_in_W];
    else Wnz = W->nz[col_in_W];

    if (LD->pinv != NULL)
    {
        ladel_int_vector_copy(W->i + W->p[col_in_W], Wnz, work->array_int_ncol3);
        ladel_double_vector_copy(W->x + W->p[col_in_W], Wnz, work->array_double_ncol1);
        ladel_permute_sparse_vector(W, col_in_W, LD->pinv, work);
        row_in_L = LD->pinv[row_in_L];
    }
    
    /* 1. Solve lower triangular system L11*D11*l12 = W12 */
    
    
    for (index = W->p[col_in_W]; index < W->p[col_in_W] + Wnz; index++)
    {
        row = W->i[index];
        l12[row] = W->x[index];
        if (row > row_in_L)
        {
            l31_pattern->set[l31_pattern->size_set] = row;
            l31_pattern->size_set++;
        }
    }

    start = ladel_etree_dfs(W, sym, col_in_W, row_in_L);
    for (index_in_pattern = start; index_in_pattern < ncol; index_in_pattern++)
    {
        row = sym->pattern[index_in_pattern];
        temp = l12[row];
        /* 2. d22 = c22 - l12^T D11*l12 */
        d22 -= temp*temp*Dinv[row];
        
        l12[row] *= Dinv[row];
        /* Gaussian elimination */  
        for (index = L->p[row]; index < (L->p[row] + L->nz[row]) && L->i[index] < row_in_L; index++)
        {
            l12[L->i[index]] -= L->x[index]*temp;
        }
        /* 3. l32 = (c32 - L31*D11*l12)/d22 */
        ladel_set_set(set_L31, L->i + index, (L->p[row] + L->nz[row])-index, ncol);
        ladel_set_union(l31_pattern, set_L31, difference, offset, insertions, row_in_L);
        for (index2 = (L->p[row] + L->nz[row] - 1); index2 >= index; index2--)
        {
            l12[L->i[index2]] -= L->x[index2]*temp;
            
            /* Shift the columns down by one to make room for l12*/
            L->i[index2+1] = L->i[index2];
            L->x[index2+1] = L->x[index2]; 
        }
        /* Insert l12[row] */
        L->i[index] = row_in_L;
        L->x[index] = l12[row];
        l12[row] = 0;
        L->nz[row]++;
        if (etree[row] == NONE || etree[row] > row_in_L) etree[row] = row_in_L;
    }

    /* Insert l31 */
    d22 = Dinv[row_in_L] = 1/d22;
    L->nz[row_in_L] = l31_pattern->size_set;
    LADEL_FOR(index, L, row_in_L)
    {
        L->i[index] = row = l31_pattern->set[index-L->p[row_in_L]];
        L->x[index] = d22*l12[row];
        l12[row] = 0;
    }
    if (l31_pattern->size_set > 0) 
    {
        etree[row_in_L] = L->i[L->p[row_in_L]];
    }
    
    /* Reset the diagonal element from W (this was not used and so not reset yet).*/
    l12[row_in_L] = 0;

    /* 4. w = l32*sqrt(abs(d22)) */
    /* 5. Update or downdate L33*D33*L33^T = L33*D33*L33^T - sign(d22)*w*w^T */
    status = ladel_rank1_update(LD, sym, L, row_in_L, 1/sqrt(LADEL_ABS(d22)), d22 < 0, work);


    /* Restore W if permuted*/
    if (LD->pinv != NULL)
    {
        ladel_int_vector_copy(work->array_int_ncol3, Wnz, W->i + W->p[col_in_W]);
        ladel_double_vector_copy(work->array_double_ncol1, Wnz, W->x + W->p[col_in_W]);
    }

    return status;
}

ladel_int ladel_row_del(ladel_factor *LD, ladel_symbolics *sym, ladel_int row_in_L, ladel_work *work)
{
    if (!LD || !sym || !work) return FAIL;
    ladel_int status, index_of_row_in_L, index, found, col,
               top, bottom, top_row, bottom_row, middle, middle_row;
    ladel_double d22_old;
    ladel_sparse_matrix *L = LD->L;
    ladel_int *etree = sym->etree;

    if (LD->pinv != NULL) row_in_L = LD->pinv[row_in_L];

    /* 1. Delete row_in_L l12 */
    for (col = 0; col < row_in_L; col++)
    {
        if (L->nz[col] == 0) continue;
        /* Perform binary search between top and bottom of col to look for row row_in_L */
        top = L->p[col];
        top_row = L->i[top];
        bottom = L->p[col] + L->nz[col] - 1;
        bottom_row = L->i[bottom];

        found = FALSE;
        if (top_row == row_in_L)
        {
            index_of_row_in_L = top;
            found = TRUE;
        } else if (bottom_row == row_in_L)
        {
            index_of_row_in_L = bottom;
            found = TRUE;
        } else if (bottom_row < row_in_L || top_row > row_in_L)
        {
            continue;
        } else
        {
            while (top < bottom)
            {
                middle = (top + bottom)/2;
                middle_row = L->i[middle];
                if (middle_row == row_in_L)
                {
                    index_of_row_in_L = middle;
                    found = TRUE;
                    break;
                } else if (middle_row > row_in_L)
                {
                    bottom = middle - 1;
                } else
                {
                    top = middle + 1;
                }
            }

            middle = (top + bottom)/2;
            middle_row = L->i[middle];
            if (middle_row == row_in_L)
            {
                index_of_row_in_L = middle;
                found = TRUE;
            }
        }
        if (found == TRUE)
        {
            /* Delete the entry at index_of_row_in_L and move the rest of the column one position up. */
            for (index = index_of_row_in_L; index < L->p[col] + L->nz[col] - 1; index++)
            {
                L->i[index] = L->i[index+1];
                L->x[index] = L->x[index+1];
            }
            L->nz[col]--;

            /* Adjust the etree accordingly. */
            if (etree[col] == row_in_L)
            {
                if (index_of_row_in_L < L->p[col] + L->nz[col]) etree[col] = L->i[index_of_row_in_L];
                else etree[col] = NONE;
            }     
        } 
    }
    
    /* 2. d22_new = 1 */
    d22_old = 1.0/LD->Dinv[row_in_L];
    LD->Dinv[row_in_L] = 1;

    /* 4. w = l32_old*sqrt(d22_old) */
    /* 5. Update or downdate L33*D33*L33^T = L33*D33*L33^T + sign(d22_old)*w*w^T */
    status = ladel_rank1_update(LD, sym, L, row_in_L, sqrt(LADEL_ABS(d22_old)), d22_old > 0, work);
    /* 3. Delete column l32 */
    L->nz[row_in_L] = 0;
    etree[row_in_L] = NONE;

    return status;
}