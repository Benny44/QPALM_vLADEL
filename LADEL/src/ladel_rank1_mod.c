#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel_rank1_mod.h"
#include "ladel_debug_print.h"

ladel_int ladel_add_nonzero_pattern_to_col_of_L(ladel_sparse_matrix *L, ladel_int col, ladel_set *col_set, ladel_set *set, ladel_set *difference, ladel_int* offset, ladel_int* insertions)
{
    ladel_int start = L->p[col], status;
    ladel_set_set(col_set, L->i + start, L->nz[col], L->p[col+1] - L->p[col]); 
    status = ladel_set_union(col_set, set, difference, offset, insertions, col);
    
    /* For now it is assumed the user has allocated enough space. If not, the error is passed on.
    Technically, some reallocation code could go here to include cases where one doesn't know the
    maximum number of nonzeros in each column of L. */
    if (status == MAX_SET_SIZE_EXCEEDED) return MAX_SET_SIZE_EXCEEDED;
    if (status == SET_HAS_NOT_CHANGED) return SET_HAS_NOT_CHANGED;

    /* Move old nz values in L to correct positions and initialize new nz values in L to zero */
    ladel_int index;
    for (index = L->nz[col]-1; index >= 0; index--) 
        L->x[start+index+offset[index]] = L->x[start+index];
        
    /* Initialize new nz values in L to zero */
    for (index = 0; index < difference->size_set; index++)
        L->x[start+insertions[index]] = 0;

    L->nz[col] = col_set->size_set;

    return status;
}

ladel_int ladel_set_union(ladel_set *first_set, ladel_set *second_set, ladel_set *difference, 
                            ladel_int *offset, ladel_int *insertions, ladel_int threshold)
{
    ladel_int *set1 = first_set->set; 
    ladel_int size_set1 = first_set->size_set; 
    ladel_int max_size_set1 = first_set->max_size_set;
    ladel_int *set2 = second_set->set;
    ladel_int size_set2 = second_set->size_set;

    ladel_int *dif = difference->set;
    difference->size_set = 0;

    ladel_int index1 = 0, index2, row1, row2, index_dif = 0;

    /* Special case: second set is empty => do nothing */
    if (size_set2 == 0) return SET_HAS_NOT_CHANGED;

    /* Special case: first set is empty => just copy */
    if (size_set1 == 0) 
    {
        for (index2 = 0; index2 < size_set2; index2++)
        {
            row2 = set2[index2];
            if (row2 <= threshold) continue;
            insertions[index1] = index1;
            set1[index1] = dif[index1] = row2; 
            index1++;
        }
        first_set->size_set = difference->size_set = index1;
        if (index1 == 0) return SET_HAS_NOT_CHANGED;
        else return SET_HAS_CHANGED;
    }

    row1 = set1[0];
/* Construct difference set and offsets -----------------------------------*/ 
    for (index2 = 0; index2 < size_set2; index2++)
    {
        row2 = set2[index2];
        if (row2 <= threshold) continue;

        for (; index1 < first_set->size_set && row1 < row2; index1++) 
        {
            row1 = set1[index1];
            offset[index1] = index_dif;
            if (row1 >= row2) break; 
        }
        if (row1 > row2)/*add elem of set2 in place index1 in set1*/
        {
            dif[index_dif] = row2;
            index_dif++;

            size_set1++;
            if (size_set1 > max_size_set1) return MAX_SET_SIZE_EXCEEDED;
        }
        else if (row1 < row2) /*append the rest of set2 to the end of set1*/
        {
            for (; index2 < size_set2; index2++, size_set1++, index_dif++)
            {
                if (size_set1 == max_size_set1) return MAX_SET_SIZE_EXCEEDED;
                dif[index_dif] = set2[index2];
                insertions[index_dif] = index1+index_dif;
            }
        }
    }
    
    if (index_dif == 0) return SET_HAS_NOT_CHANGED;
    
    for (; index1 < first_set->size_set; index1++) offset[index1] = index_dif;
    difference->size_set = index_dif;

/* Merge difference into the first set ------------------------------------*/

    /* Move original set1 values to the correct positions */
    for (index1 = first_set->size_set-1; index1 >= 0; index1--) set1[index1+offset[index1]] = set1[index1];

    /* Compute the positions for the new elements */
    index_dif = 0;
    for (index1 = 0; index1 < first_set->size_set; index1++)
        for (; index_dif < offset[index1]; index_dif++) 
            insertions[index_dif] = index1 + index_dif;

    /* Insert the new elements */
    for (index_dif = 0; index_dif < difference->size_set; index_dif++) set1[insertions[index_dif]] = dif[index_dif];
    
    first_set->size_set = size_set1;
    return SET_HAS_CHANGED;
}


ladel_int ladel_rank1_update(ladel_factor *LD, ladel_symbolics *sym, ladel_sparse_matrix *W, ladel_int col_in_W, ladel_double factor, ladel_int up_or_down, ladel_work* work)
{
    if (!LD || !sym || !W || !work) return FAIL;

    ladel_int *etree = sym->etree;
    ladel_sparse_matrix *L = LD->L;
    ladel_double *Dinv = LD->Dinv;

    /* TODO: account for the permutation! */
    
    ladel_int col, row, index, index_L, size_W = (W->nz == NULL) ? W->p[col_in_W+1] - W->p[col_in_W] : W->nz[col_in_W] ;
    if (size_W == 0) return SUCCESS; 
    ladel_int changed = SET_HAS_NOT_CHANGED, changed_W, changed_child;
    ladel_double sigma;
    if (up_or_down == UPDATE) sigma = 1.0;
    else if (up_or_down == DOWNDATE) sigma = -1.0;
    else return FAIL;

    ladel_set *set_W = work->set_unallocated_values1;
    ladel_set_set(set_W, W->i + W->p[col_in_W], size_W, size_W);
    ladel_set *set_L = work->set_unallocated_values2;
    ladel_set *difference = work->set_preallocated1;
    difference->size_set = 0;
    ladel_set *difference2 = work->set_preallocated2;
    difference2->size_set = 0;
    ladel_set *difference_child = work->set_preallocated3;
    difference_child->size_set = 0;
    ladel_set *set_child = work->set_unallocated_values3;

    ladel_int *offset = work->array_int_ncol1;
    ladel_int *insertions = work->array_int_ncol2;
    ladel_double *W_col = work->array_double_all_zeros_ncol1;
    
    LADEL_FOR(index, W, col_in_W)
        W_col[W->i[index]] = factor*W->x[index];

    ladel_double alpha = 1, alpha_new, gamma, w, dinv;
    ladel_int child, old_parent;

    LADEL_FOR(index, W, col_in_W)
    {
        col = W->i[index];
        changed = ladel_add_nonzero_pattern_to_col_of_L(L, col, set_L, set_W, difference_child, offset, insertions);
        if (changed == MAX_SET_SIZE_EXCEEDED) return FAIL;
        if (changed == SET_HAS_CHANGED)
        {
            child = col;
            old_parent = etree[col];
            col = etree[col] = L->i[L->p[col]]; /*new_parent, there must be one since the set has changed*/

            /* prepare the difference in the next col due to the child */
            if (col != old_parent)
                ladel_set_set(set_child, L->i+L->p[child], L->nz[child], L->p[child+1] - L->p[child]);
            
            break;
        }
    }

    if (changed == SET_HAS_CHANGED)
    {
        while (TRUE)
        {
            if (col == old_parent)
                changed_child = ladel_add_nonzero_pattern_to_col_of_L(L, col, set_L, difference_child, difference, offset, insertions);             
            else
                changed_child = ladel_add_nonzero_pattern_to_col_of_L(L, col, set_L, set_child, difference, offset, insertions);
                
            changed_W = ladel_add_nonzero_pattern_to_col_of_L(L, col, set_L, set_W, difference_child, offset, insertions);

            if (changed_child == MAX_SET_SIZE_EXCEEDED || changed_W == MAX_SET_SIZE_EXCEEDED) return FAIL;

            child = col;
            old_parent = etree[col];
            if (L->nz[col] == 0) break; /* There is no subdiag entry, so no new parent */
            col = etree[col] = L->i[L->p[col]]; /* new_parent */

            /* prepare the difference in the next col due to the child */
            if (col == old_parent)
                ladel_set_union(difference_child, difference, difference2, offset, insertions, 0);   
            else
                ladel_set_set(set_child, L->i+L->p[child], L->nz[child], L->p[child+1] - L->p[child]);                      
        }
    }

    /* numerical update */
    for (col = W->i[W->p[col_in_W]]; col != NONE; col = etree[col])
    {
        w = W_col[col];
        dinv = Dinv[col];
        alpha_new = alpha + sigma*w*w*dinv;
        gamma = w*dinv/alpha_new; /*if alpha_new == 0 then matrix not full rank */
        Dinv[col] *= alpha/alpha_new;
        alpha = alpha_new;
        for (index_L = L->p[col]; index_L < L->p[col]+L->nz[col]; index_L++)
        {
            row = L->i[index_L];
            W_col[row] -= w*L->x[index_L];
            L->x[index_L] += sigma*gamma*W_col[row];
        }
    }

    /* Return the double workspace to all zeros */
    for (index = W->i[W->p[col_in_W]]; index != NONE; index = etree[index]) W_col[index] = 0;

    return SUCCESS;

}
