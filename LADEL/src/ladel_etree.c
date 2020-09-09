#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel_debug_print.h"

ladel_int ladel_etree(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_work* work)
{
    if (!M || !sym || !work) return FAIL;

    ladel_int *etree = sym->etree;
    ladel_int index, row, col, next;
    ladel_int *ancestor = work->array_int_ncol1;

    for (col = 0; col < M->ncol; col++)
    {
        etree[col] = NONE;
        ancestor[col] = NONE;
        LADEL_FOR(index, M, col)
        {
            row = M->i[index];
            for (; row < col; row = next)
            {
                next = ancestor[row];
                ancestor[row] = col;
                if (next == NONE) 
                {
                    etree[row] = col;
                    break;
                }  
            }
        }
    }

    return SUCCESS;
}

#ifdef SIMPLE_COL_COUNTS
ladel_int ladel_etree_and_col_counts(ladel_sparse_matrix *M, ladel_symbolics *sym, ladel_work* work)
{
    if (!M || !sym || !work) return FAIL;

    ladel_int *etree = sym->etree, *col_counts = sym->col_counts;
    
    ladel_int index, row, col, next, ncol = M->ncol;
    ladel_int *touched = work->array_int_ncol1;
    
    for (col = 0; col < ncol; col++) 
    {
        col_counts[col] = 0;
        touched[col] = NONE;
    }

    for (col = 0; col < ncol; col++)
    {
        etree[col] = NONE;
        touched[col] = col;
        LADEL_FOR(index, M, col)
        {
            row = M->i[index];
            for (; row < col && touched[row] != col; row = next)
            {
                col_counts[row]++;
                touched[row] = col;
                next = etree[row];
                if (next == NONE) 
                {
                    etree[row] = col;
                    break;
                }  
            }
        }
    }

    for (col = 1; col < ncol; col++) col_counts[col] += col_counts[col-1];
    return col_counts[ncol-1];
}
#endif