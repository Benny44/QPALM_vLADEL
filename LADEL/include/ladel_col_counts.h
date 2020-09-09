/**
 * @file ladel_col_counts.h
 * @author Ben Hermans
 * @brief Computes the col counts needed for the symbolic factorization (after etree and postorder).
 * @details The routine in this file computes the column counts of the factor assuming the etree and
 * postordering have been performed prior to this. An alternative (with worse asymptotic complexity) 
 * that computes the etree and column counts in one go, is @a ladel_etree_and_col_counts in ladel_etree.c.  
 */

#ifndef LADEL_COL_COUNTS_H
#define LADEL_COL_COUNTS_H

#include "ladel_types.h"

/**
 * Computes the column counts of the factor.
 * 
 * This routine should only be called after the etree and postorder routines.
 * 
 * @param M     Sparse matrix to be analyzed
 * @param sym   Struct holding symbolic information
 * @param work  LADEL workspace
 * @return      Status
 */
ladel_int ladel_col_counts( ladel_sparse_matrix *M, 
                            ladel_symbolics     *sym, 
                            ladel_work          *work);

#endif /*LADEL_COL_COUNTS_H*/