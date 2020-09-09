/**
 * @file ladel_pattern.h
 * @author Ben Hermans
 * @brief Routines to compute the pattern of the result of a backsolve.
 * @details These routines are used to compute the pattern of the (next) row in L during the numerical factorization,
 * or more generally as a first step in performing a backsolve operation with a sparse right hand side.
 */

#ifndef LADEL_PATTERN_H
#define LADEL_PATTERN_H

#include "ladel_types.h"

/**
 * Computes the pattern of the (next) row in L.
 * 
 * @note This pattern is stored in sym->pattern[start] through sym->pattern[m-1], where start is the return
 * value of this function.
 * 
 * @param M     Matrix to be factorized
 * @param sym   Symbolics of the factorization
 * @param row   Row in L of which to compute the pattern
 * @return      Starting index for the pattern
 */
ladel_int ladel_nonzero_pattern_of_row_in_L(    ladel_sparse_matrix *M, 
                                                ladel_symbolics     *sym,
                                                ladel_int           row);

/**
 * Computes a depth-first search of the given pattern through the elimination tree.
 * 
 * This routine computes the pattern of the result of a backsolve (with a sparse right hand side).
 * 
 * @note This pattern is stored in sym->pattern[start] through sym->pattern[m-1], where start is the return
 * value of this function.
 * 
 * @param W             Right hand side vector of the backsolve
 * @param sym           Symbolics struct
 * @param col_in_W      Column in W to consider
 * @param maximum_row   Perform the Gaussian elimination only up until this row (not included)
 * @return              Starting index for the pattern
 */
ladel_int ladel_etree_dfs(  ladel_sparse_matrix *W, 
                            ladel_symbolics     *sym,
                            ladel_int           col_in_W,
                            ladel_int           maximum_row);

#endif /*LADEL_PATTERN_H*/