/**
 * @file ladel_etree.h
 * @author Ben Hermans
 * @brief Routines to compute the elimination tree of a matrix.
 * @details Computing the elimination tree is the first part of the symbolic factorization. 
 * One routine in this file simply computes the etree, the other (compile with -DSIMPLE_COL_COUNTS) 
 * computes the etree and column counts in parallel (but has worse time complexity).
 */

#ifndef LADEL_ETREE_H
#define LADEL_ETREE_H

#include "ladel_types.h"
/**
 * Computes the elimination tree of a matrix.
 * 
 * This tree is stored in sym->etree.
 * 
 * @param M     Matrix
 * @param sym   Symbolics struct for the factorization
 * @param work  LADEL workspace
 * @return      Status
 */
ladel_int ladel_etree(  ladel_sparse_matrix *M, 
                        ladel_symbolics     *sym, 
                        ladel_work          *work);

#ifdef SIMPLE_COL_COUNTS
/**
 * Computes the elimination tree and column counts of a matrix.
 * 
 * This function can be used to do the whole symbolic analysis at once, but has 
 * a worse asymptotic time complexity than using etree, postorder and col_counts
 * sequentially.
 * 
 * @param M     Matrix
 * @param sym   Symbolics struct for the factorization
 * @param work  LADEL workspace
 * @return      Status
 */
ladel_int ladel_etree_and_col_counts(   ladel_sparse_matrix *M, 
                                        ladel_symbolics     *sym, 
                                        ladel_work          *work);
#endif

#endif /*LADEL_ETREE_H*/