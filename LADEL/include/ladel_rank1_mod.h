/**
 * @file ladel_rank1_mod.h
 * @author Ben Hermans
 * @brief Helper routines for the rank1_update function.
 * @details The ladel_rank1_update function implements rank 1 updates for a sparse factorization, and is implemented
 * in ladel_rank1_mod.c.  
 */

#ifndef LADEL_RANK1_MOD_H
#define LADEL_RANK1_MOD_H

#include "ladel_types.h"

/**
 * Adds the nonzero pattern in @a set to the pattern of the @a col of @a L.
 * 
 * @param L             Sparse factor
 * @param col           Column index
 * @param col_set       Pointer to set struct that can hold the pattern of a column of L
 * @param set           Nonzero pattern to be added
 * @param difference    Output elements in set that did not occur in the column of L
 * @param offset        Output list of offset indices of original elements of the column of L
 * @param insertions    Output list of indices where new elements have been added to the column of L
 * @return              Resulting status (@a SET_HAS_CHANGED, @a SET_HAS_NOT_CHANGED, or @a MAX_SET_SIZE_EXCEEDED)
 */
ladel_int ladel_add_nonzero_pattern_to_col_of_L(ladel_sparse_matrix *L, 
                                                ladel_int           col, 
                                                ladel_set           *col_set, 
                                                ladel_set           *set, 
                                                ladel_set           *difference, 
                                                ladel_int           *offset, 
                                                ladel_int           *insertions);

/**
 * Computes the set union of two index sets.
 * 
 * The result is stored in @a first_set. This function also returns some information
 * on the set union, that is:
 * - The elements in the @a second_set which were not already in the @a first_set are stored in @a difference
 * - The difference of the indices of the original elements of @a first_set versus now is stored in @a offset
 * - The indices of the set union @a first_set where new elements have been added from @a second_set are stored in @a insertions 
 * @param first_set     First input set, on output the set union
 * @param second_set    Second input set
 * @param difference    Output elements in @a second_set that were not originally in @a first_set
 * @param offset        Output list of offset indices of original elements of the @a first_set
 * @param insertions    Output list of indices where new elements have been added to @a first_set from @a second_set
 * @param threshold     Only elements from the @a second_set with value > this threshold are considered   
 * @return              Resulting status (@a SET_HAS_CHANGED, @a SET_HAS_NOT_CHANGED, or @a MAX_SET_SIZE_EXCEEDED)
 */
ladel_int ladel_set_union(  ladel_set *first_set, 
                            ladel_set *second_set, 
                            ladel_set *difference, 
                            ladel_int *offset, 
                            ladel_int *insertions, 
                            ladel_int threshold);

#endif