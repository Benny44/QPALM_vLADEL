/**
 * @file ladel_copy.h
 * @author Ben Hermans
 * @brief Routines to copy matrices and vectors.
 */

#ifndef LADEL_COPY_H
#define LADEL_COPY_H

#include "ladel_types.h"

/**
 * Copies a matrix into a newly allocated one.
 * 
 * Note that even if there is superfluous space in M (i.e. M->nz indicates the number of nonzeros per column), 
 * the returned matrix will fit the required size exactly, and not make use of the nz array.
 * 
 * @param M     Matrix to be copied
 * @return      Copy of the matrix
 */
ladel_sparse_matrix *ladel_sparse_allocate_and_copy(ladel_sparse_matrix *M);

/** 
 * Copies a matrix (preallocated)
 * 
 * @param M         Matrix to be copied
 * @param M_copy    Resulting copy
 */
void ladel_sparse_copy( ladel_sparse_matrix *M, 
                        ladel_sparse_matrix *M_copy);

/**
 * Copies an integer vector (preallocated)
 * 
 * @param x     Vector to be copied
 * @param size  Size of the vector
 * @param y     Resulting copy
 */
void ladel_int_vector_copy( ladel_int *x, 
                            ladel_int size, 
                            ladel_int *y);

/**
 * Copies a double vector (preallocated)
 * 
 * @param x     Vector to be copied
 * @param size  Size of the vector
 * @param y     Resulting copy
 */
void ladel_double_vector_copy(  ladel_double    *x, 
                                ladel_int       size, 
                                ladel_double    *y);

#endif /*LADEL_COPY_H*/