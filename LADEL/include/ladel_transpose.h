/**
 * @file ladel_transpose.h
 * @author Ben Hermans
 * @brief Routine to compute the transpose of a matrix.
 */

#ifndef LADEL_TRANSPOSE_H
#define LADEL_TRANSPOSE_H

#include "ladel_types.h"

/**
 * Returns the transpose of a matrix, that is @f$out = M^T@f$.
 * 
 * If @a M is a pattern matrix, or if @a values==FALSE, a pattern matrix is returned
 * with the pattern of @f$M^T@f$.
 * 
 * @param M         Input matrix
 * @param values    Flag to indicate value or pattern computation
 * @param work      LADEL workspace
 * @return          @f$M^T@f$ (or its pattern)
 */ 
ladel_sparse_matrix *ladel_transpose(   ladel_sparse_matrix *M, 
                                        ladel_int           values, 
                                        ladel_work          *work);

#endif /*LADEL_TRANSPOSE_H*/