/**
 * @file ladel_add.h
 * @author Ben Hermans
 * @brief Routines to add matrices.
 * @details The add_matrices and add_matrices_pattern routines compute the scaled addition 
 * of the matrices provided or of the patterns of the matrices respectively.  
 */

#ifndef LADEL_ADD_H
#define LADEL_ADD_H

#include "ladel_types.h"

/**
 * Returns a sparse matrix @f$C = \alpha A + \beta B@f$.
 * 
 * @param alpha Scaling for the first matrix
 * @param A     The first matrix
 * @param beta  Scaling for the second matrix
 * @param B     The second matrix
 * @param work  LADEL workspace
 * @return      @f$C = \alpha A + \beta B@f$
 */ 
ladel_sparse_matrix *ladel_add_matrices(ladel_double        alpha, 
                                        ladel_sparse_matrix *A, 
                                        ladel_double        beta, 
                                        ladel_sparse_matrix *B, 
                                        ladel_work          *work);

/**
 * Returns a pattern matrix whose pattern includes the patterns of @a A and @a B.
 * 
 * @param A     The first matrix
 * @param B     The second matrix
 * @param work  LADEL workspace
 * @return      A sparse matrix whose pattern contains that of @a A and @a B
 */
ladel_sparse_matrix *ladel_add_matrices_pattern(ladel_sparse_matrix* A, 
                                                ladel_sparse_matrix *B, 
                                                ladel_work          *work);

/**
 * Returns a sparse matrix @f$C = \alpha A + \beta B@f$ if @a values==TRUE, or a pattern matrix
 * that includes the patterns of @a A and @a B if @a values==FALSE.
 * 
 * @param alpha     Scaling for the first matrix
 * @param A         The first matrix
 * @param beta      Scaling for the second matrix
 * @param B         The second matrix
 * @param values    Flag to choose between values and just pattern
 * @param work      LADEL workspace
 * @return          @f$C = \alpha A + \beta B@f$ or the sum of the patterns
 */ 
ladel_sparse_matrix *ladel_add_matrices_advanced(   ladel_double        alpha, 
                                                    ladel_sparse_matrix *A, 
                                                    ladel_double        beta, 
                                                    ladel_sparse_matrix *B, 
                                                    ladel_int           values, 
                                                    ladel_work          *work);

#endif /* LADEL_ADD_H */