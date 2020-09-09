/**
 * @file ladel_permutation.h
 * @author Ben Hermans
 * @brief Routines to permute vectors and matrices.
 */

#ifndef LADEL_PERMUTATION_H
#define LADEL_PERMUTATION_H

#include "ladel_types.h"

/**
 * Compute y, where @f$y[i] = x[p[i]]@f$.
 * 
 * @param x     Input vector
 * @param p     Permutation vector
 * @param size  Size of the vectors
 * @param y     Output vector
 */
void ladel_permute_vector(  ladel_double    *x, 
                            ladel_int       *p, 
                            ladel_int       size, 
                            ladel_double    *y);

/**
 * Compute y, where @f$y[pinv[i]] = x[i]@f$.
 * 
 * @param x     Input vector
 * @param pinv  (Inverse) permutation vector
 * @param size  Size of the vectors
 * @param y     Output vector
 */
void ladel_inverse_permute_vector(  ladel_double    *x, 
                                    ladel_int       *pinv, 
                                    ladel_int       size, 
                                    ladel_double    *y);

/**
 * Invert a permutation vector, such that @f$pinv[p[i]] = i@f$.
 * 
 * @param p     Permutation vector
 * @param pinv  Output inverse permutation vector
 * @param size  Size of the vectors
 */
void ladel_invert_permutation_vector(   ladel_int *p, 
                                        ladel_int *pinv, 
                                        ladel_int size);

/**
 * Compute @f$Mpp = M(p,p)@f$.
 * 
 * @param M     Input matrix
 * @param p     Permutation vector
 * @param Mpp   Output permuted matrix
 * @param work  LADEL workspace
 */
void ladel_permute_symmetric_matrix(ladel_sparse_matrix *M, 
                                    ladel_int           *p, 
                                    ladel_sparse_matrix *Mpp, 
                                    ladel_work          *work);

/**
 * Permute @a x(:,col) in place, such that on output @f$x(:,col) = x(p,col)@f$.
 * 
 * @param x     Input sparse vector
 * @param col   The column of x to permute
 * @param p     Permutation vector
 * @param work  LADEL workspace
 */
void ladel_permute_sparse_vector(   ladel_sparse_matrix *x, 
                                    ladel_int           col, 
                                    ladel_int           *p, 
                                    ladel_work          *work);

#endif /*LADEL_PERMUTATION_H*/