/**
 * @file ladel_scale.h
 * @author Ben Hermans
 * @brief Routines to scale the columns and rows of a matrix, or the whole matrix itself.
 * @details This also includes routines to compute the row-wise and column-wise infinity norms. 
 */

#ifndef LADEL_SCALE_H
#define LADEL_SCALE_H

#include "ladel_types.h"

/**
 * Scales the columns of @a M with the scalars in @a S, that is @f$M(:,i) = S(i)*M(:,i)@f$.
 * 
 * @param M Matrix
 * @param S Array with scale factors for the columns
 */
void ladel_scale_columns(   ladel_sparse_matrix *M, 
                            ladel_double        *S);

/**
 * Scales the rows of @a M with the scalars in @a S, that is @f$M(i,:) = S(i)*M(i,:)@f$.
 * 
 * @param M Matrix
 * @param S Array with scale factors for the rows
 */
void ladel_scale_rows(  ladel_sparse_matrix *M, 
                        ladel_double        *S);

/**
 * Scales the elements of @a M with @a s, that is @f$M = sM@f$.
 * 
 * @param M Matrix
 * @param s Scaling factor
 */
void ladel_scale_scalar(ladel_sparse_matrix *M, 
                        ladel_double        s);

/**
 * Computes the column-wise infinity norms, that is @f$norms(i) = norm(M(:,i), inf)@f$.
 *
 * @param M     Matrix
 * @param norms Output vector of infinity norms
 */ 
void ladel_infinity_norm_columns(   ladel_sparse_matrix *M, 
                                    ladel_double        *norms);

/**
 * Computes the row-wise infinity norms, that is @f$norms(i) = norm(M(i,:), inf)@f$.
 *
 * @param M     Matrix
 * @param norms Output vector of infinity norms
 */ 
void ladel_infinity_norm_rows(  ladel_sparse_matrix *M, 
                                ladel_double        *norms);


#endif /*LADEL_SCALE_H*/