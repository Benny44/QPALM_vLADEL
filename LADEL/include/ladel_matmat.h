/**
 * @file ladel_matmat.h
 * @author Ben Hermans
 * @brief Routines to compute matrix matrix products. For now only @f$MM^T@f$ and @f$M\Sigma M^T@f$, 
 * with @f$\Sigma@f$ a diagonal matrix, are supported.
 */

#ifndef LADEL_MATMAT_H
#define LADEL_MATMAT_H

#include "ladel_types.h"

/**
 * Computes @f$MM^T@f$.
 * 
 * @param M             Matrix
 * @param M_transpose   The transpose of @a M
 * @param work          LADEL workspace
 * @return              @f$MM^T@f$
 */
ladel_sparse_matrix *ladel_mat_mat_transpose(   ladel_sparse_matrix *M, 
                                                ladel_sparse_matrix *M_transpose, 
                                                ladel_work          *work);

/**
 * Computes the pattern of @f$MM^T@f$.
 * 
 * @param M             Matrix
 * @param M_transpose   The transpose of @a M
 * @param work          LADEL workspace
 * @return              The patter of @f$MM^T@f$
 */
ladel_sparse_matrix *ladel_mat_mat_transpose_pattern(   ladel_sparse_matrix *M, 
                                                        ladel_sparse_matrix *M_transpose, 
                                                        ladel_work          *work);

/**
 * Computes @f$M\Sigma M^T@f$.
 * 
 * @param M             Matrix
 * @param M_transpose   The transpose of @a M
 * @param diag          Array containing the diagonal of @f$\Sigma@f$
 * @param work          LADEL workspace
 */
ladel_sparse_matrix *ladel_mat_diag_mat_transpose(  ladel_sparse_matrix *M, 
                                                    ladel_sparse_matrix *M_transpose, 
                                                    ladel_double        *diag, 
                                                    ladel_work          *work);

/**
 * Core mat_mat_transpose function with all options (diag and value/pattern).
 * 
 * @note If @a values==FALSE, only the pattern of @f$MM^T@f$. Setting @a diag=NULL, then @f$\Sigma = I@f$
 * is assumed.
 *   
 * @param M             Matrix             
 * @param M_transpose   Transpose of M
 * @param diag          Array containing the diagonal of @f$\Sigma@f$. If NULL, @f$\Sigma = I@f$.
 * @param values        Values or pattern computation
 * @param work          LADEL workspace
 */
ladel_sparse_matrix *ladel_mat_mat_transpose_advanced(  ladel_sparse_matrix *M, 
                                                        ladel_sparse_matrix *M_transpose, 
                                                        ladel_double        *diag, 
                                                        ladel_int           values, 
                                                        ladel_work          *work);

#endif /* LADEL_MATMAT_H */