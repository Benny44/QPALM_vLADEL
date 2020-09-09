/**
 * @file ladel.h
 * @author Ben Hermans
 * @brief LADEL main include file.
 * @details The user can include this file in his projects to call any LADEL routine. 
 * This file contains the declaration of the factorization, backsolve and factorization
 * update routines provided by LADEL. 
 * @note The implementation of the rank1_update and row_add/del routines are respectively 
 * in ladel_rank1_mod.c and ladel_row_mod.c instead of in ladel.c.
 */

#ifndef LADEL_H
#define LADEL_H

#include "ladel_col_counts.h"
#include "ladel_constants.h"
#include "ladel_copy.h"
#include "ladel_debug_print.h"
#include "ladel_etree.h"
#include "ladel_global.h"
#include "ladel_ldl_numeric.h"
#include "ladel_ldl_symbolic.h"
#include "ladel_matvec.h"
#include "ladel_matmat.h"
#include "ladel_pattern.h"
#include "ladel_permutation.h"
#include "ladel_postorder.h"
#include "ladel_rank1_mod.h"
#include "ladel_row_mod.h"
#include "ladel_scale.h"
#include "ladel_transpose.h"
#include "ladel_types.h"
#include "ladel_upper_diag.h"
#include "ladel_add.h"
#include "ladel_submatrix.h"

/**
 * @name Main solver API
 * @{
 */

/**
 * Computes the @f$LDL^T@f$ factorization of @f$M@f$.
 * 
 * @param M                 Matrix to be factorized
 * @param sym               Symbolics of the factorization
 * @param ordering_method   Indicator to choose the ordering method: @a AMD, @a NO_ORDERING or @a GIVEN_ORDERING
 * @param LD                Output struct containing the @f$LDL^T@f$ factors
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_factorize(  ladel_sparse_matrix *M, 
                            ladel_symbolics     *sym, 
                            ladel_int           ordering_method, 
                            ladel_factor        **LD, 
                            ladel_work          *work);

/**
 * Computes the @f$LDL^T@f$ factorization of @f$M + \alpha \begin{bmatrix}I_{n} & \\ & 0\end{bmatrix}@f$.
 * 
 * @param M                 Matrix to be factorized
 * @param d                 Diagonal parameters @f$\alpha@f$ and @f$n@f$
 * @param sym               Symbolics of the factorization
 * @param ordering_method   Indicator to choose the ordering method: @a AMD, @a NO_ORDERING or @a GIVEN_ORDERING
 * @param LD                Output struct containing the @f$LDL^T@f$ factors
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_factorize_with_diag(ladel_sparse_matrix *M, 
                                    ladel_diag          d, 
                                    ladel_symbolics     *sym, 
                                    ladel_int           ordering_method, 
                                    ladel_factor        **LD, 
                                    ladel_work          *work);

/**
 * Computes the @f$LDL^T@f$ factorization of @f$M,@f$ but allocates based on @a Mbasis.
 * 
 * This routine is used for subsequent update routines. LADEL assumes the user knows in advance
 * the possible updates (at least the pattern), so that it can allocate the factor once. Dynamic 
 * reallocation during update routines is currently not supported.
 * 
 * @param M                 Matrix to be factorized
 * @param sym               Symbolics of the factorization
 * @param ordering_method   Indicator to choose the ordering method: @a AMD, @a NO_ORDERING or @a GIVEN_ORDERING
 * @param LD                Output struct containing the @f$LDL^T@f$ factors
 * @param Mbasis            Matrix that is used in the symbolic part of the factorization to allocate LD
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_factorize_advanced( ladel_sparse_matrix *M, 
                                    ladel_symbolics     *sym, 
                                    ladel_int           ordering_method, 
                                    ladel_factor        **LD, 
                                    ladel_sparse_matrix *Mbasis, 
                                    ladel_work          *work);

/**
 * Computes the @f$LDL^T@f$ factorization of @f$M+ \alpha \begin{bmatrix}I_{n} & \\ & 0\end{bmatrix},@f$ but allocates based on @a Mbasis.
 * 
 * This routine is used for subsequent update routines. LADEL assumes the user knows in advance
 * the possible updates (at least the pattern), so that it can allocate the factor once. Dynamic 
 * reallocation during update routines is currently not supported.
 * 
 * @param M                 Matrix to be factorized
 * @param d                 Diagonal parameters @f$\alpha@f$ and @f$n@f$
 * @param sym               Symbolics of the factorization
 * @param ordering_method   Indicator to choose the ordering method: @a AMD, @a NO_ORDERING or @a GIVEN_ORDERING
 * @param LD                Output struct containing the @f$LDL^T@f$ factors
 * @param Mbasis            Matrix that is used in the symbolic part of the factorization to allocate LD
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_factorize_advanced_with_diag(   ladel_sparse_matrix *M, 
                                                ladel_diag          d, 
                                                ladel_symbolics     *sym, 
                                                ladel_int           ordering_method, 
                                                ladel_factor        **LD, 
                                                ladel_sparse_matrix *Mbasis, 
                                                ladel_work          *work);

/**
 * Computes the @f$LDL^T@f$ factorization of @f$M,@f$ provided @a LD was allocated before.
 * 
 * @param M                 Matrix to be factorized
 * @param sym               Symbolics of the factorization
 * @param LD                Output struct containing the @f$LDL^T@f$ factors
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_factorize_with_prior_basis( ladel_sparse_matrix *M, 
                                            ladel_symbolics     *sym, 
                                            ladel_factor        *LD, 
                                            ladel_work          *work);

/**
 * Computes the @f$LDL^T@f$ factorization of @f$M+ \alpha \begin{bmatrix}I_{n} & \\ & 0\end{bmatrix},@f$ 
 * provided @a LD was allocated before.
 * 
 * @param M                 Matrix to be factorized
 * @param d                 Diagonal parameters @f$\alpha@f$ and @f$n@f$
 * @param sym               Symbolics of the factorization
 * @param LD                Output struct containing the @f$LDL^T@f$ factors
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_factorize_with_prior_basis_with_diag(   ladel_sparse_matrix *M, 
                                                        ladel_diag          d, 
                                                        ladel_symbolics     *sym, 
                                                        ladel_factor        *LD, 
                                                        ladel_work          *work);

/**
 * Computes @f$y = LDL^T \backslash rhs@f$.
 * 
 * @param LD                Factors of an @f$LDL^T@f$ factorization
 * @param rhs               Dense right-hand side 
 * @param y                 Output vector
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_dense_solve(const ladel_factor  *LD, 
                            const ladel_double  *rhs, 
                            ladel_double        *y, 
                            ladel_work          *work);

/**
 * Updates an @f$LDL^T@f$ factorization.
 * 
 * If LD contains the factors of @f$M@f$ on input, LD will contain the factors of @f$M + up\textunderscore or\textunderscore down*factor^2*ww^T,@f$
 * where @f$w = W(:, col\textunderscore in\textunderscore W).@f$
 * 
 * @param LD                Factors of an @f$LDL^T@f$ factorization
 * @param sym               Associated symbolic information
 * @param W                 Sparse matrix containing the vector for the update
 * @param col_in_W          Column of W that is used for the update
 * @param factor            Scaling factor for the update vector
 * @param up_or_down        Flag indicating @a UPDATE or @a DOWNDATE 
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_rank1_update(   ladel_factor        *LD, 
                                ladel_symbolics     *sym, 
                                ladel_sparse_matrix *W, 
                                ladel_int           col_in_W, 
                                ladel_double        factor, 
                                ladel_int           up_or_down, 
                                ladel_work          *work);

/**
 * Updates an @f$LDL^T@f$ factorization.
 * 
 * If LD contains the factors of 
 * @f$M = \begin{bmatrix}
 * M_{11}   & 0 & M_{13} \\
 * 0        & 1 & 0      \\
 * M_{13}^T & 0 & M_{33}
 * \end{bmatrix} @f$ on input, LD will contain the factors of 
 * @f$\begin{bmatrix}
 * M_{11}   & m_{12} & M_{13}   \\
 * m_{12}^T & m_{22} & m_{32}^T \\
 * M_{13}^T & m_{32} & M_{33}
 * \end{bmatrix} @f$ on output, where @f$m = W(:, col\textunderscore in\textunderscore W),@f$ and 
 * @f$row\textunderscore in\textunderscore L = n + 1@f$ with @f$[n,n] = size(M_{11}).@f$
 * 
 * @param LD                Factors of an @f$LDL^T@f$ factorization
 * @param sym               Associated symbolic information
 * @param row_in_L          Position of the added row/column
 * @param W                 Sparse matrix containing the vector for the update
 * @param col_in_W          Column of W that is used for the update
 * @param diag              Diagonal element @f$m_{22}@f$
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_row_add(ladel_factor        *LD, 
                        ladel_symbolics     *sym, 
                        ladel_int           row_in_L, 
                        ladel_sparse_matrix *W, 
                        ladel_int           col_in_W, 
                        ladel_double        diag, 
                        ladel_work          *work);

/**
 * Updates an @f$LDL^T@f$ factorization.
 * 
 * If LD contains the factors of 
 * @f$\begin{bmatrix}
 * M_{11}   & m_{12} & M_{13}   \\
 * m_{12}^T & m_{22} & m_{32}^T \\
 * M_{13}^T & m_{32} & M_{33}
 * \end{bmatrix} @f$ on input, LD will contain the factors of
 * @f$M = \begin{bmatrix}
 * M_{11}   & 0 & M_{13} \\
 * 0        & 1 & 0      \\
 * M_{13}^T & 0 & M_{33}
 * \end{bmatrix} @f$ on output, where 
 * @f$row\textunderscore in\textunderscore L = n + 1@f$ with @f$[n,n] = size(M_{11}).@f$
 * 
 * @param LD                Factors of an @f$LDL^T@f$ factorization
 * @param sym               Associated symbolic information
 * @param row_in_L          Position of the deleted row/column
 * @param work              LADEL workspace
 * @return                  Status
 */
ladel_int ladel_row_del(ladel_factor    *LD, 
                        ladel_symbolics *sym, 
                        ladel_int       row_in_L, 
                        ladel_work      *work);

/**
 * @}
 */

#endif /*LADEL_H*/