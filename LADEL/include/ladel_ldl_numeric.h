/**
 * @file ladel_ldl_numeric.h
 * @author Ben Hermans
 * @brief The numerical part of the factorization.
 */

#ifndef LDL_NUMERIC_H
#define LDL_NUMERIC_H

#include "ladel_types.h"

/**
 * Numerical part of the factorization of @f$M + \alpha \begin{bmatrix}I_{n} & \\ & 0\end{bmatrix}@f$.
 * 
 * @param Mpp   The matrix M or its symmetric permutation
 * @param d     Diagonal parameters @f$\alpha@f$ and @f$n@f$
 * @param sym   Symbolics of the factorization
 * @param LD    Output factor struct
 * @param work  LADEL workspace
 * @return      Status
 */
ladel_int ladel_ldl_numeric_with_diag(  ladel_sparse_matrix *Mpp, 
                                        ladel_diag          d, 
                                        ladel_symbolics     *sym, 
                                        ladel_factor        *LD, 
                                        ladel_work          *work);

/**
 * Numerical part of the factorization of @f$M@f$.
 * 
 * @param Mpp   The matrix M or its symmetric permutation
 * @param sym   Symbolics of the factorization
 * @param LD    Output factor struct
 * @param work  LADEL workspace
 * @return      Status
 */
ladel_int ladel_ldl_numeric(ladel_sparse_matrix *Mpp, 
                            ladel_symbolics     *sym, 
                            ladel_factor        *LD, 
                            ladel_work          *work);

#endif /*LDL_NUMERIC_H*/