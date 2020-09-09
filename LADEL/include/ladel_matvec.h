/**
 * @file ladel_matvec.h
 * @author Ben Hermans
 * @brief Routines to compute matrix vector products. 
 */

#ifndef LADEL_MATVEC_H
#define LADEL_MATVEC_H

#include "ladel_global.h"
#include "ladel_types.h"

/**
 * Computes @f$y = Mx@f$ or @f$y += Mx@f$ (for @a reset==TRUE of @a FALSE respectively).
 * 
 * @param M     Matrix
 * @param x     Vector
 * @param y     Resulting matrix vector product
 * @param reset If @a TRUE, @a y is first set to zero before adding @f$Mx@f$
 */
void ladel_matvec(  const ladel_sparse_matrix   *M, 
                    const ladel_double          *x, 
                    ladel_double                *y, 
                    ladel_int                   reset);

/**
 * Computes @f$y = M^Tx@f$ or @f$y += M^Tx@f$ (for @a reset==TRUE of @a FALSE respectively).
 * 
 * @param M     Matrix
 * @param x     Vector
 * @param y     Resulting matrix vector product
 * @param reset If @a TRUE, @a y is first set to zero before adding @f$M^Tx@f$
 */
void ladel_tpose_matvec(const ladel_sparse_matrix   *M, 
                        const ladel_double          *x, 
                        ladel_double                *y, 
                        ladel_int                   reset);

/**
 * Computes (for a symmetric matrix M) @f$y = Mx@f$ or @f$y += Mx@f$ (for @a reset==TRUE of @a FALSE respectively).
 * 
 * @param M     Symmetric Matrix (@a M->symmetry==UPPER or @a LOWER)
 * @param x     Vector
 * @param y     Resulting matrix vector product
 * @param reset If @a TRUE, @a y is first set to zero before adding @f$Mx@f$
 */
void ladel_symmetric_matvec(const ladel_sparse_matrix   *M, 
                            const ladel_double          *x, 
                            ladel_double                *y,
                            ladel_int                   reset);

#endif /* LADEL_MATVEC_H */