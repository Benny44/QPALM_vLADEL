/**
 * @file ladel_upper_diag.h
 * @author Ben Hermans
 * @brief Routine to keep only the upper diagonal elements of a (symmetric) matrix.
 */

#ifndef LADEL_UPPER_DIAG
#define LADEL_UPPER_DIAG

#include "ladel_types.h"

/**
 * Return @f$triu(M)@f$, that is the upper diagonal elements of @a M.
 * 
 * @param M Input matrix, contains only upper diagonal elements on output.
 */
void ladel_to_upper_diag(ladel_sparse_matrix *M);

#endif /* LADEL_UPPER_DIAG */