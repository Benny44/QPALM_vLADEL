/**
 * @file ladel_postorder.h
 * @author Ben Hermans
 * @brief Routine to compute the postordering of the elimination tree.
 */


#ifndef LADEL_POSTORDER_H
#define LADEL_POSTORDER_H

#include "ladel_types.h"

/**
 * Compute the postordering of the elimination tree.
 * 
 * @note The result is in sym->post.
 * 
 * @param M     Sparse matrix to be analyzed
 * @param sym   Symbolics of the factorization (sym->post contains the output)
 * @param work  LADEL workspace
 * @return      Status
 */
ladel_int ladel_postorder(  ladel_sparse_matrix *M, 
                            ladel_symbolics     *sym, 
                            ladel_work          *work);

#endif /*LADEL_POSTORDER_H*/