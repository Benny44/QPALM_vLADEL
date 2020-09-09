/**
 * @file ladel_debug_print.h
 * @author Ben Hermans
 * @brief Routines to print out matrices and vectors.
 * @details The (majority of) routines here print output that can be entered into Matlab for further inspection of the matrix/vector.
 */

#ifndef LADEL_DEBUG_PRINT_H
#define LADEL_DEBUG_PRINT_H

#include "ladel_global.h"
#include "ladel_types.h"

/**
 * Prints output that can be copied into Matlab to retrieve the given matrix.
 * 
 * @param M Matrix to be printed
 */
void ladel_print_sparse_matrix_matlab(ladel_sparse_matrix *M);

/**
 * Prints the contents of all entries in a matrix.
 * 
 * The matrix may have unused entries if the nz field is active. These entries are also printed
 * by this function for inspection.
 * 
 * @param M Matrix to be printed
 */
void ladel_print_sparse_matrix_entries(ladel_sparse_matrix *M);

/**
 * Prints output that can be copied into Matlab to retrieve the given double vector.
 * 
 * @param x     Vector to be printed
 * @param len   Length of the vector to be printed
 */
void ladel_print_dense_vector_matlab(   ladel_double    *x, 
                                        size_t          len);

/**
 * Prints output that can be copied into Matlab to retrieve the given double vector.
 * 
 * @param x     Vector to be printed
 * @param len   Length of the vector to be printed
 */
void ladel_print_dense_int_vector_matlab(   ladel_int   *x, 
                                            size_t      len);

/**
 * Prints output that can be copied into Matlab to retrieve the given factor.
 * 
 * This function basically prints L and D_inv.
 * 
 * @param LD    Factor to be printed
 */
void ladel_print_factor_matlab(ladel_factor *LD);

/**
 * Print the content of an index set.
 * 
 * @param set   Set to be printed
 */
void ladel_print_set(ladel_set *set);

#endif