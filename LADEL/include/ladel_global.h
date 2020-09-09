/**
 * @file ladel_global.h
 * @author Ben Hermans
 * @brief Memory allocation routines.
 */

#ifndef LADEL_GLOBAL_H
#define LADEL_GLOBAL_H

#include "ladel_types.h"
#include "ladel_constants.h"
#include "ladel_copy.h"
#include <stdlib.h>

/**
 * Version of malloc (for mex or for regular C).
 * 
 * If the malloc fails, this function will return NULL.
 * 
 * @param n     Number of blocks
 * @param size  Size of block
 * @return      Pointer to the allocated memory
 */
void *ladel_malloc( ladel_int   n, 
                    size_t      size);

/**
 * Version of calloc (for mex or for regular C).
 * 
 * If the calloc fails, this function will return NULL.
 * 
 * @param n     Number of blocks
 * @param size  Size of block
 * @return      Pointer to the allocated memory
 */
void *ladel_calloc( ladel_int   n, 
                    size_t      size);

/**
 * Version of free (for mex or for regular C).
 * 
 * @param p  Pointer to the memory to be freed
 * @return   NULL
 */
void *ladel_free(void* p);

/**
 * Version of realloc (for mex or for regular C).
 * 
 * If the realloc fails, this function will return the original pointer.
 * 
 * @param p         Pointer to the memory
 * @param n         Number of blocks
 * @param size      Size of block
 * @param status    Status to indicate success
 * @return          Pointer to the reallocated memory (or the same in case )
 */
void *ladel_realloc(void        *p, 
                    ladel_int   n, 
                    size_t      size, 
                    ladel_int   *status);

#ifdef MATLAB
#include "mex.h"
#define ladel_print mexPrintf
#else
#define ladel_print printf /**< Print function (for mex or for regular C) */
#endif

/**
 * Free a sparse matrix (and return NULL).
 * 
 * @param M     Matrix to be freed
 * @return      NULL
 */
ladel_sparse_matrix *ladel_sparse_free(ladel_sparse_matrix *M);

/**
 * Allocate a sparse matrix.
 * 
 * @param nrow      Number of rows
 * @param ncol      Number of columns
 * @param nzmax     Maximum number of nonzeros
 * @param symmetry  UNSYMMETRIC, or store only UPPER part of a matrix
 * @param values    Indicate value or pattern matrix
 * @param nz        Indicate the use of the nz field to list nonzeros per column
 * @return          Allocated matrix
 */
ladel_sparse_matrix *ladel_sparse_alloc(ladel_int nrow, 
                                        ladel_int ncol,
                                        ladel_int nzmax, 
                                        ladel_int symmetry,
                                        ladel_int values, 
                                        ladel_int nz);
/**
 * Allocate a sparse empty matrix (used in special cases).
 * 
 * @param nrow      Number of rows
 * @param ncol      Number of columns
 * @param symmetry  UNSYMMETRIC, or store only UPPER part of a matrix
 * @param values    Indicate value or pattern matrix
 * @param nz        Indicate the use of the nz field to list nonzeros per column
 * @return          Allocated empty matrix
 */
ladel_sparse_matrix *ladel_sparse_alloc_empty(  ladel_int nrow, 
                                                ladel_int ncol, 
                                                ladel_int symmetry, 
                                                ladel_int values, 
                                                ladel_int nz);

/**
 * Reallocate a sparse matrix with a new size.
 * 
 * @note If nzmax <= 0, the matrix is shrinked to fit all the current elements.
 * 
 * @param M         Sparse matrix
 * @param nzmax     New maximum number of nonzeros 
 * @return          Status 
 */
ladel_int ladel_sparse_realloc( ladel_sparse_matrix *M, 
                                ladel_int           nzmax);

/**
 * Free a symbolics struct (and return NULL).
 *  
 * @param sym   Symbolics struct 
 * @return      NULL
 */
ladel_symbolics *ladel_symbolics_free(ladel_symbolics *sym);

/**
 * Allocate a symbolics struct.
 *   
 * @param ncol  Number of columns of matrix to be analyzed 
 * @return      Allocated symbolics struct
 */
ladel_symbolics *ladel_symbolics_alloc(ladel_int ncol);

/**
 * Free a factor.
 *  
 * @param LD    Factors of an @f$LDL^T@f$ factorization
 * @return      NULL
 */
ladel_factor *ladel_factor_free(ladel_factor *LD);

/**
 * Allocate a factors struct.
 *  
 * @param sym   Symbolics struct 
 * @return      Factors struct
 */
ladel_factor *ladel_factor_allocate(ladel_symbolics *sym);

/**
 * Free a set.
 *  
 * @param set   Set to be freed
 * @return      NULL
 */
ladel_set *ladel_set_free(ladel_set *set);

/**
 * Allocate a set struct.
 *   
 * @param max_size  Maximum size of the set 
 * @return          Allocated set
 */
ladel_set *ladel_set_allocate(ladel_int max_size);

/**
 * Fill in the fields of the given set.
 *   
 * @param set           Set to be filled in 
 * @param set_vals      Array that represents the set
 * @param size_set      Current size of the set
 * @param max_size_set  Maximum size of the set
 */
void ladel_set_set( ladel_set *set, 
                    ladel_int *set_vals, 
                    ladel_int size_set, 
                    ladel_int max_size_set);

/**
 * Free a LADEL workspace.
 *  
 * @param work  LADEL workspace to be freed
 * @return      NULL
 */
ladel_work *ladel_workspace_free(ladel_work* work);

/**
 * Allocate a LADEL workspace.
 *  
 * @param ncol  Size of the structures in the workspace (= ncol of matrix to be factorized)
 * @return      Allocated LADEL workspace
 */
ladel_work *ladel_workspace_allocate(ladel_int ncol);


#endif /*LADEL_GLOBAL_H*/