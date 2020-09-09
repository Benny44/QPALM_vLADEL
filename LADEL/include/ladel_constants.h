/**
 * @file ladel_constants.h
 * @author Ben Hermans
 * @brief Constants and macros used in LADEL.
 */

#ifndef LADEL_CONSTANTS_H
#define LADEL_CONSTANTS_H

/**
 * @name Constants
 * @{
 */


#define TRUE 1 /**< For booleans */
#define FALSE 0 /**< For booleans */

#define SUCCESS 1 /**< For status returns */
#define FAIL -1 /**< For status returns */

#define NONE -1 /**< Indicates a root (a node without parent), or an untouched node */

#define UNSYMMETRIC 0 /**< No symmetry is assumed in the matrix */
#define UPPER 1 /**< Use only upper part of matrix */
#define LOWER -1 /**< Use only lower part of matrix */

#define AMD 1 /**< Ordering method during the symbolic part of the factorization */
#define NO_ORDERING 0 /**< No ordering is performed during the symbolic part of the factorization */
#define GIVEN_ORDERING 2 /**< The ordering was computed previously and is already stored in sym->p */

#define MARKED 1 /**< Indicate whether a node is marked */
#define UNMARKED 0 /**< Indicate whether a node is not marked */

#define SET_HAS_CHANGED TRUE /**< Possible return value of ladel_set_union indicating the set has changed */
#define SET_HAS_NOT_CHANGED FALSE /**< Possible return value of ladel_set_union indicating the set has not changed */
#define MAX_SET_SIZE_EXCEEDED NONE /**< Possible return value of ladel_set_union indicating the set has grown beyond the maximum size */

#define UPDATE TRUE /**< Flag in rank1_update to perform an update */
#define DOWNDATE FALSE /**< Flag in rank1_update to perform a downdate */

/**
 * @}
 */

/**
 * @name Macros
 * @{
 */

#define LADEL_MAX(a, b) ((a) > (b) ? (a) : (b)) /**< Return the maximum of two numbers */
#define LADEL_MIN(a, b) ((a) > (b) ? (b) : (a)) /**< Return the minimum of two numbers */
#define LADEL_ABS(a) ((a) < 0 ? -(a) : (a))     /**< Return the absolute value a number */

#define LADEL_FOR(index, M, col) for(index = (M)->p[(col)]; index < (((M)->nz) ? (M)->p[(col)] + (M)->nz[(col)] : (M)->p[(col)+1]); index++) /**< Loop through a column of a sparse matrix */

#define IS_ROOT(col, etree) ((etree)[(col)] == NONE) /**< Check whether a node is a root (i.e. has no parent) */

#define MARK(nodes, k) (nodes[(k)] = MARKED) /**< Mark the k-th node */
#define UNMARK(nodes, k) (nodes[(k)] = UNMARKED) /**< Unmark the k-th node */
#define IS_MARKED(nodes, k) (nodes[(k)] == MARKED) /**< Check whether the k-th node is marked */

/**
 * @}
 */

#endif /*LADEL_CONSTANTS_H*/