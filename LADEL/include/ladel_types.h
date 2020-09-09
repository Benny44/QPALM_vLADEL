/**
 * @file ladel_types.h
 * @author Ben Hermans
 * @brief Structures and types used in LADEL routines
 */

#ifndef LADEL_TYPES_H
#define LADEL_TYPES_H

#ifdef DFLOAT
    typedef float ladel_double; /**< Type for floating point numbers (default: double) */
#else
    typedef double ladel_double; /**< Type for floating point numbers (default: double) */
#endif

#ifdef DLONG
    typedef long ladel_int; /**< Type for integer numbers (default: long) */
#else
    typedef int ladel_int; /**< Type for integer numbers (default: long) */
#endif

/**
 * @brief Sparse matrix in compressed column storage.
 */
typedef struct compressed_column_sparse_matrix 
{
    ladel_int nzmax;        /**< number of nonzeros */
    ladel_int nrow;         /**< number of rows */
    ladel_int ncol;         /**< number of columns */

    ladel_int *p;           /**< column pointers (size ncol+1) */
    ladel_int *i;           /**< row pointers (size nzmax) */
    ladel_double *x;        /**< numerical values (size nzmax) */
    ladel_int *nz;          /**< (optional) number of elements in each column (size ncol) */

    ladel_int values;       /**< has numerical values */
    ladel_int symmetry;     /**< type of symmetry (@a UNSYMMETRIC, @a UPPER or @a LOWER) */    

} ladel_sparse_matrix;

/**
 * @brief Structure capturing symbolic information used for and during the factorization.
 */
typedef struct symbolic_cholesky_information
{
    ladel_int ncol;         /**< number of columns in the analyzed matrix */
    ladel_int *etree;       /**< eliminations tree*/
    ladel_int *postorder;   /**< postordiring of the elimination tree */
    ladel_int *col_counts;  /**< column counts, stored as column pointers */
    ladel_int *p;           /**< fill-reducing ordering (possibly AMD) */
    ladel_int *pinv;        /**< inverse permutation vector */
    ladel_int *pattern;     /**< stores the nonzero pattern of a row of L */ 
    ladel_int *nodes;       /**< keeps track of which nodes have been marked */
} ladel_symbolics;

/**
 * @brief Factors of an @f$LDL^T@f$ factorization.
 */
typedef struct ldl_factors
{
    ladel_int ncol;         /**< number of columns in the analyzed matrix */
    ladel_sparse_matrix *L; /**< L in LDL' factorization */
    ladel_double *D;        /**< D in LDL' factorization (stored as vector), not used but is useful for returning */
    ladel_double *Dinv;     /**< D^-1 in LDL' factorization (stored as vector) */
    ladel_int *p;           /**< permutation vector */
    ladel_int *pinv;        /**< inverse permutation vector */
} ladel_factor;

/**
 * @brief Set of integers
 */
typedef struct ladel_set_struct {
    ladel_int *set;         /**< List of integers representing the set */
    ladel_int size_set;     /**< Size of the list */
    ladel_int max_size_set; /**< Maximum (allocated) size of the list */
} ladel_set;

/**
 * @brief Column of a sparse matrix
 */ 
typedef struct ladel_col_struct {
    ladel_int *i;       /**< List of row indices */
    ladel_double *x;    /**< List of values */
    ladel_int nz;       /**< Number of elements in the column */
    ladel_int nzmax;    /**< Maximum number of elements in the column (allocated space) */
} ladel_col;

/**
 * @brief Structure representing a multiple of the identity matrix
 */
typedef struct ladel_diag_struct {
    ladel_double diag_elem; /**< Scalar */
    ladel_int diag_size;    /**< Size of the matrix */
} ladel_diag;

/**
 * @brief Workspace required for various routines in LADEL 
 */
typedef struct workspace
{
    ladel_set *set_preallocated1;               /**< A preallocated set structure */
    ladel_set *set_preallocated2;               /**< A preallocated set structure */
    ladel_set *set_preallocated3;               /**< A preallocated set structure */
    ladel_set *set_unallocated_values1;         /**< An unallocated set structure */
    ladel_set *set_unallocated_values2;         /**< An unallocated set structure */
    ladel_set *set_unallocated_values3;         /**< An unallocated set structure */
    ladel_int *array_int_ncol1;                 /**< An array of @a ncol integers */
    ladel_int *array_int_ncol2;                 /**< An array of @a ncol integers */
    ladel_int *array_int_ncol3;                 /**< An array of @a ncol integers */
    ladel_int *array_int_ncol4;                 /**< An array of @a ncol integers */
    ladel_int *array_int_ncol_flag;             /**< An array of @a ncol integers, assumed to be < @a flag */
    ladel_int flag;                             /**< Flag used to mark nodes, used in conjunction with @a array_int_ncol_flag */
    ladel_double *array_double_all_zeros_ncol1; /**< An array of @a ncol doubles, on input and output this should be all zeros */
    ladel_double *array_double_ncol1;           /**< An array of @a ncol doubles */
} ladel_work;

#endif /*LADEL_TYPES_H*/