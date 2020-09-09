#ifndef LADEL_MEX_UTIL_H
#define LADEL_MEX_UTIL_H

#include "mex.h"
#include "ladel_types.h"

ladel_sparse_matrix *ladel_get_sparse_from_matlab(const mxArray *M_mex, ladel_sparse_matrix *M, ladel_int symmetry);

mxArray *ladel_put_matlab_from_sparse(ladel_sparse_matrix *M);

ladel_sparse_matrix *ladel_convert_factor_to_sparse(ladel_sparse_matrix *L);


#endif /* LADEL_MEX_UTIL_H */