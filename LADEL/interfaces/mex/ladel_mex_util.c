#include "mex.h"
#include "ladel_types.h"
#include "ladel_global.h"
#include "ladel_constants.h"
#include "ladel_mex_util.h"

ladel_sparse_matrix *ladel_get_sparse_from_matlab(const mxArray *M_mex, ladel_sparse_matrix *M, ladel_int symmetry)
{
    M->nrow = mxGetM(M_mex);
    M->ncol = mxGetN(M_mex);
    M->p = (ladel_int *) mxGetJc(M_mex);
    M->i = (ladel_int *) mxGetIr(M_mex);
    M->x = (ladel_double *) mxGetPr(M_mex);
    M->values = TRUE;
    M->symmetry = symmetry;
    M->nzmax = M->p[M->ncol];
    M->nz = NULL;
    return M;
}

mxArray *ladel_put_matlab_from_sparse(ladel_sparse_matrix *M)
{
    if (!M || M->x == NULL) 
        mexErrMsgTxt ("Tried converting an unitialized sparse matrix to matlab.") ;

    mxArray *M_matlab ;
    
    M_matlab = mxCreateSparse(M->nrow, M->ncol, M->nzmax, mxREAL);
    mxSetJc(M_matlab, M->p);
    mxSetIr(M_matlab, M->i);
    mxSetPr(M_matlab, M->x);

    return (M_matlab);
}

ladel_sparse_matrix *ladel_convert_factor_to_sparse(ladel_sparse_matrix *L)
{
    /* The difference is that the number of nonzeros in L is indicated by L->nz,
    *  whereas in L_out they are given by L_out->p */
   if (!L) return NULL;
   ladel_int index_L, index_L_out = 0, col, ncol = L->ncol, nrow = L->nrow, Lnz = 0;
   for (col = 0; col < ncol; col++) Lnz += L->nz[col];
   ladel_sparse_matrix *L_out = ladel_sparse_alloc(nrow, ncol, Lnz, L->symmetry, L->values, FALSE);
   if (!L_out) return NULL;
   L_out->p[0] = 0;
   for (col = 0; col < ncol; col++)
   {
       for (index_L = L->p[col]; index_L < L->p[col]+L->nz[col]; index_L++)
       {
           L_out->i[index_L_out] = L->i[index_L];
           L_out->x[index_L_out] = L->x[index_L];
           index_L_out++;
       }
       L_out->p[col+1] = index_L_out; 
   }
   return L_out;
}