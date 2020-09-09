#include "ladel_global.h"
#include "ladel_types.h"
#include "ladel_matvec.h"

void ladel_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset)
{
    ladel_int index, col;

    if (reset) for (index = 0; index < M->nrow; index++) y[index] = 0;

    for (col = 0; col < M->ncol; col++)
        LADEL_FOR(index, M, col)
            y[M->i[index]] += M->x[index] * x[col];
}

void ladel_tpose_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset)
{
    ladel_int index, col;

    if (reset) for (index = 0; index < M->ncol; index++) y[index] = 0; 
    
    for (col = 0; col < M->ncol; col++)
        LADEL_FOR(index, M, col)
            y[col] += M->x[index] * x[M->i[index]];
}

void ladel_symmetric_matvec(const ladel_sparse_matrix *M, const ladel_double *x, ladel_double *y, ladel_int reset)
{
    ladel_int index, col;

    if (reset) for (index = 0; index < M->ncol; index++) y[index] = 0;

    for (col = 0; col < M->ncol; col++)
        LADEL_FOR(index, M, col)
            y[M->i[index]] += (M->i[index] == col) ? 0.0 : M->x[index] * x[col];

    ladel_tpose_matvec(M, x, y, FALSE);
}