
#include "qpalm.h"

QPALMData *qpalm_julia_set_data(c_int           n, 
                                c_int           m, 
                                solver_sparse   *Q, 
                                solver_sparse   *A, 
                                c_float         *q, 
                                c_float         c, 
                                c_float         *bmin, 
                                c_float         *bmax) 
{
  QPALMData *data   = c_calloc(1, sizeof(QPALMData));
  data->n           = n;           
  data->m           = m;                   
  data->bmin        = bmin;      
  data->bmax        = bmax;       
  data->q           = q;          
  data->c           = c;

  data->A           = ladel_sparse_allocate_and_copy(A); 
  data->Q           = ladel_sparse_allocate_and_copy(Q);

  /* Conversion of indices */
  c_int index;
  for (index = 0; index < data->A->nzmax; index++) data->A->i[index]--;
  for (index = 0; index < data->Q->nzmax; index++) data->Q->i[index]--;
  for (index = 0; index < data->n+1; index++)
  {
      data->A->p[index]--;
      data->Q->p[index]--;
  } 
  return data;
}
