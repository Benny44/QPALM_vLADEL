# LADEL
Quasidefinite sparse LDL factorization package with rank 1 updates and (symmetric) row/column adds and deletes. 

## Documentation

The documentation for LADEL can be found [online](https://benny44.github.io/LADEL/).

## License

LADEL is licensed under LGPL v3.0. This library relies on the following external module:
* AMD, licensed under [BSD 3](https://github.com/Benny44/LADEL/blob/master/amd/License.txt), Copyright (c), 1996-2015, Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff, see (Amestoy et al., 2004).

LADEL contains an implementation of several sparse matrix routines outlined in (Davis, 2006), separate from the reference implementation in [CSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse/tree/master/CSparse), licensed under [LGPL v2.1](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/CSparse/Doc/License.txt).

## References
* Amestoy, Patrick R., Timothy A. Davis, and Iain S. Duff. "Algorithm 837: AMD, an approximate minimum degree ordering algorithm." ACM Transactions on Mathematical Software (TOMS) 30.3 (2004): 381-388.
* Davis, Timothy A. Direct Methods for Sparse Linear Systems. Philadelphia: SIAM, 2006. Print. Fundamentals and Algorithms.