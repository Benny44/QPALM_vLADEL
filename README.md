# Proximal Augmented Lagrangian method for Quadratic Programs

| **QPALM is now maintained at [kul-optec/QPALM](https://github.com/kul-optec/QPALM)** |
|:---:|

---

<details>
<summary>View the original README</summary>

[![Coverage Status](https://coveralls.io/repos/github/Benny44/QPALM_vLADEL/badge.svg?branch=master)](https://coveralls.io/github/Benny44/QPALM?branch=master)

Platform | CI Status
---------|:---------
Linux    | [![Linux Build Status](https://travis-ci.org/Benny44/QPALM_vLADEL.svg?env=BADGE=linux&branch=master)](https://travis-ci.com/Benny44/QPALM)
OSX      | [![OSX Build Status](https://travis-ci.org/Benny44/QPALM_vLADEL.svg?env=BADGE=osx&branch=master)](https://travis-ci.com/Benny44/QPALM)

QPALM is a numerical optimization package that finds stationary points of (possibly **nonconvex**) quadratic programs, that is 
```
minimize        0.5 x' Q x + q' x

subject to      l <= A x <= u
```






<!-- ## Installation

First of all, clone this repo with all the submodules! For this, run
```
git clone https://github.com/Benny44/QPALM_vLADEL.git
git submodule update --init --recursive
```

### **Matlab**
* To install the mex interface of QPALM, add QPALM and its subfolders to the matlab path. Then run qpalm_make.m. You can test whether QPALM is working using 
the examples/qpalm_mex_demo.m and examples/qpalm_mex_nonconvex_demo.m.
### **C**
* To install a C-callable library, check [Bintray](https://bintray.com/benny44/generic/QPALM) for the binaries. These were compiled against [lapack](https://anaconda.org/conda-forge/lapack) from conda. So install miniconda and run the following commands
```
conda install -c conda-forge lapack
export LD_LIBRARY_PATH=path-to-miniconda/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=path-to-qpalm-binaries/lib/:$LD_LIBRARY_PATH
```
### **Python**
The python interface has been compiled for python version 3.7. If you want to use a different version, do your own install with the instructions on custom compilation below.

Follow the instructions for installing the C-library above. Then in an open terminal, do
```
export LD_LIBRARY_PATH=path-to-qpalm-binaries/interfaces/python/build/lib/:$LD_LIBRARY_PATH
python3 path-to-qpalm-binaries/interfaces/python/qpalm_python_demo.py
```

### **Julia**
See [QPALM.jl](https://github.com/kul-forbes/QPALM.jl/tree/856c70d2be99a24e5d9a6391be45cf40c48947d4) for the instructions on installing the Julia interface.

## Custom Compilation
If you wish to do a custom compilation of the shared libraries, take a look at buildCustom.sh. First install the dependencies
```
conda install -c conda-forge lapack
```

Then change the following lines near the top of the script
```
export MINICONDA_LIB=path-to-miniconda/lib
export MINICONDA_INCLUDE=path-to-miniconda/include
```

Furthermore, change the cmake line to have whatever flags you want. To build the release version (with tests), use
```
cmake $curdir -DCMAKE_BUILD_TYPE=release -DCOVERAGE=ON
```
To build the python interface, use instead
```
cmake path-to-QPALM -DCMAKE_BUILD_TYPE=release -DINTERFACES=OFF -DUNITTESTS=OFF -DPYTHON=ON
```

Finally, run the buildCustom.sh script
```
chmod 755 buildCustom.sh
./buildCustom.sh
```

## Code Examples

Basic demos are available for the different ways to call the solver:
* For the mex interface of QPALM, check out examples/qpalm_mex_demo.m and examples/qpalm_mex_nonconvex_demo.m.
* For the C-version of QPALM, check out examples/qpalm_demo.c.
* For the python interface of QPALM, check out interfaces/python/qpalm_python_demo.py.
* For the Julia interface of QPALM, check out any of the files in interfaces/QPALM.jl/test/. -->

## Documentation

You can now find the the documentation [here](https://benny44.github.io/QPALM_vLADEL/). This includes all information you need to get started using QPALM.

## Benchmarks

Check out the paper below for detailed benchmark tests comparing QPALM with state-of-the-art solvers.

 * [QPALM: A Newton-type Proximal Augmented Lagrangian Method for Quadratic Programs](https://arxiv.org/abs/1911.02934)
 * [QPALM: A Proximal Augmented Lagrangian Method for Nonconvex Quadratic Programs](https://arxiv.org/abs/2010.02653)

## Citing

If you use QPALM in your research, please cite the following paper
```
@inproceedings{hermans2019qpalm,
	author      = {Hermans, B. and Themelis, A. and Patrinos, P.},
	booktitle   = {58th IEEE Conference on Decision and Control},
	title       = {{QPALM}: {A} {N}ewton-type {P}roximal {A}ugmented {L}agrangian {M}ethod for {Q}uadratic {P}rograms},
	year        = {2019},
	volume      = {},
	number      = {},
	pages       = {},
	doi         = {},
	issn        = {},
	month       = {Dec.},
}
```

## License

QPALM is licensed under LGPL v3.0. Some modules are used in this software: 
* LADEL: authored by Ben Hermans and licensed under [LGPL-v3](https://github.com/Benny44/LADEL/blob/master/LICENSE).

* LOBPCG: the version of LOBPCG used here was written by Ben Hermans and licensed under the GNU Lesser General Public License v3.0, see [LOBPCG/LICENSE](https://github.com/Benny44/LOBPCG/blob/master/LICENSE).
* LAPACK: authored by The University of Tennessee and The University of Tennessee Research Foundation, The University of California Berkeley, and The University of Colorado Denver, and licensed under BSD-3, see [here](https://github.com/Reference-LAPACK/lapack/blob/master/LICENSE).
* Minunit: a minimal unit testing framework for C, modified from the version by David Siñuela Pastor and licensed under MIT, see [here](https://github.com/siu/minunit/blob/master/MIT-LICENSE.txt).

</details>
