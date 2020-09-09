# QPALM.jl

This repository provides a Julia interface to the [QPALM](https://github.com/Benny44/QPALM)
QP solver.

## Installation

In order to use QPALM.jl, you will need the QPALM shared library, which you can find [here](https://bintray.com/benny44/generic/QPALM/1.0#files/). You will find a lib folder which you need to copy to path-to-QPALM.jl/deps. Also, install miniconda if you haven't yet, and download the dependencies
```
conda install -c conda-forge lapack
conda install -c intel mkl
```
Then, in a terminal where you wish to run QPALM, add the relevant paths to the LD_LIBRARY_PATH environment variable:
```
export LD_LIBRARY_PATH=path-to-miniconda/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=path-to-QPALM.jl/deps/lib/:$LD_LIBRARY_PATH
```

Then, from the Julia (>= 1.0) REPL, do

```julia
] dev path-to-QPALM.jl
```

To test the correct installation of the Julia interface, you can run the unit tests:
```julia
include("test/runtests.jl")
```

## Usage

Given the QP

```
minimize    (1/2) x' Q x + q' x
subject to  bmin <= A x <= bmax
```

this is solved by

```julia
using QPALM
model = QPALM.Model()
QPALM.setup!(model, Q=Q, q=q, A=A, bmin=bmin, bmax=bmax, settings...)
results = QPALM.solve!(model)
```

where `settings...` are keyword arguments specifying the solver options
to use. They have the same name and type as the underlying C API,
so please refer to [QPALM's documentation](https://benny44.github.io/QPALM/structQPALMSettings.html)
on how to set these and their semantics. Leaving the settings unspecified
will run the solver with default options.
