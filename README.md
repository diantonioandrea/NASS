# NASS

_**N**eon **A**ccelerated **S**ketched **S**olver_

## Introduction

An implementation of a [Randomized Sketched GMRES with Truncated Arnoldi Orthogonalization](https://doi.org/10.48550/arXiv.2111.00113) (_k-sGMRES_) for real linear systems utilizing [_Neon_](https://developer.arm.com/Architectures/Neon) intrinsics, written in `C++23`.

## Table of Contents

- [Introduction](#introduction)
- [Table of Contents](#table-of-contents)
- [Setup](#setup)
    - [Cloning the Repository](#cloning-the-repository)
    - [Compiling Tests](#compiling-tests)
    - [Flags](#flags)
- [Usage](#usage)
    - [Running Tests](#running-tests)
        - [`Test_sGMRES`](#test_sgmres)

## Setup

### Cloning the Repository

Clone the repository from [here](https://github.com/diantonioandrea/NASS):

```bash
git clone git@github.com:diantonioandrea/NASS.git
```

### Compiling Tests

Compile the tests located under `src/`:

```bash
make
```

### Flags

Some flags in the [Makefile](./Makefile) modify the behaviour of the code.

- `-DMEMORY_PRIORITY`: The algorithm prioritizes memory by building and discarding the basis and storing only the needed vectors at any given time.
- `-DNDEBUG`: Disables debugging.
- `-DNVERBOSE`: Disables verbosity.
- `-DNEON32`: Enables 32-bit Neon instructions instead of (default) 64-bits instructions.

## Usage

Every method developed in **NASS**, similarly to [**NAQRA**](https://github.com/diantonioandrea/NAQRA), follows a structured naming convention with three parts, separated by underscores:

1. **Function**: Describes the primary operation or behavior of the method.
2. **Input(s)**: Specifies the type or nature of the inputs.
3. **Output**: Indicates the type or nature of the resulting output.

All methods are documented directly within the code for clarity and ease of use.

### Running Tests

Tests executables are located under `executables/`.

#### `Test_sGMRES`

`Test_sGMRES` evaluates the performance of the `sGMRES` algorithm for solving a linear system using a square sparse matrix and a randomly generated right-hand side, requiring the following inputs:

1. **Path to a matrix file**: The matrix must be stored in column-major order and in the [`.mtx`](https://math.nist.gov/MatrixMarket/formats.html#MMformat) format.
2. **Dimension of the Krylov subspace**: An integer specifying the number of basis vectors to use.
3. **Arnoldi truncation level** *(optional)*: An integer defining the truncation level for the Arnoldi process used to construct the Krylov subspace basis, defaults to `4`.

The following command demonstrates how to run `Test_sGMRES`:

```bash
./executables/Test_sGMRES.out data/5M.mtx 100