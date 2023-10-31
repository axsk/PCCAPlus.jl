# PCCAPlus

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://axsk.github.io/PCCAPlus.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://axsk.github.io/PCCAPlus.jl/dev)
[![Build Status](https://github.com/axsk/PCCAPlus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/axsk/PCCAPlus.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/axsk/PCCAPlus.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/axsk/PCCAPlus.jl)

A [KISS](https://en.wikipedia.org/wiki/KISS_principle) style implementation of PCCA+ (Robust Perron Cluster Analysis) [1,2] with support for non-reversible systems [3].
For a similar python implementation see also the [cmdtools](https://github.com/zib-cmd/cmdtools/) package.

## Basic usage

```julia
using PCCAPlus

P=rand(10,10)
P = P ./ sum(P, dims=2) # row stochastic matrix

# basic PCCA+ clustering with 2 clusters (using no weighting and the ISA initial guess only)
chi = pcca(P, 2)        # 

using KrylovKit
using SparseArrays
P = sprand(100,100, 0.1)
P = P ./ sum(P, dims=2) # sparse row stochastic matrix

# solve the PCCA+ problem weighted with the stationary density 
# and optimize for crispness, using the KrylovKit.jl eigensolver
chi = pcca(P, 2; pi=:auto, optimize=true, solver=KrylovSolver())
```

For sparse matrix support, add either the `ArnoldiMethod.jl` or `KrylovKit.jl` and pass the corresponding `ArnoldiSolver()` or `KrylovSolver()` as a solver.

## References
1. 2006, M. Weber: Meshless Methods in Conformation Dynamics
2. 2013, S. Röblitz, M. Weber: Fuzzy Spectral Clustering by PCCA+
3. 2018, K. Fackeldey, A. Sikorski, M. Weber: Spectral Clustering for Non-Reversible Markov Chains