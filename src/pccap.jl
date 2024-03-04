import Arpack: eigs
using LinearAlgebra


""" pcca(T, n; pi=nothing, optimize=false, solver=BaseSolver())

performs the pcca clustering on the transition matrix `T` with `n` clusters.
Here `T` can be either a stochastic propagator with row sum 1 or a rate matrix with row sum 0.

- `pi` is a density for weighting the result. `pi=:stationary` uses the stationary density
- `optimize` uses the optimization from Roeblitz (2013) to improve crispness
- `solver` is the solver to use for computing the schur decomposition. 
   The BaseSolver() is built in. ArnoldiSolver() and KrylovSolver() require ArnoldiMethod.jl and KrylovKit.jl respectively and provide support for sparse matrices.

Returns the membership matrix `chi` assigning to each state a membership to one of the `n` clusters.
"""

function pcca(T::AbstractMatrix, n::Integer; pi=nothing, optimize=false, solver=BaseSolver())
    israte = isratematrix(T)
    if pi == :stationary
        pi = stationarydensity(T, israte)
    end
    X = schurvectors(T, pi, n, israte, solver)
    #assertstructure(X, pi)
    A = basistransform(X, optimize)
    chi = X * A
    return (; chi, X, A)
end

function isratematrix(T::AbstractMatrix)
    s = sum(T[1, :])
    isapprox(s, 0, atol=1e-8) && return true
    isapprox(s, 1, atol=1e-8) && return false
    error("given matrix is neither a rate nor a probability matrix")
end

function stationarydensity(T, israte=isratematrix(T))
    which = israte ? :SM : :LM
    pi = eigs(T', nev=1, which=which)[2] |> vec
    #pi = eigvecs(T, sortby=real)[:,end]
    @assert isreal(pi)
    pi = abs.(pi)
    pi = pi / sum(pi)
end

function basistransform(X::AbstractMatrix, optimize::Bool)
    n = size(X, 2)
    A = innersimplexalgorithm(X)
    if n > 2 && optimize
        A = opt(A, X)
    end
    return A
end

function crispassignments(chi)
    assignments = mapslices(argmax, chi, dims=2) |> vec
end

""" This checks whether the resulting subspace basis satisfies
a) orthonormality wrt. `pi`
b) that the first column is 1 which is required by `feasiblize!`
"""
function assertstructure(X, pi)
    @assert X' * Diagonal(pi) * X ≈ I
    @assert all(isapprox.(X[:, 1], 1, atol=1e-8))
end

function assertstructure(X, pi::Nothing)
    @assert X' * X ./ size(X, 1) ≈ I
    @assert all(isapprox.(X[:, 1], 1, atol=1e-8))
end

# compute initial guess based on indexmap
innersimplexalgorithm(X) = feasiblize!(inv(X[indexmap(X), :]), X)

function indexmap(X)
    # get indices of rows of X to span the largest simplex
    rnorm(x) = sqrt.(sum(abs2.(x), dims=2)) |> vec
    ind = zeros(Int, size(X, 2))
    for j in 1:length(ind)
        rownorm = rnorm(X)
        # store largest row index
        ind[j] = argmax(rownorm)
        if j == 1
            # translate to origin
            X = X .- X[ind[1], :]'
        else
            # remove subspace
            X = X / rownorm[ind[j]]
            vt = X[ind[j], :]'
            X = X - X * vt' * vt
        end
    end
    return ind
end

# Algorithm 3.10 from 2006, Weber, Meshless Methods in Conformation dynamcis
# this ensurses that A leads to a feasible clustering χ
# Note: it only works if X[:,1] .== 1
function feasiblize!(A, X)
    A[:, 1] = -sum(A[:, 2:end], dims=2)
    A[1, :] = -minimum(X[:, 2:end] * A[2:end, :], dims=1)
    A / sum(A[1, :])
end

# crispness criterion, cf. Roeblitz (2013)
# only applies if X is normalized (eq. 8)
# TODO: check does this apply?
function roeblitzcrit(A)
    n = size(A, 1)
    trace = 0
    for i = 1:n, j = 1:n
        trace += A[i, j]^2 / A[1, j]
    end
    return trace
end

using Optim
function opt(A0, X)
    A = copy(A0)
    Av = view(A, 2:size(A, 1), 2:size(A, 2)) # view on the variable part

    function obj(a)
        Av[:] = a
        -roeblitzcrit(feasiblize!(A, X))
    end

    result = optimize(obj, Av[:], NelderMead())
    Av[:] = result.minimizer
    return feasiblize!(A, X)
end



function isreversible(P, pi=stationarydensity(P))
    db = [pi[i] * P[i, j] - pi[j] * P[j, i] for i = 1:5, j = 1:5]
    all(isapprox.(db, 0, atol=1e-8))
end
