### Schur Eigenspace solvers

# TODO: there should be at least warnings when cutting conjugate eigenvalues

""" Uses Julia's built in LinearAlgebra.schur solver. This has no support for sparse matrices """
struct BaseSolver end
""" Wrapper around the ArnoldiSolver.jl schur solver """
struct ArnoldiSolver end
""" Wrapper around the KrylovKit.jl schur solvers """
struct KrylovSolver end

function schurvectors(T, pi, n, israte, solver)
    israte && warn("Current implementation only uses :LR as selection criterion, even for Q matrices.")
    D = Diagonal(sqrt.(pi))
    Tp = D * T * D^-1
    Qp = schurvecs(Tp, n, israte, solver)
    X = D^-1 * Qp
    X .*= sign(X[1, 1])  # fix the sign for the constant eigenvector
    return X
end

function schurvectors(T, pi::Nothing, n, israte, solver)
    X = schurvecs(T, n, israte, solver)
    X ./= X[1, 1] # renormalize so that first column becomes 1
    return X
end

using SparseArrays
function schurvecs(T, n, israte, ::BaseSolver)
    issparse(T) && error("The `BaseSolver` does not support sparse matrices")
    S = schur(T)
    Q, λ = selclusters!(S, n, israte)
    return Q
end

# select the schurvectors corresponding to the n largest real part of eigenvalues
function selclusters!(S, n, israte)
    sortby = israte ? real.(S.values) : abs.(S.values)
    ind = sortperm(sortby, rev=true)          # get indices for dominant eigenvalues
    select = zeros(Bool, size(ind))           # create selection vector
    select[ind[1:n]] .= true
    S = ordschur!(S, select)                  # reorder selected vectors to the left
    if !isapprox(S.T[n+1, n], 0)                       # check if we are cutting along a schur block
        @error("conjugated eigenvector missing")
        display(S.T)
    end
    S.vectors[:, 1:n], S.values[1:n]       # select first n vectors
end
