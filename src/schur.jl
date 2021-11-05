### Schur Eigenspace solvers

struct BaseSolver end
struct ArnoldiSolver end
struct KrylovSolver end



import ArnoldiMethod
function schurvectors(T, pi, n, israte, ::ArnoldiSolver)
	D = Diagonal(sqrt.(pi))
	T̃ = D * T * D^-1
	s, hist = ArnoldiMethod.partialschur(T̃; nev=n, which=israte ? ArnoldiMethod.LR() : ArnoldiMethod.LM())
	X̃ = collect(s.Q)
	X = D^-1 * X̃
	X ./= X[1,1]
	X, s.eigenvalues
end

import KrylovKit
function schurvectors(T, pi, n, israte, ::KrylovSolver)
	D = Diagonal(sqrt.(pi))
	T̃ = D * T * D^-1
	R, Q, v, info = KrylovKit.schursolve(T̃, pi, n, israte ? :LR : :LM, KrylovKit.Arnoldi())
	X̃ = reduce(hcat, Q)
	X = D^-1 * X̃
	X = X ./ X[1,1]
	X, v
end

function schurvectors(T, pi, n, israte, ::BaseSolver)
    Tw = Diagonal(sqrt.(pi))*T*Diagonal(1 ./ sqrt.(pi)) # rescale to keep markov property
    Sw = schur!(Tw)                       # returns orthonormal vecs by def
    Xw, λ = selclusters!(Sw, n, israte)
    X  = Diagonal(1 ./sqrt.(pi)) * Xw              # scale back
    X  = X[1,1]>0 ? X : -X
    X, λ
end

# select the schurvectors corresponding to the n abs-largest eigenvalues
# if reverse==true select highest abs value, otherwise select lowest (for rate matrices)
function selclusters!(S, n, ratematrix)
    ind = sortperm(abs.(S.values), rev=!ratematrix) # get indices for dominant eigenvalues
    select = zeros(Bool, size(ind))           # create selection vector
    select[ind[1:n]] .= true
    S = ordschur!(S, select)                  # reorder selected vectors to the left
    if !isapprox(S.T[n+1, n], 0)                       # check if we are cutting along a schur block
        @error("conjugated eigenvector missing")
        display(S.T)
    end
    S.vectors[:,1:n], S.values[1:n]       # select first n vectors
end
