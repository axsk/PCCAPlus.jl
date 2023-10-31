module KrylovExt
import PCCAPlus
import KrylovKit




function PCCAPlus.schurvecs(T, n, israte, ::PCCAPlus.KrylovSolver)
    which = israte ? :LR : :LM
    R, Qs, = KrylovKit.schursolve(T, rand(size(T, 1)), n, which, KrylovKit.Arnoldi())
    Q = reduce(hcat, Qs)
    Q[:, 1:n]
end

end # module
