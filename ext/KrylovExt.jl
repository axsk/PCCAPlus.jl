module KrylovExt
import PCCA
import KrylovKit




function PCCA.schurvecs(T, n, israte, ::PCCA.KrylovSolver)
    which = israte ? :LR : :LM
    R, Qs, = KrylovKit.schursolve(T, rand(size(T, 1)), n, which, KrylovKit.Arnoldi())
    Q = reduce(hcat, Qs)
    Q[:, 1:n]
end

end # module
