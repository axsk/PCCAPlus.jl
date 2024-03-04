module KrylovExt

import PCCAPlus
import KrylovKit

function PCCAPlus.schurvecs(T, n, israte, k::PCCAPlus.KrylovSolver)
    which = israte ? :LR : :LM
    R, Qs, vals, info = KrylovKit.schursolve(T, rand(size(T, 1)), n, which, KrylovKit.Arnoldi(; k.kwargs...))
    @show vals, info
    Q = reduce(hcat, Qs)
    Q[:, 1:n]
end

end # module
