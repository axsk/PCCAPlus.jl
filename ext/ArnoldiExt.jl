module ArnoldiExt

import PCCA
import ArnoldiMethod

function PCCA.schurvecs(T, n, israte, ::PCCA.ArnoldiSolver)
    which = israte ? ArnoldiMethod.LR() : ArnoldiMethod.LM()
    Q = ArnoldiMethod.partialschur(T; nev=n, which)[1].Q
    Q[:, 1:n]
end

end # module