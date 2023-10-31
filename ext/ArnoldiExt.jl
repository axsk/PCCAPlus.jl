module ArnoldiExt

import PCCAPlus
import ArnoldiMethod

function PCCAPlus.schurvecs(T, n, israte, ::PCCAPlus.ArnoldiSolver)
    which = israte ? ArnoldiMethod.LR() : ArnoldiMethod.LM()
    Q = ArnoldiMethod.partialschur(T; nev=n, which)[1].Q
    Q[:, 1:n]
end

end # module