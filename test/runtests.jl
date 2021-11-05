using PCCA
using Test

@testset "PCCA.jl" begin
    # Write your tests here.
	x = rand(10, 10)
	x = x ./ sum(x, dims=2)
	c = pcca(x, 2)

	function randomstochasticmatrix(n, reversible=true)
		P = rand(n,n)
		if reversible
			P = (P + P') / 2
		end
		P ./= sum(P, dims=2)
	end



	@testset "Reversible = $rev" for rev in [true, false]
		P = [randomstochasticmatrix(3+mod(i,12), rev) for i in 1:10]
		@testset "Method $method" for method in [PCCA.BaseSolver, PCCA.ArnoldiSolver, PCCA.KrylovSolver]
			for P in P, n in 2:size(P,1)-1
				χ = pcca(P, n, solver=method(), optimize=true)
				a = PCCA.crispassignments(χ)

				#@show χ
				@test all(sum(χ, dims=2) .≈ 1)
			end
		end
	end

end
