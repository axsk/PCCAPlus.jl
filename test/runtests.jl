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

	P = randomstochasticmatrix(8)

	for method in [PCCA.BaseSolver, PCCA.ArnoldiSolver, PCCA.KrylovSolver]

		@testset "$method" begin
				χ = pcca(P, 2, solver=method())
				a = PCCA.crispassignments(χ)

				#@show χ
				@test all(sum(χ, dims=2) .≈ 1)
		end
	end
end
