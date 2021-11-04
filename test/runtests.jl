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

	function test()
		χ = pcca(randomstochasticmatrix(8), 2)
		a = PCCA.crispassignments(χ)

		@test all(sum(χ, dims=2) .≈ 1)
	end

	test()
end
