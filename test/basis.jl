using ..TrihpFEM.Poly
using Test
using StaticArrays
using Polynomials

@test_throws ArgumentError("Degrees does not satisfy p conformity.") StandardBasis((1,2,4))

B = StandardBasis((2,3,4))
p₀ = Polynomial((1.,))
p₁ = Polynomial((0.,1.))
p₂ = Polynomial((-1.,0.,3))/2
p₃ = Polynomial((0.,-3,0,5))/2
l₂ = [p₀,p₁,p₂]
l₃ = [p₀,p₁,p₂,p₃]

k = 0
for b in B
    (i,j) = degs(b)
    @test l₂[i+1](0.)*l₃[j+1](0.) ≈ b([0.,0.])
    @test l₂[i+1](0.23)*l₃[j+1](0.77) ≈ b(SVector(0.23,0.77))
    global k = k+1
end

@test length(B) == k