using ..TrihpFEM.Poly
using ..TrihpFEM.Integration
using Test
using LinearAlgebra

p = BiPoly((0.,2.,5.),(1,-1.,0,-3))
q = BiPoly((1,0.,-3.),(2,1,-1.))

@test ref_integrate(p) ≈ 142/21
@test ref_integrate(q) ≈ -4/15
@test ref_integrate(p+q) ≈ ref_integrate(p)+ref_integrate(q)

A = rand(2,2)
B₁ = StandardBasis((2,3,4))
B₂ = StandardBasis((2,3,4))

sch = gmquadrature(Val(2),15)

f(x) = sin(π*x[1])*cos(π*x[2])
@test ref_integrate(f,sch) ≈ 1/π

for φ ∈ B₁
    for ψ ∈ B₂
        i₁ = ref_integrate((A*∇(φ))⋅∇(ψ))
        i₂ = ref_integrate((A*∇(φ))⋅∇(ψ),sch)
        @test isapprox(i₁,i₂,atol=1e-12,rtol=1e-10)
    end
end




