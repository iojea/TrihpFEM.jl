using ..TrihpFEM.Meshes
using ..TrihpFEM.Spaces
using ..TrihpFEM.Forms
using ..TrihpFEM.Problems

using Test


vert = [0.0 0.0;1.0 0.0;1.0 1.0;0.0 1.0]'
Ω  = hpmesh(vert, sqrt(2) / 2)
dΩ = Measure(Ω)

@test typeof(dΩ)<:Measure
@test dΩ.sch.degree == 3

DΩ = Measure(Ω,13)
@test typeof(DΩ)<:Measure
@test DΩ.sch.degree == 13
@test typeof(DΩ.sch.degree) == UInt8

A = rand(2,2)
g(x) = x[1]*x[2]
@form a(u,v) = ∫((A*∇(u))⋅∇(v))*dΩ + ∫(2(u*v))*dΩ
@form b(v) = ∫(g*v)*dΩ
c = Form(IntegrationTerm((u,v)->∇(u)⋅∇(v),nothing,dΩ))

@test typeof(a)<:Form{2}
@test typeof(b)<:Form{1}
@test typeof(c)<:Form{2}

@test a.terms[1].factor == A
@test b.terms[1].factor == g
@test Forms.coefftype(a.terms[1]) isa ConstantCoeff
@test Forms.coefftype(b.terms[1]) isa VariableCoeff

@term b₀(v) = ∫(2*v)*dΩ
@term b₁(v) = ∫([1,2]⋅∇(v))*dΩ
b = b₀ + b₁
@test b₀ isa IntegrationTerm{ConstantCoeff,1}
@test b isa Form{1}

S = StdScalarSpace()
h(x) = x[1]
prob = FEProblem(a,b,S,h)

@test typeof(prob)<:FEProblem{StdScalarSpace}
