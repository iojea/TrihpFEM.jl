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
@test typeof(DΩ.sch.degree) == Int32

A = rand(2,2)
f(x) = x[1]*x[2]
@form a(u,v) = ∫((A*∇(u))⋅∇(v))*dΩ
@form b(v) = ∫(f*v)*dΩ

S = StdScalarSpace()
g(x) = x[1]
prob = FEProblem(a,b,S,g)

@test typeof(prob)<:FEProblem{StdScalarSpace}
