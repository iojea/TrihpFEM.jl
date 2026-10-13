using ..TrihpFEM.Meshes
using ..TrihpFEM.Tensors
using ..TrihpFEM.Forms
using ..TrihpFEM.Problems

using Test


vert = [0.0 0.0;1.0 0.0;1.0 1.0;0.0 1.0]'
Ω  = hpmesh(vert, sqrt(2) / 2)
dΩ = Measure(Ω)

@test typeof(dΩ)<:Measure
@test dΩ.sch.degree == 3

DΩ = Measure(Ω,17)
@test typeof(DΩ)<:Measure
@test DΩ.sch.degree == 17
@test typeof(DΩ.sch.degree) == Int8

A = Tensor{2,2}(rand(),rand(),rand(),rand())
g(x) = x[1]*x[2]
u = ShapeFunction(Trial())
v = ShapeFunction(Test())
a = ∫((A*∇(u))⋅∇(v))*dΩ + ∫(2(u*v))*dΩ
b = ∫(g*v)*dΩ

c =  ∫(∇(u)⋅∇(v))*dΩ

@test typeof(a)<:Form{2}
@test typeof(b)<:Term{1}
@test typeof(c)<:Term{2}


b₀ = ∫(2*v)*dΩ
b₁ = ∫(Tensor{1,2}(1,2)⋅∇(v))*dΩ
b = b₀ + b₁
@test b₀ isa Term{1}
@test b isa Form{1}

h(x) = x[1]
prob = FEProblem(a,b,h)

@test typeof(prob)<:FEProblem
