using TrihpFEM
using LinearAlgebra

u = ShapeFunction(Trial())
v = ShapeFunction(Test())
Ω = circmesh(0.1)
dΩ = Measure(Ω)

a = ∫(∇(u)⋅∇(v))*dΩ
b = ∫(2v)*dΩ

prob = FEProblem(a,b)

U = solve(prob)

plotsol(U)

