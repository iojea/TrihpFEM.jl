using ..TrihpFEM.Spaces
using Test

S = StdScalarSpace()
a = rand()
A = rand(2,2)
f(x) = x[1]+x[2]
@test order(a) == Order{0}()
@test order(A) == Order{0}()
@test order(f) == Order{0}()
@test order(S) == Order{(0,)}()
@test order(∇(S)) == Order{(1,)}()
@test order(∇(S)⋅∇(S)) == Order{(1,1)}()
@test order(divergence(∇(S)))== Order{(2,)}()
@test order(Δ(S)) == Order{(2,)}()

p = BiPoly((1.,2,5.),(2.,-1))
@test coefftype(p)==Constant()
@test coefftype(a)==Constant()
@test coefftype(A)==Constant()
@test coefftype(f)==Variable()
@test coefftype(A*∇(S))==Constant()
@test coefftype(f*S)==Variable()

@test basis(S,(1,2,3)) == StandardBasis((1,2,3))



