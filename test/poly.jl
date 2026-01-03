using ..TrihpFEM.Poly
using Test
using StaticArrays
using LinearAlgebra


p = BiPoly((1.,2,5.),(2.,-1))
@test indeterminates(p) == (:x,:y)
@test indeterminate(p.px) == :x
a = rand()
b = rand()
x = SVector(a,b)
y = [a,b]
@test p(x) == p(y) == p(a,b)
@test p(2,-1) ≈ 24
@test (b*p)(x) ≈ b*p(x)
q = BiPoly((1.,2.),(-2,0.,1.))
@test q(2,-1) ≈ -3.
pq = p*q
@test pq(x) ≈ p(x)*q(x)
qxp = q.px*p
@test qxp(y) == q.px(a)*p(a,b)
pyq = q*p.py
@test pyq(x) == p.py(b)*q(a,b)
ζ = BiPoly((0.,0.,1.,1.),(1.,2.),:x,:α)
@test_throws ArgumentError("Indeterminates does not match.") ζ*p
@test_throws ArgumentError("Indeterminates does not match.") p*ζ.py
@test zero(p)==zero(q)==zero(typeof(p))
@test q⋅p==q*p

v = PolyVectorField([p,q])
@test v(x)[1] == p(x)
@test v[x][2] == q(x)
@test length(v) == 2length(p)
@test [vv[x] for vv in v] == [p(x),q(x)]
@test (a*v)(a,b) == a*(v(a,b)) == (v*a)(a,b)
@test v(a,b)[1] == v[1](a,b)
A = rand(2,2)
w = A*v
@test typeof(w)::VectorTensor
@test size(w) == (2,)
row = [3 -2]
s = BiPoly(Tuple(rand() for _ in 1:rand(1:8)),Tuple(rand() for _ in 1:rand(1:8)),:x,:y)
@test typeof(row*s)::BiPoly
@test_throws DimensionMismatch() row'*s
sv = v*s
@test vs(a,b) == v(x)*s(x)

r = p+q
@test typeof(r)::PolySum
@test r(x) ≈ p(x)+q(x)
@test zero(p)==zero(p)
rs = r*s
@test rs(a,b) == r(x)*s(x)
@test (r+s)(x) == r(x)+s(x)
@test (r+rs)(x) == r(x)+rs(x)

aff = AffineToRef{Float64}()
v₁ = SVector(2rand()-1 for _ in 1:2)
v₂ = SVector(2rand()-1 for _ in 1:2)
v₃ = SVector(2rand()-1 for _ in 1:2)
affine!(aff,[v₁,v₂,v₃])
@test aff(-1,1) == v₁
@test aff(-1,-1) == v₂
@test aff(1,-1) == v₃

aff₂ = AffineToRef([√2/4 -√2/4;√2/4 √2/4],[0,√2/2])
area(aff₂) ≈ √2/4


f(x) = 2exp(x[1]+x[2])

wf2 = w*f*2
z = SVector(0.5,0.)
@test wf2(z,aff₂) ≈ 4*w(z)*exp(√2/2)

