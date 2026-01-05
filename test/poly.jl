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
@test p(2,-1) ≈ 75
@test (b*p)(x) ≈ b*p(x)
q = BiPoly((1.,2.),(-2,0.,1.))
@test q(2,-1) ≈ -5.
pq = p*q
@test pq(x) ≈ p(x)*q(x)
qxp = q.px*p
@test qxp(y) ≈ q.px(a)*p(a,b)
pyq = q*p.py
@test pyq(x) ≈ p.py(b)*q(a,b)
ζ = BiPoly((0.,0.,1.,1.),(1.,2.),:x,:α)
@test_throws ArgumentError("Indeterminates does not match.") ζ*p
@test_throws ArgumentError("Indeterminates does not match.") p*ζ.py
@test zero(p)==zero(q)==zero(typeof(p))
@test q⋅p==q*p

v = PolyVectorField([p,q])
@test v(x)[1] == p(x)
@test v(x)[2] == q(x)
@test length(v) == 2length(p)
@test v(x) == [p(x),q(x)]
@test (a*v)(a,b)≈ a*(v(a,b)) ≈ (v*a)(a,b)
@test v(a,b)[1] == v[1](a,b)
A = rand(2,2)
w = A*v
@test typeof(w)<:PolyVectorField
@test size(w) == (2,)
row = [3 -2]
@test typeof(row*w)<:PolyScalarField
s = BiPoly(Tuple(rand() for _ in 1:rand(1:8)),Tuple(rand() for _ in 1:rand(1:8)),:x,:y)
@test_throws DimensionMismatch() row'*w
sv = v*s
@test sv(a,b) ≈ v(x)*s(x)

r = p+q
@test typeof(r)<:PolySum
@test r(x) ≈ p(x)+q(x)
@test zero(p)==zero(p)
rs = r*s
@test rs(a,b) ≈ r(x)*s(x)
@test (r+s)(x) ≈ r(x)+s(x)
@test (r+rs)(x) ≈ r(x)+rs(x)

aff = AffineToRef{Float64}()
v₁ = SVector(rand()+2,rand()+1)
v₂ = SVector(rand(),rand()+2)
v₃ = SVector(rand()+1,rand())
affine!(aff,[v₁,v₂,v₃])
@test aff(SVector(-1,1)) ≈ v₁
@test aff(SVector(-1,-1)) ≈ v₂
@test aff(SVector(1,-1)) ≈ v₃
@test area(aff) ≈ 0.5abs(v₁[1]*(v₂[2]-v₃[2])+v₂[1]*(v₃[2]-v₁[2])+v₃[1]*(v₁[2]-v₂[2]))

aff₂ = AffineToRef([√2/4 -√2/4;√2/4 √2/4],[0,√2/2])
@test area(aff₂) ≈ 1/2


f(x) = 2exp(x[1]+x[2])

wf2 = (w*f)*2
z = SVector(-0.5,0.)
@test wf2(z,aff₂) ≈ 4*w(z)*exp(√2/4)
z = SVector(0.,0.5)
wf2(z,aff₂) ≈ 4*w(z)*exp(√2/2)

@test repr(p) == "(1.0 + 2.0*x + 5.0*x^2)(2.0 - 1.0*y)"
@test repr(zero(p)) == "(0.0,)"

px = derivative(p,:x)
gp = ∇(p)
@test px(a,b) ≈ (2+10a)*p.py(b) 
@test gp[1](x) ≈ px(x)
@test gp[2](x) ≈ -p.px(x[1])
@test typeof(gp)<:PolyVectorField

@test_throws ArgumentError("z is not an indeterminate of the polynomial") derivative(p,:z)

Lp = Δ(p)
@test typeof(Lp) == typeof(p)
@test Lp(a,b) ≈ 10p.py(b)

Lq = Δ(q)
Lr = Δ(r)
@test Lr(x) ≈ Lp(x) + Lq(x)
@test Δ(Lr) == zero(p)
@test repr(r) == "(1.0 + 2.0*x + 5.0*x^2)(2.0 - 1.0*y) + (1.0 + 2.0*x)(-2.0 + 1.0*y^2)"

