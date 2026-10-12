using ..TrihpFEM.Meshes
using ..TrihpFEM.Tensors
using ..TrihpFEM.PolyFields
using Test
using StaticArrays
using LinearAlgebra


p = ProductPoly((1.0, 2, 5.0), (2.0, -1))
a = rand()
b = rand()
x = Point2D(a, b)
y = [a, b]
@test p(x) == p(y)
@test p(Point2D(2, -1)) ≈ 75
@test (b * p)(x) ≈ b * p(x)
q = ProductPoly((1.0, 2.0), (-2, 0.0, 1.0))
@test q([2, -1]) ≈ -5.0
pq = p * q
@test pq(x) ≈ p(x) * q(x)
qxp = q.polys[1] * p
@test qxp(y) ≈ q.polys[1](a)*p(x)
pyq = q * p.polys[2]
@test pyq(x) ≈ p.polys[2](b) * q(x)
@test zero(p) == zero(q) == zero(typeof(p))
# @test q ⋅ p == q * p

v = Tensor{1,2}(p, q)
@test v(x)[1] == p(x)
@test v(x)[2] == q(x)
@test length(v) == 2length(p)
@test v(x) == [p(x), q(x)]
@test (a * v)(x) ≈ a * v(x) ≈ (v * a)(x)
@test v(x)[1] == v[1](x)
A = rand(2, 2)
w = A * v
@test typeof(w) <: PolyVectorField
@test size(w) == (2,)
row = [3 -2]
@test typeof(row * w) <: PolyScalarField
s = ProductPoly(Tuple(rand() for _ in 1:rand(1:8)), Tuple(rand() for _ in 1:rand(1:8)), :x, :y)
@test_throws DimensionMismatch() row' * w
sv = v * s
@test sv(x) ≈ v(x) * s(x)

r = p + q
@test typeof(r) <: PolySum
@test r(x) ≈ p(x) + q(x)
@test zero(p) == zero(p)
rs = r * s
@test rs(x) ≈ r(x) * s(x)
@test (r + s)(x) ≈ r(x) + s(x)
@test (r + rs)(x) ≈ r(x) + rs(x)

aff = AffineToRef{Float64}()
v₁ = SVector(rand() + 2, rand() + 1)
v₂ = SVector(rand(), rand() + 2)
v₃ = SVector(rand() + 1, rand())
affine!(aff, [v₁, v₂, v₃])
@test aff(SVector(-1, 1)) ≈ v₁
@test aff(SVector(-1, -1)) ≈ v₂
@test aff(SVector(1, -1)) ≈ v₃
@test area(aff) ≈ 0.5abs(v₁[1] * (v₂[2] - v₃[2]) + v₂[1] * (v₃[2] - v₁[2]) + v₃[1] * (v₁[2] - v₂[2]))

aff₂ = AffineToRef([√2 / 4 -√2 / 4;√2 / 4 √2 / 4], [0, √2 / 2])
@test area(aff₂) ≈ 1 / 2


f(x) = 2exp(x[1] + x[2])

wf2 = (w * f) * 2
z = SVector(-0.5, 0.0)
@test wf2(z, aff₂) ≈ 4 * w(z) * exp(√2 / 4)
z = SVector(0.0, 0.5)
wf2(z, aff₂) ≈ 4 * w(z) * exp(√2 / 2)

@test repr(p) == "(1.0 + 2.0*x + 5.0*x^2)(2.0 - 1.0*y)"
@test repr(zero(p)) == "(0.0,)"
@test sprint(show, p) == "(1.0 + 2.0*x + 5.0*x^2)(2.0 - 1.0*y)"

px = derivative(p, :x)
gp = ∇(p)
@test px(x) ≈ (2 + 10a) * p.polys[2](b)
@test gp[1](x) ≈ px(x)
@test gp[2](x) ≈ -p.polys[1](x[1])
@test typeof(gp) <: PolyVectorField
@test_throws ArgumentError("z is not an indeterminate of the polynomial") derivative(p, :z)
@test_throws ArgumentError("z is not an indeterminate of the polynomial") derivative(r, :z)

gpp = PolyVectorField([g for g in gp])
gpp(x) ≈ gp(x)
gpgp = gp + gp
gp2 = 2gp
pv = p * v
rsq = r * r
@test gpgp(x) ≈ gp2(x)
@test pv[1](x) ≈ p(x)^2
@test rsq(x) ≈ p(x)^2 + 2p(x) * q(x) + q(x)^2
@test repr(gp) == "ProductPoly{Float64, :x, :y}[(2.0 + 10.0*x)(2.0 - 1.0*y), (1.0 + 2.0*x + 5.0*x^2)(-1.0)]"

Lp = Δ(p)
@test typeof(Lp) == typeof(p)
@test Lp(a, b) ≈ 10p.polys[2](b)

Lq = Δ(q)
Lr = Δ(r)
@test Lr(x) ≈ Lp(x) + Lq(x)
@test Δ(Lr) == zero(p)
@test repr(r) == "(1.0 + 2.0*x + 5.0*x^2)(2.0 - 1.0*y) + (1.0 + 2.0*x)(-2.0 + 1.0*y^2)"

voutgp = v ⊗ gp
@test typeof(voutgp) <: PolyMatrixField
@test size(voutgp) == (2, 2)
@test voutgp[1, 2](x) ≈ p(x) * gp[2](x)

dvx = derivative(v,:x)
@test dvx(x)[1] ≈ derivative(p,:x)(x)
