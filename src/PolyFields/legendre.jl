abstract type AbstractBasis{D} end


"""
```
    LegendreIterator(N)
```
Builds an iterator over de Legendre polynomials.
```
    julia> for p in LegendreIterator(4)
           println(p)
           end
    1.0
    1.0*x
    -0.5 + 1.5*x^2
    -1.5*x + 2.5*x^3
```
"""
struct LegendreIterator{P <: Integer, F <: Number, X} <: AbstractBasis{1}
    N::P
    LegendreIterator{P, F, X}(N) where {P, F, X} = new{P, F, X}(P(N))
end
LegendreIterator(N::P) where {P <: Integer} = LegendreIterator{P, Float64, :x}(N)
LegendreIterator(t::NTuple{1,P}) where {P<:Integer} = LegendreIterator(t[1])

Base.IteratorSize(::Type{<:LegendreIterator}) = Base.HasLength()
Base.length(l::LegendreIterator{P, F, X}) where {P, F, X} = l.N + 1

function Base.iterate(::LegendreIterator{P, F, X}) where {P, F, X}
    q = one(ImmutablePolynomial{F, X})
    z = zero(ImmutablePolynomial{F, X})
    return q, (0, q, z)
end

function Base.iterate(l::LegendreIterator{P, F, X}, state) where {P, F, X}
    if state[1] == l.N
        return nothing
    else
        _iterate(l, state)
    end
end

function _iterate(::LegendreIterator{P, F, X}, state) where {P, F, X}
    n, p, pm = state
    if n == 0
        q = ImmutablePolynomial((zero(F), one(F)), X)
        return q, (1, q, p)
    else
        q = ((2n + 1)ImmutablePolynomial((zero(F), one(F)), X) * p - n * pm) / (n + 1)
        return q, (n + 1, q, p)
    end
end

###########################
#     STANDARD BASIS
###########################


struct StandardBasis{P <: Integer, F <: Number} <: AbstractBasis{2}
    degs::NTuple{3, P}
    Lx::LegendreIterator{P, F, VARIABLE_NAMES[1]}
    Ly::LegendreIterator{P, F, VARIABLE_NAMES[2]}
    function StandardBasis{P, F}(p) where {P, F}
        p = P.(p)
        p[1] + p[2] >= p[3] || throw(ArgumentError("Degrees does not satisfy p conformity."))
        Lx = LegendreIterator{P, F, VARIABLE_NAMES[1]}(p[1])
        Ly = LegendreIterator{P, F, VARIABLE_NAMES[2]}(p[2])
        return new{P, F}(p, Lx, Ly)
    end
end
StandardBasis(p::NTuple{3, P}) where P = StandardBasis{P, Float64}(p)
StandardBasis(p₁::P, p₂::P, p₃::P) where P = StandardBasis((p₁, p₂, p₃))
Base.IteratorSize(::Type{<:StandardBasis}) = Base.HasLength()
Base.length(sb::StandardBasis) = sum(min(sb.degs[2], sb.degs[3] - j) for j in 0:sb.degs[1]) + sb.degs[1] + 1

function Base.iterate(sb::StandardBasis{P, F}) where {P, F}
    (; Lx, Ly) = sb
    px, stx = iterate(Lx)
    py, sty = iterate(Ly)
    return ProductPoly(px, py), (stx, sty)
end

function Base.iterate(sb::StandardBasis{P, F}, state) where {P, F}
    (; degs, Lx, Ly) = sb
    (p₁, p₂, p₃) = degs
    stx, sty = state
    if (stx[1] == p₁) && (sty[1] == min(p₃ - p₁, p₂))
        return nothing
    elseif sty[1] < min(p₃ - stx[1], p₂)
        return _iteratey(Ly, sty, stx)
    else
        return _iteratex(Lx, Ly, stx)
    end
end

function _iteratex(Lx, Ly, stx)
    py, stynew = iterate(Ly)
    px, stxnew = iterate(Lx, stx)
    return ProductPoly(px, py), (stxnew, stynew)
end

function _iteratey(Ly, sty, stx)
    py, stynew = iterate(Ly, sty)
    return ProductPoly(stx[2], py), (stx, stynew)
end
