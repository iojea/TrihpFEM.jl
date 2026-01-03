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
struct LegendreIterator{I<:Integer,F<:Number,X}
    N::I
    LegendreIterator{I,F,X}(N) where {I,F,X} = new{I,F,X}(I(N))
end
LegendreIterator(N::I) where {I<:Integer} = LegendreIterator{I,Float64,:x}(N)
Base.IteratorSize(::Type{<:LegendreIterator}) = Base.HasLength()
Base.length(l::LegendreIterator{I,F,X}) where {I,F,X} = l.N+1

function Base.iterate(::LegendreIterator{I,F,X}) where {I,F,X}
    q = one(ImmutablePolynomial{F,X})
    z = zero(ImmutablePolynomial{F,X})
    q,(0,q,z)
end

function Base.iterate(l::LegendreIterator{I,F,X},state) where {I,F,X}
    if state[1] == l.N
        return nothing
    else
        _iterate(l,state)
    end
end

function _iterate(::LegendreIterator{I,F,X},state) where {I,F,X}
    n,p,pm = state
    if n==0
        q = ImmutablePolynomial((zero(F),one(F)),X)
        return q,(1,q,p)
    else
        q = ((2n+1)ImmutablePolynomial((zero(F),one(F)),X)*p - n*pm)/(n+1)
        return q,(n+1,q,p)
    end
end

###########################
#     STANDARD BASIS
###########################


struct StandardBasis{P<:Integer,F<:Number,X,Y}
    degs::NTuple{3,P}
    Lx::LegendreIterator{P,F,X}
    Ly::LegendreIterator{P,F,Y}
    function StandardBasis{P,F,X,Y}(p) where {P,F,X,Y}
        p = P.(p)
        p[1]+p[2] >= p[3] || throw(ArgumentError("Degrees does not satisfy p conformity."))
        Lx = LegendreIterator{P,F,X}(p[1])
        Ly = LegendreIterator{P,F,Y}(p[2])
        new{P,F,X,Y}(p,Lx,Ly)
    end
end
StandardBasis(p::NTuple{3,P}) where P = StandardBasis{P,Float64,:x,:y}(p)
StandardBasis(p₁::P,p₂::P,p₃::P) where P = StandardBasis((p₁,p₂,p₃))
Base.IteratorSize(::Type{<:StandardBasis}) = Base.HasLength()
Base.length(sb::StandardBasis) = sum(min(sb.degs[2],sb.degs[3]-j) for j in 0:sb.degs[1]) +sb.degs[1] + 1

function Base.iterate(sb::StandardBasis{P,F,X,Y}) where {P,F,X,Y}
    (;Lx,Ly) = sb
    px,stx = iterate(Lx)
    py,sty = iterate(Ly)
    BiPoly(px,py),(stx,sty)
end

function Base.iterate(sb::StandardBasis{P,F,X,Y},state) where {P,F,X,Y}
    (;degs,Lx,Ly) = sb
    (p₁,p₂,p₃) = degs
    stx,sty = state
    if (stx[1] == p₁) && (sty[1] == min(p₃-p₁,p₂))
        return nothing
    elseif sty[1] < min(p₃-stx[1],p₂)
        return _iteratey(Ly,sty,stx)
    else
        return _iteratex(Lx,Ly,stx)
    end
end

function _iteratex(Lx,Ly,stx)
    py,stynew = iterate(Ly) 
    px,stxnew = iterate(Lx,stx)
    return BiPoly(px,py),(stxnew,stynew)
end

function _iteratey(Ly,sty,stx)
    py,stynew = iterate(Ly,sty)
    return BiPoly(stx[2],py),(stx,stynew)
end

### The goal of DoubleStandardBasis is to have a double iterator in order to have only one function for building matrices or vectors. 
struct DoubleStandardBasis{P<:Integer,F<:Number,X,Y}
    degs::NTuple{3,P}
    Bᵢ::StandardBasis{P,F,X,Y}
    Bⱼ::StandardBasis{P,F,X,Y}
    function DoubleStandardBasis{P,F,X,Y}(p) where {P,F,X,Y}
        p = P.(p)
        Bᵢ = StandardBasis{P,F,X,Y}(p)
        Bⱼ = StandardBasis{P,F,X,Y}(p)
        new{P,F,X,Y}(p,Bᵢ,Bⱼ)
    end
end
DoubleStandardBasis(p::NTuple{3,P}) where P = DoubleStandardBasis{P,Float64,:x,:y}(p)
DoubleStandardBasis(p₁::P,p₂::P,p₃::P) where P = DoubleStandardBasis((p₁,p₂,p₃))

Base.IteratorSize(::Type{<:DoubleStandardBasis}) = Base.HasShape()
Base.size(dsb::DoubleStandardBasis) = (length(dsb.Bᵢ),length(dsb.Bⱼ))

function Base.iterate(dsb::DoubleStandardBasis{P,F,X,Y}) where {P,F,X,Y}
    (;Bᵢ,Bⱼ) = dsb
    bᵢ,sbᵢ = iterate(Bᵢ)
    bⱼ,sbⱼ = iterate(Bⱼ)
    (bᵢ,bⱼ),(sbᵢ,sbⱼ)
end

recover_actual_poly(statesb) = BiPoly(statesb[1][2],statesb[2][2])
    

function Base.iterate(dsb::DoubleStandardBasis{P,F,X,Y},state) where {P,F,X,Y}
    (;Bᵢ,Bⱼ) = dsb
    sbᵢ,sbⱼ = state
    outᵢ = iterate(Bᵢ,sbᵢ)
    if isnothing(outᵢ)
        bᵢ,sbᵢ = iterate(Bᵢ)
        outⱼ = iterate(Bⱼ,sbⱼ)
        if isnothing(outⱼ)
            return nothing
        else
            bⱼ,sbⱼ = outⱼ
            return (bᵢ,bⱼ),(sbᵢ,sbⱼ)
        end
    else
        bᵢ,sbᵢ = outᵢ
        bⱼ = recover_actual_poly(sbⱼ)
        return (bᵢ,bⱼ),(sbᵢ,sbⱼ)
    end
end


