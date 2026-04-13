"""

    Quadrature

A numerical quadrature scheme. 2D-schemes can be created by calling `gmquadrature`. However,`TrihpFEM` uses precomputed schemes, for performance reasons. 1D-schemes are only available in precomputed version. 
"""
struct Quadrature{D,R,V<:AbstractVector,P<:Integer}
    weights::R
    points::V #check dimensions
    degree::P
end

"""
    _jump(i)
    _jump(i,dim)
A very simple function that computes recursively the jump between the number of points needed for a quadrature of degree `i`, with respect to the quadrature of degree `i-1`. Note that `i` is not the actual degree but `i=deg÷2`.
"""
@inline _jump(i) = i
@inline _jump(dim,i) = dim == 1 ? _jump(i) : sum(1:_jump(dim-1,i))


"""
   numberofpoints(deg,dim)
Computes the number of points for a quadrature of odd degree `deg`. It works for `dim=2,3`.
"""
numberofpoints(dim,deg)= _numberofpoints(dim,1+deg÷2)
@inline _numberofpoints(dim,deg) =  deg == 1 ? 1 : _numberofpoints(dim,deg-1)+_jump(dim,deg)


"""
    gmquadrature(::Val{D}, degree,Tref)

Builds a `Quadrature` for dimension `D` over the simplex `Tref`. If the last argument is omitted, a reference triangle with vertices `(-1,-1),(-1,1),(1,1)` is used.

# Arguments
- `D`: dimension
- `degree`: desired polynomial degree of accuracy (must be odd)
- `Tref` (optional): simplex where the quadrature is located. A matrix of `D⨱(D+1)` with the vertices.

# Output
- A `Quadrature` of the desired degree.
"""
function tref(::Type{F}) where {F}
    return SMatrix{2,3,F}([-1 1 1;-1 -1 1])
end

function gmquadrature(::Val{D},degree::P,Tref) where {D,P<:Integer}
    T = eltype(Tref)
    L = numberofpoints(D,degree)
    _gmquadrature(T,Val(D),Val(L),degree,Tref)
end
function gmquadrature(::Type{F},d::Val{D},degree::P) where{F,D,P}
    Tref = tref(F)
    gmquadrature(d,degree,Tref)
end
gmquadrature(d::Val{D},degree::P) where {D,P} = gmquadrature(Float64,d,degree)
    



_szero(::Val{D},::Type{T}) where {D,T} = MVector{D,T}(zero(T) for _ in 1:D)
_szero(::Val{1},::Type{T}) where T = zero(T)
function _sszero(::Val{D},::Val{L},::Type{T}) where {D,L,T}
    collect_as(FixedSizeVector,(_szero(Val(D),T) for _ in 1:L))
end

"""
    DegCombination{L,D}

An iterator for building all combinations of integers in tuples of length `L` and sum `D`. For example `DegCombination{3,2}()` iterates over `(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)`. These combinations are necessary for building two dimensional Grundmann-Moeller schemes.
The iterator allows lazy creation which avoids the allocations incurred in `GrundmannMoeller.jl`.
"""
struct DegCombination{L,D} end

Base.iterate(::DegCombination{1,0}) = nothing
Base.iterate(::DegCombination{1,0},st) = nothing

# Base.iterate(::DegCombination{L,0}) where L= SVector{L,Int}(zeros(Int,L)) 

function Base.iterate(::DegCombination{L,D}) where {L,D}
    ini = D
    tail = SVector{L-1,Int}(zeros(Int,L-1))
    t = SVector{L,Int}(ini,tail...)
    return t,(ini,tail)
end

function Base.iterate(e::DegCombination{L,D},st) where {L,D}
    (_,tail) = st
    if tail[end]==D
        return nothing
    else
        _iterate(e,st)
    end
end

function _iterate(::DegCombination{L,D},st) where {L,D}
    (ini,tail) = st
    if length(tail)==1
        ini -= 1
        tail = SVector{L-1,Int}(tail[1]+1,tail[2:end]...)
        t = SVector{L,Int}(ini,tail...)
        return t,(ini,tail)
    else
        tailite = iterate(DegCombination{L-1,D-ini}(),(tail[1],tail[2:end]))
        if isnothing(tailite)
            ini = ini - 1
            tail = SVector{L-1,Int}(D-ini,zeros(Int,L-2)...)
            t = SVector{L,Int}(ini,tail...)
            return t,(ini,tail)
       else
            t = SVector{L,Int}(ini,tailite[1]...)
            return t,(ini,tailite[1])
       end
    end
end
    
function Base.length(::DegCombination{L,D}) where {L,D}
    (D==0 || L==1) ? 1 : sum(length(DegCombination{L-1,i}()) for i in 0:D)
end

function _gmquadrature(::Type{T}, ::Val{D}, ::Val{L}, degree::P,Tref) where {T,D,L,P}
    D::Int
    @assert degree ≥ 0
    @assert isodd(degree)
    order = (degree - 1) ÷ 2
    weights = FixedSizeVector{Float64}(undef,L)
    points = _sszero(Val(D),Val(L),T)
    j = 1
    for i in 0:order
        exponents = DegCombination{D+1,order-i}()
        w = T((-1)^i) * big(degree + D - 2 * i)^degree / (big(2)^(2 * order) *
             factorial(big(i)) *
             factorial(big(degree + D - i)))
        weights[j:j+length(exponents)-1] .= w
        for part in exponents
            points[j] .= Tref*(map(p -> T(2 * p + 1) / (degree + D - 2 * i), part))
            j += 1
        end
    end
    weights /= sum(weights)
    rw = RepVec(weights)
    F = eltype(first(points))
    return Quadrature{D,typeof(rw),typeof(points),P}(rw, points, degree)
end

