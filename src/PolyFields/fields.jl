###########################################################################
#################              SCALAR FIELDS              #################
###########################################################################
abstract type PolyScalarField{D,F<:Number}  end

coefficienttype(p::PolyScalarField{D,F}) where {D,F} = F
# indeterminate(::AbstractPolynomial{T, X}) where {T, X} = X


Base.length(::PolyScalarField) = 1
Base.iterate(t::PolyScalarField) = (t, nothing)
Base.iterate(::PolyScalarField, st) = nothing


###########################################################################
#################              ProductPoly                #################
###########################################################################

"""

    ProductPoly{F,X,Y} <: PolySacalarField{F,X,Y}
    function PolyTensorField{D}(polys...) where D
        N = round(Int,log(length(polys)/))
    end

A bivariate polynomial with coefficients of type `F` defined by the product of two univariate polynomials, one on variable `X` and the other on variable `Y`.

For construction it is preferred to pass two tuples of coefficients.
# Examples
```
    julia> ϕ = ProductPoly((2,1.,3.),(3.,2.))
        (2.0 + 1.0*x + 3.0*x^2)(3.0 + 2.0*y)
    julia> ϕ(1,-1)
        6.0
    julia> ϕ([1,-1])
        6.0
```
If necessary, the indeterminates can be specified:
```
    julia> ϕ = ProductPoly((2,1.,3.),(3.,2.),:z,:ξ)
        (2.0 + 1.0*z + 3.0*z^2)(3.0 + 2.0*ξ)
```
"""
_coeffs(p::AbstractPolynomial) = coeffs(p)
_coeffs(p::Tuple) = p

struct ProductPoly{D,F} <: PolyScalarField{D,F}
    polys::NTuple{D,ImmutablePolynomial{F}}
    function ProductPoly{F}(pols...) where F
        D = length(pols)
        data = tuple((ImmutablePolynomial{F,VARIABLE_NAMES[i]}(F.(_coeffs(pols[i]))) for i in 1:D)...)
        if any(p==zero(p) for p in data)
            ps = tuple((ImmutablePolynomial{F,x}(zero(F)) for x in VARIABLE_NAMES[1:D])...)
            return new{D,F}(ps)
        else
            return new{D,F}(data)
        end
    end
end
function ProductPoly(pols...)
    pols = tuple((promote(p...) for p in pols)...)
    F = promote_type((eltype(p) for p in pols)...)
    ProductPoly{F}(pols...)
end

ProductPoly{D,F}(p::ProductPoly{D}) where {D,F} = ProductPoly{F}(p.polys...)
_apply(s,v) = s(v)

(s::ProductPoly)(x::AbstractVector) = mapreduce(_apply,*,s.polys,x)

degs(p::ImmutablePolynomial) = tuple(length(p.coeffs)-1,)
degs(p::ProductPoly) = tuple((length(r) - 1 for r in p.polys)...)

function Base.promote(t::ProductPoly{D,F}, s::ProductPoly{D,T}) where {D,F,T}
    R = promote_type(F,T)
    ProductPoly{D,R}(t), ProductPoly{D,R}(s)
end

Base.:*(a::Number, p::ProductPoly) = ProductPoly((a*q for q in p.polys)...)
Base.:*(p::ProductPoly, a::Number) = a * p

Tensors.otimes(n::Number,p::PolyScalarField) = n*p
Tensors.otimes(p::PolyScalarField,n::Number) = n*p

function Base.:*(p::AbstractPolynomial{T,X}, q::ProductPoly{F}) where {T,F,X}
    i = findfirst(X.==VARIABLE_NAMES)
    data = tuple((j==i ? qq*p : qq for (j,qq) in enumerate(q.polys)))
end


Base.:*(p::ProductPoly{D}, q::ProductPoly{D}) where D = ProductPoly(map(*,p.polys,q.polys)...)

Base.zero(p::ProductPoly{D,F}) where {D,F}= ProductPoly{F}((zero(pp) for pp in p.polys)...)

Base.zero(::Type{ProductPoly{D,F}}) where {D,F} = ProductPoly((zero(F),) for _ in 1:D)

Base.one(p::ProductPoly{D,F}) where {D,F} = ProductPoly{F}((one(pp) for pp in p.polys)...)

Base.one(::Type{ProductPoly{D,F}}) where {D,F} = ProductPoly{F}(((one(F),) for _ in VARIABLE_NAMES[1:D])...)

Base.convert(::Type{T}, x::N) where {T <: PolyScalarField, N <: Number} = x * one(T)

###########################################################################
#################                PolySum                  #################
###########################################################################

"""
```
   PolySum{F,X,Y} <: PolyScalarField{F,X,Y}
```
A struct for storing the sum of two `PolyScalarField`s. A typical application is the divergence of a `PolyVectorField`, usually formed by the sum of two `ProductPoly`s. Note that the terms of the sum can be any type of `PolyScalarField`, including `PolySum`s.
"""
struct PolySum{D,F} <: PolyScalarField{D,F}
    left::PolyScalarField{D,F}
    right::PolyScalarField{D,F}
end

(p::PolySum)(x...) = p.left(x)+p.right(x)

function Base.:+(p::P, q::Q) where {D,F,P <: PolyScalarField{D,F}, Q <: PolyScalarField{D,F}}
    return if p == zero(p)
        q
    elseif q == zero(q)
        p
    else
        PolySum(p, q)
    end
end

Base.zero(::PolySum{D,F}) where {D,F} = zero(ProductPoly{D,F})


#LinearAlgebra.dot(p::PolyVectorField{F,X,Y},q::PolyVectorField{F,X,Y})  where {F,X,Y} = p.s1*q.s1 + p.s2*q.s2

Base.:*(n::Number, ps::PolySum) = n * ps.left + n * ps.right
Base.:*(ps::PolySum, n::Number) = n * ps
Base.:*(ps::PolySum, p::ProductPoly) = ps.left * p + ps.right * p
Base.:*(p::ProductPoly, ps::PolySum) = ps * p


Base.promote_rule(::Type{ProductPoly{D,F}},::Type{PolySum{D,F}}) where {D,F} = PolySum{D,F}
Base.convert(::Type{PolySum{D,F}},t::ProductPoly{2,T}) where {D,F,T} = PolySum(t,zero(t))


### These functions are defined in order to broadcast properly. Essentially: a PolyScalarField shoul behave like a number under broadcasting: f.(p) when p isa PolyScalarField is equivalent to f(p)
Broadcast.broadcastable(p::PolyScalarField) = p
Base.ndims(::Type{T}) where T<:PolyScalarField = 0
Base.ndims(::PolyScalarField) = 0
Base.size(::PolyScalarField) = ()
Base.getindex(p::PolyScalarField,::CartesianIndex{0}) = p
#------VER ESTO

# Base.:*(ps::PolySum, p::PolyVectorField) = PolyVectorField([ps * p.s1, ps * p.s2])
# Base.:*(p::PolyVectorField, ps::PolySum) = ps * p
# Base.:*(ps::PolySum, qs::PolySum) = ps.left * qs.left + ps.right * qs.left + ps.right * qs.left + ps.right * qs.right

# Base.:+(p::PolyTensorField, q::PolyTensorField) = PolyTensorField(p .+ q)

# (o::GeneralField)(x, t::AffineToRef) = o.op((evaluate(evaltype(arg), arg, t, x) for arg in o.args)...)
