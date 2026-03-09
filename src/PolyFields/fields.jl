abstract type PolyField{F, X, Y} end

indeterminate(::AbstractPolynomial{T, X}) where {T, X} = X
"""

    indeterminates(p::PolyField)
    
returns the indeterminates of `p`, for example `:x` and `:y`.
"""
indeterminates(::PolyField{F, X, Y}) where {F, X, Y} = (X, Y)

###########################################################################
#################              SCALAR FIELDS              #################
###########################################################################
#
abstract type PolyScalarField{F, X, Y} <: PolyField{F, X, Y} end
Base.length(::PolyScalarField) = 1

Base.iterate(t::PolyScalarField) = (t, nothing)
Base.iterate(::PolyScalarField, st) = nothing
LinearAlgebra.:⋅(p::PolyScalarField, q::PolyScalarField) = p * q
#
# _zerofill(t::NTuple{K,T},N) where {K,T} = K>=N ? t : (t...,zeros(T,N-K)...)

######################################
####      BiPoly     ####
######################################
"""

    BiPoly{F,X,Y} <: PolySacalarField{F,X,Y}

A bivariate polynomial with coefficients of type `F` defined by the product of two univariate polynomials, one on variable `X` and the other on variable `Y`.

For construction it is preferred to pass two tuples of coefficients.
# Examples
```
    julia> ϕ = BiPoly((2,1.,3.),(3.,2.))
        (2.0 + 1.0*x + 3.0*x^2)(3.0 + 2.0*y)
    julia> ϕ(1,-1)
        6.0
    julia> ϕ([1,-1])
        6.0
```
If necessary, the indeterminates can be specified:
```
    julia> ϕ = BiPoly((2,1.,3.),(3.,2.),:z,:ξ)
        (2.0 + 1.0*z + 3.0*z^2)(3.0 + 2.0*ξ)
```
"""
struct BiPoly{F, X, Y} <: PolyScalarField{F, X, Y}
    px::ImmutablePolynomial{F, X, N} where {N}
    py::ImmutablePolynomial{F, Y, M} where {M}
    function BiPoly(
            px::ImmutablePolynomial{F},
            py::ImmutablePolynomial{F}, X, Y
        ) where {F}
        X == Y && throw(ArgumentError("Indeterminates must be different"))
        if px == zero(px) || py == zero(py)
            return new{F, X, Y}(zero(px), zero(py))
        else
            return new{F, X, Y}(px, py)
        end
    end
end

function BiPoly(t1::Tuple, t2::Tuple, X = :x, Y = :y)
    t1 = promote(t1...)
    t2 = promote(t2...)
    F = promote_type(eltype(t1), eltype(t2))
    N = length(t1)
    M = length(t2)
    if N == 0 || M == 0
        return zero(BiPoly{F, X, Y})
    else
        p1 = ImmutablePolynomial(NTuple{N, F}(convert.(F, t1)), X)
        p2 = ImmutablePolynomial(NTuple{M, F}(convert.(F, t2)), Y)
        return BiPoly(p1, p2, X, Y)
    end
end

(s::BiPoly)(x, y) = s.px(x) * s.py(y)
(s::BiPoly)(x::T) where {T <: AbstractVector} = s.px(x[1]) * s.py(x[2])

degs(p::BiPoly) = (length(p.px) - 1, length(p.py) - 1)

Base.promote(t::BiPoly{F, X, Y}, s::BiPoly{F, X, Y}) where {F, X, Y} = (t, s)
Base.:*(a::Number, p::BiPoly) = BiPoly(a * p.px, p.py, indeterminates(p)...)
Base.:*(p::BiPoly, a::Number) = a * p


function Base.:*(p::AbstractPolynomial, q::BiPoly)
    return if indeterminate(p) === indeterminate(q.px)
        BiPoly(p * q.px, q.py, indeterminates(q)...)
    elseif indeterminate(p) === indeterminate(q.py)
        BiPoly(q.px, p * q.py, indeterminates(q)...)
    else
        throw(ArgumentError("Indeterminates does not match."))
    end
end
Base.:*(q::BiPoly, p::AbstractPolynomial) = p * q

function Base.:*(p::BiPoly{F}, q::BiPoly{F}) where {F}
    indeterminates(p) == indeterminates(q) || throw(ArgumentError("Indeterminates does not match."))
    return BiPoly(p.px * q.px, p.py * q.py, indeterminates(p)...)
end

function Base.zero(p::BiPoly{F, X, Y}) where {F, X, Y}
    return BiPoly(zero(p.px), zero(p.py), X, Y)
end
function Base.zero(::Type{BiPoly{F, X, Y}}) where {F, X, Y}
    return BiPoly(zero(ImmutablePolynomial{F, X}), zero(ImmutablePolynomial{F, Y}), X, Y)
end
function Base.one(p::BiPoly{F, X, Y}) where {F, X, Y}
    return BiPoly(one(p.px), one(p.py), X, Y)
end
function Base.one(::Type{BiPoly{F, X, Y}}) where {F, X, Y}
    return BiPoly(one(ImmutablePolynomial{F, X}), one(ImmutablePolynomial{F, Y}), X, Y)
end

Base.convert(::Type{T}, x::N) where {F, X, Y, T <: PolyField{F, X, Y}, N <: Number} = x * one(T)
###############################
#       TENSOR FIELDS
###############################
"""

     PolyTensorField{F,X,Y,T<:PolyScalarField{F,X,Y},N} <: PolyField{F,X,Y}

  A struct for storing a tensor-field  of `N` dimensions formed by a `PolyScalarField` of type `T` with coeffitients of type `F` on the indeterminates `X` and `Y`.   There are useful aliases:
  - `PolyVectorField`, for vector fields, with `N=1` and
  - `PolyMatrixField`, for matrix fields, with `N=2`.

# Examples

```
  julia> p = BiPoly((1,3.,4.),(1,2.))
  julia> q = BiPoly((0,1.),(1.,))
  julia> v = PolyTensorField([p,q])  
```

`v` is `PolyVectorField` with `T` given by `BiPoly(Float64,:x,:y)`.

If the type of the component fields is not uniform, they are promoted to a common type.

```
  julia> p = BiPoly((1,3.,4.),(1,2.))
  julia> q = BiPoly((0,1.),(1.))
  julia> r = p+q # PolySum
  julia> v = PolyTensorField([p,r])
    2-element Vector{BiPoly{Float64, :x, :y}}:
      (1.0 + 3.0*x + 4.0*x^2)(1.0 + 2.0*y)
      (1.0*x)(1.0)

  julia> v(1,-2)
    2-element Vector{Float64}:
      -24.0
        1.0  
```
"""
struct PolyTensorField{F, X, Y, T <: PolyScalarField{F, X, Y}, N} <: PolyField{F, X, Y}
    tensor::FixedSizeArrayDefault{T, N}
    function PolyTensorField(arr::AbstractArray{T, N}) where {F, X, Y, T <: PolyScalarField{F, X, Y}, N}
        tensor = FixedSizeArrayDefault{T, N}(arr)
        return new{F, X, Y, T, N}(tensor)
    end
end
function PolyTensorField{T, N}() where {N, T <: PolyScalarField}
    a = FixedSizeArrayDefault{T, N}(undef, [2 for _ in 1:N]...)
    return PolyTensorField(a)
end

# PolyVectorField
const PolyVectorField{F, X, Y, T} = PolyTensorField{F, X, Y, T, 1}
PolyVectorField(x::AbstractArray{T, 1}) where {T} = PolyTensorField(x)

#PolyMatrixField
const PolyMatrixField{F, X, Y, T} = PolyTensorField{F, X, Y, T, 2}
PolyMatrixField(x::AbstractArray{T, 2}) where {T} = PolyTensorField(x)
function (v::PolyTensorField{F, X, Y, T, N})(x, y) where {F, X, Y, T, N}
    z = FixedSizeArrayDefault{F, N}(undef, size(v.tensor)...)
    for i in eachindex(v.tensor)
        z[i] = v.tensor[i](x, y)
    end
    return z
end

function (v::PolyTensorField{F, X, Y, T, N})(x) where {F, X, Y, T, N}
    return v(x[1], x[2])
end

Base.IteratorSize(::PolyTensorField{F, X, Y, T, N}) where {F, X, Y, T, N} = Base.HasShape{N}()
Base.length(p::T) where {T <: PolyTensorField} = length(p.tensor)
Base.size(p::T) where {T <: PolyTensorField} = size(p.tensor)
function Base.iterate(p::T, st = nothing) where {T <: PolyTensorField}
    return isnothing(st) ? iterate(p.tensor) : iterate(p.tensor, st)
end
Base.getindex(p::T, i...) where {T <: PolyTensorField} = getindex(p.tensor, i...)
Base.zero(p::T) where {T <: PolyTensorField} = PolyTensorField(zero(p.tensor))

Base.:*(a::N, p::T) where {N <: Number, T <: PolyTensorField} = PolyTensorField(a * p.tensor)
Base.:*(p::T, a::N) where {N <: Number, T <: PolyTensorField} = a * p

# This function needs improvements to avoid the type instability.
function Base.:*(A::AbstractArray, p::PolyTensorField{F, X, Y, T, N}) where {F, X, Y, T, N}
    sA = size(A); sp = size(p)
    sA[end] == sp[1] || throw(DimensionMismatch())
    M = A * p.tensor
    if length(M) > 1
        c = Array{PolyScalarField{F, X, Y}, length(size(M))}(undef, size(M)...)
        c .= A * p.tensor
        return PolyTensorField(c)
    else
        return M[1]
    end
end

function Base.:*(p::PolyScalarField, v::T) where {T <: PolyTensorField}
    issubset(indeterminates(p), indeterminates(v)) || throw(ArgumentError("Fields have different indeterminates"))
    w = similar(v.tensor)
    for i in eachindex(v.tensor)
        w[i] = p * v.tensor[i]
    end
    return T(w)
end
Base.:*(v::T, p::PolyScalarField) where {T <: PolyTensorField} = p * v


Base.promote_type(T::Type{<:Number}, P::Type{<:PolyField}) = P

# """
#    outerprod(v...)
#    outerproduct(v::PolyField,w::PolyField)
# Outer product of tensors. Can be used as a binary operator with ⊗. 
# """
# function outerprod(v...)
#     dims = tuple(Iterators.flatten(size.(v))...)
#     T = promote_type(eltype.(v)...)
#     z = FixedSizeArrayDefault{T, length(dims)}(undef, dims...)
#     for (i, k) in enumerate(Iterators.product(v...))
#         z[i] = prod(k)
#     end
#     return z
# end
# const ⊗ = outerprod

# contract(x, y) = sum(x .* y)

# outerprod(v::PolyTensorField, w::PolyTensorField) = PolyTensorField(v.tensor ⊗ w.tensor)
# outerprod(v, w::PolyScalarField) = v * w
# outerprod(v::AbstractArray, w::PolyTensorField) = PolyTensorField(v ⊗ w.tensor)


# function _outer(p::PolyVectorField, q::PolyVectorField)
#     @cast a[i, j] := p.tensor[i] * q.tensor[j]
#     return PolyTensorField(a)
# end

LinearAlgebra.dot(p::PolyVectorField, q::PolyVectorField) = outerprod(p.tensor, q.tensor)


###############################
#           POLYSUM
###############################

"""
```
   PolySum{F,X,Y} <: PolyScalarField{F,X,Y}
```
A struct for storing the sum of two `PolyScalarField`s. A typical application is the divergence of a `PolyVectorField`, usually formed by the sum of two `BiPoly`s. Note that the terms of the sum can be any type of `PolyScalarField`, including `PolySum`s.
"""
struct PolySum{F, X, Y} <: PolyScalarField{F, X, Y}
    left::PolyScalarField{F, X, Y}
    right::PolyScalarField{F, X, Y}
end

(p::PolySum)(x, y) = p.left(x, y) + p.right(x, y)
(p::PolySum)(x) = p.left(x) + p.right(x)


function Base.:+(p::P, q::Q) where {F, X, Y, P <: PolyField{F, X, Y}, Q <: PolyField{F, X, Y}}
    return if p == zero(p)
        q
    elseif q == zero(q)
        p
    else
        PolySum(p, q)
    end
end

function Base.zero(::PolySum{F, X, Y}) where {F, X, Y}
    return zero(BiPoly{F, X, Y})
end
#LinearAlgebra.dot(p::PolyVectorField{F,X,Y},q::PolyVectorField{F,X,Y})  where {F,X,Y} = p.s1*q.s1 + p.s2*q.s2

Base.:*(n::Number, ps::PolySum) = n * ps.left + n * ps.right
Base.:*(ps::PolySum, n::Number) = n * ps
Base.:*(ps::PolySum, p::BiPoly) = ps.left * p + ps.right * p
Base.:*(p::BiPoly, ps::PolySum) = ps * p
Base.:*(ps::PolySum, p::PolyVectorField) = PolyVectorField([ps * p.s1, ps * p.s2])
Base.:*(p::PolyVectorField, ps::PolySum) = ps * p
Base.:*(ps::PolySum, qs::PolySum) = ps.left * qs.left + ps.right * qs.left + ps.right * qs.left + ps.right * qs.right

Base.:+(p::PolyTensorField, q::PolyTensorField) = PolyTensorField(p .+ q)

###########################################################################################
#########################           AffineToRef        ###########################
###########################################################################################
"""

    AffineToRef{F}

A struct for defining and updating an affine transformation from the reference triangle to some other triangle. 
"""
struct AffineToRef{F <: Number}
    A::Tensor{2, 2, F}
    b::Tensor{1, 2, F}
    function AffineToRef{F}(vert) where {F <: Number}
        A = affinetoref_matrix(F, vert)
        b = affinetoref_vec(F, vert)
        return new{F}(A, b)
    end
end
function affinetoref_matrix(::Type{F}, vert) where {F}
    return Tensor{2, 2, F}((i, j) -> vert[mod1(j, 2)][i] - vert[3][i])
end
function affinetoref_vec(::Type{F}, vert) where {F}
    return Tensor{1, 2, F}(i -> (vert[1][i] + vert[3][i]) / 2)
end

(aff::AffineToRef{F})(x) where {F} = aff.A * x + aff.b


jac(aff::AffineToRef) = det(aff.A)
# area(x,y,z) = 0.5abs(x[1]*(y[2]-z[2])+y[1]*(z[2]-x[2])+z[1]*(x[2]-y[2]))
# area(v::Vector) = area(v...)
area(t::AffineToRef) = 2jac(t)


##########################################################################################
###################                GENERALFIELD                 ##########################
##########################################################################################
"""

    GeneralField
    
A `struct` for defining fields that mix user defined functions with `PolyField`s.

`GeneralField`s are necessary for properly distinguish the method of integration.

Consider, for example, the form `a(u,v) = ∫(∇(u)⋅∇(v)*dΩ`. In order to assembly the associated matrix, local terms can be fully pre-computed on the reference triangle, and then just transformed into each mesh triangle. Note that only **one** local tensor is computed for each combination of degrees, so this precomputation is very fast when certain combination of degrees are repeated along the mesh. When applied to elements of a basis, `∇(u)⋅∇(v)` produces a `PolyScalarField` that can be directly integrated. 

On the other hand, the precomputation cannot be carried out in the same way for `b(v) = ∫(f*v)*dΩ` with `f` a function, since `f` takes different values on each mesh triangle. Hence `f*v` produces a `GeneralField` which integrates differently. In particular, it evaluates differently for integration.

"""
struct GeneralField <: Function
    op
    args::Tuple
end


Base.:*(f::Function, p::PolyField) = GeneralField(*, (f, p))
Base.:*(p::PolyField, f::Function) = f * p
LinearAlgebra.dot(f::Function, p::PolyField) = GeneralField(dot, (f, p))
LinearAlgebra.dot(p::PolyField, f::Function) = dot(f, p)

function Base.:*(n::Number, g::GeneralField)
    ind = findfirst(isa.(g.args, PolyField))
    if !isnothing(ind)
        newargs = (g.args[1:(ind - 1)]..., n * g.args[ind], g.args[(ind + 1):end]...)
        return GeneralField(g.op, newargs)
    else
        return GeneralField(*, (n, g))
    end
end
Base.:*(g::GeneralField, n::Number) = n * g

# A trait for evaluation of Field
abstract type EvalType end
struct Eval <: EvalType end
struct Compose <: EvalType end
struct Pass <: EvalType end

evaltype(_) = Pass()
evaltype(::Function) = Compose()
evaltype(::PolyField) = Eval()

evaluate(::Eval, f, t::AffineToRef, x) = f(x)
evaluate(::Compose, f, t::AffineToRef, x) = f(t(x))
evaluate(::Pass, f, t, x) = f

(o::GeneralField)(x, t::AffineToRef) = o.op((evaluate(evaltype(arg), arg, t, x) for arg in o.args)...)
