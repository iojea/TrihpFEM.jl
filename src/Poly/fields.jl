abstract type PolyField{F,X,Y} end

indeterminate(::AbstractPolynomial{T,X}) where {T,X} = X
"""

    indeterminates(p::PolyField)
    
returns the indeterminates of `p`, for example `:x` and `:y`.
"""
indeterminates(::PolyField{F,X,Y}) where {F,X,Y} = (X,Y)

###########################################################################
#################              SCALAR FIELDS              #################
###########################################################################
#
abstract type PolyScalarField{F,X,Y} <: PolyField{F,X,Y} end
Base.length(::PolyScalarField) = 1
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
struct BiPoly{F,X,Y} <: PolyScalarField{F,X,Y}
    px::ImmutablePolynomial{F,X,N} where N
    py::ImmutablePolynomial{F,Y,M} where M
    function BiPoly(px::ImmutablePolynomial{F},
                    py::ImmutablePolynomial{F},X,Y) where F
        X==Y && throw(ArgumentError("Indeterminates must be different"))
        if px == zero(px) || py == zero(py)
            return new{F,X,Y}(zero(px),zero(py))
        else
            return new{F,X,Y}(px,py)
        end
    end
end

function BiPoly(t1::Tuple,t2::Tuple,X=:x,Y=:y)
    t1 = promote(t1...)
    t2 = promote(t2...)
    F = promote_type(eltype(t1),eltype(t2))
    N = length(t1)
    M = length(t2)
    if N == 0 || M == 0
        return zero(BiPoly{F,X,Y})
    else
        p1 = ImmutablePolynomial(NTuple{N,F}(convert.(F,t1)),X)
        p2 = ImmutablePolynomial(NTuple{M,F}(convert.(F,t2)),Y)
        return BiPoly(p1,p2,X,Y)
    end
end

(s::BiPoly)(x,y) = s.px(x)*s.py(y)
(s::BiPoly)(x::T) where T<:AbstractVector = s.px(x[1])*s.py(x[2])

degs(p::BiPoly) = (length(p.px)-1,length(p.py)-1)

Base.promote(t::BiPoly{F,X,Y},s::BiPoly{F,X,Y}) where {F,X,Y}= (t,s)
Base.:*(a::Number,p::BiPoly) = BiPoly(a*p.px,p.py)
Base.:*(p::BiPoly,a::Number) = a*p

Base.iterate(t::BiPoly) = t
Base.iterate(::BiPoly,st) = nothing

function Base.:*(p::AbstractPolynomial,q::BiPoly)
    if indeterminate(p)===indeterminate(q.px)
        BiPoly(p*q.px,q.py)
    elseif indeterminate(p)===indeterminate(q.py)
        BiPoly(q.px,p*q.py)
    else
        throw(ArgumentError("Indeterminates does not match."))
    end
end
Base.:*(q::BiPoly,p::AbstractPolynomial) = p*q

function Base.:*(p::BiPoly{F},q::BiPoly{F}) where {F}
    indeterminates(p)==indeterminates(q) || throw(ArgumentError("Indeterminates does not match."))
    BiPoly(p.px*q.px,p.py*q.py)
end

LinearAlgebra.:⋅(p::BiPoly,q::BiPoly) = p*q
function Base.zero(p::BiPoly{F,X,Y}) where {F,X,Y}
    BiPoly(zero(p.px),zero(p.py),X,Y)
end
function Base.zero(::Type{BiPoly{F,X,Y}}) where {F,X,Y}
    BiPoly(zero(ImmutablePolynomial{F,X}),zero(ImmutablePolynomial{F,Y}),X,Y)
end

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
struct PolyTensorField{F,X,Y,T<:PolyScalarField{F,X,Y},N} <: PolyField{F,X,Y}
    tensor::Array{T,N}
end
function PolyTensorField{T,N}() where {N,T<:PolyScalarField}
    a = Array{T,N}(undef,[2 for _ in 1:N]...)
    PolyTensorField(a)
end

# PolyVectorField
const PolyVectorField{F,X,Y,T} = PolyTensorField{F,X,Y,T,1}
PolyVectorField(x::AbstractArray{T,1}) where T = PolyTensorField(x)

#PolyMatrixField
const PolyMatrixField{F,X,Y,T} = PolyTensorField{F,X,Y,T,2}
PolyMatrixField(x::AbstractArray{T,2}) where T = PolyTensorField(x)
function (v::PolyTensorField{F,X,Y,T,N})(x,y) where {F,X,Y,T,N}
    z = Array{F,N}(undef,size(v.tensor)...)
    for i in eachindex(v.tensor)
        z[i] = v.tensor[i](x,y)
    end
    z
end


function (v::PolyTensorField{F,X,Y,T,N})(x) where {F,X,Y,T,N}
    v(x[1],x[2])
end

Base.IteratorSize(::PolyTensorField{F,X,Y,T,N}) where {F,X,Y,T,N} = Base.HasShape{N}()
Base.length(p::PolyTensorField) = length(p.tensor)
Base.size(p::PolyTensorField) = size(p.tensor)
Base.iterate(p::PolyTensorField,st=nothing) = iterate(p.tensor,st)
Base.getindex(p::PolyTensorField,i) = getindex(p.tensor,i)

Base.:*(a::Number,p::PolyTensorField) = PolyTensorField(a*p.tensor)
Base.:*(p::PolyTensorField,a::Number) = a*p

# This function needs improvements to avoid the type instability. 
function Base.:*(A::AbstractArray,p::PolyTensorField{F,X,Y,T,N}) where {F,X,Y,T,N}
    sA = size(A); sp = size(p)
    sA[end] == sp[1] || throw(DimensionMismatch())
    M = A*p.tensor
    if length(M)>1
        c = Array{PolyScalarField{F,X,Y},length(size(M))}(undef,size(M)...)
        c .= A*p.tensor
        return PolyTensorField(c)
    else
        return M[1]
    end
end

function Base.:*(p::PolyScalarField,v::T) where T<:PolyTensorField
    issubset(indeterminates(p),indeterminates(v)) || throw(ArgumentError("Fields have different indeterminates"))
    w = similar(v.tensor)
    for i in eachindex(v.tensor)
        w[i] = p*v.tensor[i]
    end
    T(w)
end
Base.:*(v::T,p::PolyScalarField) where T<:PolyTensorField = p*v


function _outer(p::PolyVectorField,q::PolyVectorField)
    @cast a[i,j] := p.tensor[i]*q.tensor[j]
    PolyTensorField(a)
end

LinearAlgebra.dot(p::PolyVectorField,q::PolyVectorField) = dot(p.tensor,q.tensor)


###############################
#           POLYSUM
###############################

"""
```
   PolySum{F,X,Y} <: PolyScalarField{F,X,Y}
```
A struct for storing the sum of two `PolyScalarField`s. A typical application is the divergence of a `PolyVectorField`, usually formed by the sum of two `BiPoly`s. Note that the terms of the sum can be any type of `PolyScalarField`, including `PolySum`s.
"""
struct PolySum{F,X,Y} <: PolyScalarField{F,X,Y}
    left::PolyScalarField{F,X,Y}
    right::PolyScalarField{F,X,Y}
end

(p::PolySum)(x,y) = p.left(x,y)+p.right(x,y)
(p::PolySum)(x) = p.left(x) + p.right(x)

function Base.:+(p::P,q::Q) where {F,X,Y,P<:PolyField{F,X,Y},Q<:PolyField{F,X,Y}}
    if p==zero(p)
        q
    elseif q == zero(q)
        p
    else
        PolySum(p,q)
    end
end

function Base.zero(::PolySum{F,X,Y}) where {F,X,Y}
    zero(BiPoly{F,X,Y})
end
#LinearAlgebra.dot(p::PolyVectorField{F,X,Y},q::PolyVectorField{F,X,Y})  where {F,X,Y} = p.s1*q.s1 + p.s2*q.s2

Base.:*(n::Number,ps::PolySum) = n*ps.left + n*ps.right
Base.:*(ps::PolySum,n::Number) = n*ps
Base.:*(ps::PolySum,p::BiPoly) = ps.left*p + ps.right*p
Base.:*(p::BiPoly,ps::PolySum) = ps*p
Base.:*(ps::PolySum,p::PolyVectorField) = PolyVectorField(ps*p.s1,ps*p.s2)
Base.:*(p::PolyVectorField,ps::PolySum) = ps*p
Base.:*(ps::PolySum,qs::PolySum) = ps.left*qs.left + ps.right*qs.left + ps.right*qs.left + ps.right*qs.right


###########################################################################################
#########################           AffineToRef        ###########################
###########################################################################################
"""

    AffineToRef{F}

A struct for defining and updating an affine transformation from the reference triangle to some other triangle. 
"""
struct AffineToRef{F} 
    A::MMatrix{2,2,F,4}
    iA::MMatrix{2,2,F,4}
    b::MVector{2,F}
    jacobian::Base.RefValue{F}
    function AffineToRef{F}(A,b) where F<:Number
        @assert size(A)==(2,2)
        @assert size(b)==(2,)
        iA = det(A) != 0 ? inv(A) : zero(A)
        new{F}(F.(A),iA,F.(b),Base.RefValue(abs(det(A))))
    end
end
function AffineToRef(A,b)
    TA = eltype(A)
    Tb = eltype(b)
    T = promote_type(TA,Tb)
    AffineToRef{T}(A,b)
end

AffineToRef{F}() where F = AffineToRef(@MMatrix(zeros(F,2,2)),@MVector(zeros(F,2)))

(aff::AffineToRef)(x) = aff.A*x+aff.b

"""
```
   affine!(aff::AffinetTransformation,vert) 
```
Updates the `AffineToRef` `aff` so that it transforms the reference triangle into the trangle of vertices given by `vert`. `vert` should be a vector of points.
"""
function affine!(aff::AffineToRef,vert)
    (;A,iA,b,jacobian) = aff
    @views A[:, 1] .= 0.5(vert[3] - vert[2])
    @views A[:, 2] .= 0.5(vert[1] - vert[2])
    @views b .= 0.5(vert[1] + vert[3])
    iA .= inv(A)
    jacobian[] = abs(det(A))
end

jac(aff::AffineToRef) = aff.jacobian[]
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


Base.:*(f::Function,p::PolyField) = GeneralField(*,(f,p))
Base.:*(p::PolyField,f::Function) = f*p
LinearAlgebra.dot(f::Function,p::PolyField) = GeneralField(dot,(f,p))
LinearAlgebra.dot(p::PolyField,f::Function) = dot(f,p)

function Base.:*(n::Number,g::GeneralField)
    ind = findfirst(isa.(g.args,PolyField))
    if !isnothing(ind)
        newargs = (g.args[1:ind-1]...,n*g.args[ind],g.args[ind+1:end]...)
        return GeneralField(g.op,newargs)
    else
        return GeneralField(*,(n,g))
    end
end
Base.:*(g::GeneralField,n::Number) = n*g

# A trait for evaluation of Field
abstract type EvalType end
struct Eval <: EvalType end
struct Compose <: EvalType end
struct Pass <: EvalType end

evaltype(_) = Pass()
evaltype(::Function) = Compose()
evaltype(::PolyField) = Eval()

evaluate(::Eval,f,t::AffineToRef,x) = f(x)
evaluate(::Compose,f,t::AffineToRef,x) = f(t(x))
evaluate(::Pass,f,t,x) = f

(o::GeneralField)(x,t::AffineToRef) = o.op((evaluate(evaltype(arg),arg,t,x) for arg in o.args)...)

