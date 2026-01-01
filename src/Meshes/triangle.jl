"""

  Triangle{I}(t::NTuple{3,I}) where I<:Integer

creates a triangle for a triangulation. `I`  is the type of the indices. 
"""
struct Triangle{I} <: SetTuple{3,I}
    data::NTuple{3,I}
end

Triangle{I}(x::StaticArray) where I = Triangle(tuple(I.(x...)))
Triangle{I}(x::Base.Generator) where I = Triangle(I.(tuple(x...)))
Triangle{I}(x,y,z) where I= Triangle(I(x),I(y),I(z))
function Triangle{I}(x::T) where {I,T<:AbstractArray}
    T<:AbstractVector || throw(ArgumentError("`Triangle` can only be created from a one dimensional array."))
    length(x) == 3 || throw(DimensionMismatch("`Triangle`s store three indices."))
    Triangle(I.(tuple(x...)))
end

"""
    data(t::Triangle)
returns the tuple defining `t`.  
"""
data(t::Triangle) = t.data

"""
    _eval(t::Triangle,k)

Returns the index stored in the triangle at index `k` mod 3. 
""" 
@inline _eval(t::Triangle,k) = t[mod1(k,3)]


"""
    edges(t::Triangle)

Return a tuple of edges with type `Edge`, containing the edges of `t`.
"""
@inline edges(t::Triangle)   = Tuple(Edge(_eval(t,i),_eval(t,i+1)) for i in 1:3)


""" 
    longestedge(t::Triangle)

Returns an `Edge` with the longest edge of `t`.
"""
 @inline longestedge(t::Triangle) = Edge(_eval(t,1),_eval(t,2))


"""

  triangle(t,p::Vector)
  triangle(t,p::AbstractMatrix)

constructs an `Triangle` from a list of indices `t` and a list of points `p`. `p` can be a vector of vectors or a matrix with point coordinates on each column. The output has the verices sorted so that the first edge is the longest.
This is the preferred constructor for a `Triangle`.
"""
function triangle(::Type{I},t,p::Vector) where {I}
    maxi = argmax(norm.(p[t[SVector(1,2,3)]].- p[t[SVector(2,3,1)]]))
    Triangle(t[one(I)*mod1.(maxi:maxi+2,3)])
end
function triangle(::Type{I},t,p::AbstractMatrix) where {I}
    maxi = argmax(sum(abs2,p[:,t[SVector(1,2,3)]] - p[:,t[SVector(2,3,1)]],dims=1))[2]
    Triangle(t[one(I)*mod1.(maxi:maxi+2,3)])
end
triangle(t,p) = triangle(eltype(t),t,p)

# Attributes

"""
    TriangleAttributes{P,F}(refine,η,ηₚ) where {P<:Integer,F<:AbstractFloat}

constructs a `struct` for storing attributes of a triangle. These attributes are:
+ `refine`: 
    - `0`: not marked for refinement. 
    - `1`: marked for refinement of _green_ type.
    - `2`: marked for refinement of _blue_ type.
    - `3`: marked for refinement of _red_ type. 
+ `η`: estimate for the local error. 
+ `ηₚ`: predictor of local error, based on previos estimations.  
The types can be inferred from the data:

    TriangleAttributes(refine,η,ηₚ)

If only the `refine` argument is passed, `η` and `ηₚ` are initialized as `0.`
If no arguments are passed, `refine` is initialized as `Int8(0)`.
"""
struct TriangleAttributes{P<:Integer,F<:AbstractFloat} 
    refine::Base.RefValue{P}
    η::Base.RefValue{F}
    ηₚ::Base.RefValue{F}
    TriangleAttributes{P,F}(val,η,ηₚ) where{P,F} = new{P,F}(Ref(P(val)),Ref(F(η)),Ref(F(ηₚ)))
    TriangleAttributes{P,F}() where {P,F} = new{P,F}(Ref(zero(P)),Ref(zero(F)),Ref(zero(F)))
end
function TriangleAttributes(r,η,ηₚ)
    z = promote(η,ηₚ)
    TriangleAttributes{typeof(r),eltype(z)}(r,z...)
end
TriangleAttributes() = TriangleAttributes(zero(UInt8),0.0,0.0)

"""
    ismarked(t::TriangleAttributes)
checks if `t` is marked for refinement.  
"""
@inline ismarked(t::TriangleAttributes)  = t.refine[] > 0

"""
    ismarked(t::TriangleAttributes)
checks if `t` is marked for refinement as `:green`.  
"""
@inline isgreen(t::TriangleAttributes)   = t.refine[] == 1

"""
    ismarked(t::TriangleAttributes)
checks if `t` is marked for refinement as `:blue`  
"""
@inline isblue(t::TriangleAttributes)    = t.refine[] ==2

"""
    ismarked(t::TriangleAttributes)
checks if `t` is marked for refinement as `:ref`  
"""
@inline isred(t::TriangleAttributes)     = t.refine[] == 3

"""
    mark!(t::TriangleAttributes,k::Integer=3)
marks `t` for refinemente according to `k`. `k=1` is `:green`, `k=2` is `:blue`, `k=3` is `:red`. 
"""
@inline mark!(t::TriangleAttributes,k::Integer=3)   = t.refine[] = k

"""
    setη!(t::TriangleAttributes,η)
sets the local estimator of the error.  
"""
@inline setη!(t::TriangleAttributes,η)   = t.η[] = η


"""
    setη!(t::TriangleAttributes,η)
sets the _previous_ local estimator of the error.  
"""
@inline setηₚ!(t::TriangleAttributes,ηₚ) = t.ηₚ[] = ηₚ

