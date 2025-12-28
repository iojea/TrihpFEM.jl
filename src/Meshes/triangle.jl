"""

  Triangle{I}(t::NTuple{3,I}) where I<:Integer

creates a triangle for a triangulation. `I`  is the type of the indices. 
"""
struct Triangle{I} <: SetTuple{3,I}
    data::NTuple{3,I}
end


"""
    _eval(t::Triangle,k)

Returns the index stored in the triangle at index `k` mod 3. 
""" 
@inline _eval(t::Triangle,k) = t[mod1(k,3)]


"""
    edges(t::Triangle)

Return a tuple of edges with type `Edge`, containing the edges of `t`.
"""
@inline edges(t::Triangle)   = tuple(Edge(_eval(t,i),_eval(t,i+1)) for i in 1:3)


""" 
    longestedge(t::Triangle)

Returns an `Edge` with the longest edge of `t`.
"""
 @inline longestedge(t::Triangle) = Edge(_eval(t,1),_eval(t,2))


"""
  triangle(t::T,p::AbstractMatrix)

constructs an `Triangle` where the vertices are given by the columns of `p[:,t]`. Hence, the triangle is defined by the indices stored in `t`, but sorted in such a way that the first edge is the longest. 
"""
function triangle(::Type{I},t,p::AbstractMatrix) where {I}
    maxi = argmax(sum(abs2,p[:,t[SVector(1,2,3)]] - p[:,t[SVector(2,3,1)]],dims=1))[2]
    Triangle(t[one(I)*mod1.(maxi:maxi+2,3)])
end



# Properties

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
    TriangleAttributes(val::P,η::F,ηₚ::F) where {P,F} = new(Ref(val),Ref(η),Ref(ηₚ))
    TriangleAttributes{P,F}(val,η,ηₚ) where{P,F} = new(P(val),F(η),F(ηₚ))
end
TriangleAttributes{P,F}() where {P,F} = TriangleAttributes(zero(P),zero(F),zero(F))

@inline ismarked(t::TriangleAttributes)  = t.refine[] > 0
@inline isgreen(t::TriangleAttributes)   = t.refine[] == 1
@inline isblue(t::TriangleAttributes)    = t.refine[] ==2
@inline isred(t::TriangleAttributes)     = t.refine[] == 3
@inline mark!(t::TriangleAttributes,k)   = t.refine[] = k
@inline mark!(t::TriangleAttributes)     = mark!(t,3)
@inline setη!(t::TriangleAttributes,η)   = t.η[] = η
@inline setηₚ!(t::TriangleAttributes,ηₚ) = t.ηₚ[] = ηₚ

