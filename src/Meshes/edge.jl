"""
    Edge{I} <: SetTuple
a struct for storing an edge of a triangulation. It contains an `NTuple{2,I}` that stores the indices of the two vertices. 
"""
struct Edge{I} <: SetTuple{2, I}
    data::NTuple{2, I}
end
Edge{I}(x::StaticArray) where {I} = Edge(I.(tuple(x...)))
Edge{I}(x::Base.Generator) where {I} = Edge(I.(tuple(x...)))
Edge{I}(x, y) where {I} = Edge((I(x), I(y)))
function Edge{I}(x::T) where {I, T <: AbstractArray}
    T <: AbstractVector || throw(ArgumentError("`Edge` can only be created from a one dimensional array."))
    length(x) == 2 || throw(DimensionMismatch("`Edge`s store two indices."))
    return Edge(I.(tuple(x...)))
end

"""
   data(e::Edge)
returns the tuple defining `e`. 
"""
data(e::Edge) = getproperty(e, :data)


"""
    EdgeAttributes(deg::P,tag::P,refine::Bool,v::MVector{I,P}) where {I<:Integer,P<:Integer}

constructs a `struct` for storing attributes of an edge. These attributes are:
+ `degree`: degree of the polynomial approximator on the edge.
+ `tag`: a tag indicating if the edge belongs to the boundary of the domain, or to an interface or to the interior. 
+ `refine`: `true` if the edge is marked for refinement.
+ `adjacent`: an `MVector{2,I}` storing the two vertices adjacent to the edge. `0` is used in the absence of adjacent vertex (boundary edge).
"""
struct EdgeAttributes{I<:Integer,P <: Integer}
    degree::Base.RefValue{P}
    tag::Base.RefValue{P}
    refine::Base.RefValue{Bool}
    seen::Base.RefValue{I}
    adjacent::MVector{2,I}
    EdgeAttributes{I,P}(d, m, r,s,a) where {I,P} = new{I,P}(Ref(P(d)), Ref(P(m)), Ref(r),Ref(s),a)
end
EdgeAttributes(::Type{I},d::P, m::P, r::Bool) where {I,P} = EdgeAttributes{I,P}(d, m, r, zero(I),MVector{2,I}(zero(I),zero(I)))
# EdgeAttributes(I,d::P, m::P, r::Bool) where {P} = EdgeAttributes(d, m, r, zero(I),MVector{2,I}(zero(I),zero(I)))

"""
    ismarked(e::EdgeAttributes)
returs `true` if `e` is marked for refinement. 
"""
@inline ismarked(e::EdgeAttributes) = e.refine[]
"""
    degree(e::EdgeAttributes)
returs the degree of `e`. 
"""
@inline degree(e::EdgeAttributes) = e.degree[]
"""
    tag(e::EdgeAttributes)
returs the tag of `e`. The tag indicates if `e` is a boundary edge with Dirichlet or Neumann condition, an interior boundary, etc. 
"""
@inline tag(e::EdgeAttributes) = e.tag[]

"""
    istagged(e::EdgeAttributes,i::Integer)
check if `e` has tag `i`.

It can be used passing only `i` to create a function that checks for tag `i`:
    istagged(i) 
"""
@inline istagged(e::EdgeAttributes, i::Integer) = tag(e) == i
@inline istagged(i::Integer) = Base.Fix{2}(istagged, i)

"""
    settag!(e::EdgeAttributes,i::Integer)
sets the tag of `e` to `i`.  
"""
@inline settag!(e::EdgeAttributes{P}, i::Integer) where {P} = e.tag[] = P.(i)

"""
    mark!(e::EdgeAttributes)
marks `e` for refinement.  
"""
@inline mark!(e::EdgeAttributes) = e.refine[] = true

"""
    setdegree!(e::EdgeAttributes,deg)
sets the degree of `e`.  
"""
@inline setdegree!(e::EdgeAttributes, deg) = e.degree[] = deg

"""
    isinterior(e::EdgeAttributes)
returns `true` if the edge is an interior one.  
"""
@inline isinterior(e::EdgeAttributes) = e.tag[] == 0

"""
    isboundary(e::EdgeAttributes)
returns `true` if the edge is a boundary one.  
"""
@inline isboundary(e::EdgeAttributes) = e.tag[] > 0

"""
   isseen(e::EdgeAttributes)
returns `true` if the edge has been seen in the past (in a refinement process)
"""
@inline isseen(e::EdgeAttributes) = e.seen[]>0

"""
   seen(e::EdgeAttributes)
returns the value of the `seen` field. This field is used in the refinement process to indicate the index of the new vertex obtained by bisecting the edge.
"""
@inline seen(e::EdgeAttributes) = e.seen[]
"""
   seen!(e::EdgeAttributes,i)
sets the `seen` field of the edge to `i`. This is done for indicating the index of the new vertex obtained by bisecting the edge.
"""
@inline seen!(e::EdgeAttributes,i) = e.seen[]=i

"""
    adjacents(e::EdgeAttributes)
returns the two vertices adjacent to the edge.
"""
@inline adjacents(e::EdgeAttributes) = e.adjacent

"""
    setadjacents!(e::EdgeAttributes{I,P},v) where {I,P}
set the two vertices adjacent to the edge.
"""
@inline setadjacents!(e::EdgeAttributes,v) = e.adjacent[:] .= v


"""
    setadjacent!(e::EdgeAttributes{I,P},i,v) where {I,P}
set the two vertices adjacent to the edge.
"""
@inline setadjacent!(e::EdgeAttributes,i,a) = e.adjacent[i] = a

