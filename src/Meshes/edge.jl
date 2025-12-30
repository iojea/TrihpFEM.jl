"""
    Edge{I} <: SetTuple
a struct for storing an edge of a triangulation. It contains an `NTuple{2,I}` that stores the indices of the two vertices. 
"""
struct Edge{I} <: SetTuple{2,I}
    data::NTuple{2,I}
end
Edge{I}(x::StaticArray) where I = Edge(tuple(I.(x)))
Edge{I}(x::Base.Generator) where I = Edge(I.(tuple(x)))
Edge{I}(x...) where I = Edge(I.(x))
Edge(x,y) = Edge(promote(x,y))
function Edge{I}(x::AbstractArray) where I
    typeof(x)<:AbstractVector || throw(ArgumentError("Edge can only be created from a one dimensional array."))
    Edge{I}(I.(x))
end
    
data(e::Edge) = getproperty(e,:data)


"""
    EdgeAttributes(degree::P,marker::P,refine::Bool) where P<:Integer

constructs a `struct` for storing attributes of an edge. These attributes are:
+ `degree`: degree of the polynomial approximator on the edge.
+ `marker`: a marker indicating if the edge belongs to the boundary of the domain, or to an interface or to the interior. 
+ `refine`: `true` if the edge is marked for refinement. 
"""
struct EdgeAttributes{P<:Integer} 
    degree::Base.RefValue{P}
    marker::Base.RefValue{P}
    refine::Base.RefValue{Bool}
    #adjacent::SVector{2,I}
end
EdgeAttributes(d::P,m::P,r::Bool) where P<:Integer = EdgeAttributes(Ref(d),Ref(m),Ref(r))
EdgeAttributes(e::EdgeAttributes)  = EdgeAttributes(e.degree,e.marker,e.refine)


"""
    ismarked(e::EdgeAttributes)
returs `true` if `e` is marked for refinement. 
"""
@inline ismarked(e::EdgeAttributes)        = e.refine[]
"""
    degree(e::EdgeAttributes)
returs the degree of `e`. 
"""
@inline degree(e::EdgeAttributes)          = e.degree[]
"""
    marker(e::EdgeAttributes)
returs the marker of `e`. The marker indicates if `e` is a boundary edge with Dirichlet or Neumann condition, an interior boundary, etc. 
"""
@inline marker(e::EdgeAttributes)          = e.marker[]
@inline ismarked(e::EdgeAttributes,i::Integer) = marker(e)==i
@inline ismarked(i::Integer) = Base.Fix{2}(ismarked,i)
@inline setmarker!(e::EdgeAttributes{P},i::Integer)where P = e.marker[] = P.(i)
"""
    mark!(e::EdgeAttributes)
marks `e` for refinement.  
"""
@inline mark!(e::EdgeAttributes)           = e.refine[] = true
"""
    setdegree!(e::EdgeAttributes,deg)
sets the degree of `e`.  
"""
@inline setdegree!(e::EdgeAttributes,deg)  = e.degree[] = deg
"""
    isinterior(e::EdgeAttributes)
returns `true` if the edge is an interior one.  
"""
@inline isinterior(e::EdgeAttributes)      = e.marker[] == 0


            
