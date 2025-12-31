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

"""
   data(e::Edge)
returns the tuple defining `e`. 
"""
data(e::Edge) = getproperty(e,:data)


"""
    EdgeAttributes(degree::P,tag::P,refine::Bool) where P<:Integer

constructs a `struct` for storing attributes of an edge. These attributes are:
+ `degree`: degree of the polynomial approximator on the edge.
+ `tag`: a tag indicating if the edge belongs to the boundary of the domain, or to an interface or to the interior. 
+ `refine`: `true` if the edge is marked for refinement. 
"""
struct EdgeAttributes{P<:Integer} 
    degree::Base.RefValue{P}
    tag::Base.RefValue{P}
    refine::Base.RefValue{Bool}
    #adjacent::SVector{2,I}
    EdgeAttributes{P}(d,m,r) where P = new{P}(Ref(P(d)),Ref(P(m)),Ref(r))
end
EdgeAttributes(d::P,m::P,r::Bool) where P<:Integer = EdgeAttributes{P}(d,m,r)
# EdgeAttributes(e::EdgeAttributes)  = EdgeAttributes(e.degree,e.tag,e.refine)


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
@inline istagged(e::EdgeAttributes,i::Integer) = tag(e)==i
@inline istagged(i::Integer) = Base.Fix{2}(istagged,i)

"""
    settag!(e::EdgeAttributes,i::Integer)
sets the tag of `e` to `i`.  
"""
@inline settag!(e::EdgeAttributes{P},i::Integer) where P = e.tag[] = P.(i)

"""
    mark!(e::EdgeAttributes)
marks `e` for refinement.  
"""
@inline mark!(e::EdgeAttributes) = e.refine[] = true

"""
    setdegree!(e::EdgeAttributes,deg)
sets the degree of `e`.  
"""
@inline setdegree!(e::EdgeAttributes,deg) = e.degree[] = deg

"""
    isinterior(e::EdgeAttributes)
returns `true` if the edge is an interior one.  
"""
@inline isinterior(e::EdgeAttributes) = e.tag[] == 0

"""
    isboundary(e::EdgeAttributes)
returns `true` if the edge is a boundary one.  
"""
@inline isboundary(e::EdgeAttributes) = e.tag[]>0


            
