const TriangleList{I,P,F} = Dictionary{Triangle{I},TriangleAttributes{P,F}} where {I,P,F}
const EdgeList{I,P} = Dictionary{Edge{I},EdgeAttributes{P}} where {I,P}

"""
    HPTriangulation
abstract type for triangulations.
"""
abstract type HPTriangulation end


############################################
####         Degrees of Freedom         ####
############################################
"""

    DOF{I<:Integer}
an `struct` for storing the degrees of freedom of a mesh. It can be initialized as an empty structure with
```
    julia> DOF(UInt8);
```
In practice, a `DOF` is created empty an then filled using `degrees_of_freedom!(mesh)`.
"""
struct DOF{I<:Integer}
    n::Base.RefValue{I}
    by_edge::Dictionary{Edge{I}, SArray{S,I, 1} where S<:Tuple}
    by_tri::Dictionary{Triangle{I}, MArray{S,I, 1} where S<:Tuple} 
end

function DOF{I}() where I<:Integer
    by_edge = Dictionary{Edge{I},SArray{S,I,1} where S<:Tuple}()
    by_tri= Dictionary{Triangle{I}, MArray{S,I, 1} where S<:Tuple}()
    DOF{I}(Base.RefValue{I}(zero(I)),by_edge,by_tri)
end

"""
   empty!(d::DOF)

empty both dictionaries stored in `d` (degrees of freedom by edge and by triangle). 
"""
function Base.empty!(d::DOF{I}) where I
    empty!(d.by_edge)
    empty!(d.by_tri)
    d.n[] = zero(I)
end
Base.isempty(d::DOF{I}) where I = d.n[] == zero(I)

############################################
####      Domain Mesh construction      ####
############################################
"""
    $(SIGNATURES)

A mesh for `HP` finite element methods. Its fields are: 
    + `points::Vector{SVector{2,F}}`: a vector of points
    + `trilist::TriangleList{I,P,F}`: set of trilist
    + `edgelist::EdgeList{I,P}`: set of edgelist
    + `dofs::DOF{I}`: auxiliary data for integrating the 

    HPMesh(tri::TriangulationIO)
builds an `HPMesh` from a `Triangulate.TriangulatioIO` struct.

    Meshes of type `HPMesh` can also be constructed using the helper functions such us `circmesh`, `rectmesh`, `squaremesh`. For a more general constructor, check the docs for `hpmesh`.
"""
struct HPMesh{F<:Real,I<:Integer,P<:Integer} <: HPTriangulation
    points::Vector{SVector{2,F}} 
    trilist::TriangleList{I,P,F}
    edgelist::EdgeList{I,P}
    dofs::DOF{I} 
end

function HPMesh(v::Vector{SVector{2,F}},
                t::TriangleList{I,P,F},
                e::EdgeList{I,P}) where {I,P,F}
     HPMesh(v,t,e,DOF{I}())
end

function HPMesh(mat::AbstractMatrix,tris::TriangleList,edgs::EdgeList)
    points = SVector{2,eltype(mat)}.(eachcol(mat))
    HPMesh(points,tris,edgs)
end

function HPMesh{F,I,P}(tri::TriangulateIO) where {F,I,P}
    (;pointlist,trianglelist,edgelist,edgemarkerlist) = tri
    points = SVector{2,F}.(eachcol(pointlist))
    edgelist = maybeconvert(I,edgelist)
    triangles = dictionary([triangle(I,t,pointlist) => TriangleAttributes{P,F}() for t in eachcol(trianglelist)])
    edges = dictionary([Edge(e) => EdgeAttributes(one(P),P(edgemarkerlist[i]),false) for (i,e) in enumerate(eachcol(edgelist))] )
    HPMesh(points,triangles,edges)
end

"""
    $(SIGNATURES)

A helper function for creating `hp`-meshes from boundary data.
For a simple polygonal mesh, it is enough to provide a matrix of `2×N` with the vertices of the polygon and a mesh size `h`.

```
julia> vertices = [-1. 0.;1. 0.;1.5 1.;0. 1.5;-1. 1.]';

julia> m = hpmesh(vertices,1.5)
HPMesh{Float64, Int32, UInt8}
  + 6 nodes.
  + 4 triangles.
  + 9 edges.
6-element Vector{StaticArraysCore.SVector{2, Float64}}:
 [-1.0, 0.0]
 [1.0, 0.0]
 [1.5, 1.0]
 [0.0, 1.5]
 [-1.0, 1.0]
 [0.0, 0.0]

4-element Dictionaries.Dictionary{Triangle{Int32}, TrihpFEM.Meshes.TriangleAttributes{UInt8, Float64}}
 Int32[6, 5, 1] │ :noref
 Int32[4, 2, 3] │ :noref
 Int32[2, 4, 6] │ :noref
 Int32[6, 4, 5] │ :noref


9-element Dictionaries.Dictionary{Edge{Int32}, TrihpFEM.Meshes.EdgeAttributes{UInt8}}
 Int32[1, 6] │ (0x01, :∂𝔇, :noref)
 Int32[6, 5] │ (0x01, :Ω°, :noref)
 Int32[5, 1] │ (0x01, :∂𝔇, :noref)
 Int32[4, 2] │ (0x01, :Ω°, :noref)
 Int32[2, 3] │ (0x01, :∂𝔇, :noref)
 Int32[3, 4] │ (0x01, :∂𝔇, :noref)
 Int32[4, 6] │ (0x01, :Ω°, :noref)
 Int32[6, 2] │ (0x01, :∂𝔇, :noref)
 Int32[4, 5] │ (0x01, :∂𝔇, :noref)
```

For more complex meshes, a matrix of `segments` and a vector of `tags` can be passed. `segments` is a matrix of integers with size `2×S` indicating how vertices should be joined. `tags` is a vector of integers that impose a tag on each segment. The primary goal of tags is to indicate if a piece of boundary will hold Dirichlet (`tag==1`) or Neumann (`tag==2`) conditions. If ommited, Dirichlet conditions will be assumed. In the following example we create a mesh of a square where Neumann conditions are imposed on the upper half.

```
julia> vert = [0. 0.;1. 0.;1. 0.5;1. 1.;0. 1.;0. 0.5]';

julia> segs = [1 2;2 3;3 4;4 5;5 6;6 1]';

julia> tags = [1,1,2,2,2,1];

julia> m = hpmesh(vert,1.5;segments=segs,tags=tags)
HPMesh{Float64, Int32, UInt8}
  + 8 nodes.
  + 6 triangles.
  + 13 edges.
8-element Vector{StaticArraysCore.SVector{2, Float64}}:
 [0.0, 0.0]
 [1.0, 0.0]
 [1.0, 0.5]
 [1.0, 1.0]
 [0.0, 1.0]
 [0.0, 0.5]
 [0.5, 0.0]
 [0.5, 1.0]

6-element Dictionaries.Dictionary{Triangle{Int32}, TrihpFEM.Meshes.TriangleAttributes{UInt8, Float64}}
 Int32[7, 6, 1] │ :noref
 Int32[3, 7, 2] │ :noref
 Int32[6, 3, 8] │ :noref
 Int32[8, 3, 4] │ :noref
 Int32[3, 6, 7] │ :noref
 Int32[6, 8, 5] │ :noref


13-element Dictionaries.Dictionary{Edge{Int32}, TrihpFEM.Meshes.EdgeAttributes{UInt8}}
 Int32[6, 1] │ (0x01, :∂𝔇, :noref)
 Int32[1, 7] │ (0x01, :∂𝔇, :noref)
 Int32[7, 6] │ (0x01, :Ω°, :noref)
 Int32[7, 2] │ (0x01, :∂𝔇, :noref)
 Int32[2, 3] │ (0x01, :∂𝔇, :noref)
 Int32[3, 7] │ (0x01, :Ω°, :noref)
 Int32[6, 3] │ (0x01, :Ω°, :noref)
 Int32[3, 8] │ (0x01, :Ω°, :noref)
 Int32[8, 6] │ (0x01, :Ω°, :noref)
 Int32[3, 4] │ (0x01, :∂𝔑, :noref)
 Int32[4, 8] │ (0x01, :∂𝔑, :noref)
 Int32[8, 5] │ (0x01, :∂𝔑, :noref)
 Int32[5, 6] │ (0x01, :∂𝔑, :noref)
```
Notice that some edges are marked as `∂𝔑`, i.e.: Neumann boundary.

For the markers, `1` indicates _Dirichlet boundary_, whereas `2` stands for _Neumann boundary_. If preferred, a vector of symbols (`:dirichlet` or `:neumann`) can be used:

```
julia> vert = [0. 0.;1. 0.;1. 0.5;1. 1.;0. 1.;0. 0.5]';

julia> segs = [1 2;2 3;3 4;4 5;5 6;6 1]';

julia> tags = [:dirichlet,:dirichlet,:neumann,:neumann,:neumann,:dirichlet];

julia> m = hpmesh(vert,1.5;segments=segs,tags=tags)
HPMesh{Float64, Int32, UInt8}
  + 8 nodes.
  + 6 triangles.
  + 13 edges.
8-element Vector{StaticArraysCore.SVector{2, Float64}}:
 [0.0, 0.0]
 [1.0, 0.0]
 [1.0, 0.5]
 [1.0, 1.0]
 [0.0, 1.0]
 [0.0, 0.5]
 [0.5, 0.0]
 [0.5, 1.0]

6-element Dictionaries.Dictionary{Triangle{Int32}, TrihpFEM.Meshes.TriangleAttributes{UInt8, Float64}}
 Int32[7, 6, 1] │ :noref
 Int32[3, 7, 2] │ :noref
 Int32[6, 3, 8] │ :noref
 Int32[8, 3, 4] │ :noref
 Int32[3, 6, 7] │ :noref
 Int32[6, 8, 5] │ :noref


13-element Dictionaries.Dictionary{Edge{Int32}, TrihpFEM.Meshes.EdgeAttributes{UInt8}}
 Int32[6, 1] │ (0x01, :∂𝔇, :noref)
 Int32[1, 7] │ (0x01, :∂𝔇, :noref)
 Int32[7, 6] │ (0x01, :Ω°, :noref)
 Int32[7, 2] │ (0x01, :∂𝔇, :noref)
 Int32[2, 3] │ (0x01, :∂𝔇, :noref)
 Int32[3, 7] │ (0x01, :Ω°, :noref)
 Int32[6, 3] │ (0x01, :Ω°, :noref)
 Int32[3, 8] │ (0x01, :Ω°, :noref)
 Int32[8, 6] │ (0x01, :Ω°, :noref)
 Int32[3, 4] │ (0x01, :∂𝔑, :noref)
 Int32[4, 8] │ (0x01, :∂𝔑, :noref)
 Int32[8, 5] │ (0x01, :∂𝔑, :noref)
 Int32[5, 6] │ (0x01, :∂𝔑, :noref)
```


Furthermore, meshes with holes can also be constructed. In this case, the `segments` should include the segments corresponding to the interior boudaries (always in a positively oriented order), and a parameter `holes` should also be passed. `holes` is a matrix where each column contains a point included in one hole. Only one point per hole is necessary. 

```
julia> vert = [0. 0.;1. 0.;1. 0.5;1. 1.;0. 1.;0. 0.5;0.25 0.25;0.5 0.25;0.5 0.5;0.25 0.5]';

julia> segs = [1 2;2 3;3 4;4 5;5 6;6 1;7 8;8 9;9 10;10 7]';

julia> tags = [1,1,2,2,2,1,1,1,2,2];

julia> hole = [0.3 0.3]';

julia> m = hpmesh(vert,0.1;segments=segs,tags=tags,holes=hole)
HPMesh{Float64, Int32, UInt8}
  + 18 nodes.
  + 21 triangles.
  + 39 edges.
18-element Vector{StaticArraysCore.SVector{2, Float64}}:
 [0.0, 0.0]
 [1.0, 0.0]
 [1.0, 0.5]
 [1.0, 1.0]
 [0.0, 1.0]
 [0.0, 0.5]
 [0.25, 0.25]
 [0.5, 0.25]
 [0.5, 0.5]
 [0.25, 0.5]
 [0.5, 0.0]
 [0.0, 0.75]
 [0.5, 1.0]
 [0.75, 0.0]
 [0.75, 0.375]
 [0.25, 1.0]
 [0.375, 0.75]
 [0.75, 0.6875]

21-element Dictionaries.Dictionary{Triangle{Int32}, TrihpFEM.Meshes.TriangleAttributes{UInt8, Float64}}
    Int32[6, 1, 7] │ :noref
   Int32[1, 11, 7] │ :noref
  Int32[10, 12, 6] │ :noref
   Int32[6, 7, 10] │ :noref
  Int32[9, 17, 10] │ :noref
 Int32[18, 13, 17] │ :noref
 Int32[12, 17, 16] │ :noref
  Int32[14, 8, 11] │ :noref
  Int32[2, 15, 14] │ :noref
   Int32[2, 3, 15] │ :noref
  Int32[15, 18, 9] │ :noref
  Int32[18, 15, 3] │ :noref
   Int32[8, 15, 9] │ :noref
  Int32[12, 16, 5] │ :noref
 Int32[17, 12, 10] │ :noref
   Int32[7, 11, 8] │ :noref
 Int32[16, 17, 13] │ :noref
  Int32[14, 15, 8] │ :noref
  Int32[18, 17, 9] │ :noref
   Int32[3, 4, 18] │ :noref
  Int32[4, 13, 18] │ :noref


39-element Dictionaries.Dictionary{Edge{Int32}, TrihpFEM.Meshes.EdgeAttributes{UInt8}}
   Int32[1, 7] │ (0x01, :Ω°, :noref)
   Int32[7, 6] │ (0x01, :Ω°, :noref)
   Int32[6, 1] │ (0x01, :∂𝔇, :noref)
  Int32[1, 11] │ (0x01, :∂𝔇, :noref)
  Int32[11, 7] │ (0x01, :Ω°, :noref)
  Int32[6, 10] │ (0x01, :Ω°, :noref)
 Int32[10, 12] │ (0x01, :Ω°, :noref)
  Int32[12, 6] │ (0x01, :∂𝔑, :noref)
  Int32[7, 10] │ (0x01, :∂𝔑, :noref)
  Int32[10, 9] │ (0x01, :∂𝔑, :noref)
  Int32[9, 17] │ (0x01, :Ω°, :noref)
 Int32[17, 10] │ (0x01, :Ω°, :noref)
 Int32[13, 17] │ (0x01, :Ω°, :noref)
 Int32[17, 18] │ (0x01, :Ω°, :noref)
 Int32[18, 13] │ (0x01, :Ω°, :noref)
             ⋮ │ ⋮
  Int32[3, 15] │ (0x01, :Ω°, :noref)
 Int32[15, 18] │ (0x01, :Ω°, :noref)
  Int32[18, 9] │ (0x01, :Ω°, :noref)
  Int32[9, 15] │ (0x01, :Ω°, :noref)
  Int32[3, 18] │ (0x01, :Ω°, :noref)
   Int32[9, 8] │ (0x01, :∂𝔇, :noref)
  Int32[8, 15] │ (0x01, :Ω°, :noref)
  Int32[5, 12] │ (0x01, :∂𝔑, :noref)
  Int32[16, 5] │ (0x01, :∂𝔑, :noref)
   Int32[8, 7] │ (0x01, :∂𝔇, :noref)
 Int32[13, 16] │ (0x01, :∂𝔑, :noref)
   Int32[3, 4] │ (0x01, :∂𝔑, :noref)
  Int32[4, 18] │ (0x01, :Ω°, :noref)
  Int32[4, 13] │ (0x01, :∂𝔑, :noref)
```
"""
function hpmesh(vertices,h;
                segments=nothing,
                tags=nothing,
                holes=nothing)
    F = eltype(vertices)
    P = UInt8
    tri = TriangulateIO()
    tri.pointlist = F.(vertices)
    if isnothing(segments)
        segments = _boundary_segments(size(vertices,2))
    end
    if isnothing(tags)
        tags = ones(P,size(vertices,2))
    elseif eltype(tags)==Symbol
        tags = [BOUNDARY_DICT[m] for m in tags]
    end
    tri.segmentlist = segments
    tri.segmentmarkerlist = tags
    if !isnothing(holes)
        tri.holelist = holes
    end
    maxarea  = Printf.@sprintf "%0.15f" h^2/2
    minangle = Printf.@sprintf "%0.15f" 30.
    (tri,_) = triangulate("pea$(maxarea)q$(minangle)Q",tri)
    I = eltype(tri.edgelist)
    HPMesh{F,I,P}(tri)
end

"""
    setboundary!(s,condition,m::HPMesh)
sets the boundary edges of mesh `m` which satisfy `condition` to kind `s`, where `s` can be `:dirichlet`, `:neumann` (or alternatively, `1` for Dirichlet or `2` for Neumann). Notice that `condition` is evaluated at the vertices of each edge. An edge is tagged when both vertices satisfy `condition`.

See also [`setdirichlet!`](@ref) and [`setneumann!`](@ref).
"""
function setboundary!(i::Integer,
                       condition::Function,
                       m::HPMesh{F,I,P}) where {F,I,P}
    (;points,edgelist) = m
    for ea in pairs(edgelist)
        e = first(ea)
        if condition(points[e[1]]) && condition(points[e[2]])
            settag!(last(ea),P(i))
        end
    end
end
function setboundary!(s::Symbol,
                       condition::Function,
                       m::HPTriangulation)
    i = BOUNDARY_DICT[s]
    setboundary!(i,condition,m)    
end
"""
    setdirichlet!(condition,m::HPMesh)
Marks as `:dirichlet` all boundary edges of `m` satisfying `condition`. Notice that `condition` is evaluated at the vertices of each edge. An edge is tagged when both vertices satisfy `condition`.
See also [`setneumann!`](@ref).
"""
setdirichlet!(condition::Function,m::HPMesh) = setboundary!(1,condition,m)
"""
    setneumann!(condition,m::HPMesh)
Marks as `:neumann` all boundary edges of `m` satisfying `condition`. Notice that `condition` is evaluated at the vertices of each edge. An edge is tagged when both vertices satisfy `condition`.
See also [`setdirichlet!`](@ref).
"""
setneumann!(condition::Function,m::HPMesh) = setboundary!(2,condition,m)


############################################
####        Auxiliary functions         ####
############################################
"""
    maybeconvert(::Type{T},arr) where T
If necessary converts the elements of `arr` to type `T`.  
"""
function maybeconvert(::Type{T},arr::AbstractArray) where T
    eltype(arr)==T ? arr : convert.(T,arr)
end

"""
    degtype(::HPMesh{F,I,P}) wehere {F,I,P}
returns `P`  
"""
degtype(::HPMesh{F,I,P}) where {F,I,P} = P
"""
    inttype(::HPMesh{F,I,P}) wehere {F,I,P}
returns `I`  
"""
inttype(::HPMesh{F,I,P}) where {F,I,P} = I
"""
    floattype(::HPMesh{F,I,P}) wehere {F,I,P}
returns `F`  
"""
floattype(::HPMesh{F,I,P}) where {F,I,P} = F


function Base.copy(mesh::HPMesh)
    HPMesh(deepcopy(mesh.points),deepcopy(mesh.trilist),deepcopy(mesh.edgelist))
end

"""
  edges(m)
retuns the edges of `m`, being `m` an `HPMesh` or an `EdgeList`.  
"""
@inline edges(list::T) where T<:EdgeList = keys(list)
@inline edges(mesh::T) where T<:HPMesh = keys(mesh.edgelist)
"""
  triangles(m) 
retuns the triangles of `m`, being `m` and `HPMesh` or a `TriList`.
"""
@inline triangles(list::T) where T<:TriangleList = keys(list)
@inline triangles(mesh::T) where T<:HPMesh = keys(mesh.trilist)

"""
  dof(m)
retuns the degrees of freedom of the mesh `m`.
"""
@inline dof(mesh) = mesh.dofs

"""
    isgreen(t::Triangle,m::HPMesh
retuns `true` if triangle `t` in `m` is marked green.
"""
@inline isgreen(t::Triangle,m::HPMesh) = m.trilist[t].refine[] == 1
"""
    isblue(t::Triangle,m::HPMesh
retuns `true` if triangle `t` in `m` is marked blue.
"""
@inline isblue(t::Triangle,m::HPMesh)  = m.trilist[t].refine[] == 2
"""
    isred(t::Triangle,m::HPMesh
retuns `true` if triangle `t` in `m` is marked red.
"""
@inline isred(t::Triangle,m::HPMesh)  = m.trilist[t].refine[] == 3

"""
    $(SIGNATURES)

Checks if 'e' is stored as presented or in reverse order. 
"""
function _same_order(e::Edge{I},elist::EdgeList{I,P}) where {I,P}
    _,(_,k) = gettoken(elist,e)
    oe      = gettokenvalue(keys(elist),k)
    oe == e
end

"""
    $(SIGNATURES)

returns a permutation of contiguous indices `ind` such that `v[ind]` is ordered. 
"""
function psortperm(v)
    i    = argmin(v)
    if v[i] ≤ v[mod1(i+1,3)] ≤ v[mod1(i+2,3)]
        ind = mod1.(i:i+2,3)
    else
        ind = mod1.(i:-1:i-2,3)
    end
    return ind
end

"""
  $(SIGNATURES)  
Given a triangle  `t`, of  a mesh `mesh` returs the degrees of the edges of `t`
"""
function degrees(t::Triangle{I},mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    eds = edges(t)
    p = Tuple(degree.(getindices(edgelist,eds)))
    sort(p)
end

"""
  $(SIGNATURES)

Given a triangle `t` belonging to a `mesh`, it returns the (sorted) degrees
`p₁<=p₂<=p₃` of its edges, and the edges also sorted according to the degrees.
"""
function psortededges(t::Triangle{I},mesh::HPMesh{F,I,P}) where {F,I,P}
    p,nod = psortednodes(t,mesh)
    eds   = [Edge(nod[SVector(i,mod1(i+1,3))]) for i in 1:3]
    p,eds
end

"""
   $(SIGNATURES)

Given a triangle `t` belonging to a `mesh`, it returns the degrees `p₁<=p₂<=p₃` of the edges (sorted) and the nodes of `t` also sorted according to the degrees. 
"""
function psortednodes(t::Triangle{I},mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    eds = edges(t)
    p   = degree.(getindices(edgelist,eds))
    pind = psortperm(p)
    p[pind],first.(eds[pind])
end

"""
    $(SIGNATURES)

returns a list containing the indices of the vertices of `mesh` tha lie in its boundary. 
"""
function boundarynodes(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    v = zeros(I,sum(>(0),tag.(edgelist)))
    i = 1
    for e in pairs(filter(isboundary,edgelist))
        if (last(e))>0
            for j in first(e)
                if j ∉ v
                    v[i] = j
                    i   += 1
                end
            end
        end
    end
    v
end


            
"""
    _boundary_segments(n)
creates pairs of indices for building a sequence of segments from a list of points.  
"""
_boundary_segments(n) = reduce(hcat,[i,mod1(i+1,n)] for i in 1:n)


"""

    compute_dimension(p₁,p₂,p₃)
    compute_dimension(p₁,p₂)   
    compute_dimension(p₁)
    compute_dimension(t::Tuple)

Computes the dimension of the space ℓp₁p₂p₃. 
"""
compute_dimension(p₁::Integer, p₂::Integer, p₃::Integer) = sum(min(p₂, p₃ - j) + 1 for j = 0:p₁);
compute_dimension(p₁::Integer, p₂::Integer) = compute_dimension(p₁, p₂, p₂)
compute_dimension(p₁::Integer) = compute_dimension(p₁, p₁)
compute_dimension(t::T) where {T<:AbstractArray} = compute_dimension(t...)
compute_dimension(t::Tuple) = compute_dimension(t...)

"""
    degrees_of_freedom_by_edge!(mesh::HPMesh{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}

Creates a dictionary (from `Dictionaries.jl`) where the keys are the edges of `mesh` and the values are vectors with indices corresponding to the nodal degrees of freedom. 
"""
function degrees_of_freedom_by_edge!(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;points,edgelist,dofs) = mesh 
    (;by_edge)  = dofs
    i        = size(points,2)+1
    for edge in edges(edgelist)
        med  = collect(I(i):I(i+degree(edgelist[edge])-2))
        v = SVector{length(med)+2,I}(edge[1],med...,edge[2])
        set!(by_edge,edge,v)
        i   += degree(edgelist[edge])-1
    end
end


"""
    degrees_of_freedom!(mesh::HPMesh{F,I,P}) where {F,I,P}

Creates a dictionary (from `Dictionaries.jl`) where the keys are the triangles of the mesh and the values are vectores storing the indices of the corresponding degrees of freedom. 

Internally, `degrees_of_freedom_by_edge` is called in order to obtain the nodal degrees of freedom, and then the degrees of freedom corresponding to bubble functions are computed. 

If a dictionary of degrees of freedom by edge has already been computed, it is recommended to run: 

    degrees_of_freedom!(mesh::HPMesh{F,I,P},by_edge::Dictionary{Edge{I},Vector{I}}) where {F,I,P}
"""
function degrees_of_freedom!(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist,trilist,dofs) = mesh
    (;n,by_edge,by_tri) = dofs
    if isempty(by_edge)
        degrees_of_freedom_by_edge!(mesh)
    end
    k       = maximum(maximum.(by_edge))+1 #first non-edge dof
    for t in triangles(trilist)
        p,t_edges = psortededges(t,mesh)
        newdofs = @MVector zeros(I,compute_dimension(p))
        set!(by_tri,t,newdofs)
        j = 1 #counter of dof in current triangle
        @inbounds for i in 1:3
            newdof = by_edge[t_edges[i]]
            if  _same_order(t_edges[i],edgelist)
                newdofs[j:j+length(newdof)-2] .= newdof[1:end-1]
            else
                newdofs[j:j+length(newdof)-2] .= reverse(newdof[2:end])
            end 
            j += length(newdof)-1
        end
        newdofs[j:end] = k:k+(length(newdofs)-j)
        k += length(newdofs)-j + 1
    end
    n[] = k-1
end


"""
    tagged_dof(mesh::HPMesh{F,I,P},tag) where {F,I,P}

Returs a list of indices corresponding to the degrees of freedom labeled with `tag`. If a vector of labels is passed, it returs degrees of freedom marked with any of them. 
"""
function tagged_dof(mesh::HPMesh{F,I,P},tag::N) where {F,I,P,N<:Integer}
    tagged_dof(mesh,[tag])
end


"""
    tagged_dof(mesh::HPMesh,tags)
returns a list of all degrees of freedom tagged with a tag in `tag`.
"""
function tagged_dof(mesh::HPMesh{F,I,P},tags) where {F,I,P}
    (;dofs) = mesh
    msg = "The degrees of freedom of the mesh have not been computed."
    isempty(dofs) && throw(ArgumentError(msg))
    _tagged_dof(mesh,tags)
end

"""
    _tagged_dof(mesh::HPMesh,tags)
internal function. Returns a list of all degrees of freedom tagged with a tag in `tag`.  
"""
function _tagged_dof(mesh::HPMesh{F,I,P},tags) where {F,I,P}
    (;edgelist,dofs) = mesh
    (;by_edge) = dofs
    v = Vector{I}()
    for e in pairs(edgelist)
        if tag(last(e)) in tags
            push!(v,by_edge[first(e)]...)
        end
    end
    return unique(v)
end

"""
    boundary_dof(mesh::HPMesh{F,I,P}) where {F,I,P}

Returns a list of the degrees of freedom lying at the boundary of the mesh.
Internally, it calls `degrees_of_freedom_by_edge` in order to obtain the nodal degrees of freedom. If a dof by edge dictionary has already been computed, it is recommended to run: 

    boundary_dof(mesh,by_edge) 
"""
function boudary_dof(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;dofs) = mesh; (;by_edge) = dofs;
    msg = "The degrees of freedom of the mesh have not been computed."
    isempty(dofs) && throw(ArgumentError(msg))
    tagged_dof(mesh,(1,2))
end


########
