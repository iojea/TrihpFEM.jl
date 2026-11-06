const TriangleList{F, I, P} = Dictionary{Triangle{I}, TriangleAttributes{F, P}} where {F, I, P}
const EdgeList{I, P} = Dictionary{Edge{I}, EdgeAttributes{I, P}} where {I, P}

abstract type Point{D, T} <: FieldVector{D, T} end

"""
   Point2D(x,y)

creates a 2D point with coordinates `x` and `y`. They can be accessed with number indexing or by field. 
"""
struct Point2D{T <: Real} <: Point{2, T}
    x::T
    y::T
end

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
a struct for storing the degrees of freedom of a mesh. It can be initialized as an empty structure with
```
    julia> DOF{Int64}();
```
In practice, a `DOF` is created empty an then filled using `degrees_of_freedom!(mesh)`.
"""
struct DOF{I <: Integer}
    n::Base.RefValue{I}
    by_edge::Dictionary{Edge{I}, SArray{S, I, 1} where {S <: Tuple}}
    by_tri::Dictionary{Triangle{I}, MArray{S, I, 1} where {S <: Tuple}}
    gen::Base.RefValue{I}
end

function DOF{I}() where {I <: Integer}
    by_edge = Dictionary{Edge{I}, SArray{S, I, 1} where {S <: Tuple}}()
    by_tri = Dictionary{Triangle{I}, MArray{S, I, 1} where {S <: Tuple}}()
    return DOF{I}(Base.RefValue{I}(zero(I)), by_edge, by_tri, Ref{I}(zero(I)))
end

"""
   empty!(d::DOF)

empty both dictionaries stored in `d` (degrees of freedom by edge and by triangle). 
"""
function Base.empty!(d::DOF{I}) where {I}
    empty!(d.by_edge)
    empty!(d.by_tri)
    d.n[] = zero(I)
    d.gen[] = zero(I)
    return nothing
end
Base.isempty(d::DOF{I}) where {I} = d.n[] == zero(I)
Base.length(d::DOF) = d.n[]
Base.copy(dof::DOF{I}) where {I} = DOF{I}(
    Ref{I}(length(dof)),
    deepcopy(dof.by_edge),
    deepcopy(dof.by_tri),
    Ref{I}(gen(dof))
)
gen(dof::DOF{I}) where {I} = dof.gen[]
setgen!(dof::DOF{I}, g::I) where {I} = dof.gen[] = g

############################################
####      Domain Mesh construction      ####
############################################
"""
    $(SIGNATURES)

A mesh for `HP` finite element methods. Its fields are: 
    + `points::Vector{SVector{2,F}}`: a vector of points
    + `trilist::TriangleList{I,P,F}`: set of trilist
    + `edgelist::EdgeList{I,P}`: set of edgelist
    + `dofs::DOF{I}`: degrees of freedom.
    + `cmatrix::CMatrix{P,F}`: memoization for C matrix (change from std to mixed basis).
    + `gen::RefValue{I}`: generation, keeps track of modifications to avoid inconsistencies.
    + `lock::RefValue{Bool}`: lock to avoid direct usage of the mesh, once it is passed to a Problem.

HPMesh(tri::TriangulationIO)
builds an `HPMesh` from a `Triangulate.TriangulatioIO` struct.

    Meshes of type `HPMesh` can also be constructed using the helper functions such us `circmesh`, `rectmesh`, `squaremesh`. For a more general constructor, check the docs for `hpmesh`.
"""
struct HPMesh{F <: Real, I <: Integer, P <: Integer} <: HPTriangulation
    points::Vector{Point2D{F}}
    trilist::TriangleList{F, I, P}
    edgelist::EdgeList{I, P}
    dofs::DOF{I}
    cmatrix::CMatrix{F, P}
    gen::Base.RefValue{I}
    lock::Base.RefValue{Bool}
    HPMesh{F, I, P}(v, t, e, d, c, g, l) where {F, I, P} = new{F, I, P}(v, t, e, d, c, Ref{I}(g), Ref{Bool}(l))
end

function HPMesh(
        v::Vector{Point2D{F}},
        t::TriangleList{F, I, P},
        e::EdgeList{I, P}
    ) where {F, I, P}
    return HPMesh{F, I, P}(v, t, e, DOF{I}(), CMatrix{F, P}(), one(I), false)
end

function HPMesh(mat::AbstractMatrix{F}, tris::TriangleList, edgs::EdgeList) where {F}
    points = Point2D{F}.(eachcol(mat))
    println(typeof(points))
    return HPMesh(points, tris, edgs)
end


function HPMesh{F, I, P}(tri::TriangulateIO) where {F, I, P}
    (; pointlist, trianglelist, edgelist, edgemarkerlist) = tri
    points = Point2D{F}.(eachcol(pointlist))
    edgelist = maybeconvert(I, edgelist)
    trianglelist = maybeconvert(I, trianglelist)
    triangles = dictionary([triangle(I, t, pointlist) => TriangleAttributes{F, P}() for t in eachcol(trianglelist)])
    edges = dictionary([Edge{I}(e) => EdgeAttributes(I, one(P), P(edgemarkerlist[i]), false) for (i, e) in enumerate(eachcol(edgelist))])
    computeneighbors!(edges, triangles)
    return HPMesh(points, triangles, edges)
end


"""
   computeneighbors!(tri::TriangulateIO())
computes neighboring  
"""
function computeneighbors!(edgelist::EdgeList, triangles::TriangleList)
    for t in keys(triangles)
        for (i, e) in enumerate(edges(t))
            ep = edgelist[e]
            j = 1 + (ep.adjacent[1] > 0)
            setadjacent!(ep, j, _eval(t, i + 2))
        end
    end
    return
end

"""jldoctest

    $(SIGNATURES)

A helper function for creating `hp`-meshes from boundary data.
For a simple polygonal mesh, it is enough to provide a matrix of `2×N` with the vertices of the polygon and a mesh size `h`.

```
julia> vertices = [-1. 0.;1. 0.;1.5 1.;0. 1.5;-1. 1.]';

julia> m = hpmesh(vertices,1.5)
julia> m = hpmesh(vertices,1.5)
HPMesh{Float64, Int64, Int8}
  + 6 nodes.
  + 4 triangles.
  + 9 edges.
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
```


Furthermore, meshes with holes can also be constructed. In this case, the `segments` should include the segments corresponding to the interior boudaries (always in a positively oriented order), and a parameter `holes` should also be passed. `holes` is a matrix where each column contains a point included in one hole. Only one point per hole is necessary. 

```
julia> vert = [0. 0.;1. 0.;1. 0.5;1. 1.;0. 1.;0. 0.5;0.25 0.25;0.5 0.25;0.5 0.5;0.25 0.5]';

julia> segs = [1 2;2 3;3 4;4 5;5 6;6 1;7 8;8 9;9 10;10 7]';

julia> tags = [1,1,2,2,2,1,1,1,2,2];

julia> hole = [0.3 0.3]';

julia> m = hpmesh(vert,0.1;segments=segs,tags=tags,holes=hole)
HPMesh{Float64, Int64, Int8}
  + 171 nodes.
  + 294 triangles.
  + 465 edges.
```

Neumann boundary nodes can be set after creation: see [`setneumann!`](@ref).

There are functions for creating simple meshes: see [`circmesh`](@ref), [`rectmesh`](@ref)

There are functions for creating simple meshes: see [`circmesh`](@ref), [`rectmesh`](@ref).
"""
function hpmesh(
        vertices, h;
        segments = nothing,
        tags = nothing,
        holes = nothing,
        I = isnothing(segments) ? Int64 : eltype(segments),
        P = Int8,
        F = eltype(vertices)
    )
    typemax(I) < typemax(Int64) && @warn "You are using $I for indices. Take into account that this imposes hard limits on the mesh: a combination of a large number of triangles and  relatively large polynomial degrees, produced during `hp`-refinement, can give indices exceding `typemax($I)==$(typemax(I))`."
    tri = TriangulateIO()
    tri.pointlist = F.(vertices)
    if isnothing(segments)
        segments = _boundary_segments(size(vertices, 2))
    end
    if isnothing(tags)
        tags = ones(P, size(vertices, 2))
    elseif eltype(tags) == Symbol
        tags = [BOUNDARY_DICT[m] for m in tags]
    end
    tri.segmentlist = segments
    tri.segmentmarkerlist = tags
    if !isnothing(holes)
        tri.holelist = holes
    end
    maxarea = Printf.@sprintf "%0.15f" h^2 / 2
    minangle = Printf.@sprintf "%0.15f" 30.0
    (tri, _) = triangulate("pena$(maxarea)q$(minangle)Q", tri)
    return HPMesh{F, I, P}(tri)
end

"""
    setboundary!(s,condition,m::HPMesh)
sets the boundary edges of mesh `m` which satisfy `condition` to kind `s`, where `s` can be `:dirichlet`, `:neumann` (or alternatively, `1` for Dirichlet or `2` for Neumann). Notice that `condition` is evaluated at the vertices of each edge. An edge is tagged when both vertices satisfy `condition`.

See also [`setdirichlet!`](@ref) and [`setneumann!`](@ref).
"""
function setboundary!(
        i::Integer,
        condition::Function,
        m::HPMesh{F, I, P}
    ) where {F, I, P}
    (; points, edgelist) = m
    for ea in pairs(edgelist)
        e = first(ea)
        if tag(last(ea)) > 0 && condition(points[e[1]]) && condition(points[e[2]])
            settag!(last(ea), P(i))
        end
    end
    return
end
function setboundary!(
        s::Symbol,
        condition::Function,
        m::HPTriangulation
    )
    i = BOUNDARY_DICT[s]
    return setboundary!(i, condition, m)
end
"""
    setdirichlet!(condition,m::HPMesh)
Marks as `:dirichlet` all boundary edges of `m` satisfying `condition`. Notice that `condition` is evaluated at the vertices of each edge. An edge is tagged when both vertices satisfy `condition`.
See also [`setneumann!`](@ref).
"""
setdirichlet!(condition::Function, m::HPMesh) = setboundary!(1, condition, m)
"""
    setneumann!(condition,m::HPMesh)
Marks as `:neumann` all boundary edges of `m` satisfying `condition`. Notice that `condition` is evaluated at the vertices of each edge. An edge is tagged when both vertices satisfy `condition`.
See also [`setdirichlet!`](@ref).
"""
setneumann!(condition::Function, m::HPMesh) = setboundary!(2, condition, m)


############################################
####        Auxiliary functions         ####
############################################
"""
    maybeconvert(::Type{T},arr) where T
If necessary converts the elements of `arr` to type `T`.  
"""
function maybeconvert(::Type{T}, arr::AbstractArray) where {T}
    return eltype(arr) == T ? arr : convert.(T, arr)
end

"""
    degtype(::HPMesh{F,I,P}) wehere {F,I,P}
returns `P`  
"""
degtype(::HPMesh{F, I, P}) where {F, I, P} = P
"""
    inttype(::HPMesh{F,I,P}) wehere {F,I,P}
returns `I`  
"""
inttype(::HPMesh{F, I, P}) where {F, I, P} = I
"""
    floattype(::HPMesh{F,I,P}) wehere {F,I,P}
returns `F`  
"""
floattype(::HPMesh{F, I, P}) where {F, I, P} = F

"""
    gen(m::HPMesh)
returns the generation number of `m`. The generation number is bumped after every refinement process.     
"""
gen(m::HPMesh) = m.gen[]
"""
    upgen!(m::HPMesh)
Bumps the generation of a mesh. This is done after refinement.
"""
upgen!(m::HPMesh{F, I, P}) where {F, I, P} = m.gen[] += one(I)

"""
    islocked(m::HPMesh)
Indicates if the mesh is locked to prevent direct modification.
"""
islocked(m::HPMesh) = m.lock[]
"""
    lock!(m::HPMesh)
locks the mesh to prevent direct modification.
"""
lock!(m::HPMesh) = m.lock[] = true

function Base.copy(mesh::HPMesh{F, I, P}) where {F, I, P}
    return HPMesh{F, I, P}(
        deepcopy(mesh.points),
        deepcopy(mesh.trilist),
        deepcopy(mesh.edgelist),
        copy(mesh.dofs),
        gen(mesh),
        islocked(mesh)
    )
end

"""
  edges(m)
retuns the edges of `m`, being `m` an `HPMesh` or an `EdgeList`.  
"""
@inline edges(list::T) where {T <: EdgeList} = keys(list)
@inline edges(mesh::T) where {T <: HPMesh} = keys(mesh.edgelist)
@inline edgedict(mesh::HPMesh) = mesh.edgelist
"""
  triangles(m) 
retuns the triangles of `m`, being `m` and `HPMesh` or a `TriList`.
"""
@inline triangles(list::T) where {T <: TriangleList} = keys(list)
@inline triangles(mesh::T) where {T <: HPMesh} = keys(mesh.trilist)

@inline tridict(mesh::HPMesh) = mesh.trilist

"""
  dof(m)
retuns the degrees of freedom of the mesh `m`.
"""
@inline dof(mesh) = mesh.dofs


"""
    $(SIGNATURES)

Checks if 'e' is stored as presented or in reverse order. 
"""
function _same_order(e::Edge{I}, elist::EdgeList{I, P}) where {I, P}
    _, t = gettoken(elist, e)
    oe = gettokenvalue(keys(elist), t)
    return oe == e
end

"""
    $(SIGNATURES)

returns a permutation of contiguous indices `ind` such that `v[ind]` is ordered. 
"""
function psortperm(v)
    i = argmin(v)
    if v[i] ≤ v[mod1(i + 1, 3)] ≤ v[mod1(i + 2, 3)]
        ind = mod1.(i:(i + 2), 3)
    else
        ind = mod1.(i:-1:(i - 2), 3)
    end
    return ind
end

"""
  $(SIGNATURES)  
Given a triangle  `t`, of  a mesh `mesh` returs the degrees of the edges of `t`
"""
function degrees(t::Triangle{I}, mesh::HPMesh{F, I, P}) where {F, I, P}
    (; edgelist) = mesh
    eds = edges(t)
    p = Tuple(degree.(getindices(edgelist, eds)))
    return sort(p)
end

"""
  $(SIGNATURES)

Given a triangle `t` belonging to a `mesh`, it returns the (sorted) degrees
`p₁<=p₂<=p₃` of its edges, and the edges also sorted according to the degrees.
"""
function psortededges(t::Triangle{I}, mesh::HPMesh{F, I, P}) where {F, I, P}
    p, nod = psortednodes(t, mesh)
    eds = [Edge(nod[SVector(i, mod1(i + 1, 3))]) for i in 1:3]
    return p, eds
end

"""
   $(SIGNATURES)

Given a triangle `t` belonging to a `mesh`, it returns the degrees `p₁<=p₂<=p₃` of the edges (sorted) and the nodes of `t` also sorted according to the degrees. 
"""
function psortednodes(t::Triangle{I}, mesh::HPMesh{F, I, P}) where {F, I, P}
    (; edgelist) = mesh
    eds = edges(t)
    p = degree.(getindices(edgelist, eds))
    pind = psortperm(p)
    return p[pind], first.(eds[pind])
end
psortednodes(e::Edge{I}, mesh::HPMesh{F, I, P}) where {F, I, P} = degree(mesh.edgelist[e]), e


"""
    _boundary_segments(n)
creates pairs of indices for building a sequence of segments from a list of points.  
"""
_boundary_segments(n) = mapreduce(i -> [i, mod1(i + 1, n)], hcat, 1:n)


"""
    degrees_of_freedom_by_edge!(mesh::HPMesh{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}

Creates a dictionary (from `Dictionaries.jl`) where the keys are the edges of `mesh` and the values are vectors with indices corresponding to the nodal degrees of freedom. 
"""
function degrees_of_freedom_by_edge!(mesh::HPMesh{F, I, P}) where {F, I, P}
    (; points, edgelist, dofs) = mesh
    (; by_edge) = dofs
    i = size(points, 2) + 1
    for edge in edges(edgelist)
        med = I(i):I(i + degree(edgelist[edge]) - 2)
        v = SVector{length(med) + 2, I}(edge[1], med..., edge[2])
        set!(by_edge, edge, v)
        i += degree(edgelist[edge]) - 1
    end
    return i
end


"""
    degrees_of_freedom!(mesh::HPMesh{F,I,P}) where {F,I,P}

Creates a dictionary (from `Dictionaries.jl`) where the keys are the triangles of the mesh and the values are vectores storing the indices of the corresponding degrees of freedom. 

Internally, `degrees_of_freedom_by_edge` is called in order to obtain the nodal degrees of freedom, and then the degrees of freedom corresponding to bubble functions are computed. 

If a dictionary of degrees of freedom by edge has already been computed, it is recommended to run: 

    degrees_of_freedom!(mesh::HPMesh{F,I,P},by_edge::Dictionary{Edge{I},Vector{I}}) where {F,I,P}
"""
function degrees_of_freedom!(mesh::HPMesh{F, I, P}) where {F, I, P}
    (; edgelist, trilist, dofs) = mesh
    (; n, by_edge, by_tri) = dofs
    if gen(mesh) != gen(dofs)
        empty!(dofs)
        degrees_of_freedom_by_edge!(mesh)
        k = maximum(maximum.(by_edge)) + 1 #first non-edge dof
        for t in triangles(trilist)
            p, t_edges = psortededges(t, mesh)
            newdofs = @MVector zeros(I, compute_dimension(p))
            j = 1 #counter of dof in current triangle
            @inbounds for i in 1:3
                newdof = by_edge[t_edges[i]]
                if _same_order(t_edges[i], edgelist)
                    newdofs[j:(j + length(newdof) - 2)] .= newdof[1:(end - 1)]
                else
                    newdofs[j:(j + length(newdof) - 2)] .= reverse(newdof[2:end])
                end
                j += length(newdof) - 1
            end
            newdofs[j:end] = k:(k + (length(newdofs) - j))
            set!(by_tri, t, newdofs)
            k += length(newdofs) - j + 1
            setgen!(dofs, gen(mesh))
        end
        n[] = k - 1
    end
    return n[]
end


"""
    tagged_dof(mesh::HPMesh{F,I,P},tag) where {F,I,P}

Returs a list of indices corresponding to the degrees of freedom labeled with `tag`. If a vector of labels is passed, it returs degrees of freedom marked with any of them. 
"""
function tagged_dof(mesh::HPMesh{F, I, P}, tag::N) where {F, I, P, N <: Integer}
    return tagged_dof(mesh, [tag])
end


"""
    tagged_dof(mesh::HPMesh,tags)
returns a list of all degrees of freedom tagged with a tag in `tag`.
"""
function tagged_dof(mesh::HPMesh{F, I, P}, tags) where {F, I, P}
    (; dofs) = mesh
    msg = "The degrees of freedom of the mesh have not been computed."
    isempty(dofs) && throw(ArgumentError(msg))
    return _tagged_dof(mesh, tags)
end
function tagged_dof(cond::Function, mesh::HPMesh{F, I, P}) where {F, I, P}
    (; dofs) = mesh
    msg = "The degrees of freedom of the mesh have not been computed."
    isempty(dofs) && throw(ArgumentError(msg))
    return _tagged_dof(cond, mesh)
end

"""
    _tagged_dof(mesh::HPMesh,tags)
internal function. Returns a list of all degrees of freedom tagged with a tag in `tag`.  
"""
_tagged_dof(mesh::HPMesh, tags) = _tagged_dof(in(tags), mesh)
function _tagged_dof(cond::Function, mesh::HPMesh{F, I, P}) where {F, I, P}
    (; edgelist, dofs) = mesh
    (; by_edge) = dofs
    v = Vector{I}()
    for e in pairs(edgelist)
        if cond(tag(last(e)))
            push!(v, by_edge[first(e)]...)
        end
    end
    return unique(v)
end

dof(t::Triangle, mesh::HPMesh) = mesh.dofs.by_tri[t]
dof(e::Edge, mesh::HPMesh) = mesh.dofs.by_edge[e]
ndof(mesh::HPMesh) = mesh.dofs.n[]
elements(mesh::HPMesh) = keys(mesh.trilist)

########
