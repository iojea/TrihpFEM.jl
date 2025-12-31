"""
    BoundaryHPMesh(mesh::HPMesh,kind)

See also [`dirichletboundary`](@ref), [`neumannboundary`](@ref).

creates a lazy container for accessing the part of the boundary of a domain `mesh` with boundary type `kind`. 

# Examples:

```@docs
julia> Ω = circmesh(0.1);

julia> ∂Ω = BoundaryHPMesh(Ω,:dirichlet);
```

Returns the Dirichlet boundary of `Ω`. This is equivalent to:

```@docs
julia> Ω = circmesh(0.1);

julia> ∂Ω = BoundaryHPMesh(Ω,1);
```

A more complex situation is:

```@docs
julia> vert = [0. 0.;1. 0.;1. 1.;0. 1.]';

julia> segs = [1 2;2 3;3 4;4 1]';

julia> tags = [1,2,1,2];

julia> Ω = hpmesh(vert,0.1;segments=segs,tags=tags);

julia> ΓD = BoundaryHPMesh(Ω,:dirichlet);

julia> ΓN = BoundaryHPMesh(Ω,:neumann);
```
Alternatively, `1` can be used for `:dirichlet` and `2` for `:neumann`.

"""
struct BoundaryHPMesh{F,I,P} <:HPTriangulation
    mesh::HPMesh{F,I,P}
    kind::P
    function BoundaryHPMesh(m::HPMesh{F,I,P},i::Integer) where {F,I,P}
        new{F,I,P}(m,P(i))
    end
end

function BoundaryHPMesh(m::HPMesh{F,I,P},s::Symbol) where {F,I,P}
    BoundaryHPMesh(m,BOUNDARY_DICT[s])
end

"""
   dirichletboundary(m::HPMesh)
returns the Dirichlet boundary of `m` as a `BoundaryHPMesh`. 
"""
dirichletboundary(m::HPMesh) = BoundaryHPMesh(m,1)
"""
   dirichletboundary(m::HPMesh)
returns the Neumann boundary of `m` as a `BoundaryHPMesh`. 
"""
neumannboundary(m::HPMesh) = BoundaryHPMesh(m,2)



"""
    edges(bm::BoundaryMesh)
returns the edges of `bm` as an `EdgeList` (a `Dictionary{Edge,EdgeAttributes}`). 
"""
edges(bm::BoundaryHPMesh) = filter(istagged(bm.kind),bm.mesh.edgelist)


"""

    domainmesh(m::Triangulation)

returns the underlying mesh in an `Triangulation`. If `m` is an `HPMesh`, then `domainmesh(m)===m`. If `m` is a `BoundaryHPMesh`, it returns the mesh of the domain part of whose boundary is represented by `m`.
"""
domainmesh(m::BoundaryHPMesh) = m.mesh
domainmesh(m::HPMesh) = m



