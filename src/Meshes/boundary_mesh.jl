"""

    BoundaryHPMesh{F,I,P} <: HPTriangulation

A `struct` for accessing a part of the boundary of a domain.

```julia
  julia> Ω = circmesh(0.1);
  julia> ∂Ω = BoundaryHPMesh(mesh,1)
```
"""
struct BoundaryHPMesh{F,I,P} <:HPTriangulation
    mesh::HPMesh{F,I,P}
    marker::P
end

"""

    domainmesh(m::Triangulation)

returns the underlying mesh in an `Triangulation`. If `m` is an `HPMesh`, then `domainmesh(m)===m`. If `m` is a `BoundaryHPMesh`, it returns the mesh of the domain part of whose boundary is represented by `m`.
"""
domainmesh(m::BoundaryHPMesh) = m.mesh
domainmesh(m::HPMesh) = m



