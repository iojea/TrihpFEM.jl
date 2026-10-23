abstract type AbstractMeasure{F <: Real, I <: Integer, P <: Integer} end

struct Measure{F <: Real, I <: Integer, P <: Integer, Q <: Quadrature{2}} <: AbstractMeasure{F, I, P}
    mesh::HPMesh{F, I, P}
    sch::Q
end

function Measure(mesh::HPMesh{F, I, P}, degsch) where {F <: Real, I <: Integer, P <: Integer}
    sch = quadrature(F, P, Val(2), Val(degsch))
    return Measure{F, I, P, typeof(sch)}(mesh, sch)
end

function Measure(mesh::HPMesh{F, I, P}) where {F, I, P}
    (; edgelist) = mesh
    maxdeg = maximum(Meshes.degree.(edgelist))
    return Measure(mesh, 2maxdeg + 1)
end


struct BoundaryMeasure{F <: Real, I <: Integer, P <: Integer, Q <: Quadrature{1}} <: AbstractMeasure{F, I, P}
    kind::P
    mesh::HPMesh{F, I, P}
    sch::Q
end

function BoundaryMeasure(mesh::HPMesh{F, I, P}, kind = Integer) where {F, I, P}
    (; edgelist) = mesh
    maxdeg = maximum(Meshes.degree.(edgelist))
    return BoundaryMeasure(kind, mesh, 2maxdeg + 1)
end


mesh(m::AbstractMeasure) = m.mesh

elements(m::Measure) = m.mesh.trilist
elements(bm::BoundaryMeasure) = filter(e -> tag(e) == bm.kind, bm.mesh.edgelist)

dim(::Measure) = 2
dim(::BoundaryMeasure) = 1

basis(t::Triangle, m::AbstractMeasure) = StandardBasis(first(psortednodes(t, m.mesh)))
function basis(e::Edge, m::BoundaryMeasure)
    els = elements(m)
    isin, tok = gettoken(els, e)
    isin || throw(ArgumentError("The `Edge` does not belong to the given `BoundaryMeasure`."))
    att = gettokenvalue(els, tok)
    adj = adjacents(att)
    last(adj) == 0 || throw(ArgumentError("`Edge` seems to be interior to the domain."))
    t = Triangle(e[1], e[2], adj[1])
    p = first(psortednodes(t, m.mesh))
    return StandardBasis(p)
end

function dof(ed::Edge, m::Measure)
    indofs, token = gettoken(m.mesh.dofs.by_edge, ed)
    indofs || throw(ArgumentError("No degrees of freedom are associated to the edge $ed.l Maybe you need to compute the degrees of freedom of the mesh, using `degrees_of_freedom!`."))
    return gettokenvalue(m.mesh.dofs.by_edge, token)
end

function dof(t::Triangle, m::Measure)
    indofs, token = gettoken(m.mesh.dofs.by_tri, t)
    indofs || throw(ArgumentError("No degrees of freedom are associated to the triangle $t.l Maybe you need to compute the degrees of freedom of the mesh, using `degrees_of_freedom!`."))
    return gettokenvalue(m.mesh.dofs.by_tri, token)
end
