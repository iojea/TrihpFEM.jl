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
    mesh::HPMesh{F, I, P}
    kind::P
    sch::Q
end

BoundaryMeasure(m::HPMesh{F,I,P},kind::Integer,sch::Q) where {F,I,P,Q<:Quadrature{1}} = BoundaryMeasure{F,I,P,Q}(m,kind,sch)

function BoundaryMeasure(mesh::HPMesh{F, I, P}, kind:: Integer) where {F, I, P}
    (; edgelist) = mesh
    maxdeg = maximum(Meshes.degree.(edgelist))
    return BoundaryMeasure(mesh, kind,2maxdeg + 1)
end

function BoundaryMeasure(mesh::HPMesh{F, I, P}, kind:: Integer,deg::Integer) where {F, I, P}
    (; edgelist) = mesh
    sch = quadrature(F,P,Val(1),Val(deg))
    return BoundaryMeasure(mesh,kind,sch)
end


mesh(m::AbstractMeasure) = m.mesh

elements(m::Measure) = pairs(m.mesh.trilist)
elements(bm::BoundaryMeasure) = pairs(filter(e -> tag(e) == bm.kind, bm.mesh.edgelist))

elementtype(::Measure{F,I,P}) where {F,I,P} = Triangle{I}
elementtype(::BoundaryMeasure{F,I,P}) where {F,I,P} = Edge{I}

dim(::Measure) = 2
dim(::BoundaryMeasure) = 1



function _extract(t::Triangle,m::Measure{F,I,P}) where {F,I,P}
    p,_ = psortednodes(t,mesh(m))
    return ((zero(P),p...),StandardBasis(p))
end
function _extract(e::Edge,m::BoundaryMeasure{F,I,P}) where {F,I,P}
    els = elements(m)
    isin,tok = gettoken(els,e)
    isin || throw(ArgumentError("The `Edge` does not belong to the given `BoundaryMeasure`"))
    _,att = gettokenvalue(els,tok)
    adj = adjacents(att)
    last(adj) == 0 || throw(ArgumentError("`Edge` seems to be interior to the domain."))
    t = Triangle(e[1], e[2], adj[1])
    p,t_eds = psortededges(t, mesh(m))
    j = findfirst(isequal(e),t_eds)
    return ((P(j),p...),StandardBasis(p))
end

function dof(ed::Edge, m::Measure)
    indofs, token = gettoken(m.mesh.dofs.by_edge, ed)
    indofs || throw(ArgumentError("No degrees of freedom are associated to the edge $ed. gngn Maybe you need to compute the degrees of freedom of the mesh, using `degrees_of_freedom!`."))
    return gettokenvalue(m.mesh.dofs.by_edge, token)
end

function dof(t::Triangle, m::Measure)
    indofs, token = gettoken(m.mesh.dofs.by_tri, t)
    indofs || throw(ArgumentError("No degrees of freedom are associated to the triangle $t.l Maybe you need to compute the degrees of freedom of the mesh, using `degrees_of_freedom!`."))
    return gettokenvalue(m.mesh.dofs.by_tri, token)
end
