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


mesh(m::AbstractMeasure) = m.mesh

elements(m::Measure) = m.mesh.trilist
elementtype(::Measure{F, I, P}) where {F, I, P} = Triangle{I}

dim(::Measure) = 2


function Meshes.dof(ed::Edge, m::AbstractMeasure)
    indofs, token = gettoken(m.mesh.dofs.by_edge, ed)
    indofs || throw(ArgumentError("No degrees of freedom are associated to the edge $ed. gngn Maybe you need to compute the degrees of freedom of the mesh, using `degrees_of_freedom!`."))
    return gettokenvalue(m.mesh.dofs.by_edge, token)
end

function Meshes.dof(t::Triangle, m::AbstractMeasure)
    indofs, token = gettoken(m.mesh.dofs.by_tri, t)
    indofs || throw(ArgumentError("No degrees of freedom are associated to the triangle $t.l Maybe you need to compute the degrees of freedom of the mesh, using `degrees_of_freedom!`."))
    return gettokenvalue(m.mesh.dofs.by_tri, token)
end

Meshes.dof(am::T) where {T <: AbstractMeasure} = dof(am.mesh)
Meshes.ndof(am::T) where {T <: AbstractMeasure} = ndof(am.mesh)
Meshes.degrees_of_freedom!(am::AbstractMeasure) = degrees_of_freedom!(am.mesh)
Meshes.floattype(::AbstractMeasure{F, I, P}) where {F, I, P} = F
Meshes.inttype(::AbstractMeasure{F, I, P}) where {F, I, P} = I
Meshes.degtype(::AbstractMeasure{F, I, P}) where {F, I, P} = P

tridof(m::T) where {T <: AbstractMeasure} = m.mesh.dofs.by_tri
