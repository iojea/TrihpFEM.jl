struct BoundaryMeasure{F <: Real, I <: Integer, P <: Integer, Q <: QuadratureRule{1}} <: AbstractMeasure{F, I, P}
    mesh::HPMesh{F, I, P}
    kind::P
    sch::Q
end

BoundaryMeasure(m::HPMesh{F, I, P}, kind::Integer, sch::Q) where {F, I, P, Q <: QuadratureRule{1}} = BoundaryMeasure{F, I, P, Q}(m, kind, sch)

function BoundaryMeasure(mesh::HPMesh{F, I, P}, kind::Integer) where {F, I, P}
    (; edgelist) = mesh
    maxdeg = maximum(Meshes.degree.(edgelist))
    return BoundaryMeasure(mesh, kind, 2maxdeg + 1)
end

function BoundaryMeasure(mesh::HPMesh{F, I, P}, kind::Integer, deg::Integer) where {F, I, P}
    (; edgelist) = mesh
    sch = rule(Interval();degree=deg)
        #quadrature(F, P, Val(1), Val(deg))
    return BoundaryMeasure(mesh, kind, sch)
end
elements(bm::BoundaryMeasure) = filter(e -> tag(e) == bm.kind, bm.mesh.edgelist)
elementtype(::BoundaryMeasure{F, I, P}) where {F, I, P} = Edge{I}
dim(::BoundaryMeasure) = 1
