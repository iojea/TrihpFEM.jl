struct SkeletonMeasure{F <: Real, I <: Integer, P <: Integer, Q <: QuadratureRule{1}} <: AbstractMeasure{F, I, P}
    mesh::HPMesh{F, I, P}
    kind::P
    sch::Q
end

SkeletonMeasure(m::HPMesh{F, I, P}, kind::Integer, sch::Q) where {F, I, P, Q <: Quadrature{1}} = SkeletonMeasure{F, I, P, Q}(m, kind, sch)

function SkeletonMeasure(mesh::HPMesh{F, I, P}, kind::Integer) where {F, I, P}
    (; edgelist) = mesh
    maxdeg = maximum(Meshes.degree.(edgelist))
    return SkeletonMeasure(mesh, kind, 2maxdeg + 1)
end

function SkeletonMeasure(mesh::HPMesh{F, I, P}, kind::Integer, deg::Integer) where {F, I, P}
    (; edgelist) = mesh
    sch = rule(Interval();degree=deg)
        #quadrature(F, P, Val(1), Val(deg))
    return SkeletonMeasure(mesh, kind, sch)
end
elements(bm::SkeletonMeasure) = filter(e -> tag(e) == bm.kind, bm.mesh.edgelist)
elementtype(::SkeletonMeasure{F, I, P}) where {F, I, P} = Edge{I}
dim(::SkeletonMeasure) = 1
