abstract type Solution end

struct FESolution{O <: DiffOperator, M <: HPMesh, V <: AbstractVector} <: Solution
    mesh::M
    vals::V
end

FESolution(m::M, v::V) where {M, V} = FESolution{Identity, M, V}(m, v)
function (::Gradient)(s::FESolution{Identity, M, V}) where {M, V}
    return FESolution{Gradient, M, V}(s.mesh, s.vals)
end

underlyingmesh(s::FESolution) = s.mesh
values(s::FESolution) = s.vals
