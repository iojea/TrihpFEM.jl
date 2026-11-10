abstract type Solution end

struct FESolution{O <: DiffOperator, M <: HPMesh, V <: AbstractVector} <: Solution
    mesh::M
    vals::V
end

FESolution(m::M, v::V) where {M, V} = FESolution{Identity, M, V}(m, v)
function (::Gradient)(s::FESolution{Identity, M, V}) where {M, V}
    return FESolution{Gradient, M, V}(s.mesh, s.vals)
end

function (::Laplacian)(s::FESolution{Identity, M, V}) where {M, V}
    return FESolution{Laplacian, M, V}(s.mesh, s.vals)
end


Meshes.mesh(s::FESolution) = s.mesh
values(s::FESolution) = s.vals
Forms.operator(::FESolution{O}) where O = O()
Forms.chain_operator(::FESolution{Gradient}) = DiffMatrix()
Forms.chain_operator(::FESolution{AdjointGradient}) = AdjointDiffMatrix()
Forms.chain_operator(::FESolution{Identity}) = Identity()
Forms.chain_operator(::FESolution{Laplacian}) = nothing
