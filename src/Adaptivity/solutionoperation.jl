struct SolOperation{C <: CoeffType, U <: Tuple,E} <: AbstractOperation{C,U}
    parts::U
    SolOperation{C}(t...) where {C <: CoeffType} = new{C, typeof(t),1}(t)
end

const NumArray = Union{Number, AbstractArray}

Forms.isadjoint(s::FESolution) = isadjoint(operator(s))
function Forms.unroll(s::FESolution)
    return isadjoint(s) ? (s, *, chain_operator(s)) : (chain_operator(s), *, s)
end


Base.:*(a::NumArray, b::FESolution) = SolOperation{ConstantCoeff}(a, *, unroll(b)...)
Base.:*(a::FESolution, b::NumArray) = SolOperation{ConstantCoeff}(unroll(b)..., *, a)
Base.:*(a::FESolution, b::FESolution) = SolOperation{ConstantCoeff}(unroll(a)..., *, unroll(b)...)
Base.:*(a::Function, b::FESolution) = SolOperation{VariableCoeff}(a, *, unroll(b)...)
Base.:*(a::FESolution, b::Function) = SolOperation{VariableCoeff}(unroll(a)..., *, b)
Base.:*(a::SolOperation{C}, b::FESolution) where {C <: CoeffType} = SolOperation{C}(a.parts..., *, unroll(b)...)
Base.:*(a::FESolution, b::SolOperation{C}) where {C <: CoeffType} = SolOperation{C}(unroll(a)..., *, b)
Base.:*(a::NumArray, b::SolOperation{C}) where {C <: CoeffType} = SolOperation{C}(a, *, b.parts...)
Base.:*(a::SolOperation{C}, b::NumArray) where {C <: CoeffType} = SolOperation{C}(a.parts..., *, b)
function Base.:*(a::SolOperation{C₁}, b::SolOperation{C₂}) where {C₁ <: CoeffType, C₂ <: CoeffType}
    C = promote_type(C₁, C₂)
    return SolOperation{C}(a.parts..., *, b.parts...)
end

Base.:+(a::NumArray, b::FESolution) = SolOperation{ConstantCoeff}(a, +, unroll(b)...)
Base.:+(a::FESolution, b::NumArray) = SolOperation{ConstantCoeff}(unroll(b)..., +, a)
Base.:+(a::FESolution, b::FESolution) = SolOperation{ConstantCoeff}(unroll(a)..., +, unroll(b)...)
Base.:+(a::Function, b::FESolution) = SolOperation{VariableCoeff}(a, +, unroll(b)...)
Base.:+(a::FESolution, b::Function) = SolOperation{VariableCoeff}(unroll(a)..., +, b)
Base.:+(a::SolOperation{C}, b::FESolution) where {C <: CoeffType} = SolOperation{C}(a.parts..., +, unroll(b)...)
Base.:+(a::FESolution, b::SolOperation{C}) where {C <: CoeffType} = SolOperation{C}(unroll(a)..., +, b)
Base.:+(a::NumArray, b::SolOperation{C}) where {C <: CoeffType} = SolOperation{C}(a, +, b.parts...)
Base.:+(a::SolOperation{C}, b::NumArray) where {C <: CoeffType} = SolOperation{C}(a.parts..., +, b)
function Base.:+(a::SolOperation{C₁}, b::SolOperation{C₂}) where {C₁ <: CoeffType, C₂ <: CoeffType}
    C = promote_type(C₁, C₂)
    return SolOperation{C}(a.parts..., +, b.parts...)
end

Base.:-(a::NumArray, b::FESolution) = SolOperation{ConstantCoeff}(a, -, unroll(b)...)
Base.:-(a::FESolution, b::NumArray) = SolOperation{ConstantCoeff}(unroll(b)..., -, a)
Base.:-(a::FESolution, b::FESolution) = SolOperation{ConstantCoeff}(unroll(a)..., -, unroll(b)...)
Base.:-(a::Function, b::FESolution) = SolOperation{VariableCoeff}(a, -, unroll(b)...)
Base.:-(a::FESolution, b::Function) = SolOperation{VariableCoeff}(unroll(a)..., -, b)
Base.:-(a::SolOperation{C}, b::FESolution) where {C <: CoeffType} = SolOperation{C}(a.parts..., -, unroll(b)...)
Base.:-(a::FESolution, b::SolOperation{C}) where {C <: CoeffType} = SolOperation{C}(unroll(a)..., -, b)
Base.:-(a::NumArray, b::SolOperation{C}) where {C <: CoeffType} = SolOperation{C}(a, -, b.parts...)
Base.:-(a::SolOperation{C}, b::NumArray) where {C <: CoeffType} = SolOperation{C}(a.parts..., -, b)
function Base.:-(a::SolOperation{C₁}, b::SolOperation{C₂}) where {C₁ <: CoeffType, C₂ <: CoeffType}
    C = promote_type(C₁, C₂)
    return SolOperation{C}(a.parts..., -, b.parts...)
end


Base.adjoint(op::SolOperation{C}) where {C <: CoeffType} = SolOperation{C}(reverse(_adj.(op.parts))...)

LinearAlgebra.dot(a::NumArray, b::FESolution) = a' * b
LinearAlgebra.dot(a::FESolution, b::NumArray) = a' * b
LinearAlgebra.dot(a::FESolution, b::FESolution) = a' * b
LinearAlgebra.dot(a::Function, b::FESolution) = (transpose ∘ a) * b
LinearAlgebra.dot(a::FESolution, b::Function) = a' * b
LinearAlgebra.dot(a::SolOperation, b::FESolution) = a' * b
LinearAlgebra.dot(a::FESolution, b::SolOperation) = a' * b
LinearAlgebra.dot(a::SolOperation, b::SolOperation) = a' * b

Base.:-(u::FESolution) = (-1) * u

get_sol(op::SolOperation) = filter(Base.Fix2(isa, FESolution), op.parts)

Base.:^(d::SolOperation{S, F, 1}, e) where {S, F} = SolOperation{S, F, e}(d)

exponent(::SolOperation{S, F, E}) where {S, F, E} = E

# Base.:*(i::Integrand{C,T,N,SolOperation{C,T}}, m::Measure) where {C,T,N} = compute_error(i.operation, m::Measure)

# function compute_error(d::SolOperation{D, F, E}, m::DomainMeasure) where {D, F, E}
#     u = exact(d)
#     uₕ = numerical(d)
#     (; mesh, sch) = m
#     underlyingmesh(uₕ) === mesh || throw(ArgumentSolOperation("The numerical solution and the `Measure` are defined on different meshes."))
#     err = zero(floattype(mesh))
#     for el in Measures.elements(mesh)
#         degs, _ = psortednodes(el, mesh)
#         aff = AffineToRef(mesh.points[el])
#         locdof = dof(el, mesh)
#         dim = length(locdof)
#         C = aux[degs].C
#         U(x) = (uₕ.vals[locdof] ⋅ (C' * [φ(x) for φ in StandardBasis(degs)]) - u(aff(x)))^E
#         err += jac(aff) * ref_integrate(U, sch)
#     end
#     return err
# end

"""

    estimate_order(problem,u;boundary_projection=nothing,E=2,iterations=3,deg=11)
solves `problem` iteratively refining all triangles in the mesh. After each iteration the `E` norm of the error is computed, with respect to the exact solution `u`.  The optional argument `boundary_projection` is used to project boundary nodes to the real boundary (when working on non-polygonal domains). The estimated order of convergence is returned. 

"""
# function estimate_order(problem, u; boundary_projection = nothing, E = 2, iterations = 3, deg = 11)
#     uₕ = solve(problem)
#     Ω = domainmesh(problem)
#     dΩ = Measure(Ω, deg)
#     e = zeros(iterations)
#     e[1] = (∫(u - uₕ)^E * dΩ)^(1 / E)
#     for k in 2:iterations
#         Meshes.mark!.(Ω.edgelist)
#         Meshes._h_conformity!(Ω)
#         Meshes.refine!(Ω)
#         if !isnothing(boundary_projection)
#             boundary_projection(Ω)
#         end
#         update!(problem)
#         uₕ = solve(problem)
#         e[k] = (∫(u - uₕ)^E * dΩ)^(1 / E)
#     end
#     ord = (log.(e[2:end]) - log.(e[1:(end - 1)])) ./ (-log(2))
#     return e, ord
# end
