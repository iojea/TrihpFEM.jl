struct Error{S <: Solution, F <: Function, E}
    uₕ::S
    u::F
end
Error(d, m) = Error{typeof(d), typeof(m), 1}(d, m)
Error(d, m, p) = Error{typeof(d), typeof(m), p}(d, m)

Base.:-(uₕ::Solution, u::Function) = Error(uₕ, u)
Base.:-(u::Function, uₕ::Solution) = Error(uₕ, u)
Base.:^(d::Error{S, F, 1}, e) where {S, F} = Error{S, F, e}(d.uₕ, d.u)

exponent(s::Error{S, F, E}) where {S, F, E} = E
exact(s::Error) = s.u
numerical(s::Error) = s.uₕ
Forms.∫(d::Error) = d

Base.:*(d::Error, m::Measure) = compute_error(d::Error, m::Measure)
function compute_error(d::Error{D, F, E}, m::Measure) where {D, F, E}
    u = exact(d)
    uₕ = numerical(d)
    (; mesh, aux, sch) = m
    underlyingmesh(uₕ) === mesh || throw(ArgumentError("The numerical solution and the `Measure` are defined on different meshes."))
    err = zero(floattype(mesh))
    for el in Measures.elements(mesh)
        degs, _ = psortednodes(el, mesh)
        aff = AffineToRef(mesh.points[el])
        locdof = dof(el, mesh)
        dim = length(locdof)
        C = aux[degs].C
        U(x) = (uₕ.vals[locdof] ⋅ (C' * [φ(x) for φ in StandardBasis(degs)]) - u(aff(x)))^E
        err += jac(aff) * ref_integrate(U, sch)
    end
    return err
end

"""

    estimate_order(problem,u;boundary_projection=nothing,E=2,iterations=3,deg=11)
solves `problem` iteratively refining all triangles in the mesh. After each iteration the `E` norm of the error is computed, with respect to the exact solution `u`.  The optional argument `boundary_projection` is used to project boundary nodes to the real boundary (when working on non-polygonal domains). The estimated order of convergence is returned. 

"""
function estimate_order(problem, u; boundary_projection = nothing, E = 2, iterations = 3, deg = 11)
    uₕ = solve(problem)
    Ω = domainmesh(problem)
    dΩ = Measure(Ω, deg)
    e = zeros(iterations)
    e[1] = (∫(u - uₕ)^E * dΩ)^(1 / E)
    for k in 2:iterations
        Meshes.mark!.(Ω.edgelist)
        Meshes._h_conformity!(Ω)
        Meshes.refine!(Ω)
        if !isnothing(boundary_projection)
            boundary_projection(Ω)
        end
        update!(problem)
        uₕ = solve(problem)
        e[k] = (∫(u - uₕ)^E * dΩ)^(1 / E)
    end
    ord = (log.(e[2:end]) - log.(e[1:(end - 1)])) ./ (-log(2))
    return e, ord
end
