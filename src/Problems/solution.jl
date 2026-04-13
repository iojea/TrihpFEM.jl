abstract type Solution end

struct FESolution{O<:DiffOperator,M<:HPMesh,V<:AbstractVector} <: Solution
    mesh::M
    vals::V
end

FESolution(m::M,v::V) where {M,V} = FESolution{Identity,M,V}(m,v)
function (::Gradient)(s::FESolution{Identity,M,V}) where {M,V}
    FESolution{Gradient,M,V}(s.mesh,s.vals)
end

underlyingmesh(s::FESolution) = s.mesh
values(s::FESolution) = s.vals


struct Error{S<:Solution,F<:Function,E}
    uₕ::S
    u::F
end
Error(d,m) = Error{typeof(d),typeof(m),1}(d,m)
Error(d,m,p) = Error{typeof(d),typeof(m),p}(d,m)

Base.:-(uₕ::Solution,u::Function) = Error(uₕ,u)
Base.:-(u::Function,uₕ::Solution) = Error(uₕ,u)
Base.:^(d::Error{S,F,1},e) where {S,F} = Error{S,F,e}(d.uₕ,d.u)

exponent(s::Error{S,F,E}) where {S,F,E} = E
exact(s::Error) = s.u
numerical(s::Error) = s.uₕ

Forms.∫(d::Error) = d

Base.:*(d::Error,m::Measure) = compute_error(d::Error,m::Measure)
function compute_error(d::Error{D,F,E},m::Measure) where {D,F,E}
    u = exact(d)
    uₕ = numerical(d)
    (;mesh,aux,sch) = m
    underlyingmesh(uₕ) === mesh || throw(ArgumentError("The numerical solution and the `Measure` are defined on different meshes."))
    err = zero(floattype(mesh))
    for el in Measures.elements(mesh)
        degs, _ = psortednodes(el, mesh)
        aff = AffineToRef(mesh.points[el])
        locdof = dof(el,mesh)
        dim = length(locdof)
        C = aux[degs].C
        U(x)   = (uₕ.vals[locdof]⋅(C'*[φ(x) for φ in StandardBasis(degs)])-u(aff(x)))^E
        err += jac(aff)*ref_integrate(U,sch)
    end
    err
end 

