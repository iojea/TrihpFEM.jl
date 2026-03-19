"""
   collapse(aff:AffineToRef,t::Tensor)
performs an appropriate contraction of the tensor `t` with the change of variables `aff`. 
"""
function collapse(aff::AffineToRef, t::Tensor{4, 2})
    iA = inv(aff.A)
    return iA'⊡t⊡iA
end

collapse(aff::AffineToRef, t::Tensor{2, 2}) = inv(aff.A) ⊡ t
collapse(aff::AffineToRef, t::Tensor{1, 1}) = one(SymmetricTensor{1, 1}) ⋅ t


"""
    _initvectors(I,F,ℓ)
creates two vectors of type `I` for indices, and a vector of type `F` for values, all of them with size `ℓ`.  
"""
function _initvectors(::HPMesh{F, I, P}, ℓ) where {F, I, P}
    ivec = FixedSizeArray{I, 1}(undef, ℓ)
    fill!(ivec, zero(I))
    jvec = FixedSizeArray{I, 1}(undef, ℓ)
    fill!(jvec, zero(I))
    vals = FixedSizeArray{F, 1}(undef, ℓ)
    fill!(vals, zero(F))
    return ivec, jvec, vals
end
function _init_rhs(::HPMesh{F,I,P},ℓ) where {F,I,P}
    vec = FixedSizeArray{F,1}(undef,ℓ)
    fill!(vec,zero(F))
end

"""
   tensorize(f)
When appropriate, converts `f` to a tensor or to a function returning a tensor. This is done in order to perform an outer product using `f`.

**Important:** `tensorize` transposes `f` if `f` is a matrix. 
"""
tensorize(factor::Number) = factor
function tensorize(factor::AbstractArray)
    ord = ndims(factor)
    return Tensor{ord, 2}(factor')
end

"""
    outer(operand,factor)
performs the appropriate operation for the tensor product between the `operand` that comes from integrating basis functions and `factor`, that can be a matrix, vector or number.   
"""
outer(operand::Tensor{2, 2}, factor::Tensor{2, 2}) = otimesl(operand, factor)
outer(operand::Tensor{1, 2}, factor::Tensor{1, 2}) = factor ⊗ operand
outer(operand::Tensor{1, 1}, factor) = factor * operand


"""
   ref_tensors(inte::Integrand,degs)
Compute the local tensors produced by `inte` over a reference element with degrees given by `degs`. These tensors are later contracted wih the corresponding `affineToRef` change of variables. 
"""
function ref_tensors(inte::Integrand{ConstantCoeff, T, Order{B}, 2}, degs) where {T, B}
    dim = sum(B) > 0 ? 2 : 1
    ord = max(1, sum(B))
    (; factor, funs) = inte
    b₁ = basis(funs[1], degs); b₂ = basis(funs[2], degs)
    o₁, o₂ = operator.(funs)
    f = tensorize(factor)
    return collect_as(FixedSizeArrayDefault, (outer(_tensor(o₁(φ), o₂(ψ), Val(sum(B)),f) for φ in b₁, ψ in b₂))
end

_tensor(f, g, ::Val{2}) = Tensor{2, 2}((i, j) -> Integration.ref_integrate(f[i] * g[j]))
_tensor(f, g, ::Val{1}) = Tensor{1, 2}(i -> Integration.ref_integrate(f[i] * g))
_tensor(f, g, ::Val{0}) = Tensor{1, 1}((Integration.ref_integrate(f * g),))

"""
   assembly_matrix(form::Form{2})
assembles the matrix corresponding to the bilinear form  `form`.
"""
function assembly_matrix(term::Term{C, O, T, 2, M}) where {C, O, T, 2, M}
    return assembly_matrix(Form{2}((term,)))
end
function assembly_matrix(form::Form{2})
    (; terms) = form
    mesh = domainmesh(first(terms))
    ℓ = degrees_of_freedom!(mesh)
    N = sum(map(length, mesh.dofs.by_tri) .^ 2)
    ivec, jvec, vals = _initvectors(mesh, N)
    for t in terms
        add_to_matrix!(ivec, jvec, vals, t)
    end
    return sparse(ivec, jvec, vals, ℓ, ℓ)
end

"""
    add_to_matrix!(ivec,jvec,vals,t::Term)
integrates `t` adding the results to `ivec`,`jvec` and `vals`, for later building a sparse matrix.
"""
function add_to_matrix!(ivec, jvec, vals, t::Term{ConstantCoeff, O, T, 2, M}) where {O, T, M}
    (; integrand, measure) = t
    (; aux, mesh) = measure
    (; points, trilist, dofs) = mesh
    (; by_tri) = dofs
    r = 1
    tensordict = Dictionary()
    for tri in keys(trilist)
        degs, _ = psortednodes(tri, mesh)
        isin, token = gettoken(tensordict, degs)
        if isin
            loctensor = gettokenvalue(tensordict, token)
        else
            loctensor = ref_tensors(integrand, degs)
            set!(tensordict, degs, loctensor)
        end
        aff = AffineToRef(points[tri])
        doft = by_tri[tri]
        dim = length(doft)
        C = aux[degs].C
        v = collect_as(FixedSizeArrayDefault, (collapse(aff, loctensor[i, j]) for i in 1:dim, j in 1:dim))
        v .= jac(aff) * C' * v * C
        println("dim:", dim)
        println("ℓ:", length(ivec))
        println("r:", r)
        ivec[r:(r + dim^2 - 1)] .+= repeat(doft, dim)
        jvec[r:(r + dim^2 - 1)] .+= repeat(doft, inner = dim)
        vals[r:(r + dim^2 - 1)] .+= v[:]
        r += dim^2
    end
    return
end


"""
   assembly_rhs(form::Form{1})
assembles the rhs corresponding to the bilinear form  `form`.
"""
function assembly_rhs(term::Term{C, O, T, 1, M}) where {C, O, T, 1, M}
    return assembly_rhs(Form{1}((term,)))
end
function assembly_rhs(form::Form{1})
    (; terms) = form
    mesh = domainmesh(first(terms))
    ℓ    = ndofs(mesh)
    vals = _init_rhs(mesh,N)
    for t in terms
        add_to_rhs!(vals, t)
    end
    return vals
end

"""
    add_to_rhs!(vals,t::Term)
integrates `t` adding the results to `vals`.
"""
function add_to_rhs!(vals, t::Term{ConstantCoeff, O, T, 1, M}) where {O, T, M}
    (; integrand, measure) = t
    (; aux, mesh) = measure
    (; points, trilist, dofs) = mesh
    (; by_tri) = dofs
    r = 1
    tensordict = Dictionary()
    for tri in keys(trilist)
        degs, _ = psortednodes(tri, mesh)
        isin, token = gettoken(tensordict, degs)
        if isin
            loctensor = gettokenvalue(tensordict, token)
        else
            loctensor = ref_tensors(integrand, degs)
            set!(tensordict, degs, loctensor)
        end
        aff = AffineToRef(points[tri])
        doft = by_tri[tri]
        dim = length(doft)
        C = aux[degs].C
        v = collect_as(FixedSizeArrayDefault, (collapse(aff, loctensor[i, j]) for i in 1:dim, j in 1:dim))
        v .= jac(aff) * C' * v * C
        println("dim:", dim)
        println("ℓ:", length(ivec))
        println("r:", r)
        ivec[r:(r + dim^2 - 1)] .+= repeat(doft, dim)
        jvec[r:(r + dim^2 - 1)] .+= repeat(doft, inner = dim)
        vals[r:(r + dim^2 - 1)] .+= v[:]
        r += dim^2
    end
    return
end


# oprod(t₁::T, t::S) where {
# T <: Tensor, S <: Tensor} = t₁ ⊗ t₂
# oprod(::Nothing, t) = t

# function build_local_tensor(term::Term{C,Order{B},N,M}) where {C<:Union{ConstantCoeff,NoCoeff},B,M}
#     (;integrand,measure) = term
#     (;factor,funs) = integrand
#     ord = 2^sum(B)
#     dim = sum(B)>0 ? 2 : 1
#     b = basis.(funs)
#     for (i,φs) in Iterators.product(b...)
#         K[i] = oprod(factor,ref_integrate(φs...))
#     end
# end


# function build_local_tensor(term::IntegrationTerm, ::Order{B}, base) where {B}
#     F = floattype(domainmesh(term))
#     n = length(base)
#     M = FizedSizeArrayDefault{Tensor{4, 2, F}, 2}(undef, (n, n))
#     for φ in base
#         for ψ in base
#             M[i, j] = term.constant ⊗ term.polyfun(φ, ψ)
#         end
#     end
#     return M
# end


# """
#   integrate(form::Form{2},space::AbstractSpace)

# Integrates the `Form` `form` using the basis of the space `space`.
# """
# function integrate(form::Form{2}, space::Spaces.AbstractSpace)
#     mesh = domainmesh(form)
#     N = length(dof(mesh))
#     ℓ = sum(Base.Fix{2}(^, 2) ∘ length, dof(mesh).by_tri)
#     ivec, jvec, vals = _initvectors(mesh, ℓ)
#     for term in terms(form)
#         mock = polyfun(space, space)
#         ord = order(mock)
#         buildmatrix!(term, ord, ivec, jvec, vals, space)
#     end
#     return sparse(ivec, jvec, vals, N, N)
# end


# # Integration with constant coefficients
# function buildmatrix!(term::IntegrationTerm{ConstantCoeff, 2}, ord::Order{B}, ivec, jvec, vals, space) where {B}
#     (; measure) = term
#     (; aux, mesh) = measure
#     (; points, trilist, dofs) = mesh
#     (; by_tri) = dofs
#     F = eltype(vals)
#     Itype = eltype(ivec)
#     tensordict = Dictionary{NTuple{3, Itype}, FixedSizeArray{F, 2 + sum(B)}}()
#     aff = AffineToRef{F}()
#     r = 1
#     for tri in keys(trilist)
#         p, _ = psortednodes(tri, mesh)
#         isin, token = gettoken(tensordict, p)
#         if isin
#             loctensor = gettokenvalue(tensordict, token)
#         else
#             loctensor = build_local_tensor(term, ord, basis(space, p))
#             set!(tensordict, p, loctensor)
#         end
#         affine!(aff, points[tri])
#         COl = collapser(ord, aff)
#         doft = by_tri[tri]
#         dim = length(doft)
#         C = aux[p].C
#         v .= jac(aff) * C' * v * C
#         ivec[r:(r + dim^2 - 1)] .+= repeat(doft, dim)
#         jvec[r:(r + dim^2 - 1)] .+= repeat(doft, inner = dim)
#         vals[r:(r + dim^2 - 1)] .+= v[:]
#         r += dim^2
#     end
#     return
# end


# function integrate(form::Form{1}, space::Spaces.AbstractSpace)
#     (; integrands, measures) = form
#     mesh = domainmesh(first(measures))
#     F = floattype(mesh)
#     N = degrees_of_freedom!(mesh)
#     ℓ = sum(Base.Fix{2}(^, 2) ∘ length, tridofs(mesh))
#     rhs = FixedSizeArray{F, 1}(undef, ℓ)
#     fill!(zero(F), rhs)
#     for (fun, measure) in zip(integrands, measures)
#         mock = fun(space, space)
#         CT = coefftype(mock)
#         ord = order(mock)
#         rhs .+= buildvector(CT, ord::Order{B}, fun, measure, space)
#     end
#     return
# end

# # function buildvector(::ConstantCoeff,ord::Order{B},fun,measure,space)
# #     (;aux,mesh) = measure
# #     (;points,trilist,dofs) = mesh
# #     (;by_tri) = dofs

# # end
# # Integration with variable coefficients
# # function buildmatrix!(::Variable,ord::Order{B},ivec,jvec,vals,fun,measure,base) where B
# #     (;aux,mesh,sch) = measure
# #     (;points,trilist,dofs) = mesh
# #     (;by_tri) = dofs
# #     (;points,weights) = sch
# #     F = eltype(vals)
# #     aff = AffineToRef{F}()
# #     Ascale = collapser(ord,aff);
# #     r = 1
# #     for tri in keys(trilist)
# #         p,_ = psortednodes(tri,mesh)
# #         affine!(aff,points[tri])
# #         Ascale = collapser(ord,aff);
# #         doft = by_tri[tri]
# #         dim = length(doft)
# #         v = build_local_tensor(Variable(),ord::Order{B},fun,sch,base)
# #         C = aux[p].C
# #         @tensor v[i,j] := (loctensor[i,j,w,k,l]*Ascale[k,l])*sch.weights[w]
# #         v .= jac(aff)*C'*v*C

# #         ivec[r:r+dim^2-1] .+= repeat(doft,dim)
# #         jvec[r:r+dim^2-1] .+= repeat(doft,inner=dim)
# #         vals[r:r+dim^2-1] .+= v[:]
# #         r += dim^2
# #     end
# # end


# # function build_local_tensor(::Variable,::Val{2},::Order{B},fun,sch,base) where B
# #     N = length(base)
# #     Nw = length(sch.points)
# #     local_tensor = FixedSizeArray{Float64}(undef,N,N)
# #     for (j,φⱼ) in enumerate(base)
# #         for (i,φᵢ) in enumerate(base)
# #             for (k,point) in enumerate(sch.points)
# #                 local_tensor[i,j,k,..] .= fun(φᵢ,φⱼ)(point)
# #             end
# #         end
# #     end
# #     return local_tensor
# # end


# ##########################################################
# ##########################################################
# ##########################################################
# ##########################################################
# # function integrate(::Type{Spaces.Order{B}},fun::PolyField,m::Measure{M}) where {B,F,I,P,M<:HPMesh{F,I,P}}
# #     (;mesh,aux) = m
# #     degrees_of_freedom!(mesh)
# #     dims = len(B)+sum(B)
# #     if sum(B)>0
# #         T = MArray{Tuple{2,2},F,2}()
# #     end
# #     tensors = Dict{NTuple{3,P},Array{F,dims}}()
# #     (;points,trilist,DOFs) = mesh
# #     (;by_tri) = DOFs
# #     ℓ = sum(x->length(x)^2,by_tri)
# #     J = Vector{Int32}(undef,ℓ)
# #     K = Vector{Int32}(undef,ℓ)
# #     V = Vector{Float64}(undef,ℓ)
# #     affₜ = AffineToRef{F}()
# #     r = 1

# #     @inbounds for t in triangles(trilist)
# #         dofT = dof[t]
# #         p, pnod = psortednodes(t, mesh)
# #         update!(affₜ, view(points, pnod))
# #         if p in keys(tensors)
# #             local_tensor = tensors[p]
# #         else
# #             local_tensor = build_local_tensor(fun,p)
# #             tensors[p] = local_tensor
# #         end

# #         v = zeros(dim, dim)
# #         for j = 1:dim, i = 1:j
# #             v[i, j] = S[i, j, :] ⋅ z
# #         end
# #         v = dAₜ * C' * Symmetric(v) * C
# #         i = repeat(dofT, dim)
# #         j = repeat(dofT, inner = dim)
# #         J[r:r+dim^2-1] = i
# #         K[r:r+dim^2-1] = j
# #         V[r:r+dim^2-1] = v
# #         r += dim^2
# #     end
# #     sparse(J, K, V)
# # end


# # @inbounds function ref_integrate(fun, scheme::QScheme{N,T}) where {N,T}
# #     @assert N > 0

# #     ws = scheme.weights
# #     ps = scheme.points
# #     @assert length(ws) == length(ps)

# #     p1 = ps[1]
# #     R = typeof(ws[1] * fun(p1))

# #     s = zero(R)
# #     @simd for i in 1:length(ws)
# #         w = ws[i]
# #         p = ps[i]
# #         s += w * fun(p)
# #     end

# #     return 2s / factorial(N - 1)
# # end


# ### This function has a BUG! calc_vol is not defined. But for now we do not need this function.
# # @inbounds function integrate(fun, scheme::QScheme{N,T},
# #                              vertices::SMatrix{D,N,U}) where {N,T,D,U}
# #     @assert N > 0
# #     @assert N >= D + 1

# #     ws = scheme.weights
# #     ps = scheme.points
# #     @assert length(ws) == length(ps)

# #     p1 = ps[1]
# #     x1 = (vertices * p1)::SVector{D}
# #     X = typeof(x1)
# #     R = typeof(ws[1] * fun(x1))

# #     s = zero(R)
# #     @simd for i in 1:length(ws)
# #         w = ws[i]
# #         p = ps[i]
# #         x = vertices * p
# #         s += w * fun(x)
# #     end

# #     # If `U` is an integer type, then Linearalgebra.det converts to
# #     # floating-point values; we might want a different type
# #     vol = R(calc_vol(vertices)) / factorial(N - 1)
# #     return vol * s
# # end
# # function integrate(fun, scheme::QScheme, vertices::SMatrix)
# #     return error("Wrong dimension for vertices matrix")
# # end
# # @inbounds function integrate(fun, scheme::QScheme{N,T},
# #                              vertices::SVector{N,SVector{D,U}}) where {N,T,D,U}
# #     return integrate(fun, scheme,
# #                      SMatrix{D,N,U}(vertices[n][d] for d in 1:D, n in 1:N))
# # end
# # function integrate(fun, scheme, vertices::SVector{N,<:SVector}) where {N}
# #     return error("Wrong dimension for vertices array")
# # end

# # function integrate(kernel, scheme::QScheme{N},
# #                    vertices::AbstractVector) where {N}
# #     @assert length(vertices) == N
# #     @assert N > 0
# #     D = length(vertices[1])
# #     @assert N >= D + 1
# #     vertices′ = SVector{N}(map(SVector{D}, vertices))
# #     return integrate(kernel, scheme, vertices′)
# # end

# # function integrate(kernel, scheme::QScheme{N},
# #                    vertices::AbstractMatrix) where {N}
# #     @assert size(vertices, 1) == N
# #     @assert N > 0
# #     D = size(vertices, 2)
# #     @assert N >= D + 1
# #     vertices′ = SMatrix{N,D}(vertices)'
# #     return integrate(kernel, scheme, vertices′)
# # end


# # function integrate(::Type{Spaces.ConstantCoeff},::Type{Spaces.Order{B}},op,m::Measure{M}) where {B,F,I,P,M<:HPMesh{F,I,P}}
# #     (;mesh,aux) = m
# #     degrees_of_freedom!(mesh)
# #     dims = len(B)+sum(B)
# #     tensors = Dict{NTuple{3,P},Array{F,dims}}()
# #     (;trilist,DOFs) = mesh
# #     (;by_tri) = DOFs
# #     ℓ = sum(x->length(x)^2,by_tri)
# #     J = Vector{Int32}(undef,ℓ)
# #     K = Vector{Int32}(undef,ℓ)
# #     V = Vector{Float64}(undef,ℓ)
# #     Aₜ = MMatrix{2,2}(zeros(2,2))
# #     iAₜ = MMatrix{2,2}(zeros(2,2))

# #     r = 1
# #     @inbounds for t in triangles(trilist)
# #         dofT = dof[t]
# #         p, pnod = psortednodes(t, mesh)
# #         transform_matrix!(Aₜ, view(points, :, pnod))
# #         iAₜ .= inv(Aₜ)
# #         iAₜ .= iAₜ * iAₜ'
# #         z = vec(iAₜ)
# #         dAₜ = abs(det(Aₜ))
# #         v = zeros(dim, dim)
# #         for j = 1:dim, i = 1:j
# #             v[i, j] = S[i, j, :] ⋅ z
# #         end
# #         v = dAₜ * C' * Symmetric(v) * C
# #         i = repeat(dofT, dim)
# #         j = repeat(dofT, inner = dim)
# #         J[r:r+dim^2-1] = i
# #         K[r:r+dim^2-1] = j
# #         V[r:r+dim^2-1] = v
# #         r += dim^2
# #     end
# #     sparse(J, K, V)
# # end
