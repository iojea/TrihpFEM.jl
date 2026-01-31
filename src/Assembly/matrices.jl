"""
  collapser(::Order,aff::AffineToRef)

computes de matrix used for tensor collapsing depending on the order of the operator.   
"""
collapser(::Order{(0,0)},aff::AffineToRef) = I(2)
collapser(::Order{(1,0)},aff::AffineToRef) = aff.iA
collapser(::Order{(1,1)},aff::AffineToRef) = aff.iA'*aff.iA

"""
    _initvectors(I,F,ℓ)
creates two vectors of type `I` for indices, and a vector of type `F` for values, all of them with size `ℓ`.  
"""
function _initvectors(I,F,ℓ)
    ivec = FixedSizeArray{I,1}(undef,ℓ)
    fill!(ivec,zero(I))
    jvec = FixedSizeArray{I,1}(undef,ℓ)
    fill!(jvec,zero(I))
    vals = FixedSizeArray{F,1}(undef,ℓ)
    fill!(vals,zero(F))
    ivec,jvec,vals
end


"""
  integrate(form::Form{2},space::AbstractSpace)

Integrates the `Form` `form` using the basis of the space `space`.    
"""
function integrate(form::Form{2},space::Spaces.AbstractSpace)
    (;integrands,measures) = form
    mesh = domainmesh(first(measures))
    F = floattype(mesh)
    Itype = inttype(mesh)
    N = degrees_of_freedom!(mesh)
    ℓ = sum(Base.Fix{2}(^,2)∘length,dof(mesh).by_tri)
    ivec,jvec,vals = _initvectors(Itype,F,ℓ)
    for (fun,measure) in zip(integrands,measures)
        mock = fun(space,space)
        CT = coefftype(mock)
        ord = order(mock)
        buildmatrix!(CT,ord,ivec,jvec,vals,fun,measure,space)
    end
    sparse(ivec,jvec,vals,N,N)
end


# Integration with constant coefficients
function buildmatrix!(::ConstantCoeff,ord::Order{B},ivec,jvec,vals,fun,measure,space) where B
    (;aux,mesh) = measure
    (;points,trilist,dofs) = mesh
    (;by_tri) = dofs
    F = eltype(vals)
    Itype = eltype(ivec)
    tensordict = Dictionary{NTuple{3,Itype},FixedSizeArray{F,2+sum(B)}}()
    aff = AffineToRef{F}()
    r = 1
    for tri in keys(trilist)
        p,_ = psortednodes(tri,mesh)
        isin,token = gettoken(tensordict,p)
        if isin
            loctensor = gettokenvalue(tensordict,token)
        else
            loctensor = build_local_tensor(ConstantCoeff(),ord,fun,basis(space,p))
            set!(tensordict,p,loctensor)
        end
        affine!(aff,points[tri])
        Ascale = collapser(ord,aff);
        doft = by_tri[tri]
        dim = length(doft)
        C = aux[p].C
        @tensor v[i,j] := loctensor[i,j,k,l]*Ascale[k,l]
        v .= jac(aff)*C'*v*C
        ivec[r:r+dim^2-1] .+= repeat(doft,dim)
        jvec[r:r+dim^2-1] .+= repeat(doft,inner=dim)
        vals[r:r+dim^2-1] .+= v[:]
        r += dim^2
    end
end

function build_local_tensor(::ConstantCoeff,::Order{B},fun,base) where B
    inner_dim  = length(base)
    outer_dims = sum(B)
    dims = (inner_dim,inner_dim,(2 for _ in 1:outer_dims)...)
    local_tensor = FixedSizeArray{Float64}(undef,dims...)
    for (j,ψ) in enumerate(base)
        for (i,φ) in enumerate(base)
            local_tensor[i,j,..] .= ref_integrate(fun(φ,ψ))
        end
    end
    return local_tensor
end

function integrate(form::Form{1},space::Spaces.AbstractSpace)
    (;integrands,measures) = form
    mesh = domainmesh(first(measures))
    F = floattype(mesh)
    N = degrees_of_freedom!(mesh)
    ℓ = sum(Base.Fix{2}(^,2)∘length,tridofs(mesh))
    rhs = FixedSizeArray{F,1}(undef,ℓ)
    fill!(zero(F),rhs)
    for (fun,measure) in zip(integrands,measures)
        mock = fun(space,space)
        CT = coefftype(mock)
        ord = order(mock)
        rhs .+= buildvector(CT,ord::Order{B},fun,measure,space)
    end
end

# function buildvector(::ConstantCoeff,ord::Order{B},fun,measure,space)
#     (;aux,mesh) = measure
#     (;points,trilist,dofs) = mesh
#     (;by_tri) = dofs
    
# end
# Integration with variable coefficients
# function buildmatrix!(::Variable,ord::Order{B},ivec,jvec,vals,fun,measure,base) where B
#     (;aux,mesh,sch) = measure
#     (;points,trilist,dofs) = mesh
#     (;by_tri) = dofs
#     (;points,weights) = sch
#     F = eltype(vals)
#     aff = AffineToRef{F}()
#     Ascale = collapser(ord,aff);
#     r = 1
#     for tri in keys(trilist)
#         p,_ = psortednodes(tri,mesh)
#         affine!(aff,points[tri])
#         Ascale = collapser(ord,aff);
#         doft = by_tri[tri]
#         dim = length(doft)
#         v = build_local_tensor(Variable(),ord::Order{B},fun,sch,base)
#         C = aux[p].C
#         @tensor v[i,j] := (loctensor[i,j,w,k,l]*Ascale[k,l])*sch.weights[w]
#         v .= jac(aff)*C'*v*C

#         ivec[r:r+dim^2-1] .+= repeat(doft,dim)
#         jvec[r:r+dim^2-1] .+= repeat(doft,inner=dim)
#         vals[r:r+dim^2-1] .+= v[:]
#         r += dim^2
#     end
# end


# function build_local_tensor(::Variable,::Val{2},::Order{B},fun,sch,base) where B
#     N = length(base)
#     Nw = length(sch.points)
#     local_tensor = FixedSizeArray{Float64}(undef,N,N)
#     for (j,φⱼ) in enumerate(base)
#         for (i,φᵢ) in enumerate(base)
#             for (k,point) in enumerate(sch.points)
#                 local_tensor[i,j,k,..] .= fun(φᵢ,φⱼ)(point)
#             end
#         end
#     end
#     return local_tensor
# end


##########################################################
##########################################################
##########################################################
##########################################################
# function integrate(::Type{Spaces.Order{B}},fun::PolyField,m::Measure{M}) where {B,F,I,P,M<:HPMesh{F,I,P}}
#     (;mesh,aux) = m
#     degrees_of_freedom!(mesh)
#     dims = len(B)+sum(B)
#     if sum(B)>0
#         T = MArray{Tuple{2,2},F,2}()
#     end
#     tensors = Dict{NTuple{3,P},Array{F,dims}}()
#     (;points,trilist,DOFs) = mesh
#     (;by_tri) = DOFs
#     ℓ = sum(x->length(x)^2,by_tri)
#     J = Vector{Int32}(undef,ℓ)
#     K = Vector{Int32}(undef,ℓ)
#     V = Vector{Float64}(undef,ℓ)
#     affₜ = AffineToRef{F}()
#     r = 1
    
#     @inbounds for t in triangles(trilist)
#         dofT = dof[t]
#         p, pnod = psortednodes(t, mesh)
#         update!(affₜ, view(points, pnod))
#         if p in keys(tensors)
#             local_tensor = tensors[p]
#         else
#             local_tensor = build_local_tensor(fun,p)
#             tensors[p] = local_tensor
#         end
            
#         v = zeros(dim, dim)
#         for j = 1:dim, i = 1:j
#             v[i, j] = S[i, j, :] ⋅ z
#         end
#         v = dAₜ * C' * Symmetric(v) * C
#         i = repeat(dofT, dim)
#         j = repeat(dofT, inner = dim)
#         J[r:r+dim^2-1] = i
#         K[r:r+dim^2-1] = j
#         V[r:r+dim^2-1] = v
#         r += dim^2
#     end
#     sparse(J, K, V)
# end


# @inbounds function ref_integrate(fun, scheme::QScheme{N,T}) where {N,T}
#     @assert N > 0

#     ws = scheme.weights
#     ps = scheme.points
#     @assert length(ws) == length(ps)

#     p1 = ps[1]
#     R = typeof(ws[1] * fun(p1))

#     s = zero(R)
#     @simd for i in 1:length(ws)
#         w = ws[i]
#         p = ps[i]
#         s += w * fun(p)
#     end

#     return 2s / factorial(N - 1)
# end


### This function has a BUG! calc_vol is not defined. But for now we do not need this function.
# @inbounds function integrate(fun, scheme::QScheme{N,T},
#                              vertices::SMatrix{D,N,U}) where {N,T,D,U}
#     @assert N > 0
#     @assert N >= D + 1

#     ws = scheme.weights
#     ps = scheme.points
#     @assert length(ws) == length(ps)

#     p1 = ps[1]
#     x1 = (vertices * p1)::SVector{D}
#     X = typeof(x1)
#     R = typeof(ws[1] * fun(x1))

#     s = zero(R)
#     @simd for i in 1:length(ws)
#         w = ws[i]
#         p = ps[i]
#         x = vertices * p
#         s += w * fun(x)
#     end

#     # If `U` is an integer type, then Linearalgebra.det converts to
#     # floating-point values; we might want a different type
#     vol = R(calc_vol(vertices)) / factorial(N - 1)
#     return vol * s
# end
# function integrate(fun, scheme::QScheme, vertices::SMatrix)
#     return error("Wrong dimension for vertices matrix")
# end
# @inbounds function integrate(fun, scheme::QScheme{N,T},
#                              vertices::SVector{N,SVector{D,U}}) where {N,T,D,U}
#     return integrate(fun, scheme,
#                      SMatrix{D,N,U}(vertices[n][d] for d in 1:D, n in 1:N))
# end
# function integrate(fun, scheme, vertices::SVector{N,<:SVector}) where {N}
#     return error("Wrong dimension for vertices array")
# end

# function integrate(kernel, scheme::QScheme{N},
#                    vertices::AbstractVector) where {N}
#     @assert length(vertices) == N
#     @assert N > 0
#     D = length(vertices[1])
#     @assert N >= D + 1
#     vertices′ = SVector{N}(map(SVector{D}, vertices))
#     return integrate(kernel, scheme, vertices′)
# end

# function integrate(kernel, scheme::QScheme{N},
#                    vertices::AbstractMatrix) where {N}
#     @assert size(vertices, 1) == N
#     @assert N > 0
#     D = size(vertices, 2)
#     @assert N >= D + 1
#     vertices′ = SMatrix{N,D}(vertices)'
#     return integrate(kernel, scheme, vertices′)
# end


# function integrate(::Type{Spaces.ConstantCoeff},::Type{Spaces.Order{B}},op,m::Measure{M}) where {B,F,I,P,M<:HPMesh{F,I,P}}
#     (;mesh,aux) = m
#     degrees_of_freedom!(mesh)
#     dims = len(B)+sum(B)
#     tensors = Dict{NTuple{3,P},Array{F,dims}}()
#     (;trilist,DOFs) = mesh
#     (;by_tri) = DOFs
#     ℓ = sum(x->length(x)^2,by_tri)
#     J = Vector{Int32}(undef,ℓ)
#     K = Vector{Int32}(undef,ℓ)
#     V = Vector{Float64}(undef,ℓ)
#     Aₜ = MMatrix{2,2}(zeros(2,2))
#     iAₜ = MMatrix{2,2}(zeros(2,2))
    
#     r = 1
#     @inbounds for t in triangles(trilist)
#         dofT = dof[t]
#         p, pnod = psortednodes(t, mesh)
#         transform_matrix!(Aₜ, view(points, :, pnod))
#         iAₜ .= inv(Aₜ)
#         iAₜ .= iAₜ * iAₜ'
#         z = vec(iAₜ)
#         dAₜ = abs(det(Aₜ))
#         v = zeros(dim, dim)
#         for j = 1:dim, i = 1:j
#             v[i, j] = S[i, j, :] ⋅ z
#         end
#         v = dAₜ * C' * Symmetric(v) * C
#         i = repeat(dofT, dim)
#         j = repeat(dofT, inner = dim)
#         J[r:r+dim^2-1] = i
#         K[r:r+dim^2-1] = j
#         V[r:r+dim^2-1] = v
#         r += dim^2
#     end
#     sparse(J, K, V)
# end
    

