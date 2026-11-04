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

"""

    _init_rhs(m::HPMesh{F,I,P},ℓ)
creates a rhs vector of type  `I` and length `ℓ`
"""
function _init_rhs(::HPMesh{F, I, P}, ℓ) where {F, I, P}
    vec = FixedSizeArray{F, 1}(undef, ℓ)
    fill!(vec, zero(F))
    return vec
end


"""

    assembly_matrix(form::Form{2}) 
assembles the matrix corresponding to the bilinear form  `form`.
"""
function assembly_matrix(term::Term{2})
    return assembly_matrix(Form{2}((term,)))
end
function assembly_matrix(form::Form{2})
    (; terms) = form
    mesh = getmesh(first(terms))
    ℓ = degrees_of_freedom!(mesh)
    N = sum(map(length, mesh.dofs.by_tri) .^ 2)
    ivec, jvec, vals = _initvectors(mesh, N)
    for t in terms
        add_to_matrix!(ivec, jvec, vals, t)
    end
    return sparse(ivec, jvec, vals, ℓ, ℓ)
end


PolyFields.AffineToRef(tri::Triangle,m::HPMesh) = AffineToRef(m.points[tri])

function PolyFields.AffineToRef(e::Edge,m::HPMesh)
    adj = adjacents(m.edgelist[e])
    adj[2]==0 || throw(ArgumentError("Is not possible to create an `AffineToRef` transformation from an interior edge."))
    tri = Triangle(e[1],e[2],adj[1])
    AffineToRef(tri,m)
end


"""
   add_to_matrix!(ivec,jvec,vals,t::Term{2}) 
"""
function add_to_matrix!(ivec,jvec,vals,t::Term{2})
    meas = getmeasure(t)
    mesh = getmesh(t)
    reftensor = RefTensor(t)
    collapser = Collapser(t)
    nmax = maximum(length.(mesh.dofs.by_tri))
    v = zeros(floattype(mesh), nmax, nmax)
    r = 1
    for (el,_) in elements(meas)
        aff = AffineToRef(el,getmesh(t))
        degs,b₁ = _extract(el,meas)
        b₂ = StandardBasis(degs[2:end])
        for (j,ψ) in enumerate(b₂)
            for (i,φ) in enumerate(b₁)
                v[i,j] = _collapse(reftensor,collapser,el,aff,φ,ψ)
            end
        end
        C = changetostd(t,degs)
        locdof = dof(el, mesh)
        dim = length(locdof)
        @views v[1:dim, 1:dim] .= jac(aff) * C' * v[1:dim,1:dim] * C
        ivec[r:(r + dim^2 - 1)] .= repeat(locdof, dim)
        jvec[r:(r + dim^2 - 1)] .= repeat(locdof, inner = dim)
        @views vals[r:(r + dim^2 - 1)] .+= v[1:dim, 1:dim][:]
        r += dim^2
    end
end


function _collapse(rt::RefTensor{K,T,Nothing},
                   c::Collapser,
                   el::MeshElement,
                   aff::AffineToRef,
                   bs...) where {K,T}
    rt(el,bs...)⋅c(aff)
end

function _collapse(rt::RefTensor{K,T,R},
                   c::Collapser,
                   el::MeshElement,
                   aff::AffineToRef,
                   bs...) where {K,T,R<:Quadrature}
    rt(el,bs...)⋅c(aff,rt.sch)
end

_dot(a::Number,b::Number) = a*b
_dot(a,b) = dot(a,b)

"""

    assembly_rhs(form::Form{1})
assembles the rhs corresponding to the bilinear form  `form`.
"""
function assembly_rhs(term::Term{1})
    return assembly_rhs(Form{1}((term,)))
end
function assembly_rhs(form::Form{1})
    (; terms) = form
    mesh = getmesh(first(terms))
    ℓ = ndof(mesh)
    vals = _init_rhs(mesh, ℓ)
    for t in terms
        add_to_rhs!(vals, t)
    end
    return vals
end

# ConstantCoeff version
"""

    add_to_rhs!(vals,t::Term)
integrates `t` adding the results to `vals`.
"""
function add_to_rhs!(vals, t::Term{1})
    mesh = getmesh(t)
    meas = getmeasure(t)
    reftensor = RefTensor(t)
    collapser = Collapser(t)
    nmax = maximum(length.(mesh.dofs.by_tri))
    v = zeros(floattype(mesh), nmax)
    r = 1
    for (el,_) in elements(meas)
        degs,b = _extract(el,meas)
        aff = AffineToRef(el,mesh)
        for (i,ψ) in enumerate(b)
                v[i] = _collapse(reftensor,collapser,el,aff,ψ)
        end
        C = changetostd(t,degs)
        
        locdof = dof(el, mesh)
        dim = length(locdof)
        mul!(view(vals,locdof),C',view(v,1:dim))
        vals[locdof] .*= jac(aff)
        # @views vals[locdof] .+= jac(aff) * C' * v[1:dim]
    end
end


