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
# function assembly_matrix(term::Term{C, O, T, 2, M}) where {C, O, T, M}
#     return assembly_matrix(Form{2}((term,)))
# end
# function assembly_matrix(form::Form{2})
#     (; terms) = form
#     mesh = domainmesh(first(terms))
#     ℓ = degrees_of_freedom!(mesh)
#     N = sum(map(length, mesh.dofs.by_tri) .^ 2)
#     ivec, jvec, vals = _initvectors(mesh, N)
#     for t in terms
#         add_to_matrix!(ivec, jvec, vals, t)
#     end
#     return sparse(ivec, jvec, vals, ℓ, ℓ)
# end

#### Constant Coefficients
"""
   add_to_matrix!(ivec,jvec,vals,lt::LocalTensor) 
# """
# function add_to_matrix!(ivec,jvec,vals,lt::LocalTensor{K,T}) where {K,T<:Tensor}
#     for el in elements(measure)
#         B = collapser(el,lt)
#     end
# end

# function build_loc_tensor(lt::LocalTensor,element)
# end



function collapser(aff::AffineToRef,lt::LocalTensor{K,T}) where {K,T}
    sfs = get_shape_functions(lt.term.integrand.operation,lt.term.measure.mesh)
    ops = chain_operator.(sfs)
    mats = (o(aff) for o in ops)
    target_dims = ndims(T)
    _collapser(Val(target_dims),mats...)
end


function _collapser(::Val{L},mat) where L
    ndims(mat) == L || throw(ArgumentError("Collapser has `ndims` $(ndims(mat)), but $L are expected"))
    return mat
end

function _collapser(::Val{L},mat1,mat2) where L
    l1 = ndims(mat1)
    l2 = ndims(mat2)
    __collapser(Val(L),Val(l1),Val(l1+l2),mat1,mat2)
end

__collapser(::Val{L},::Val,::Val{L},m₁,m₂) where {L} = m₁⊗m₂
__collapser(::Val{L},::Val{L},::Val,m₁,m₂) where {L} = m₁*m₂ 
__collapser(::Val{L},::Val{L},::Val{L},m₁,m₂) where {L} = m₁*m₂
