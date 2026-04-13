"""

    collapse(aff:AffineToRef,t::Tensor)
   
performs an appropriate contraction of the tensor `t` with the change of variables `aff`. `t` is obtained by integrating a form, including a constant factor, on the reference triangle. The contraction depends on the order and dimension of the tensor `t` which in turn is determined by the differential operators involved in the bilinear form, and the dimensions of the constant factor.

    collapse(factors,aff::AffineToRef,tensors)

performs a similar operation when the factor is not constant. In this case the factor is not included in the tensor. Morever: a sequence of factors and tensors is obtained by evaluating the factor and the tensor on the nodes of a quadrature rule. 
"""
function collapse(aff::AffineToRef, t::Tensor{4, 2})
    iA = inv(aff.A)
    return iA⊡t⊡iA'
end
collapse(aff::AffineToRef, t::Tensor{2, 2}) = inv(aff.A) ⊡ t
collapse(aff::AffineToRef, t::Tensor{1, 1}) = Tensor{1,1}((1,)) ⋅ t
function collapse(factors,aff::AffineToRef,tensors)
    iA = inv(aff.A)
    fac = (_collapser(f,aff) for f in factors)
    sum(_contraction(f,t) for (f,t) in zip(fac,tensors))
end

"""
   _contraction(f,t)
performs a single or double contraction depending on the arguments.  
"""
_contraction(f,t) = f⋅t
_contraction(f::Tensor{2,2},t::Tensor{2,2}) = f⊡t


"""

#     _collapser(::Order,aff::AffineToRef)
# computes de matrix used for tensor collapsing depending on the order of the operator.

    _collapser(aff:AffineToRef,t::Tensor)
performs an appropriate contraction of the tensor `t` with the change of variables `aff`.
"""
function _collapser(f::Tensor{2,2},aff::AffineToRef)
    iA = inv(aff.A)
    Tensor{2,2}((i,j) -> iA[i,:]⋅f'⋅iA[:,j])
end
function _collapser(f::Tensor{1,2},aff::AffineToRef)
    iA = inv(aff.A)
    f⋅iA
end
_collapser(f::Tensor{1,1},aff::AffineToRef) = f
# collapser(::Order{(0,)}, ::AffineToRef) = one(SymmetricTensor{1, 1})
# collapser(::Order{(0, 0)}, ::AffineToRef) = one(SymmetricTensor{1, 1})
# collapser(::Order{(1, 0)}, aff::AffineToRef) = inv(aff.A)
# collapser(::Order{(0, 1)}, aff::AffineToRef) = inv(aff.A)
# function collapser(::Order{(1, 1)}, aff::AffineToRef)
#     iA = inv(aff.A)
#     return otimesl(iA', iA)
# end
# collapser(f::T, aff::AffineToRef) where {T <: Integrand} = collapser(order(f), aff)

"""

    variablefactors(f,aff::AffineToRef,sch::Quadrature)

Receives a non-constant factor `f` (a `Function`), an affine transformation `aff` and a quadrature scheme, `sch` and computes a sequence of `Tensor`s produced by evaluating `f` on the nodes of `sch`. 
"""
function variablefactors(f,aff::AffineToRef,sch::Quadrature{D,F,P}) where {D,F,P}
    s = size(f(zeros(F,D)))
    dim = length(s)
    ord = dim > 0 ? 2 : 1
    if dim > 0
        return (Tensor{ord,dim}(f(aff(x))) for x in sch.points)
    else
        return (Tensor{ord,1}((f(aff(x)),)) for x in sch.points)
    end
end
    
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
    return vec
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
function outer(operand::Tensor{2, 2}, factor::Number)
    f = SymmetricTensor{2,2}((factor,zero(factor),factor))
    otimesl(operand,f)
end
outer(operand::Tensor{1, 2}, factor::Tensor{1, 2}) = factor ⊗ operand
outer(operand::Tensor{1, 1}, factor::Number) = factor * operand


"""

    ref_tensors(inte::Integrand{ConstantCoeff,T,Order{B},N},degs) where {T,B,N}
Computes the local tensors produced by `inte` over a reference element with degrees given by `degs`. These tensors are later contracted with the corresponding `AffineToRef` change of variables.

    ref_tensors(inte::Integrand{VariableCoeff,T,Order{B},N},degs,sch) where {T,B,N}
does the same for a `VariableCoeff` integrand. In this case instead of a `Tensor` a sequence of `Tensor`s is given where each `Tensor` corresponds to the evaluation of the integrand in a node of the quadrature scheme `sch`.
    """
function ref_tensors(inte::Integrand{ConstantCoeff, T, Order{B}, 2}, degs) where {T, B}
    (; factor, funs) = inte
    b₁ = basis(funs[1], degs); b₂ = basis(funs[2], degs)
    o₁, o₂ = operator.(funs)
    f = tensorize(factor)
    return collect_as(FixedSizeArrayDefault, (outer(_tensor(o₁(φ), o₂(ψ), Val(sum(B))),f) for φ in b₁, ψ in b₂))
end

function ref_tensors(inte::Integrand{ConstantCoeff,T,Order{B},1},degs) where {T,B}
    (;factor,funs) = inte
    b = basis(funs[1],degs)
    o = operator(funs[1])
    f = tensorize(factor)
    return collect_as(FixedSizeArrayDefault,(outer(_tensor(o(φ),Val(sum(B))),f) for φ in b))
end
function ref_tensors(inte::Integrand{VariableCoeff,T,Order{B},2},degs,sch) where {T,B}
    (;funs) = inte
    (;points,weights)=sch
    b₁ = basis(funs[1],degs); b₂ = basis(funs[2],degs)
    o₁,o₂ = operator.(funs)
    return collect_as(FixedSizeArrayDefault,((w*_quadtensor(o₁(φ)(x),o₂(ψ)(x),Val(sum(B))) for (x,w) in zip(points,weights)) for φ in b₁, ψ in b₂))
end
function ref_tensors(inte::Integrand{VariableCoeff,T,Order{B},1},degs,sch) where {T,B}
    (;funs) = inte
    (;points,weights)=sch
    b = basis(funs[1],degs);
    o = operator(funs[1])
    return collect_as(FixedSizeArrayDefault,((w*_quadtensor(o(φ)(x),Val(sum(B))) for (x,w) in zip(points,weights)) for φ in b))
end

"""

    _tensor(f,g,::Val{N}) where N
    _tensor(f,::Val{N})
an auxiliary function that builds a tensor by integrating `f*g` in the reference triangle. `N` is given to indicate the order of the `Tensor`. `N==2` corresponds to two derivatives (∇⋅∇), `N==1` to one derivative 
"""
_tensor(f, g, ::Val{2}) = Tensor{2, 2}((i, j) -> ref_integrate(f[i] * g[j]))
_tensor(f, g, ::Val{1}) = Tensor{1, 2}(i -> ref_integrate(f[i] * g))
_tensor(f, g, ::Val{0}) = Tensor{1, 1}((ref_integrate(f * g),))
_tensor(f,::Val{0}) = Tensor{1,1}((ref_integrate(f),))



"""

    _quadtensor(f,g,::Val{N}) where N
    _quadtensor(f,::Val{N})
similar to `_tensor` but for evaluating the functions on quadrature nodes instead of integrating them directly. This is the version used for non-constant coefficients.
"""
_quadtensor(f,g,::Val{2}) = Tensor{2,2}((i,j) -> f[i]*g[j])
_quadtensor(f,g,::Val{1}) = Tensor{1,2}(i->f[i]*g)
_quadtensor(f,g,::Val{0}) = Tensor{1,1}((f*g,))
_quadtensor(f,::Val{0}) = Tensor{1,1}((f,))



"""

    assembly_matrix(form::Form{2})
assembles the matrix corresponding to the bilinear form  `form`.
"""
function assembly_matrix(term::Term{C, O, T, 2, M}) where {C, O, T, M}
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
integrates the `Term` `t` adding the results to `ivec`,`jvec` and `vals`, for later building a sparse matrix.
"""
# ConstantCoeff version
function add_to_matrix!(ivec, jvec, vals, t::Term{ConstantCoeff, O, T, 2, M}) where {O, T, M}
    (; integrand, measure) = t
    (; mesh,aux) = measure
    r = 1
    tensordict = Dictionary()
    nmax = maximum(length.(mesh.dofs.by_tri))
    v = zeros(floattype(mesh),nmax,nmax)
    for el in Measures.elements(measure)
        degs, _ = psortednodes(el, mesh)
        isin, token = gettoken(tensordict, degs)
        if isin
            loctensor = gettokenvalue(tensordict, token)
        else
            loctensor = ref_tensors(integrand, degs)
            set!(tensordict, degs, loctensor)
        end
        aff = AffineToRef(mesh.points[el])
        locdof = dof(el,mesh)
        dim = length(locdof)
        C = aux[degs].C
        for i in 1:dim, j in 1:dim
            v[i,j] = collapse(aff, loctensor[i,j])
        end
       
        @views v[1:dim,1:dim] .= jac(aff) * C' * v[1:dim,1:dim] * C
        ivec[r:(r + dim^2 - 1)] .= repeat(locdof, dim)
        jvec[r:(r + dim^2 - 1)] .= repeat(locdof, inner = dim)
        @views vals[r:(r + dim^2 - 1)] .+= v[1:dim,1:dim][:]
        r += dim^2
    end
end

# VariableCoeff version
function add_to_matrix!(ivec, jvec, vals, t::Term{VariableCoeff, O, T, 2, M}) where {O, T, M}
    (; integrand, measure) = t
    (; mesh,aux,sch) = measure
    r = 1
    tensordict = Dictionary()
    nmax = maximum(length.(mesh.dofs.by_tri))
    v = zeros(floattype(mesh),nmax,nmax)
    for el in elements(measure)
        degs, _ = psortednodes(el, mesh)
        isin, token = gettoken(tensordict, degs)
        if isin
            loctensor = gettokenvalue(tensordict, token)
        else
            loctensor = ref_tensors(integrand, degs,sch)
            set!(tensordict, degs, loctensor)
        end
        aff = AffineToRef(mesh.points[el])
        locdof = dof(el,mesh)
        dim = length(locdof)
        C = aux[degs].C
        factors = variablefactors(integrand.factor,aff,sch)
        for i in 1:dim, j in 1:dim
            v[i,j] = collapse(aff, loctensor[i,j])
││      end
        # v = collect_as(FixedSizeArrayDefault, (collapse(factors,aff, loctensor[i,j]) for i in 1:dim, j in 1:dim))
        copyto!(v[1:dim,1:dim],jac(aff) * C' * v[1:dim,1:dim] * C)
        # v .= jac(aff) * C' * v * C
        ivec[r:(r + dim^2 - 1)] .= repeat(locdof, dim)
        jvec[r:(r + dim^2 - 1)] .= repeat(locdof, inner = dim)
        vals[r:(r + dim^2 - 1)] .+= v[1:dim,1:dim]
        r += dim^2
    end
end

"""

    assembly_rhs(form::Form{1})
assembles the rhs corresponding to the bilinear form  `form`.
"""
function assembly_rhs(term::Term{C, O, T, 1, M}) where {C, O, T, M}
    return assembly_rhs(Form{1}((term,)))
end
function assembly_rhs(form::Form{1})
    (; terms) = form
    mesh = domainmesh(first(terms))
    ℓ    = ndof(mesh)
    vals = _init_rhs(mesh,ℓ)
    for t in terms
        add_to_rhs!(vals, t)
    end
    return vals
end

"""

    add_to_rhs!(vals,t::Term)
integrates `t` adding the results to `vals`.
"""
# ConstantCoeff version
function add_to_rhs!(vals, t::Term{ConstantCoeff, O, T, 1, M}) where {O, T, M}
    (; integrand, measure) = t
    (; mesh,aux) = measure
    r = 1
    tensordict = Dictionary()
    nmax = maximum(length.(mesh.dofs.by_tri))
    v = zeros(floattype(mesh),nmax)
    for el in Measures.elements(measure)
        degs, _ = psortednodes(el, mesh)
        isin, token = gettoken(tensordict, degs)
        if isin
            loctensor = gettokenvalue(tensordict, token)
        else
            loctensor = ref_tensors(integrand, degs)
            set!(tensordict, degs, loctensor)
        end
        aff = AffineToRef(mesh.points[el])
        locdof = dof(el,mesh)
        dim = length(locdof)
        C = aux[degs].C
        v[1:dim] .= (collapse(aff, loctensor[i]) for i in 1:dim)
        @views vals[locdof] .+= jac(aff) * C' * v[1:dim]
        r += dim
    end
end

# VariableCoeff version
function add_to_rhs!(vals, t::Term{VariableCoeff, O, T, 1, M}) where {O, T, M}
    (; integrand, measure) = t
    (; mesh,aux,sch) = measure
    r = 1
    tensordict = Dictionary()
    for el in Measures.elements(measure)
        degs, _ = psortednodes(el, mesh)
        isin, token = gettoken(tensordict, degs)
        if isin
            loctensor = gettokenvalue(tensordict, token)
        else
            loctensor = ref_tensors(integrand, degs, sch)
            set!(tensordict, degs, loctensor)
        end
        aff = AffineToRef(mesh.points[el])
        locdof = dof(el,mesh)
        dim = length(locdof)
        C = aux[degs].C
        factors = variablefactors(integrand.factor,aff,sch)
        v = collect_as(FixedSizeArrayDefault, (collapse(factors,aff, loctensor[i]) for i in 1:dim))
        v .= jac(aff) * C' * v
        vals[locdof] .+= v[:]
        r += dim
    end
end


