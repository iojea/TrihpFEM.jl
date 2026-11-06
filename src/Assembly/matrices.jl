"""
    _initvectors(::AbstractMeasure,ℓ)

creates two vectors of type `I` for indices, and a vector of type `F` for values, all of them with size `ℓ`.  
"""
function _initvectors(::AbstractMeasure{F, I, P}, ℓ) where {F, I, P}
    ivec = FixedSizeArray{I, 1}(undef, ℓ)
    fill!(ivec, zero(I))
    jvec = FixedSizeArray{I, 1}(undef, ℓ)
    fill!(jvec, zero(I))
    vals = FixedSizeArray{F, 1}(undef, ℓ)
    fill!(vals, zero(F))
    return ivec, jvec, vals
end

"""

    _init_rhs(m::AbstractMeasure{F,I,P},ℓ)
creates a rhs vector of type  `I` and length `ℓ`
"""
function _init_rhs(::AbstractMeasure{F, I, P}, ℓ) where {F, I, P}
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
    meas = measure(first(terms))
    ℓ = degrees_of_freedom!(meas)
    N = sum(map(length, tridof(meas)) .^ 2)
    ivec, jvec, vals = _initvectors(meas, N)
    for t in terms
        add_to_matrix!(ivec, jvec, vals, t)
    end
    return sparse(ivec, jvec, vals, ℓ, ℓ)
end


"""
   add_to_matrix!(ivec,jvec,vals,t::Term{2}) 
"""
function add_to_matrix!(ivec, jvec, vals, t::Term{2})
    meas = measure(t)
    reftensor = RefTensor(t)
    collapser = Collapser(t)
    nmax = maximum(length.(tridof(meas)))
    v = zeros(floattype(meas), nmax, nmax)
    r = 1
    for el in keys(elements(meas))
        aff = affine(el, meas)
        degs, b₁ = keybasis(el, meas)
        b₂ = StandardBasis(degs[2:end])
        for (j, ψ) in enumerate(b₂)
            for (i, φ) in enumerate(b₁)
                v[i, j] = _collapse(reftensor, collapser, el, aff, φ, ψ)
            end
        end
        C = changetostd(t, degs)
        locdof = dof(el, meas)
        dim = length(locdof)
        @views v[1:dim, 1:dim] .= jac(aff) * C' * v[1:dim, 1:dim] * C
        ivec[r:(r + dim^2 - 1)] .= repeat(locdof, dim)
        jvec[r:(r + dim^2 - 1)] .= repeat(locdof, inner = dim)
        @views vals[r:(r + dim^2 - 1)] .+= v[1:dim, 1:dim][:]
        r += dim^2
    end
    return
end


function _collapse(
        rt::RefTensor{K, T, Nothing},
        c::Collapser,
        el::MeshElement,
        aff::AffineToRef,
        bs...
    ) where {K, T}
    return rt(el, bs...) ⋅ c(aff)
end

function _collapse(
        rt::RefTensor{K, T, R},
        c::Collapser,
        el::MeshElement,
        aff::AffineToRef,
        bs...
    ) where {K, T, R <: Quadrature}
    return rt(el, bs...) ⋅ c(aff, rt.sch)
end

_dot(a::Number, b::Number) = a * b
_dot(a, b) = dot(a, b)

"""

    assembly_rhs(form::Form{1})
assembles the rhs corresponding to the bilinear form  `form`.
"""
function assembly_rhs(term::Term{1})
    return assembly_rhs(Form{1}((term,)))
end
function assembly_rhs(form::Form{1})
    (; terms) = form
    meas = measure(first(terms))
    ℓ = Measures.ndof(meas)
    vals = _init_rhs(meas, ℓ)
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
    meas = measure(t)
    reftensor = RefTensor(t)
    collapser = Collapser(t)
    nmax = maximum(length.(tridof(meas)))
    v = zeros(floattype(meas), nmax)
    r = 1
    for el in keys(elements(meas))
        degs, b = keybasis(el, meas)
        aff = affine(el, meas)
        for (i, ψ) in enumerate(b)
            v[i] = _collapse(reftensor, collapser, el, aff, ψ)
        end
        C = changetostd(t, degs)

        locdof = dof(el, meas)
        dim = length(locdof)
        mul!(view(vals, locdof), C', view(v, 1:dim))
        vals[locdof] .*= jac(aff)
        # @views vals[locdof] .+= jac(aff) * C' * v[1:dim]
    end
    return
end
