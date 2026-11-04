#####################################################################################
# TERMs and FORMs
#####################################################################################
abstract type AbstractForm end

"""
   Term 
"""
struct Term{N, C <: CoeffType, T <: Tuple, M <: AbstractMeasure} <: AbstractForm
    integrand::Integrand{C, T, N}
    measure::M
    Term(i::Integrand{C, T, N}, m::M) where {C, T, N, M} = new{N, C, T, M}(i, m)
end
coefftype(::Term{N, C}) where {N, C} = C()
order(::Term{C, T, N}) where {C, T, N} = T


struct Form{N} <: AbstractForm
    terms::Tuple{Vararg{Term{N}}}
end

function changetostd(t::Term, degs::NTuple{4,P}) where P<:Integer
    _changetostd(t,Val(Int(degs[1])),degs[2:end])
end

_changetostd(t::Term,::Val{0},degs) = t.measure.mesh.cmatrix(degs)
function _changetostd(t::Term,v::Val{N},degs) where N
    C = t.measure.mesh.cmatrix(degs)
    inds = _indices(Val(Int(N)),degs)
    view(C,inds,inds)
end

_indices(::Val{0},degs::NTuple{3,P}) where P<:Integer = 1:length(StandardBasis(degs))
_indices(::Val{1},degs::NTuple{3,P}) where P<:Integer = 1:(degs[1]+1)
_indices(::Val{2},degs::NTuple{3,P}) where P<:Integer = (degs[1]+1):(degs[1]+degs[2]+1)
_indices(::Val{3},degs::NTuple{3,P}) where P<:Integer = mod1.((degs[1]+degs[2]+1):(sum(degs)+1),sum(degs))



getmesh(t::Term) = t.measure.mesh
getmeasure(t::Term) = t.measure

function Base.:*(inte::Integrand, meas::AbstractMeasure)
    return Term(inte, meas)
end
function Base.:-(term::Term)
    (; integrand, measure) = term
    (; operation) = integrand
    return Term(Integrand((-1) * operation), measure)
end

Base.:+(t₁::Term{N}, t₂::Term{N}) where {N} = Form((t₁, t₂))
Base.:-(t₁::Term{N}, t₂::Term{N}) where {N} = Form((t₁, -t₂))
Base.:+(form::Form{N}, term::Term{N}) where {N} = Form((form.terms..., term))
Base.:+(term::Term{N}, form::Form{N}) where {N} = Form((term, form.terms...))
Base.:-(form::Form{N}, term::Term{N}) where {N} = Form((form.terms..., -term))
Base.:-(term::Term{N}, form::Form{N}) where {N} = Form((term, (-t for t in form.terms)...))
Base.:+(f₁::Form{N}, f₂::Form{N}) where {N} = Form(((t for t in f₁.terms)..., (t for t in f₂.terms)...))
Base.:-(f₁::Form{N}, f₂::Form{N}) where {N} = Form(((t for t in f₁.terms)..., (-t for t in f₂.terms)...))
