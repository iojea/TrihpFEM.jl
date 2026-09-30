#####################################################################################
# TERMs and FORMs
#####################################################################################
abstract type AbstractForm end

"""
   Term 
"""
struct Term{N,C <: CoeffType, T<:Tuple,M<:Measure} <: AbstractForm
    integrand::Integrand{C, T,N}
    measure::M
    Term(i::Integrand{C,T,N},m::M) where {C,T,N,M} = new{N,C,T,M}(i,m)
end
coefftype(::Term{N,C}) where {N,C} = C()
order(::Term{C,T,N}) where {C,T,N} = T


struct Form{N} <: AbstractForm
    terms::Tuple{Vararg{Term{N}}}
end

Meshes.domainmesh(t::Term) = domainmesh(t.measure)
function Base.:*(inte::Integrand, meas::Measure)
    return Term(inte, meas)
end
function Base.:-(term::Term)
    (; integrand, measure) = term
    (; operation) = integrand
    Term(Integrand((-1)*operation),measure)
end

Base.:+(t₁::Term{N},t₂::Term{N}) where N = Form((t₁,t₂))
Base.:-(t₁::Term{N}, t₂::Term{N}) where N = Form((t₁, -t₂))
Base.:+(form::Form{N},term::Term{N}) where N = Form((form.terms...,term))
Base.:+(term::Term{N},form::Form{N}) where N = Form((term,form.terms...))
Base.:-(form::Form{N},term::Term{N}) where N = Form((form.terms...,-term))
Base.:-(term::Term{N},form::Form{N}) where N = Form((term,(-t for t in form.terms)...))
Base.:+(f₁::Form{N},f₂::Form{N}) where N = Form(((t for t in f₁.terms)...,(t for t in f₂.terms)...))
Base.:-(f₁::Form{N},f₂::Form{N}) where N = Form(((t for t in f₁.terms)...,(-t for t in f₂.terms)...))
