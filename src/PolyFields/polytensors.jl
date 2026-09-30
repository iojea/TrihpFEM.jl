###########################################################################
##################              Tensor Fields              ################
###########################################################################

const PolyTensorField{N,D,T} = Tensor{N,D,T} where {N,D,T<:PolyScalarField{D}}
const PolyVectorField{D,T} = PolyTensorField{1,D,T} where {D,T<:PolyScalarField{D}}
const PolyMatrixField{D,T} = PolyTensorField{D,D,T} where {D,T<:PolyScalarField{D}}

(t::PolyTensorField{N,D,T})(x) where {N,D,T<:PolyScalarField} = Tensor{N,D}((p(x) for p in t)...)

# in order to be able to transpose tensor fields (done recursively on the dimensions)
Base.adjoint(p::PolyScalarField) = p

otimes(t::Tensor,q::PolyScalarField) = q.*t
otimes(q::PolyScalarField,t::Tensor) = q.*t
otimes(p::PolyScalarField,q::PolyScalarField) = p*q

