module Poly

using StaticArrays
using Polynomials
using LinearAlgebra
using TensorCast
using ..Meshes

include("fields.jl")
# include("opfields.jl")
include("differentiation.jl")
include("legendre.jl")

include("show.jl")

const ∇ = gradient
const Δ = laplacian
const ⊗ = _outer
export PolyField
export BiPoly
export PolyScalarField, PolyVectorField, PolyTensorField
export GeneralField
export PolySum
export indeterminate,indeterminates,degs
export LegendreIterator,StandardBasis
export derivative,gradient, divergence, laplacian,_outer
export dot
export AffineToRef
export affine!
export jac
export area
export EvalType,Eval,Compose,Pass
export evaluate
export ∇
export Δ
export _outer,⊗

end; #module
