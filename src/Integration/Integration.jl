module Integration

using StaticArrays
using Polynomials
using FixedSizeArrays

using ..Poly

include("gmquads.jl")
include("integrate.jl")

export Quadrature
export gmquadrature
export ref_integrate


end; #module

