module Integration

using StaticArrays
using Polynomials
using FixedSizeArrays
using Collects

using ..PolyFields

include("repvec.jl")
include("gmquads.jl")
include("precomputed.jl")
include("integrate.jl")

export Quadrature
export gmquadrature
export ref_integrate
export quadrature

end; #module

