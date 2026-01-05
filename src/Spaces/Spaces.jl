module Spaces

using StaticArrays
using Dictionaries
using LinearAlgebra
using Polynomials
using ..Meshes
using ..PolyFields
include("hpspaces.jl")
# include("operation.jl")



export AbstractSpace
export StdScalarSpace, StdVectorSpace, StdTensorSpace
export ScalarSpace, VectorSpace, TensorSpace
export basis
export OperatorSpace
export gradient,divergence,laplacian
export dot
export Order
export order
export Constant, Variable
export coefftype


end; #module
