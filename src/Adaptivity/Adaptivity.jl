module Adaptivity

using ExactPredicates
using LinearAlgebra
using ..Meshes
using ..Measures
using ..Forms
using ..Problems

export Normal
export SolOperation

include("normal.jl")
include("solutionoperation.jl")

end; #module