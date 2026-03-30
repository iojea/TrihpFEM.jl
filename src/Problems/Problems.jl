module Problems

using StaticArrays
using Dictionaries
using FixedSizeArrays
using LinearAlgebra
using Makie
using ..Meshes
using ..PolyFields
using ..Integration
using ..Measures
using ..Forms
using ..Assembly

import CommonSolve: solve

include("solution.jl")
include("feproblem.jl")
include("plots.jl")

export FEProblem,FESolution,solve, plotsol

end; #module
