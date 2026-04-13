module Problems

using StaticArrays
using Dictionaries
using FixedSizeArrays
using LinearAlgebra
using Makie
using ..DifferentialOperators
using ..Meshes
using ..PolyFields
using ..Integration
using ..Measures
using ..Forms
using ..Assembly

import CommonSolve: solve

import ..DifferentialOperators: DiffOperator,Identity,Derivatex,Derivatey,Gradient,Divergence,Laplacian
import ..DifferentialOperators: ∂x,∂y,gradient,∇,divergence,laplacian,Δ

import ..Forms: ∫

include("solution.jl")
include("feproblem.jl")
# include("error.jl")
include("plots.jl")

export FEProblem,FESolution,solve, plotsol
export ∂x,∂y,gradient,∇,divergence,laplacian,Δ,∫
export error
end; #module
