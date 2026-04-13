module TrihpFEM

# In case you want to know, why the last line of the docstring below looks like it is:
# It will show the package (local) path when help on the package is invoked like     help?> TrihpFEM
# but it will interpolate to an empty string on CI server,
# preventing appearing the server local path in the documentation built there.

"""
    Package TrihpFEM v\$(pkgversion(TrihpFEM))

TrihpFEM implements an hp-adaptive Finite Element Method based on triangular meshes (2D).

\$(isnothing(get(ENV, "CI", nothing)) ? ("\n" * "Package local path: " * pathof(TrihpFEM)) : "") 
"""


using CommonSolve
using Dictionaries
using DocStringExtensions
using EllipsisNotation
using ExactPredicates
using FixedSizeArrays
using LinearAlgebra
using MacroTools
using Makie
using Markdown
using Pkg
using Polynomials
using Printf
using SparseArrays
using StaticArrays
using Test
using Triangulate

include("Meshes/Meshes.jl")
include("DifferentialOperators/DifferentialOperators.jl")
include("PolyFields/PolyFields.jl")
# include("Spaces/Spaces.jl")
include("Integration/Integration.jl")
include("Measures/Measures.jl")
include("Forms/Forms.jl")
include("Assembly/Assembly.jl")
include("Problems/Problems.jl")

using ..Meshes: Edge, Triangle, HPMesh, BoundaryHPMesh, hpmesh, plothpmesh, dirichletboundary, neumannboundary, edges, setdirichlet!, setneumann!, setdegrees!, degplot, circmesh, rectmesh, squaremesh
export Edge, Triangle, HPMesh, BoundaryHPMesh, hpmesh, plothpmesh, dirichletboundary, neumannboundary, edges, setdirichlet!, setneumann!, setdegrees!, degplot, circmesh, rectmesh, squaremesh

using ..DifferentialOperators: DiffOperator, Identity, Derivatex, Derivatey, Gradient, Divergence, Laplacian, gradient, divergence, laplacian, ∇, Δ
export DiffOperator, Identity, Derivatex, Derivatey, Gradient, Divergence, Laplacian, gradient, divergence, laplacian, ∇, Δ

using ..PolyFields: BiPoly, PolyTensorField, PolyVectorField, PolyMatrixField, AffineToRef, GeneralField, StandardBasis
export BiPoly, PolyTensorField, PolyVectorField, PolyMatrixField, AffineToRef, GeneralField, StandardBasis

# using ..Spaces: StdScalarSpace, StdVectorSpace, OperatorSpace, order, EvalType, Order, Eval, Pass, combine, basis, ∇, Δ
# export StdScalarSpace, StdVectorSpace, OperatorSpace, order, EvalType, Order, Eval, Pass, combine, basis


using ..Integration: Quadrature, gmquadrature, ref_integrate, quadrature
export Quadrature, gmquadrature, ref_integrate, quadrature

using ..Forms: basis, Form, Term, ShapeFunction, CoeffType, NoCoeff, ConstantCoeff, VariableCoeff, ∫
export basis, Term, Form, ShapeFunction, CoeffType, NoCoeff, ConstantCoeff, VariableCoeff, ∫

using ..Measures: Measure
export Measure

using ..Assembly: ref_integrate, ref_tensors, assembly_matrix
export ref_integrate, ref_tensors, assembly_matrix

using ..Problems: FEProblem, FESolution, solve, plotsol, error
export FEProblem, FESolution, solve, plotsol, error
end
