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

using StyledStrings
using CommonSolve
using Dictionaries
using DocStringExtensions
using ExactPredicates
using FixedSizeArrays
using LinearAlgebra
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
include("Tensors/Tensors.jl")
include("PolyFields/PolyFields.jl")
# include("Spaces/Spaces.jl")
include("Integration/Integration.jl")
include("Measures/Measures.jl")
include("Forms/Forms.jl")
include("Assembly/Assembly.jl")
include("Problems/Problems.jl")

using ..Meshes: Point2D, Edge, Triangle, HPMesh, BoundaryHPMesh, hpmesh, plothpmesh, plothpmesh2, dirichletboundary, neumannboundary, edges, mark!, refine!, setdirichlet!, setneumann!, setdegrees!, degplot, circmesh, circmesh_graded_center, rectmesh, squaremesh, details
export Point2D,Edge, Triangle, HPMesh, BoundaryHPMesh, hpmesh, plothpmesh, plothpmesh2, dirichletboundary, neumannboundary, edges, mark!, refine!, setdirichlet!, setneumann!, setdegrees!, degplot, circmesh, circmesh_graded_center, rectmesh, squaremesh, details

using ..DifferentialOperators: DiffOperator, Identity, Derivatex, Derivatey, Gradient, AdjointGradient, Divergence, Laplacian, gradient, divergence, laplacian, ∇, Δ,∂x,∂y
export DiffOperator, Identity, Derivatex, Derivatey, Gradient, AdjointGradient, Divergence, Laplacian, gradient, divergence, laplacian, ∇, Δ, ∂x,∂y

using ..Tensors: Tensor, otimes,⊗
export Tensor,otimes,⊗

using ..PolyFields: ProductPoly, AffineToRef, GeneralField, LegendreIterator, StandardBasis, degs
export ProductPoly, AffineToRef, GeneralField, LegendreIterator, StandardBasis,degs

# using ..Spaces: StdScalarSpace, StdVectorSpace, OperatorSpace, order, EvalType, Order, Eval, Pass, combine, basis, ∇, Δ
# export StdScalarSpace, StdVectorSpace, OperatorSpace, order, EvalType, Order, Eval, Pass, combine, basis


using ..Integration: Quadrature, gmquadrature, ref_integrate, quadrature
export Quadrature, gmquadrature, ref_integrate, quadrature

using ..Forms: basis, Form, Term, Operation, ShapeFunction, CoeffType, ConstantCoeff, VariableCoeff, ∫,Trial,Test
export basis, Term, Form, Operation, ShapeFunction, CoeffType,  ConstantCoeff, VariableCoeff,Trial,Test, ∫

using ..Measures: Measure
export Measure

using ..Assembly: LocalTensor, _eval_operation, get_shape_functions
export LocalTensor, _eval_operation, get_shape_functions

# using ..Problems: FEProblem, FESolution, solve, plotsol, error, estimate_order
# export FEProblem, FESolution, solve, plotsol, error, estimate_order
end
