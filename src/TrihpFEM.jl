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
using Makie
using Markdown
using Pkg
using Polynomials
using Printf
using SparseArrays
using StaticArrays
using TensorCast
using TensorOperations
using Test
using Triangulate

include("Meshes/Meshes.jl")
include("Poly/Poly.jl")
# Write your package code here.
#

using ..Meshes: Edge,Triangle,HPMesh,BoundaryHPMesh,hpmesh,plothpmesh,dirichletboundary,neumannboundary,edges,setdirichlet!,setneumann!
export Edge,Triangle,HPMesh,BoundaryHPMesh,hpmesh,plothpmesh,dirichletboundary,neumannboundary,edges,setdirichlet!,setneumann!

using ..Poly: BiPoly,PolyTensorField,PolyVectorField,PolyMatrixField,AffineToRef,GeneralField,StandardBasis
export BiPoly,PolyTensorField,PolyVectorField,PolyMatrixField,AffineToRef,GeneralField,StandardBasis
end
