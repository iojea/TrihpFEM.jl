module Meshes

using Dictionaries
using ExactPredicates
using LinearAlgebra
using Makie
using Markdown
using Printf
using Triangulate
using StaticArrays
using DocStringExtensions


export Edge,EdgeAttributes,Triangle,TriangleAttributes,HPMesh,DOF
export triangle,hpmesh,data
export tag,degree,dof,longestedge,tagged_dof
export ismarked,isgreen,isblue,isred,istagged,isinterior,isboundary
export inttype,floattype,degtype
export degrees_of_freedom!,isempty,empty!
export settag!,setdegree!,setboundary!,setdirichlet!,setneumann!
export degrees,psortperm,edges,triangles,psortednodes,psortededges
export plothpmesh
export BoundaryHPMesh,dirichletboundary,neumannboundary
export mark!,refine!,p_conformity!,check_p_conformity
export circmesh,circmesh_graded_center,rectmesh,squaremesh
export plothpmesh,degplot


BOUNDARY_DICT = Dict(:dirichlet=>1,:neumann=>2)

include("settuple.jl")
include("edge.jl")
include("triangle.jl")
include("mesh.jl")
include("boundary_mesh.jl")
include("refine.jl")
include("show.jl")
include("plots.jl")
include("some_meshes.jl")



const HASH_SEED = UInt === UInt64 ? 0x793bac59abf9a1da : 0xdea7f1da

end;