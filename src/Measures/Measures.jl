module Measures

    using LinearAlgebra
    using StaticArrays
    using FixedSizeArrays
    using CubatureRules
    using Dictionaries

    using ..Meshes
    using ..PolyFields
    using ..Integration

    import ..Meshes: elements


    #   NOTES
    ## Dirac measures should be implemented as subtypes of `Measure` but will probably require a specialized `Term`. At least for point Dirac measures.
    ## A basis type can be added to AbstractMeasure as a type parameter. This would probably be the best way to introduce other types of basis, if it is possible to make them compatible with the library approach.

    include("meas.jl")
    include("boundarymeas.jl")
    include("affine.jl")

    export AbstractMeasure, DomainMeasure, BoundaryMeasure, SkeletonMeasure
    export elementtype, elements, mesh, basis, dim
    export dof, ndof, tridof
    export floattype, inttype, degtype
    export affine
end; #module
