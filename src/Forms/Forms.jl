module Forms

    using LinearAlgebra
    using MacroTools
    using Tensors
    using ..Meshes
    using ..Measures
    using ..PolyFields


    import ..DifferentialOperators: DiffOperator,Identity,Derivatex,Derivatey,Gradient,Divergence,Laplacian
    import ..DifferentialOperators: ∂x,∂y,gradient,∇,divergence,laplacian,Δ

    # include("terms.jl")
    include("form.jl")

    export Term
    export Form
    export ShapeFunction
    export Integrand
    export Order,CoeffType, NoCoeff, ConstantCoeff, VariableCoeff
    export basis,order,coefftype,operator
    export ∂x,∂y,gradient,∇,divergence,laplacian,Δ,∫

end; #module
