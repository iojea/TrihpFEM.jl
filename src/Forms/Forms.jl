module Forms

    using LinearAlgebra

    using ..Meshes
    using ..Tensors
    using ..Measures
    using ..PolyFields


    import ..DifferentialOperators: DiffOperator,Identity,Derivatex,Derivatey,Gradient,AdjointGradient, Divergence,Laplacian,DiffMatrix,AdjointDiffMatrix
    import ..DifferentialOperators: ∂x,∂y,gradient,∇,divergence,laplacian,Δ

    # include("terms.jl")
    include("shapefunction.jl")
    include("operation.jl")
    include("form.jl")

    export Term
    export Form
    export ShapeFunction
    export Integrand
    export Operation
    export Trial,Test
    export Order,CoeffType, ConstantCoeff, VariableCoeff
    export basis,order,coefftype,operator,chain_operator,dim
    export ∂x,∂y,gradient,∇,divergence,laplacian,Δ,∫

end; #module
