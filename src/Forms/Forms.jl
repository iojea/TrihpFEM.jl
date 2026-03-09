module Forms

    using LinearAlgebra
    using MacroTools
    using ..Meshes
    using ..Measures
    using ..PolyFields


    import ..DifferentialOperators: DiffOperator,Identity,Derivatex,Derivatey,Gradient,Divergence,Laplacian
    import ..DifferentialOperators: ∂x,∂y,gradient,∇,divergence,laplacian,Δ

    include("terms.jl")
    include("form.jl")

    export IntegrationTerm
    export Form
    export @term
    export @form
    export CoeffType, NoCoeff, ConstantCoeff, VariableCoeff

end; #module
