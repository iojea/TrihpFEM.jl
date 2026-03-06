module Forms

    using LinearAlgebra
    using MacroTools
    using ..Meshes
    using ..Measures
    using ..Poly
    using ..Spaces

    include("terms.jl")
    include("form.jl")

    export IntegrationTerm
    export Form
    export @term
    export @form
    export CoeffType, NoCoeff, ConstantCoeff, VariableCoeff

end; #module
