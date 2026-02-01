module Forms

using LinearAlgebra
using ..Measures
using ..Poly
using ..Spaces

include("terms.jl")
include("form.jl")

export IntegrationTerm
export Form
export @term
export @form
export CoeffType,ConstantCoeff,VariableCoeff

end; #module