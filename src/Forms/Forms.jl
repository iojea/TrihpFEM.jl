module Forms

using LinearAlgebra
using ..Measures
using ..Poly

include("terms.jl")
include("form.jl")

export IntegrationTerm
export Form
export @form

export Constant, Variable
end; #module