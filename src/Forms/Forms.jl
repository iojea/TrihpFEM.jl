module Forms

using ..Measures
using ..Poly

include("terms.jl")
include("form.jl")

export IntegrationTerm
export Form
export @term
export @form

end; #module