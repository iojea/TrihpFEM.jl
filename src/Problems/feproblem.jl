"""

    FEProblem(a,b,space,g)

Defines a Finite Element Problem (no refinement) of the form `a(u,v)=b(v)` over the space `space` with Dirichlet boundary data `g`.

# Arguments
+ `a`: a bivariate form, typically defined using `@form`.
+ `b`: a univariate form, typically defined using `@form`.
+ `space`: an appropriate space. Currently, only `StdScalarSpace()` is supported.
+ `g` a function defining the Dirichlet boundary condition.
"""

struct FEProblem{S<:AbstractSpace}
    a::Form
    b::Form
    space::S
    g
end

FEProblem(a,b,space,g) = FEProblem{typeof(space)}(a,b,space,g)


# function FEProblem(a,b,space,g)
#     matrix = integrate(a,space)
#     rhs = integrate(b,space)
# end

