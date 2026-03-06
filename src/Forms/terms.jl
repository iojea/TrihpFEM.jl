"""
    abstract type CoeffType end
An abstract type for defining a treat that allow `IntegrationTerm`s to be dispatched to different methods for integration. Exact integration for constant coefficient terms and quadrature rules for variable coefficient terms.
"""
abstract type CoeffType end

struct ConstantCoeff <: CoeffType end
struct VariableCoeff <: CoeffType end
struct NoCoeff <: CoeffType end

coefftype(_) = ConstantCoeff()
coefftype(::Function) = VariableCoeff()
coefftype(::DiffOperator) = ConstantCoeff()
coefftype(::Nothing) = NoCoeff()
"""
   IntegrationTerm(pf,f,meas)
Builds an integration term of the form `term(args...) = ∫(f*pf(args...))*meas`, where `pf` is a function to be evaluated on the local basis, `f` is a factor that can be constant or variable and `meas` is a measure for integration. 
"""
struct IntegrationTerm{C <: CoeffType, N}
    polyfun
    factor
    measure
end

function IntegrationTerm(pf, f, meas)
    T = typeof(coefftype(f))
    N = first(methods(pf)).nargs - 1
    return IntegrationTerm{T, N}(pf, f, meas)
end

numargs(::IntegrationTerm{C, N}) where {C, N} = N

coefftype(::IntegrationTerm{C, N}) where {C, N} = C()

Meshes.domainmesh(t::IntegrationTerm) = domainmesh(t.measure)

function Spaces.order(term::IntegrationTerm{C, N}, space::Spaces.AbstractSpace) where {C, N}
    mock = term.polyfun((space for _ in 1:N)...)
    return order(mock)
end

"""
  A list of reserved symbols to discard when looking for the `factor` in an `IntegrationTerm`  
"""
RESERVED_SYMBS = (:*, :⋅, :-, :∂x, :∂y, :gradient, :divergence, :laplacian, :∇, :Δ, :ε)

"""
    strip_block(ex::Expr)
removes the `:block` part of an expression, returning the meaningful part.
"""
strip_block(ex::Expr) = ex.head == :block ? ex.args[2] : ex

"""
   @term ex

Builds an `IntegrationTerm` from expression `ex`. `ex` should be a one line function definition defined by an integral.

# Examples

```julia
julia> Ω = circmesh(0.1);
julia> dΩ = Measure(Ω,5);
julia> A = rand(2,2);
julia> @term a(u,v) = ∫((A*∇(u))⋅∇(v))*dΩ
```
creates an `IntegrationTerm` named `a`, with function `(u,v)->∇(u)*∇(v)`, constant `A` and measure `dΩ`
"""
macro term(ex)
    return if @capture(ex, name_(args__) = term_)
        term = strip_block(term)
        msg = "Terms does not include summations. Maybe you want to create a `Form`. See `@form`."
        occursin("+", string(term)) && throw(ArgumentError(msg))
        occursin("-", string(term)) && throw(ArgumentError(msg))
        if @capture(term, ∫(integrand_) * meas_)
            factor = get_factor(integrand, args)
            intbody = cleanfactor(integrand, factor, args)
            fun = build_polyfun(args, intbody)
            Expr(:(=), esc(name), Expr(:call, :IntegrationTerm, esc_non_params(fun, args), (esc(factor)), esc(meas)))
        else
            throw(ArgumentError("Malformed expression. The right hand side must be an integral of the form `∫(fun)*measure`"))
        end
    else
        throw(ArgumentError("Malformed expression. A term of the form `name(args...) = ∫(fun)*measure` is expected"))
    end
end

"""
   get_name(ex::Expr)
extracts the name of the function from the head of a function definition. `get_name(:(a(u)))`  returns `:a`.
"""
get_name(ex::Expr) = ex.args[1]

"""
   get_parameters(ex::Expr)
extracts a tuple of parameters from the head of a function definition. `get_parameters(:(a(u,v)))` returns `:((u,v))`.
"""
get_parameters(ex::Expr) = Expr(:tuple, ex.args[2:end]...)


"""
   process_term(ex::Expr)

Extracts the integrand and the measure from an integral, discarding the integral symbol. For example: `process_term(:(∫(u*v)*dΩ))` returns `:(u*v)` and `:(dΩ)`.
"""
function process_term(ex::Expr)
    if @capture(ex, ∫(fun_) * meas_)
        return fun, meas
    elseif @capture(ex, -∫(fun_) * meas_)
        negfun = Expr(:call, :-, fun)
        return negfun, meas
    end
end

"""
   validate_integrand(ex::Expr)
checks if a term is defined by an integral.
"""
function validate_integrand(ex)
    return ex.args[1] != :∫ && throw(ArgumentError("Terms should be given by integrals."))
end

"""
   esc_non_params(expr,par)
escape every argument in `expr` that is not a parameter defined in `par`. 
"""
function esc_non_params(expr, par)
    newargs = []
    for arg in expr.args
        if arg isa Expr
            push!(newargs, esc_non_params(arg, par))
        else
            if arg in par
                push!(newargs, arg)
            else
                push!(newargs, esc(arg))
            end
        end
    end
    return Expr(expr.head, newargs...)
end

"""
    appearsin(arg::Symbol,expr)
checks if the symbol `arg` appears in the expression. If `expr` is a constant, it returns `false`, if `expr` is a Symbol, it returns `arg===expr`. Finally, if `expr isa Expr` it looks for `arg` inside `expr`.
"""
appearsin(arg::Symbol, thing) = false
appearsin(arg::Symbol, s::Symbol) = arg === s
appearsin(arg::Symbol, ex::Expr) = any(a -> appearsin(arg, a), ex.args)


"""
    get_factor(s::Expr,par)

Extracts a factor from a expression `s`. The function is meant to operate on an expression defining an integrand. `par` is an expression containing a tuple of parameters (the arguments of the integrand). For example, if `s = :((A*∇(u)⋅∇(v)))` and `par = :((u,v))`, the function retrieves the factor `A`. 
"""
get_factor(s::Union{Number, AbstractArray}, _) = s
function get_factor(s::Symbol, par)
    return if s in par || s in RESERVED_SYMBS
        :nothing
    else
        s
    end
end
function get_factor(s::Expr, par)
    return try
        eval(s)
        s
    catch
        for arg in s.args
            f = get_factor(arg, par)
            if f != :nothing
                return f
            end
        end
        :nothing
    end
end

"""
   cleanfactor(expr,factor)
takes the expressions `expr` and `factor`, where `factor` was obtained from `expr` using `get_factor` and cleans `expr` by removing `factor`, thus returning an expression containing only an operation between the parameters. 
"""
function cleanfactor(expr, factor, par)
    return if factor != :nothing
        if expr isa Expr && length(expr.args) == length(par) + 1
            op, left, right = expr.args
            if left == factor
                return right
            elseif right == factor
                return left
            else
                newleft = cleanfactor(left, factor, par)
                newright = cleanfactor(right, factor, par)
                return Expr(:call, op, newleft, newright)
            end
        elseif expr isa Expr && length(expr.args) == length(par)
            if expr.args[1] == :-
                nofactor = cleanfactor(expr.args[2], factor, par)
                Expr(:call, :-, nofactor)
            else
                expr
            end
        end
    else
        expr
    end
end

"""
   tensorize_body(body)
changes operation in the body of a function to """
function tensorize_body(body)
    sbody = string(body)
    sbody = replace(sbody, "*" => "⊗")
    sbody = replace(sbody, "⋅" => "⊗")
    # left, right = split(sbody, '⊗')
    # ordleft = Int8(occursin('∇', left) || occursin("∂x", left) || occursin("∂y", left))
    # ordright = Int8(occursin('∇', right) || occursin("∂x", right) || occursin("∂y", right))
    return Meta.parse(sbody) #, Order{(ordleft, ordright)}()
end
"""
    build_polyfun(par,body)
builds an expression defining the function `par->body`.  
"""
function build_polyfun(par, body)
    tbody = tensorize_body(body)
    return Expr(:->, Expr(:tuple, par...), tbody)
end

# macro term(expr)
#     head,termexpr = head_and_terms(expr)
#     name = get_name(head)
#     params = get_parameters(head)
#     integrand,meas = process_term(termexpr)
#     factor = get_factor(integrand,params)
#     funbody = cleanfactor(integrand,factor,params)
#     fun = build_polyfun(params,funbody)
#     Expr(:(=),esc(name),Expr(:call,:IntegrationTerm,esc_non_params(fun,params),(esc(factor)),esc(meas)))
# end
