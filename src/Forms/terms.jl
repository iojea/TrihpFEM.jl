abstract type CoeffType end

struct ConstantCoeff <: CoeffType end
struct VariableCoeff <: CoeffType end

coefftype(_) = ConstantCoeff()
coefftype(::Function) = VariableCoeff()
coefftype(::DiffOperator) = ConstantCoeff()


struct IntegrationTerm{C<:CoeffType,N}
    polyfun
    factor
    measure
end

function IntegrationTerm(pf,f,meas)
    T = typeof(coefftype(f))
    N = first(methods(pf)).nargs-1
    IntegrationTerm{T,N}(pf,f,meas)
end
numargs(::IntegrationTerm{C,N}) where {C,N} = N


"""
  A list of reserved symbols to discard when looking for the `factor` in an `IntegrationTerm`  
"""
RESEVED_SYMBS = (:*,:⋅,:-,:∂x,:∂y,:gradient,:divergence,:laplacian,:∇,:Δ,:ε)

"""
    strip_block(ex::Expr)
removes the `:block` part of an expression, returning the meaningful part.
"""
strip_block(ex::Expr) = ex.head == :block ? ex.args[2] : ex

"""
    head_and_terms(ex::Expr)
Separates both sides of an assignment. For example: `head_and_terms(:(a(u) = 2u))` returns `:(a(u))` and `:(2*u)`.
"""
function head_and_terms(ex::Expr)
    ex.head != :(=)  && throw(ArgumentError("An assignment is needed."))
    head = ex.args[1]
    tail = strip_block(ex.args[2])
    return head,tail
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
get_parameters(ex::Expr) = Expr(:tuple,ex.args[2:end]...)


"""
   separate_terms(ex::Expr)
Given an expression `ex` containing a sum of integrals, it returns a `Vector{Expr}` containing each term (each integral).
"""
function separate_terms(ex::Expr)
    terms = Expr[]
    if ex.args[1] == :+
        for t in ex.args[2:end]
            newterms = separate_terms(t)
            push!(terms,newterms...)
        end
    elseif ex.args[1] == :-
        if length(ex.args)==2
           newterms = separate_terms(ex.args[2])
            for nt in newterms
                nnt = Expr(:call,:-,nt)
                push!(terms,nnt)
            end
        elseif length(ex.args)==3
            newterms1 = separate_terms(ex.args[2])
            push!(terms,newterms1...)
            newterms2 = separate_terms(ex.args[3])
            for nt in newterms2
                nnt = Expr(:call,:-,nt)
                push!(terms,nnt)
            end
        end
    else
        push!(terms,ex)
    end
    terms
end

"""
   process_term(ex::Expr)

Extracts the integrand and the measure from an integral, discarding the integral symbol. For example: `process_term(:(∫(u*v)*dΩ))` returns `:(u*v)` and `:(dΩ)`.
"""
function process_term(ex::Expr)
    if ex.args[1] == :*
        validate_integrand(ex.args[2])
        return ex.args[2].args[2],ex.args[3]        
    elseif ex.args[1] == :-
        integrand,measure = process_term(ex)
        negint = Expr(:call,:-,integrand)
        return negint,measure
    end
end

"""
   validate_integrand(ex::Expr)
checks if a term is defined by an integral.  
"""
function validate_integrand(ex)
    ex.args[1] != :∫ && throw(ArgumentError("Terms should be given by integrals."))
end

"""
   esc_non_params(expr,par)
escape every argument in `expr` that is not a parameter defined in `par`. 
"""
function esc_non_params(expr,par)
    newargs = []
    for arg in expr.args
        if arg isa Expr
            push!(newargs,esc_non_params(arg,par))
        else
            if arg in par.args
                push!(newargs,arg)
            else
                push!(newargs,esc(arg))
            end
        end
    end
    return Expr(expr.head,newargs...)
end

"""
    appearsin(arg::Symbol,expr)
checks if the symbol `arg` appears in the expression. If `expr` is a constant, it returns `false`, if `expr` is a Symbol, it returns `arg===expr`. Finally, if `expr isa Expr` it looks for `arg` inside `expr`.
"""
appearsin(arg::Symbol,thing) = false
appearsin(arg::Symbol,s::Symbol) = arg === s
appearsin(arg::Symbol,ex::Expr) = any(a->appearsin(arg,a),ex.args)


"""
    get_factor(s::Expr,par)

Extracts a factor from a expression `s`. The function is meant to operate on an expression defining an integrand. `par` is an expression containing a tuple of parameters (the arguments of the integrand). For example, if `s = :((A*∇(u)⋅∇(v)))` and `par = :((u,v))`, the function retrieves the factor `A`. 
"""
get_factor(s::Union{Number,AbstractArray},_) = s
function get_factor(s::Symbol,par)
    if appearsin(s,par) || s in RESEVED_SYMBS
        :nothing
    else
        s
    end
end
function get_factor(s::Expr,par)
    try eval(s)
        s
    catch
        for arg in s.args
            f = get_factor(arg,par)
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
function cleanfactor(expr,factor,par)
    if factor != :nothing
        if expr isa Expr &&length(expr.args)==length(par.args)+1
            op,left,right = expr.args
            if left == factor
                return right
            elseif right == factor
                return left
            else
                newleft = cleanfactor(left,factor,par)
                newright = cleanfactor(right,factor,par)
                return Expr(:call,op,newleft,newright)
            end
        elseif expr isa Expr && length(expr.args)==length(par.args)
            if expr.args[1] == :-
                nofactor = cleanfactor(expr.args[2],factor,par)
                Expr(:call,:-,nofactor)
            else
                expr
            end
        end
    else
        expr
    end
end

"""
    build_polyfun(par,body)
builds an expression defining the function `par->body`.  
"""
build_polyfun(par,body) = Expr(:->,par,body)

macro term(expr)
    head,termexpr = head_and_terms(expr)
    name = get_name(head)
    params = get_parameters(head)
    integrand,meas = process_term(termexpr)
    factor = get_factor(integrand,params)
    funbody = cleanfactor(integrand,factor,params)
    fun = build_polyfun(params,funbody)
    Expr(:(=),esc(name),Expr(:call,:IntegrationTerm,esc_non_params(fun,params),$(esc(factor)),esc(meas)))
end







