abstract type CoeffType end

struct Constant <: CoeffType end
struct Variable <: CoeffType end

coefftype(_) = Constant()
coefftype(::Function) = Variable()
coefftype(::DiffOperator) = Constant()


struct IntegrationTerm{C<:CoeffType}
    polyfun
    factor
    measure
end

IntegrationTerm(pf,f,meas) = IntegrationTerm{typeof(coefftype(f))}(pf,f,meas)
IntegrationTerm(pf,meas) = IntegrationTerm{Constant}(pf,I,meas)

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
    try es = Core.eval(__module__,s)
        if !(es isa DiffOperator) && !(s in (:*,:⋅,:-)) && !(appearsin(s,par))
            s
        end
    catch
        :nothing
    end
end
function get_factor(s::Expr,par)
    try Core.eval(__module__,s)
        s
    catch
        for arg in s.args
            f = get_factor(arg,par)
            if !isnothing(f)
                return f
            end
        end
        return :nothing
    end
end

"""
   cleanfactor(expr,factor)
takes the expressions `expr` and `factor`, where `factor` was obtained from `expr` using `get_factor` and cleans `expr` by removing `factor`, thus returning an expression containing only an operation between the parameters. 
"""
function cleanfactor(expr,factor)
    if length(expr.args)==n
        op,left,right = expr.args
        if left == factor
            return right
        elseif right == factor
            return left
        else
            newleft = cleanfactor(left,factor)
            newright = cleanfactor(right,factor)
            return Expr(:call,op,newleft,newright)
        end
    elseif length(expr.args)==n-1
        if expr.args[1] == :-
            nofactor = cleanfactor(expr.args[2],factor)
            Expr(:call,:-,nofactor)
        else
            expr
        end
    end
end

"""
    build_polyfun(par,body)
builds an expression defining the function `par->body`.  
"""
build_polyfun(par,body) = Expr(:->,par,body)

macro term(expr)
    head,termsexpr = head_and_terms(expr)
    name = get_name(head)
    params = get_parameters(head)
    terms = separate_terms(termsexpr)
    intterms = []
    for t in terms
        integrand,meas = process_term(t)
        factor = get_factor(integrand,params)
        funbody = cleanfactor(integrand,factor)
        fun = build_polyfun(params,funbody)
        it = Expr(:call,:IntegrationTerm,esc_non_params(fun,params),factor,esc(meas))
        push!(intterms,it)
    end
    Expr(:(=),esc(name),Expr(:call,:Form,intterms...))
end



struct Form
    terms  
end
Form(x::IntegrationTerm) = Form((x,))
Form(x...) = Form(x)






