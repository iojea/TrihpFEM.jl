
"""
See `@form`

```
   struct Form
```

A struct for storing the list of integrands and measures that define a linear of bilinear form.

It is recommended to build `Form`s using the `@form` macro. Howvever, a `Form` can be created directly by passing a function to be integrated and a `Measure`.

```
    julia> Ω = circmesh(0.1);
    julia> dΩ = Measure(Ω);
    julia> Form((u,v)->∇(u)⋅∇(v),dΩ)
```

If the `Form` has two or more terms, tuples of functions to be integrated and its corresponding  `Measures` are expected.
```
   julia> term₁(u,v) = ∇(u)⋅∇(v);
   julia> term₂(u,v) = 2u*v;
   julia> Form((term₁,term₂),(dΩ,dΩ)) 
```
But when every term is integrated with respect to the same `Measure`, it can be passed once:
```
   julia> Form((term₁,term₂),dΩ) 
```
"""
struct Form{N}
    integrands::Tuple
    measures::Tuple
end

function get_nargs(t::Tuple)
    tnargs = Tuple(first(methods(tt)).nargs-1 for tt in t)
    all(tnargs .== tnargs[1]) || throw(ArgumentError("Number of arguments vary from term to term."))
    tnargs[1]
end

function Form(f::Function,m::Measure)
    nargs = get_nargs((f,))
    Form{nargs}((f,),(m,))
end
#Form(t::Tuple,m::Measure) = Form(t,Tuple(m for _ in 1:length(t)))


"""

    _form(x...)
A convenient constructor for `Form`.
"""
function _form(x...)
    n = length(x)
    n%2==0 || throw(ArgumentError("Malformed expression. Each term must be of the form `∫(fun)*dΩ` where `fun` is some function and `dΩ` is a `Measure`."))
    k = n÷2
    nargs = get_nargs(tuple(x[1:k]...))
    all(typeof(mm)<:Measure for mm in x[k+1:end]) || throw(ArgumentError("Measures are not measures."))
    Form{nargs}(tuple(x[1:k]...),tuple(x[k+1:end]...))
end

"""

    bad_integrand(expr)
Checks for malformed integrands (e.g.: integrands with sums). 
""" 
bad_integrand(::Any) = false
function bad_integrand(expr::Expr)
    expr.head == :call && expr.args[1] == :+ && return true
    expr.head == :call && expr.args[1] == :- && length(expr.args)>2 && return true
    expr.head == :call && return any(bad_integrand.(expr.args[2:end]))
    expr.head == :block && return any(bad_integrand.(expr.args))
    return false
end

"""

    fracture(expr,arg,signed,others...)
splits expressions into integrands and corresponding measures.   
"""
function fracture(expr,arg,signed,others...)
    if expr.head == :call && expr.args[1] == :*
        measure = esc(expr.args[3])
        if expr.args[2].head == :call && expr.args[2].args[1] == :∫
            if signed 
                integrand = Expr(:call,:-,expr.args[2].args[2])
            else
                integrand = expr.args[2].args[2]
            end
            bad_integrand(integrand) && throw(error("Integrands containing sums are not supported. Please use one integral per term in your form. "))
            fun_body = Expr(:->,arg,integrand)
            return fun_body,measure,others...
        elseif expr.args[2].head == :call && expr.args[2].args[1] in (:*,:⋅)
            throw(error("Multiplication outside of integrals are not supported yet. Place the multiplication inside of the integrand instead.")) 
        end
    end
    if expr.head == :call && expr.args[1] == :+
        terms1 = fracture(expr.args[2],arg,false)
        terms2 = fracture(expr.args[3],arg,false)
        return terms1...,terms2...,others...
    elseif expr.head == :call && expr.args[1] == :-
        terms1 = fracture(expr.args[2],arg,false)
        terms2 = fracture(expr.args[3],arg,true)
        return terms1...,terms2...,others...
    end
end

"""

    get_name(expr)
extracts the name to be assigned to the `Form`  
"""
get_name(expr) = esc(expr.args[1].args[1])
"""

    get_parameters(expr)
extracts the names of the parameters of the `Form`  
"""
get_parameters(expr) = Expr(:tuple,expr.args[1].args[2:end]...)
"""

    get_terms(expr)
separate different terms of a `Form`.  
"""
get_terms(expr) = expr.args[2].args[2]

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
```
   @form form_definition
```
Creates a `Form` 
```
   julia> A = rand(2,2);
   julia> Ω = circmesh(0.1);
   julia> dΩ = Measure(Ω)
   julia> @form a(u,v) = ∫((A*∇(u))⋅∇(v))*dΩ
```
Creates a `Form` named `a`.

Linear `Forms` are built the same way
```
   julia> f(x) = x[1]^2*sin(x[2])
   julia> @form b(v) = ∫(v*f)*dΩ
```
It is also possible to create `Form`s with several terms.
```
    julia> @form a(u,v) = ∫(∇(u)*∇(v))*dΩ + ∫(3u*v)*dΩ
```
It is not allowed to have several terms under the same integral sign.
```
    julia> @form a(u,v)=∫(∇(u)⋅∇(v)+u*v)*dΩ
    ERROR: LoadError: Integrands containing sums are not supported. Please use one integral per term in your form.
```
This is mainly for performance reasons. The terms of first and second order are integrated differently, so they cannot be parts of the same function. 
""" 
macro form(expr)
       name = get_name(expr)
       par = get_parameters(expr)
       terms = get_terms(expr)
       tira = fracture(terms,par,false)
       integrands = [esc_non_params(i,par) for i in tira[1:2:end]]
       Expr(:(=),name,Expr(:call,:_form,integrands...,tira[2:2:end]...))
end
