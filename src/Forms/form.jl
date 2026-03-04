struct Form{N}
    terms
end
Form(x::IntegrationTerm{C, N}) where {C, N} = Form{N}((x,))

function Form(x...)
    all(isa.(x, IntegrationTerm)) || throw(ArgumentError("All parameters should be `IntegrationTerm`s."))
    N = numargs(x[1])
    all(numargs.(x) .== N) || throw(ArgumentError("Terms should operate on the same number of arguments."))
    return Form{N}(x)
end

terms(form::Form) = form.terms
function Base.:+(t₁::IntegrationTerm{C, N}, t₂::IntegrationTerm{D, N}) where {C, D, N}
    return Form(t₁, t₂)
end

Meshes.domainmesh(f::Form) = domainmesh(f.terms[1])

macro form(expr)
    if @capture(expr, name_(args__) = multiterms_)
        terms = separate_terms(strip_block(multiterms))
        intterms = []
        for t in terms
            integrand, meas = process_term(t)
            factor = get_factor(integrand, args)
            funbody = cleanfactor(integrand, factor, args)
            fun = build_polyfun(args, funbody)
            it = Expr(:call, :IntegrationTerm, esc_non_params(fun, args), esc(factor), esc(meas))
            push!(intterms, it)
        end
        return Expr(:(=), esc(name), Expr(:call, :Form, intterms...))
    else
        throw(ArgumentError("Malformed expression. A term of the form `name(args...) = ∫(fun)*measure` is expected"))
    end
end


"""
    head_and_terms(ex::Expr)
Separates both sides of an assignment. For example: `head_and_terms(:(a(u) = 2u))` returns `:(a(u))` and `:(2*u)`.
"""
function head_and_terms(ex::Expr)
    ex.head != :(=)  && throw(ArgumentError("An assignment is needed."))
    head = ex.args[1]
    tail = strip_block(ex.args[2])
    return head, tail
end
"""
   separate_terms(ex::Expr)
Given an expression `ex` containing a sum of integrals, it returns a `Vector{Expr}` containing each term (each integral).
"""
function separate_terms(ex::Expr)
    terms = Expr[]
    if ex.args[1] == :+
        for t in ex.args[2:end]
            newterms = separate_terms(t)
            push!(terms, newterms...)
        end
    elseif ex.args[1] == :-
        if length(ex.args) == 2
            newterms = separate_terms(ex.args[2])
            for nt in newterms
                nnt = Expr(:call, :-, nt)
                push!(terms, nnt)
            end
        elseif length(ex.args) == 3
            newterms1 = separate_terms(ex.args[2])
            push!(terms, newterms1...)
            newterms2 = separate_terms(ex.args[3])
            for nt in newterms2
                nnt = Expr(:call, :-, nt)
                push!(terms, nnt)
            end
        end
    else
        push!(terms, ex)
    end
    return terms
end
