struct Form{N}
    terms  
end
Form(x::IntegrationTerm{C,N}) where {C,N} = Form{N}((x,))

function Form(x...)
    all(isa.(x,IntegrationTerm)) || throw(ArgumentError("All parameters should be `IntegrationTerm`s."))
    N = numargs(x[1])
    all(numargs.(x).==N) || throw(ArgumentError("Terms should operate on the same number of arguments."))
    Form{N}(x)
end

function Base.:+(t₁::IntegrationTerm{C,N},t₂::IntegrationTerm{D,N}) where {C,D,N}
    Form(t₁,t₂)
end

macro form(expr)
    head,termsexpr = head_and_terms(expr)
    name = get_name(head)
    params = get_parameters(head)
    terms = separate_terms(termsexpr)
    intterms = []
    for t in terms
        integrand,meas = process_term(t)
        factor = get_factor(integrand,params)
        funbody = cleanfactor(integrand,factor,params)
        fun = build_polyfun(params,funbody)
        it = Expr(:call,:IntegrationTerm,esc_non_params(fun,params),esc(factor),esc(meas))
        push!(intterms,it)
    end
    Expr(:(=),esc(name),Expr(:call,:Form,intterms...))
end



