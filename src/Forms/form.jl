"""
    ShapeFunction(dim)
Defines a standard shape function, of dimension `dim`.
```
julia> u = ShapeFunction(1)
```
A gradient shape function can be build with `ShapeFunction{Gradient,1}()`. This means: the gradient of a shape function with `dim=1`. However the preferred way to build `ShapeFunction`s of this kind is to a apply the differential operator to the standard `ShapeFunction`:
```
julia> grad_u = ∇(u) 
```
"""
struct ShapeFunction{O<:DiffOperator,D} end
function ShapeFunction(N=1)
    (N isa Integer && 1<=N<=2) || throw(ArgumentError("Only `ShapeFunction`s of dimension `1` or `2` can be built."))
    ShapeFunction{Identity,N}()
end

(::Gradient)(::ShapeFunction{Identity,1}) = ShapeFunction{Gradient,1}()
(::Laplacian)(::ShapeFunction{Identity,1}) = ShapeFunction{Laplacian,1}()
(::Divergence)(::ShapeFunction{Gradient,1}) = ShapeFunction{Laplacian,1}()

basis(sf::ShapeFunction{O,1},degs) where O = StandardBasis(degs) 
operator(sf::ShapeFunction{O,D}) where {O,D} = O()
"""
    abstract type CoeffType end
An abstract type for defining a treat that allow `IntegrationTerm`s to be dispatched to different methods for integration. Exact integration for constant coefficient terms and quadrature rules for variable coefficient terms.
"""
abstract type CoeffType end

struct ConstantCoeff <: CoeffType end
struct VariableCoeff <: CoeffType end
struct NoCoeff <: CoeffType end

coefftype(_) = ConstantCoeff
coefftype(::Function) = VariableCoeff
coefftype(::Nothing) = NoCoeff

struct Order{B} end

order(::Type{Identity}) = 0
order(::Type{Divergence}) = 1
order(::Type{Gradient}) = 1
order(::Type{Laplacian}) = 2
order(sf::ShapeFunction{O,D}) where {O,D} = order(O)

tensorize(::Nothing) = Tensor{1,1,Float64}((1,))
function tensorize(arr::A) where A<:AbstractArray
    ord = ndims(arr)
    Tensor{ord,2,eltype(arr)}(arr)
end

    
struct Integrand{C<:CoeffType,T<:Tensor,O<:Order,N}
    factor::T
    funs::Tuple
    function Integrand(fac,tup)
        all(s isa ShapeFunction for s in tup) || throw(ArgumentError("A `Tuple` of `ShapeFunction`s is expected in order to form an `Integrand`."))
        length(tup)<=2 || throw(ArgumentError("At most two `ShapeFunction`s should be involved in an `Integrand`."))
        fac isa ShapeFunction && throw(ArgumentError("Factors should be constants (`<:Number` or `<:AbstractArray`) or functions. `ShapeFunction`s should be part of the tuple."))
        f = tensorize(fac)
        C = coefftype(fac)
        O = Order{Tuple(order(t) for t in tup)}
        N = length(tup)
        new{C,typeof(f),O,N}(f,tup)
    end
end
order(::Integrand{C,T,O,N}) where {C,T,O,N} = O()
coefftype(::Integrand{C,T,O,N}) where {C,T,O,N} = C()

Base.:*(sf₁::ShapeFunction{},sf₂::ShapeFunction) = Integrand(nothing,(sf₁,sf₂))
function Base.:*(a::T,sf::ShapeFunction) where {T<:Union{Number,AbstractArray,Function}}
    Integrand(a,(sf,))
end
Base.:*(sf::ShapeFunction,a::T) where {T<:Union{Number,AbstractArray,Function}} = a*sf
LinearAlgebra.dot(sf₁::ShapeFunction{},sf₂::ShapeFunction) = sf₁*sf₂
function LinearAlgebra.dot(a::T,sf::ShapeFunction) where {T<:Union{Number,AbstractArray,Function}}
    a*sf
end
function LinearAlgebra.dot(sf::ShapeFunction,a::T) where {T<:Union{Number,AbstractArray,Function}}
    a*sf
end

function Base.:*(inte::Integrand,sf::ShapeFunction)
    (;factor,funs) = inte
    tup = (funs...,sf)
    Integrand(factor,tup)
end
LinearAlgebra.dot(inte::Integrand,sf::ShapeFunction) = inte*sf

∫(inte::Integrand) = inte
∫(sf::ShapeFunction) = Integrand(nothing,(sf,))


abstract type AbstractForm end

struct Term{C<:CoeffType,O<:Order,N,M<:Measure} <: AbstractForm
    integrand::Integrand{C,O,N}
    measure::M
end
coefftype(t::Term{C,O,N,M}) where {C,O,N,M} = C


struct Form{N}
    terms
end

function Base.:*(inte::Integrand{C,O,N},meas::M) where {C,O,N,M}
    Term{C,O,N,M}(inte,meas)
end
function Base.:-(term::Term{C,O,N,M}) where {C,O,N,M}
    (;integrand,measure) = term
    (;factor,funs)  = integrand
    newintegrand = Integrand((-1)*factor,funs)
    Term(newintegrand,measure)
end

function Base.:+(t₁::Term{C₁,O₁,N,M₁},t₂::Term{C₂,O₂,N,M₂}) where {C₁,O₁,N,M₁,C₂,O₂,M₂}
    Form{N}((t₁,t₂))
end
function Base.:-(t₁::Term{C₁,O₁,N,M₁},t₂::Term{C₂,O₂,N,M₂}) where {C₁,O₁,N,M₁,C₂,O₂,M₂}
    Form{N}((t₁,-t₂))
end





# Form(x::IntegrationTerm{C, N}) where {C, N} = Form{N}((x,))

# function Form(x...)
#     all(isa.(x, IntegrationTerm)) || throw(ArgumentError("All parameters should be `IntegrationTerm`s."))
#     N = numargs(x[1])
#     all(numargs.(x) .== N) || throw(ArgumentError("Terms should operate on the same number of arguments."))
#     return Form{N}(x)
# end

# terms(form::Form) = form.terms
# function Base.:+(t₁::IntegrationTerm{C, N}, t₂::IntegrationTerm{D, N}) where {C, D, N}
#     return Form(t₁, t₂)
# end

# Meshes.domainmesh(f::Form) = domainmesh(f.terms[1])

# macro form(expr)
#     if @capture(expr, name_(args__) = multiterms_)
#         terms = separate_terms(strip_block(multiterms))
#         intterms = []
#         for t in terms
#             integrand, meas = process_term(t)
#             factor = get_factor(integrand, args)
#             funbody = cleanfactor(integrand, factor, args)
#             fun = build_polyfun(args, funbody)
#             it = Expr(:call, :IntegrationTerm, esc_non_params(fun, args), esc(factor), esc(meas))
#             push!(intterms, it)
#         end
#         return Expr(:(=), esc(name), Expr(:call, :Form, intterms...))
#     else
#         throw(ArgumentError("Malformed expression. A term of the form `name(args...) = ∫(fun)*measure` is expected"))
#     end
# end


# """
#     head_and_terms(ex::Expr)
# Separates both sides of an assignment. For example: `head_and_terms(:(a(u) = 2u))` returns `:(a(u))` and `:(2*u)`.
# """
# function head_and_terms(ex::Expr)
#     ex.head != :(=)  && throw(ArgumentError("An assignment is needed."))
#     head = ex.args[1]
#     tail = strip_block(ex.args[2])
#     return head, tail
# end
# """
#    separate_terms(ex::Expr)
# Given an expression `ex` containing a sum of integrals, it returns a `Vector{Expr}` containing each term (each integral).
# """
# function separate_terms(ex::Expr)
#     terms = Expr[]
#     if ex.args[1] == :+
#         for t in ex.args[2:end]
#             newterms = separate_terms(t)
#             push!(terms, newterms...)
#         end
#     elseif ex.args[1] == :-
#         if length(ex.args) == 2
#             newterms = separate_terms(ex.args[2])
#             for nt in newterms
#                 nnt = Expr(:call, :-, nt)
#                 push!(terms, nnt)
#             end
#         elseif length(ex.args) == 3
#             newterms1 = separate_terms(ex.args[2])
#             push!(terms, newterms1...)
#             newterms2 = separate_terms(ex.args[3])
#             for nt in newterms2
#                 nnt = Expr(:call, :-, nt)
#                 push!(terms, nnt)
#             end
#         end
#     else
#         push!(terms, ex)
#     end
#     return terms
# end
