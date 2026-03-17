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
struct ShapeFunction{O <: DiffOperator, D} end
function ShapeFunction(N = 1)
    (N isa Integer && 1 <= N <= 2) || throw(ArgumentError("Only `ShapeFunction`s of dimension `1` or `2` can be built."))
    return ShapeFunction{Identity, N}()
end

(::Identity)(s) = s
(::Gradient)(::ShapeFunction{Identity, 1}) = ShapeFunction{Gradient, 1}()
(::Laplacian)(::ShapeFunction{Identity, 1}) = ShapeFunction{Laplacian, 1}()
(::Divergence)(::ShapeFunction{Gradient, 1}) = ShapeFunction{Laplacian, 1}()

basis(sf::ShapeFunction{O, 1}, degs) where {O} = StandardBasis(degs)
operator(sf::ShapeFunction{O, D}) where {O, D} = O()


function (sf::ShapeFunction{O, D})(p::PolyField) where {O, D}
    D == length(p) || throw(DimensionMismatch("A `PolyField` of length $(length(p)) was used in place of a `ShapeFunction` of dimension $D."))
    op = operator(sf)
    return op(p)
end


"""
    abstract type CoeffType end
An abstract type for defining a treat that allow `IntegrationTerm`s to be dispatched to different methods for integration. Exact integration for constant coefficient terms and quadrature rules for variable coefficient terms.
"""
abstract type CoeffType end

struct ConstantCoeff <: CoeffType end
struct VariableCoeff <: CoeffType end

coefftype(_) = ConstantCoeff
coefftype(::Function) = VariableCoeff

"""
   Order{B}
A trait for storing the order of a form. It is used to deduce the contraction tensor. 
"""
struct Order{B} end

order(::Type{Identity}) = 0
order(::Type{Divergence}) = 1
order(::Type{Gradient}) = 1
order(::Type{Laplacian}) = 2
order(sf::ShapeFunction{O, D}) where {O, D} = order(O)

"""
   Integrand{C<:CoeffType,O<:Order,N}
A product of `ShapeFunction`s and a factor. `N` is the number of `ShapeFunction`s involved, `O` their order and  `C` the type of coefficient of the factor.  
"""
struct Integrand{C <: CoeffType, T, O <: Order, N}
    factor::T
    funs::Tuple
    function Integrand(fac, tup)
        all(s isa ShapeFunction for s in tup) || throw(ArgumentError("A `Tuple` of `ShapeFunction`s is expected in order to form an `Integrand`."))
        length(tup) <= 2 || throw(ArgumentError("At most two `ShapeFunction`s should be involved in an `Integrand`."))
        fac isa ShapeFunction && throw(ArgumentError("Factors should be constants (`<:Number` or `<:AbstractArray`) or functions. `ShapeFunction`s should be part of the tuple."))
        O = Order{Tuple(order(t) for t in tup)}
        N = length(tup)
        f = isnothing(fac) ? 1 : fac
        C = coefftype(f)
        T = typeof(f)
        return new{C, T, O, N}(f, tup)
    end
end
order(::Integrand{C, T, O, N}) where {C, T, O, N} = O()
coefftype(::Integrand{C, T, O, N}) where {C, T, O, N} = C()


Base.:*(sf₁::ShapeFunction{}, sf₂::ShapeFunction) = Integrand(nothing, (sf₁, sf₂))
function Base.:*(a::T, sf::ShapeFunction) where {T <: Union{Number, AbstractArray, Function}}
    return Integrand(a, (sf,))
end
Base.:*(sf::ShapeFunction, a::T) where {T <: Union{Number, AbstractArray, Function}} = a * sf
LinearAlgebra.dot(sf₁::ShapeFunction{}, sf₂::ShapeFunction) = sf₁ * sf₂
function LinearAlgebra.dot(a::T, sf::ShapeFunction) where {T <: Union{Number, AbstractArray, Function}}
    return a * sf
end
function LinearAlgebra.dot(sf::ShapeFunction, a::T) where {T <: Union{Number, AbstractArray, Function}}
    return a * sf
end

function Base.:*(inte::Integrand, sf::ShapeFunction)
    (; factor, funs) = inte
    tup = (funs..., sf)
    return Integrand(factor, tup)
end
LinearAlgebra.dot(inte::Integrand, sf::ShapeFunction) = inte * sf

function Base.:-(inte::Integrand)
    (; factor, funs) = inte
    return Integrand((-1) * factor, funs)
end
Base.:-(u::ShapeFunction) = Integrand(-1, (u,))

∫(inte::Integrand) = inte
∫(sf::ShapeFunction) = Integrand(nothing, (sf,))


#####################################################################################
#
#####################################################################################
abstract type AbstractForm end

"""
   Term 
"""
struct Term{C <: CoeffType, O <: Order, T, N, M <: Measure} <: AbstractForm
    integrand::Integrand{C, T, O, N}
    measure::M
end
coefftype(t::Term{C, T, O, N, M}) where {C, T, O, N, M} = C()
order(::Term{C, T, O, N, M}) where {C, T, O, N, M} = O()

struct Form{N} <: AbstractForm
    terms
end
Meshes.domainmesh(t::Term) = domainmesh(t.measure)
function Base.:*(inte::Integrand, meas)
    return Term(inte, meas)
end
function Base.:-(term::Term)
    (; integrand, measure) = term
    (; factor, funs) = integrand
    newintegrand = Integrand((-1) * factor, funs)
    return Term(newintegrand, measure)
end

function Base.:+(t₁::Term{C₁, O₁, N, M₁}, t₂::Term{C₂, O₂, N, M₂}) where {C₁, O₁, N, M₁, C₂, O₂, M₂}
    return Form{N}((t₁, t₂))
end
function Base.:-(t₁::Term{C₁, O₁, N, M₁}, t₂::Term{C₂, O₂, N, M₂}) where {C₁, O₁, N, M₁, C₂, O₂, M₂}
    return Form{N}((t₁, -t₂))
end
