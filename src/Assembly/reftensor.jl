const MeshElement = Union{Edge, Triangle}


"""
   RefTensor(t::Term)
   creates a memoisation struct for storing local tensors or evaluations of tensors of polinomials on quadrature points. A `RefTensor` per `Term` is necessary for integration. `ConstantCoeff` `Term`s store local numeric tensors while `VariableCoeff` `Term`s store local polynomial tensors evaluated on quadrature points.  
"""
struct RefTensor{K, T, R <: Union{Nothing, Quadrature}}
    meas::AbstractMeasure
    t::Tuple
    sch::R
    dict::Dict{K, T}
end


#### NOTICE that the constructor for `RefTensor` does not build a complete dictionary: it only initializes it with the first possible element (p = (0,0,0)). For this reason, the constructor use the `obtainkey` function, whereas the `keybasis` function is necessary when considering a new basis.

function RefTensor(t::Term{N, ConstantCoeff, T, M}) where {N, T <: Tuple, F, I, P, M <: AbstractMeasure{F, I, P}}
    (; integrand, measure) = t
    (; operation) = integrand
    D = dim(measure)
    op = _scalars_out(operation, F)
    globpart = global_part(op)
    sfs = get_shape_functions(operation)
    length(sfs) == N || throw(ArgumentError("Malformed term: a `Term{$N}` should be formed by $N `ShapeFunction`s, but $(length(sfs)) are present."))
    bs = ((first(StandardBasis(ntuple(_ -> zero(P), 3))) for _ in sfs)...,)
    key = obtainkey(measure, bs)
    poly_tensor = _topolytensor(globpart, bs...)
    tensor = ref_integrate.(poly_tensor, _val(measure))
    U = typeof(tensor)
    K = typeof(key)
    d = Dict{K, U}(key => tensor)
    return RefTensor{K, U, Nothing}(measure, globpart, nothing, d)
end

function RefTensor(t::Term{N, VariableCoeff, T, M}) where {N, T <: Tuple, F, I, P, M <: AbstractMeasure{F, I, P}}
    (; integrand, measure) = t
    (; operation) = integrand
    (; mesh, sch) = measure
    (; points, weights) = sch
    D = dim(measure)
    v = _val(measure)
    J = _jac(v)
    pts = _schpoints(points, v)
    op = _scalars_out(operation, F)
    globpart = global_part(op)
    sfs = get_shape_functions(operation)
    length(sfs) == N || throw(ArgumentError("Malformed term: a `Term{$N}` should be formed by $N `ShapeFunction`s, but $(length(sfs)) are present."))
    bs = ((first(StandardBasis(ntuple(_ -> zero(P), 3))) for _ in sfs)...,)
    key = obtainkey(measure, bs)
    poly_tensor = _topolytensor(globpart, bs...)
    # IMPORTANT: The factor 2 comes from integration on the reference triangle, with area 2.
    tensor = collect_as(FixedSizeArrayDefault, (J * w * poly_tensor(x) for (x, w) in zip(pts, weights)))
    U = typeof(tensor)
    K = typeof(key)
    d = Dict{K, U}(key => tensor)
    return RefTensor{K, U, typeof(sch)}(measure, globpart, sch, d)
end

"""
   keybasis(t::Triangle,m::Measure)
   keybasis(e::Edge,m::BoundaryMeasure)
returns the key and the basis corresponding to a element of an `AbstractMeasure`. The key is an  `NTuple{4,P} where P<:Integer` where the first number codifies the type of element:
- 0: Triangle
- 1: Edge with minimum degree in the triangle (so it matches the horizontal edge in the reference triangle)
- 2: Edge with intermediate degree in the triangle (so it matches the vertical edge in the reference triangle)
- 3: Edge with maximum degree in the triangle (so it matches the diagonal edge in the reference triangle)
The last three numbers in the key are the degrees of the triangle.

Note that this only works for _boundary edges_, where only ONE triangle is involved. 
"""
function keybasis(t::Triangle, m::Measure{F, I, P}) where {F, I, P}
    p, _ = psortednodes(t, mesh(m))
    return ((zero(P), p...), StandardBasis(p))
end
function keybasis(e::Edge, m::BoundaryMeasure{F, I, P}) where {F, I, P}
    els = elements(m)
    isin, tok = gettoken(els, e)
    isin || throw(ArgumentError("The `Edge` does not belong to the given `BoundaryMeasure`"))
    _, att = gettokenvalue(els, tok)
    adj = adjacents(att)
    last(adj) == 0 || throw(ArgumentError("`Edge` seems to be interior to the domain."))
    t = Triangle(e[1], e[2], adj[1])
    p, t_eds = psortededges(t, mesh(m))
    j = findfirst(isequal(e), t_eds)
    return ((P(j), p...), StandardBasis(p))
end


"""
    _val(::Measure)
    _val(::BoundaryMeasure)
    _val(::Triangle,::Measure)
    _val(e::Edge,m::BoundaryMeasure)

computes the first value of the key for the `Dict` stored un a `RefTensor`. Check [keybasis](@ref).
"""
_val(::Measure) = Val(0)
_val(::BoundaryMeasure) = Val(1)
_val(::Triangle, ::Measure) = Val(0)
function _val(e::Edge{I}, meas::BoundaryMeasure) where {I}
    m = mesh(meas)
    eds = elements(meas)
    adj = adjacents(last(eds[e]))
    t = Triangle{I}(e[1], e[2], adj[1])
    _, t_eds = psortededges(t, m)
    j = findfirst(isequal(e), t_eds)
    return Val(j)
end

"""
   obtainkey(::AbstractMeasure,bs)
computes the key for initiallizing the `Dict` in `RefTensor` 
"""
obtainkey(::Measure, bs) = (0, degs.(bs)...)
obtainkey(::BoundaryMeasure, bs) = (1, degs.(bs)...)

function (lt::RefTensor{K, T, Nothing})(el::MeshElement, bs...) where {K, T}
    (; meas, t, dict) = lt
    v = _val(el, meas)
    key = _tokey(v, degs.(bs)...)
    haskey(dict, key) && return dict[key]
    poly_tensor = _topolytensor(t, bs...)
    tensor = ref_integrate.(poly_tensor, v)
    dict[key] = tensor
    return tensor
end
function (lt::RefTensor{K, T, R})(el::MeshElement, bs...) where {K, T, R <: Quadrature}
    (; meas, t, sch, dict) = lt
    v = _val(el, meas)
    J = _jac(v)
    key = _tokey(v, degs.(bs)...)
    (; points, weights) = sch
    pts = _schpoints(points, v)
    haskey(lt.dict, key) && return dict[key]
    poly_tensor = _topolytensor(t, bs...)

    # IMPORTANT: The factor 2 comes from integration on the reference triangle, with area 2.
    tensor = collect_as(FixedSizeVectorDefault, J * poly_tensor(x) for (x, w) in zip(pts, weights))
    dict[key] = tensor
    return tensor
end

_tokey(::Val{N}, a...) where {N} = (N, a...)
"""
   _isglobal(a)
Decides if a factor `a` in an `Operation` belongs to its global part (`ShapeFunction`s and constants).
"""
_isglobal(::Any) = true
_isglobal(::Function) = false
_isglobal(::DiffOperator) = false


"""
   _scalars_out(op::Operation,::Type{F})
takes `Operation` `op` and moves scalars (numbers or scalar functions) to the left, returning a reorderd `Operation`. This is done because it is the easiest way to ensure that scalars operate against other factors with a simple product (*) and not an outer product (⊗). The product itself will always reduce to a simple product (because ⊗ applied to scalars is just *), but other adjacent products will produce unnessesary large tensors. E.g.: if `B` is a matrix from an `AffineToRef` and `z isa Number` the product `∇(u)'*B'*z*B*∇(v)` would be  converted into a global term `∇(u)'⊗z⊗∇(v)` and a local term `B'⊗B`. However, is the scalar `z` is extracted, the new `Operation` is `z*∇(u)'*B'*B*∇(v)` which results in the global term `z*∇u'⊗∇(v)` and the local term `B'*B`.
"""
function _scalars_out(op::Operation{C}, ::Type{F}) where {C <: CoeffType, F}
    for i in 3:2:length(op)
        if _isscalar(op[i], F)
            op = Operation{C}(op[i], *, op[1:(i - 2)]..., op[(i + 1):end]...)
        end
    end
    return op
end
_isscalar(::Any, _) = false
_isscalar(::Number, _) = true
function _isscalar(f::Function, ::Type{F}) where {F}
    return f(zeros(F, 2)) isa Number
end

"""
   global_part(op::Operation)
Extracts the global part of an `Operation` as a `Tuple`. Products are transformed (*->⊗) when needed. 
"""
function global_part(op::Operation)
    ind = findall(_isglobal, op.parts)
    vout = Vector{Any}(undef, 2length(ind) - 1)
    vout[1:2:end] .= op[ind]
    @inbounds for i in 2:length(ind)
        if ind[i] == ind[i - 1] + 2
            vout[2i - 2] = op[ind[i] - 1]
        else
            vout[2i - 2] = DICT_OP[op[ind[i] - 1]]
        end
    end
    return tuple(vout...)
end


# get_shape_functions(sf::ShapeFunction, mesh::HPTriangulation) = _adjust_to_mesh(sf, mesh)
# get_shape_functions(::Any, ::HPTriangulation) = ()
# get_shape_functions(t::Tuple, ::HPTriangulation) = filter(Base.Fix2(isa, ShapeFunction), t)
# function get_shape_functions(operation::Operation, m::HPTriangulation)
#     a, b = ((get_shape_functions(p, m) for p in operation.parts)...,)
#     a = a isa ShapeFunction ? (a,) : get_shape_functions(a, m)
#     b = b isa ShapeFunction ? (b,) : get_shape_functions(b, m)
#     return filter(Base.Fix2(isa, ShapeFunction), (a..., b...))
# end


"""
    _topolytensor(t::Tuple,trial,test)
    _topolytensor(t::Tuple,test)
    _topolytensor(t::Tuple,aff::AffineToRef)
Evaluates `t` recursively to produce a `Tensor`. `t` must be the global or local part of an `Operation`.

    _topolytensor(t::Tuple,trial,test,sch::Quadrature)
    _topolytensor(t::Tuple,test,sch::Quadrature)
    _topolytensor(t::Tuple,)
"""
_topolytensor(a::Any, _, _) = a
_topolytensor(a::Any, _) = a
_topolytensor(sf::ShapeFunction{Trial}, trial, _) = operator(sf)(trial)
_topolytensor(sf::ShapeFunction{Test}, _, test) = operator(sf)(test)
_topolytensor(t::NTuple{1}, trial, test) = _topolytensor(t[1], trial, test)
_topolytensor(t::Tuple, trial, test) = t[2](_topolytensor(t[1], trial, test), _topolytensor(t[3:end], trial, test))
_topolytensor(sf::ShapeFunction{Test}, test) = operator(sf)(test)
_topolytensor(t::NTuple{1}, test) = _topolytensor(t[1], test)
_topolytensor(t::Tuple, test) = t[2](_topolytensor(t[1], test), _topolytensor(t[3:end], test))


_schpoints(points, ::Val{0}) = points
_schpoints(points, ::Val{1}) = collect_as(FixedSizeVectorDefault, (Point2D(x, -one(x)) for x in points))
_schpoints(points, ::Val{2}) = collect_as(FixedSizeVectorDefault, (Point2D(one(x), x) for x in points))
_schpoints(points, ::Val{3}) = collect_as(FixedSizeVectorDefault, (Point2D(x, x) for x in points))

_jac(::Val{N}) where {N} = 2
_jac(::Val{3}) = 2sqrt(2)

# constantcollapser(t::Tuple) = aff::AffineToRef -> _topolytensor(t,aff)


# _topolytensor(a::Function,aff::AffineToRef) = a∘aff
# _topolytensor(t::NTuple{1},aff::AffineToRef) = _topolytensor(t[1],aff,trial,test)
# _topolytensor(t::Tuple,aff::AffineToRef) = x->t[2](_topolytensor(t[1],aff)(x),_topolytensor(t[3:end],aff))(x)


# _eval_operation(z::Any, _, _) = z
# _eval_operation(sf::ShapeFunction{Trial}, trial, _) = operator(sf)(trial)
# _eval_operation(sf::ShapeFunction{Test}, _, test) = operator(sf)(test)
# _eval_operation(z::Any, _) = z
# _eval_operation(sf::ShapeFunction, test) = operator(sf)(test)

# function _eval_operation(o::Operation{ConstantCoeff}, test)
#     efectiveop = DICT_OP[o.op]
#     return efectiveop(_eval_operation.(o.parts, test)...)
# end
# function _eval_operation(o::Operation{ConstantCoeff}, trial, test)
#     efectiveop = DICT_OP[o.op]
#     return efectiveop(_eval_operation.(o.parts, trial, test)...)
# end
