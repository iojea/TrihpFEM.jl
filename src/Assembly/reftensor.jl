"""
   RefTensor{D}
   A memoisation struct for computing the local tensor corresponding to a certain basis.
"""
struct RefTensor{K, T,R<:Union{Nothing,Quadrature}}
    t::Tuple
    sch::R
    dict::Dict{K, T}
end

function RefTensor(t::Term{N, ConstantCoeff, T, Measure{F,I,P}}) where {N,T<: Tuple,F,I,P}
    (; integrand, measure) = t
    (; operation) = integrand
    D = dim(measure)
    op = _scalars_out(operation,F,Val(D))
    globpart = global_part(op)
    sfs = get_shape_functions(operation)
    length(sfs) == N || throw(ArgumentError("Malformed term: a `Term{$N}` should be formed by $N `ShapeFunction`s, but $(length(sfs)) are present."))
    bs = ((first(StandardBasis(zeros(P,2D-1)...)) for _ in sfs)...,)
    key = degs.(bs)
    poly_tensor = _totensor(globpart, bs...)
    tensor = ref_integrate.(poly_tensor)
    U = typeof(tensor)
    K = typeof(key)
    d = Dict{K, U}(key => tensor)
    return RefTensor{K, U,Nothing}(globpart,nothing, d)
end

function RefTensor(t::Term{N, VariableCoeff, T, Measure{F,I,P}}) where {N, T <: Tuple,F,I,P}
    (; integrand, measure) = t
    (; operation) = integrand
    (;mesh,sch) = measure
    D = dim(measure)
    B = basistype(mesh)
    op = _scalars_out(operation,F,Val(D))
    globpart = global_part(op)
    sfs = get_shape_functions(operation)
    length(sfs) == N || throw(ArgumentError("Malformed term: a `Term{$N}` should be formed by $N `ShapeFunction`s, but $(length(sfs)) are present."))
    bs = ((first(B(zeros(P,2D-1)...)) for _ in sfs)...,)
    key = degs.(bs)
    poly_tensor = _totensor(globpart, bs...)
    # IMPORTANT: The factor 2 comes from integration on the reference triangle, with area 2.
    tensor = collect_as(FixedSizeArrayDefault,(2w*poly_tensor(x) for (x,w) in zip(sch.points,sch.weights)))
    U = typeof(tensor)
    K = typeof(key)
    d = Dict{K, U}(key => tensor)
    return RefTensor{K, U, typeof(sch)}(globpart, sch,d)
end




function (lt::RefTensor{K,T,Nothing})(bs...) where {K,T}
    key = degs.(bs)
    haskey(lt.dict, key) && return lt.dict[key]
    poly_tensor = _totensor(lt.t, bs...)
    tensor = ref_integrate.(poly_tensor)
    lt.dict[key] = tensor
    return tensor
end
function (lt::RefTensor{K,T,R})(bs...) where {K,T,R<:Quadrature}
    key = degs.(bs)
    (;points,weights) = lt.sch
    haskey(lt.dict, key) && return lt.dict[key]
    poly_tensor = _totensor(lt.t, bs...)
    # IMPORTANT: The factor 2 comes from integration on the reference triangle, with area 2.
    tensor = collect_as(FixedSizeVectorDefault,2w*poly_tensor(x) for (x,w) in zip(points,weights))
    lt.dict[key] = tensor
    return tensor
end


"""
   _isglobal(a)
Decides if a factor `a` in an `Operation` belongs to its global part (`ShapeFunction`s and constants).
"""
_isglobal(::Any) = true
_isglobal(::Function) = false
_isglobal(::DiffOperator) = false



"""
   _scalars_out(op::Operation,::Type{F},::Val{D})
takes `Operation` `op` and moves scalars (numbers of scalar functions) to the left, returning a reorderd `Operation`. This is done because it is the easiest way to ensure that scalars operate against other factors with a simple product (*) and not an outer product (⊗). The product itself will always reduce to a simple product (because ⊗ applied to scalars is just *), but other adjacent products will produce unnessesary large tensors. E.g.: if `B` is a matrix from an `AffineToRef` and `z isa Number` the product `∇(u)'*B'*z*B*∇(v)` would be  converted into a global term `∇(u)'⊗z⊗∇(v)` and a local term `B'⊗B`. However, is the scalar `z` is extracted, the new `Operation` is `z*∇(u)'*B'*B*∇(v)` which results in the global term `z*∇u'⊗∇(v)` and the local term `B'*B`.
"""
function _scalars_out(op::Operation{C}, ::Type{F}, ::Val{D}) where {C <: CoeffType, F, D}
    for i in 3:2:length(op)
        if _isscalar(op[i], F, Val(D))
            op = Operation{C}(op[i], *, op[1:(i - 2)]..., op[(i + 1):end]...)
        end
    end
    return op
end
_isscalar(::Any, _, _) = false
_isscalar(::Number, _, _) = true
function _isscalar(f::Function, ::Type{F}, ::Val{D}) where {F, D}
    return f(zeros(F, D)) isa Number
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
    _totensor(t::Tuple,trial,test)
    _totensor(t::Tuple,test)
    _totensor(t::Tuple,aff::AffineToRef)
Evaluates `t` recursively to produce a `Tensor`. `t` must be the global or local part of an `Operation`.

    _totensor(t::Tuple,trial,test,sch::Quadrature)
    _totensor(t::Tuple,test,sch::Quadrature)
    _totensor(t::Tuple,)
"""
_totensor(a::Any,_,_) = a
_totensor(a::Any,_) = a
_totensor(sf::ShapeFunction{Trial}, trial, _) = operator(sf)(trial)
_totensor(sf::ShapeFunction{Test}, _, test) = operator(sf)(test)
_totensor(t::NTuple{1}, trial, test) = _totensor(t[1], trial, test)
_totensor(t::Tuple, trial, test) = t[2](_totensor(t[1], trial, test), _totensor(t[3:end], trial, test))
_totensor(sf::ShapeFunction{Test},test) = operator(sf)(test)
_totensor(t::NTuple{1},test) = _totensor(t[1],test)
_totensor(t::Tuple,test) = t[2](_totensor(t[1],test),_totensor(t[3:end],test))




# constantcollapser(t::Tuple) = aff::AffineToRef -> _totensor(t,aff)


# _totensor(a::Function,aff::AffineToRef) = a∘aff
# _totensor(t::NTuple{1},aff::AffineToRef) = _totensor(t[1],aff,trial,test)
# _totensor(t::Tuple,aff::AffineToRef) = x->t[2](_totensor(t[1],aff)(x),_totensor(t[3:end],aff))(x)


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


