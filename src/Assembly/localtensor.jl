"""
   LocalTensor{D}
   A memoisation struct for computing the local tensor corresponding to a certain basis.
"""
struct LocalTensor{K,T<:Union{Tensor,Number}}
    term::Term
    dict::Dict{K,T}
end

function LocalTensor(t::Term{N,ConstantCoeff,T,M}) where {N,T<:Tuple,M<:Measure}
    (;integrand,measure) = t
    (;operation) = integrand
    P = degtype(measure.mesh)
    sfs = get_shape_functions(operation,measure.mesh)
    length(sfs) == N || throw(ArgumentError("Malformed term: a `Term{$N}` should be formed by $N `ShapeFunction`s, but $(length(sfs)) are present."))
    bs = ((first(basis(s,ntuple(_->zero(P),2dim(s)-1))) for s in sfs)...,)
    key = degs.(bs)
    poly_tensor = _eval_operation(operation,bs...)
    tensor = ref_integrate.(poly_tensor)
    U = typeof(tensor)
    K = typeof(key)
    d = Dict{K,U}(key=>tensor)
    LocalTensor{K,U}(t,d)
end

function LocalTensor(t::Term{N,VariableCoeff,T,M}) where {N,T<:Tuple,M<:Measure}
    (;integrand,measure) = t
    (;operation) = integrand
    P = degtype(measure.mesh)
    sfs = get_shape_functions(operation,measure.mesh)
    length(sfs) == N || throw(ArgumentError("Malformed term: a `Term{$N}` should be formed by $N `ShapeFunction`s, but $(length(sfs)) are present."))
    bs = ((first(basis(s,ntuple(_->zero(P),2dim(s)-1))) for s in sfs)...,)
    key = degs.(bs)
    poly_tensor = _eval_operation(operation,bs...)
    tensor = ref_integrate.(poly_tensor)
    U = typeof(tensor)
    K = typeof(key)
    d = Dict{K,U}(key=>tensor)
    LocalTensor{K,U}(t,d)
end

get_shape_functions(sf::ShapeFunction,mesh::HPTriangulation) = _adjust_to_mesh(sf,mesh)
get_shape_functions(::Any,::HPTriangulation) = ()
get_shape_functions(t::Tuple,::HPTriangulation) = filter(Base.Fix2(isa,ShapeFunction),t)
function get_shape_functions(operation::Operation,m::HPTriangulation)
    a,b = ((get_shape_functions(p,m) for p in operation.parts)...,)
    a = a isa ShapeFunction ? (a,) : get_shape_functions(a,m) 
    b = b isa ShapeFunction ? (b,) : get_shape_functions(b,m) 
    filter(Base.Fix2(isa,ShapeFunction),(a...,b...))
end


_eval_operation(z::Any,_,_) = z
_eval_operation(sf::ShapeFunction{Trial},trial,_) = operator(sf)(trial)
_eval_operation(sf::ShapeFunction{Test},_,test) = operator(sf)(test)
_eval_operation(z::Any,_) = z
_eval_operation(sf::ShapeFunction,test) = operator(sf)(test)

function _eval_operation(o::Operation{ConstantCoeff},test)
    efectiveop = DICT_OP[o.op]
    efectiveop(_eval_operation.(o.parts,test)...)
end
function _eval_operation(o::Operation{ConstantCoeff},trial,test)
    efectiveop = DICT_OP[o.op]
    efectiveop(_eval_operation.(o.parts,trial,test)...)
end

function _eval_operation(o::Operation{VariableCoeff},test)
    
end



_adjust_to_mesh(sf::ShapeFunction{T,O,N,2},::HPMesh) where {T,O,N} = sf
_adjust_to_mesh(::ShapeFunction{T,O,N,2},::BoundaryHPMesh) where {T,O,N} = ShapeFunction{T,O,N,1}()
_init_deg_tuple(P,::ShapeFunction{T,O,N,D}) where {T,O,N,D} = ntuple(_->zero(P),2D-1)

function (lt::LocalTensor)(bs...)
    key = degs.(bs)
    haskey(lt.dict,key) && return lt.dict[key]
    poly_tensor = _eval_operation(lt.term.integrand.operation,bs...)
    tensor = ref_integrate.(poly_tensor)
    lt.dict[key] = tensor
    return tensor
end
