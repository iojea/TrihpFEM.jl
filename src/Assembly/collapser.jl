"""
   _islocal(a)
Decides if a factor `a` in an `Operation` belongs to its local part: `Function`s, except products between parts of the `Operation` and `DiffOperator`s. Notice that a pure `DiffOperator` only appears as a factor in an `Operation` as a marker of the `DiffOperator` that should be applied to an `AffineToRef` as a consequence of the chain rule applied to a `ShapeFunction`. Hence, the `DiffOperator` is essentially a place holder for a change of variables that will need to be applied element-wise.
"""
_islocal(::Any) = false
_islocal(::typeof(*)) = false
_islocal(::typeof(⋅)) = false
_islocal(::Function) = true
_islocal(::DiffOperator) = true

"""
   local_part(op::Operation)
Extracts the local part of an `Operation` as a `Tuple`. Products are transformed (*->⊗) when needed. 
"""
function local_part(op::Operation)
    ind = findall(_islocal, op.parts)
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

# _totensor(diffop::DiffOperator,aff::AffineToRef) = diffop(aff)


struct Collapser{T<:Tuple,C<:CoeffType}
    lp::T
    function Collapser(t::Term{N,C}) where {N,C<:CoeffType} 
        (; integrand, measure) = t
        (; operation) = integrand
        B = basistype(measure.mesh)
        F = floattype(measure.mesh)
        op = _scalars_out(operation,F,Val(PolyFields.dim(B)))
        locpart = local_part(op)
        new{typeof(locpart),C}(locpart)
    end
end



(c::Collapser{T,ConstantCoeff})(aff::AffineToRef) where T = _totensor(c.lp,aff)

_totensor(a::Any,::AffineToRef) = a
_totensor(t::NTuple{1},aff::AffineToRef) = _totensor(t[1],aff)
_totensor(diffop::DiffOperator,aff::AffineToRef) = diffop(aff)
_totensor(t::Tuple,aff::AffineToRef) = t[2](_totensor(t[1],aff),_totensor(t[3:end],aff))


_totensor(f::Function,aff::AffineToRef,sch::Quadrature) = (f(aff(x)) for x in sch.points)
_totensor(diffop::DiffOperator,aff::AffineToRef,sch::Quadrature) = diffop(aff)
_totensor(t::Tuple,aff::AffineToRef,sch::Quadrature) = t[2].(_totensor(t[1],aff,sch),_totensor(t[3:end],aff,sch))
_totensor(t::NTuple{1},aff::AffineToRef,sch::Quadrature) = _totensor(t[1],aff,sch)

function (c::Collapser{T,VariableCoeff})(aff::AffineToRef,sch::Quadrature) where {T}
    collect_as(FixedSizeVectorDefault,_totensor(c.lp,aff,sch))
end

