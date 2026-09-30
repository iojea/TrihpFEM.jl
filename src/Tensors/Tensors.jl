module Tensors

using LinearAlgebra
using StaticArrays


"""

    Tensor{N,D,T}(x...)

builds a tensor of order `N` in dimension `D` with elements given by `x...` (converted to type `T`). The alternative constructor

    Tensor{N,D}(x...)

infers `T` from `promote(x...)`. Notice that `x` should satisfy `length(x)==D^N`. The result is of type:

    Tensor{N,D,T,L,S<:Tuple} <: StaticArray{S,T,N}
    
which represents a simple and brief implementation of tensors, relying on `StaticArrays` for Linear Algebra fast operations. 
"""
struct Tensor{N,D,T,L,S<:Tuple} <: StaticArray{S,T,N} 
    data::NTuple{L,T}
    function Tensor{N,D,T,L,S}(x::NTuple{L,T}) where {N,D,T,L,S<:Tuple}
        StaticArrays.check_array_parameters(S, T, Val{N}, Val{L})
        new{N,D,T,L,S}(x)
    end

    function Tensor{N,D,T,L,S}(x::NTuple{L,Any}) where {N,D,T,L,S<:Tuple}
        StaticArrays.check_array_parameters(S, T, Val{N}, Val{L})
        new{N,D,T,L,S}(StaticArrays.convert_ntuple(T, x))
    end
end

function Tensor{N,D,T}(x...) where {N,D,T}
    L = _compute_L(Val(N),Val(D))
    S = _compute_S(Val(N),Val(D))
    Tensor{N,D,T,L,S}(T.(x))
end

function Tensor{N,D}(x...) where {N,D}
    L = _compute_L(Val(N),Val(D))
    S = _compute_S(Val(N),Val(D))
    x = promote(x...)
    T = eltype(x)
    Tensor{N,D,T,L,S}(x)
end

function Tensor{N,D}(f::Function) where {N,D}
    Tensor{N,D}(map(f,map(idx2cartesian(Val(N),Val(D)),1:_compute_L(Val(N),Val(D))))...)
    #Tensor{N,D}(map(f,map(cartesian2idx(Val(N),Val(D)),1:_compute_L(Val(N),Val(D))))...)
end

function Tensor(a::Array{T,N}) where {T,N}
    D = length(a)%2==0 ? 2 : 3
    M = (N==2 && length(a)==D) ? 1 : N
    Tensor{M,D,T}(a...)
end

function Tensor(a::StaticArray{S,T,N}) where {S<:Tuple,T,N}
    allequal(size(a)) || error(ArgumentError("Dimensions must have equal length for constructing a `Tensor`."))
    D = length(a)%2==0 ? 2 : 3
    M = (N==2 && length(a)==D) ? 1 : N
    Tensor{M,D,T}(a...)
end

function Tensor(a::LinearAlgebra.Adjoint{T}) where {T}
    z = parent(a)
    D = length(z)%2==0 ? 2 : 3
    N = length(size(z))
    M = (N==2 && length(a)==D) ? 1 : N
    Tensor{M,D,T}(z...)'
end
    

_idx2first(i,::Val{D}) where D = mod1(i,D)
_idx2last(i,::Val{D}) where D = 1+(i-1)÷D
function idx2cartesian(i::Number,::Val{2},vd::Val{D}) where {D}
    (_idx2first(i,vd),_idx2last(i,vd))
end
function idx2cartesian(i::Number,::Val{N},vd::Val{D}) where {N,D}
    (_idx2first(i,vd),idx2cartesian(_idx2last(i,vd),Val(N-1),vd)...)
end
idx2cartesian(i::Number,::Val{1},::Val{D}) where D = i
idx2cartesian(::Val{N},::Val{D}) where {N,D} = Base.Fix2(Base.Fix2(idx2cartesian,Val(N)),Val(D))

_compute_S(::Val{N},::Val{D}) where {N,D}= Tuple{(D for _ in 1:N)...}
_compute_L(::Val{N},::Val{D}) where {N,D}= D^N


#cartesian2idx(v,::Val{N},::Val{D}) where {N,D} = 1+mapreduce(i->(v[i]-1)*D^(i-1),+,1:N)
#cartesian2idx(::Val{N},::Val{D}) where {N,D} = Base.Fix{2}(Base.Fix{3}(cartesian2idx,Val(D)),Val(N))
    
Base.getindex(t::Tensor,i::Int64) = t.data[i]

Base.zero(::Tensor{N,D,T}) where {N,D,T} = Tensor{N,D}((zero(T) for _ in 1:D^N)...)

### Adjust type inference for the result of operations (+,*) between Tensors isa Tensor
StaticArrays.similar_type(::Type{A},::Type{F},s::Size{E}) where {N,D,T,L,S<:Tuple,A<:Tensor{N,D,T,L,S},F,E} =
    Tensor{N,D,F,L,S}

    
#     StaticArrays.default_similar_type(T,s,StaticArrays.length_val(s))
# function StaticArrays.default_similar_type(::Type{T},s::Size{S},::Type{Val{E}}) where {N,D,F,T<:Tensor{N,D,F},S,E}
#     Tensor{E,D,F,prod(s),Tuple{S...}}
# end

# similar_type(::Type{A},::Type{T},s::Size{S}) where {A<:AbstractArray,T,S} = default_similar_type(T,s,length_val(s))
# default_similar_type(::Type{T}, s::Size{S}, ::Type{Val{D}}) where {T,S,D} = SArray{Tuple{S...},T,D,prod(s)}

## As a default method (fallback), Tensors of numbers evaluate to themselves
(t::Tensor{N,D,T})(_) where {N,D,T<:Number} = t 

# No-op for adjoint
Base.adjoint(t::Tensor{1,D}) where D = t


"""

   otimes(t₁::Tensor{N,D,T},t₂::Tensor{M,D,T}) where {N,M,D,T}

computes the outer product of `t₁` and `t₂`, obtaining as result a `Tensor{N+M,D,T}`. The binary operator `⊗` can be used instead of the function
 
"""
function otimes(t₁::Tensor{N,D,T},t₂::Tensor{M,D,R}) where {N,M,D,T,R}
    Tensor{N+M,D}(tuple((t₁[i]*t₂[j] for i in eachindex(t₁), j in eachindex(t₂))...)...)
end
otimes(n::Number,t::Tensor{N,D,T}) where {N,D,T} = n*t
otimes(t::Tensor{N,D,T},n::Number) where {N,D,T} = n*t
const ⊗ = otimes
export Tensor, otimes, ⊗

end; #module



# _idx2first(i,::Val{D}) where D = mod1(i,D)
# _idx2last(i,::Val{D}) where D = 1+(i-1)÷D
# function idx2cartesian(i::Number,::Val{2},vd::Val{D}) where {D}
#     (_idx2first(i,vd),_idx2last(i,vd))
# end
# function idx2cartesian(i::Number,::Val{N},vd::Val{D}) where {N,D}
#     (_idx2first(i,vd),idx2cartesian(_idx2last(i,vd),Val(N-1),vd)...)
# end

# idx2cartesian(i::Number,::Val{1},::Val{D}) where D = i


# Base.length(::Tensor{N,D,T,M}) where {N,D,T,M} = M 
# Base.eltype(::Tensor{N,D,T}) where {N,D,T} = T
# Base.ndims(::Tensor{N}) where {N} = N
# Base.size(::Tensor{N,D,T,M}) where {N,D,T,M} = NTuple{N,Int}(D for _ in 1:N)
# Base.getindex(t::Tensor,i::Integer) = t.data[i]
# Base.getindex(t::Tensor{N,D},c::CartesianIndex) where {N,D} = t.data[cartesian2idx(c,Val(N),Val(D))]
# Base.getindex(t::Tensor{N,D,T,M},i,j...) where {N,D,T,M} = t.data[cartesian2idx((i,j...),Val(N),Val(D))]
# Base.eachindex(t::Tensor) = Base.OneTo(length(t))


# LinearAlgebra.dot(t₁::Tensor{N,D,T},t₂::Tensor{N,D,T}) where {N,D,T}= mapreduce(*,+,t₁,t₂)

# function LinearAlgebra.dot(t₁::Tensor{N,D},t₂::Tensor{M,D}) where {N,M,D}
#     N<M && return _reduce_left(t₁,t₂)
#     M>N && return _reduce_right(t₁,t₂)
# end 

# function _reduce_left(t₁::Tensor{N,D},t₂::Tensor{M,D}) where {N,M,D}
#     Tensor{M-N,D}()
# end


# Base.:*(a::Number,t₁::Tensor{N,D,T,M}) where {N,D,T,M} = Tensor{N,D}((a*t for t in t₁)...)

