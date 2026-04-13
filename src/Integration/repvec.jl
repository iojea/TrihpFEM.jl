"""

    RepVec{N,I,F}
a struct for compressing vectors with many repetitions. 
"""
struct RepVec{L,N,I<:Integer,F<:Number}
    indices::NTuple{N,I}
    values::NTuple{N,F}
    function RepVec(ind,val)
        N = length(ind)
        N == length(val) || throw(ArgumentError("Indices and values must be of equal length"))
        L = last(ind)-1
        I = eltype(ind)
        F = eltype(val)
        new{L,N,I,F}(ind,val)
    end
end

function RepVec(v::T) where {F,T<:AbstractVector{F}}
    vals = tuple(unique(v)...)
    ind = zeros(Int,length(vals));
    i = 1
    j = 1
    while i<length(vals)
        if vals[i] == v[j]
            j = j+1
        else
            ind[i] = j
            i = i+1
        end
    end
    ind[i] = length(v)+1
    RepVec(tuple(ind...),vals)
end


Base.length(rv::RepVec{L,N,I,F}) where {L,N,I,F} = L

Base.iterate(rv::RepVec) = (first(rv.values),(1,1))
function Base.iterate(rv::RepVec{L,N,I,F},st) where {L,N,I,F}
    i,j = st
    j = i<rv.indices[j]-1 ? j : j+1
    j == N+1 && return nothing
    return rv.values[j],(i+1,j)
end


function Base.getindex(rv::RepVec,i)
    1<=i<= length(rv) || throw(BoundsError(rv,i))
    if i< rv.indices[1]
       return rv.values[1]
    else
       k = 1
       while rv.indices[k]<=i
           k +=1
       end
       return rv.values[k]
    end
end

Base.lastindex(rv::RepVec) = last(rv.indices)
