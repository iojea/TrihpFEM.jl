"""
see Edge, Triangle.

    SetTuple{L,I<:Integer} <: StaticArray{Tuple{L,I,1}}
    
abstract type for storing an immutable fixed size array. An `SetTuple` can be used for indexing, slicing and Linear Algebra operations. However, it behaves like a mathematical set for the purposes of hashing and equality checking via `isequal`. The goal of this implementation is to use `SetTuple`s as keys of a `Dictionary` allowing search with permutations of the `SetTuple`.

Concrete subtypes of `SetTuple` should implement the function `data` to retrieve the underlying data. 
"""

abstract type SetTuple{L,I<:Integer} <: StaticArray{Tuple{L},I,1} end

function data(v::SetTuple) end

#Indexing
Base.getindex(v::SetTuple, i::Int) = data(v)[i]


function Base.hash(s::SetTuple, h::UInt)
    hv = HASH_SEED
    for x in data(s)
        hv ⊻= hash(x)
    end
    hash(hash(hv, h),hash(typeof(s)))
end

Base.isequal(t1::T,t2::T) where T<:SetTuple = length(t1)==length(t2) && issubset(t1,t2)




