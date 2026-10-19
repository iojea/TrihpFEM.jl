struct CMatrix{T<:Real,P<:Integer}
    dict::Dict{NTuple{3,P},FixedSizeMatrixDefault{T}}
    CMatrix{T,P}() where {T<:Real,P<:Integer} = new{T,P}(Dict{NTuple{3,P},FixedSizeMatrixDefault{T}}())
end

function _matrixF!(F::FixedSizeMatrixDefault{T},p::NTuple{3,P},sb::StandardBasis) where {T<:Real,P<:Integer}
    seg₁   = range(start=-one(T),stop=one(T),length=p[1]+1)
    seg₂   = range(start=-one(T),stop=one(T),length=p[2]+1)
    seg₃   = range(start=one(T),stop=-one(T),length=p[3]+1)
    for (j,b) in enumerate(sb)
        i = 1
        @inbounds for x in seg₁[1:end-1]
            F[i,j] = b(Point2D{T}(x,-one(T)))
            i += 1
        end
        @inbounds for y in seg₂[1:end-1]
            F[i,j] = b(Point2D{T}(one(T),y))
            i += 1
        end
        @inbounds for z in seg₃[1:end-1]
            F[i,j] = b(Point2D{T}(z,z))
            i += 1
        end
    end
end

function (cm::CMatrix{T,P})(p::NTuple{3,P}) where {P<:Integer,T<:Real}
    p in keys(cm.dict) && return cm.dict[p]
    sb = StandardBasis(p)
    n = length(sb)
    nₙ = sum(p)
    F = FixedSizeMatrixDefault{T}(undef,nₙ,n)
    _matrixF!(F,p,sb)
    U,Σ,V = svd!(F,full=true)
    C = FixedSizeMatrixDefault{T}(undef,n,n)
    D = Diagonal(1 ./Σ)
    lmul!(D,U')
    mul!(view(C,:,1:nₙ),view(V,:,1:nₙ),U')
    C[:,(nₙ+1:n)] .= V[:,(nₙ+1):n]
    cm.dict[p] = C
    return C
end



