function Base.show(io::IO, p::P) where {F, P <: ProductPoly{F}}
    return if any(pp == zero(pp) for pp in p.polys)
        print(io, "$(zero(F))")
    elseif p.polys[1] == one(p.polys[1])
        printpoly(io, p.polys[2])
    elseif p.polys[2] == one(p.polys[2])
        printpoly(io, p.polys[1])
    else
        print(io, "(")
        printpoly(io, p.polys[1])
        print(io, ")(")
        printpoly(io, p.polys[2])
        print(io, ")")
    end
end

function Base.show(io::IO, p::P) where {P <: PolySum}
    print(io, p.left)
    print(io, " + ")
    return print(io, p.right)
end
