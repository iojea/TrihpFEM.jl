function Base.show(io::IO, p::P) where {F, P <: BiPoly{F}}
    return if p.px == zero(p.px) || p.py == zero(p.py)
        print(io, "($(zero(F)),)")
    elseif p.px == one(p.px)
        printpoly(io, p.py)
    elseif p.py == one(p.py)
        printpoly(io, p.px)
    else
        print(io, "(")
        printpoly(io, p.px)
        print(io, ")(")
        printpoly(io, p.py)
        print(io, ")")
    end
end

function Base.show(io::IO, p::P) where {P <: PolySum}
    print(io, p.left)
    print(io, " + ")
    return print(io, p.right)
end

function Base.show(io::IO, p::P) where {P <: PolyTensorField}
    return show(io, p.tensor)
end
