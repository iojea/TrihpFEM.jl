# printcoeffs(io::IO, pj::P, j, mimetype) where P<:AbstractPolynomial = printpoly(io,pj,mimetype)


# function showpolyterm(io::IO, ::Type{P}, pj::T, var, j, first::Bool, mimetype) where {B, T, P<:AbstractUnivariatePolynomial{B,T}}
#     if Polynomials._iszero(pj) return false end

#     pj = Polynomials.printsign(io, pj, first, mimetype)
#     if !first; print(io,"(") end
#     if Polynomials.hasone(T)
#         if !(Polynomials._isone(pj) && !(Polynomials.showone(T) || j == 0))
#             printcoeffs(io, pj, j, mimetype)
#         end
#     else
#         printcoeffs(io, pj, j, mimetype)
#     end

#     if !first; print(io,")") end
#     Polynomials.printproductsign(io, pj, j, mimetype)
#     #printbasis(io, P, j, mimetype)
#     Polynomials.printbasis(io, Polynomials.constructorof(P){T,Symbol(var)}, j, mimetype)
#     return true
# end

# function printpoly2(io::IO, p::P, mimetype=MIME"text/plain"();
#                    descending_powers=false, offset::Int=0, var=indeterminate(p),
#                    compact=false, mulsymbol="*") where {T,P<:AbstractPolynomial{T}}
#     first = true
#     printed_anything = false
#     for i in (descending_powers ? reverse(eachindex(p)) : eachindex(p))
#         ioc = IOContext(io,
#                         :compact=>get(io, :compact, compact),
#                         :multiplication_symbol => get(io, :multiplication_symbol, mulsymbol)
#                         )
#         printed = showpolyterm(ioc, P, p[i], var, i+offset, first, mimetype)
#         first &= !printed
#         printed_anything |= printed
#     end
#     printed_anything || print(io, zero(eltype(T)))
#     return nothing
# end

function Base.show(io::IO, mimetype::MIME"text/plain", p::P) where {F, P <: BiPoly{F}}
    return if p.px == zero(p.px) || p.py == zero(p.py)
        println(io, "($(zero(F)),)")
    elseif p.px == one(p.px)
        printpoly(io, p.py, mimetype)
    elseif p.py == one(p.py)
        printpoly(io, p.px, mimetype)
    else
        print(io, "(")
        printpoly(io, p.px, mimetype)
        print(io, ")(")
        printpoly(io, p.py, mimetype)
        print(io, ")")
    end
end

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

function Base.show(io::IO, mimetype::MIME"text/plain", p::P) where {P <: PolySum}
    show(io, mimetype, p.left)
    print(io, " + ")
    return show(io, mimetype, p.right)
end
function Base.show(io::IO, p::P) where {P <: PolySum}
    print(io, p.left)
    print(io, " + ")
    return print(io, p.right)
end

function Base.show(io::IO, mimetype::MIME"text/plain", p::P) where {P <: PolyTensorField}
    return show(io, mimetype, p.tensor)
end
function Base.show(io::IO, p::P) where {P <: PolyTensorField}
    return show(io, p.tensor)
end
