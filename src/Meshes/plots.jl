# Makie.convert_arguments(::Type{<:AbstractPlot}, x::Point2D) = (Point2f.(x),)


function degplot end

@recipe DegPlot (m,) begin
    linewidth = 0.5
    annotate = false
    title = ""
    colorscheme = :coolwarm
    maxdegree = nothing      # nothing => use the mesh's actual maximum degree
end

function Makie.plot!(p::DegPlot)
    msh = p[1][]
    (; points, edgelist) = msh
    nedge = length(edgelist)
    nedge == 0 && return p
    md = something(p[:maxdegree][], maximum(degree(edgelist[e]) for e in keys(edgelist)))
    cmap = to_colormap(p[:colorscheme][])        # Vector{Colorant} OR Colormap
    # sample at normalized position t in [0,1]; handle both vector and Colormap returns
    sample_color(c, t) = c isa AbstractVector ?
        c[clamp(round(Int, t * (length(c) - 1)) + 1, 1, length(c))] : c[t]
    C = typeof(sample_color(cmap, 0.5))           # concrete color type (e.g. RGBA{Float32})
    segpts = Vector{Point2D{Float32}}(undef, 2 * nedge)
    segcols = Vector{C}(undef, nedge)             # per-SEGMENT (Makie colors flat pts per-segment)
    t(d) = (clamp(d, 1, md) - 1) / max(md - 1, 1) # normalized position -> full hue span
    for (k, e) in enumerate(keys(edgelist))
        segpts[2k - 1] = points[e[1]]
        segpts[2k]     = points[e[2]]
        segcols[k]     = sample_color(cmap, t(degree(edgelist[e])))
    end
    linesegments!(p, segpts; color = segcols, linewidth = p[:linewidth][])
    return p
end



"""
    plothpmesh(m::HPMesh)

"""
function plothpmesh end

@recipe PlotHPMesh (mesh,) begin
    linewidth = 0.5
    annotate = false
    title = ""
end

function Makie.plot!(p::PlotHPMesh)
    msh = p[1][]
    (; points, trilist, edgelist) = msh
    T = eltype(points)

    # ---- Triangle faces: one poly! call, closed polygons + per-face color ----
    ntri = length(trilist)
    polys = Vector{Vector{T}}(undef, ntri)
    facecolors = Vector{Symbol}(undef, ntri)
    for (k, t) in enumerate(keys(trilist))
        a, b, c = points[t[1]], points[t[2]], points[t[3]]
        polys[k] = [a, b, c, a]  # close the loop
        attr = trilist[t]
        facecolors[k] = COLOR_DICT[attr.refine[]]
    end
    poly!(p, polys; color = facecolors, strokecolor = :white, strokewidth = 0.0)

    # ---- All edges: one linegments! call, tag-colored (marked = black) ----
    nedge = length(edgelist)
    segpts = Vector{T}(undef, 2 * nedge)
    segcols = Vector{Symbol}(undef, nedge)
    marked_idx = Int[]
    for (k, e) in enumerate(keys(edgelist))
        ea = edgelist[e]
        segpts[2k - 1] = points[e[1]]
        segpts[2k]     = points[e[2]]
        segcols[k] = ismarked(ea) ? _MARKED_EDGE_COLOR : _EDGE_TAG_COLORS[tag(ea)]
        ismarked(ea) && push!(marked_idx, k)
    end
    linesegments!(p, segpts; color = segcols, linewidth = p[:linewidth][])

    # ---- Marked edges overlay: thicker black pass on top (restores emphasis) ----
    if !isempty(marked_idx)
        mseg = Vector{T}(undef, 2 * length(marked_idx))
        for (j, k) in enumerate(marked_idx)
            mseg[2j - 1] = segpts[2k - 1]
            mseg[2j]     = segpts[2k]
        end
        linesegments!(p, mseg; color = _MARKED_EDGE_COLOR,
                      linewidth = 2 * p[:linewidth][])
    end

    # ---- Optional vertex annotations: one text! call ----
    if p[:annotate][]
        text!(p, points; text = string.(eachindex(points)))
    end

    return p
end


# @recipe(PlotSolHP, mesh,u) do scene
#     Attributes(
#                title = "",
#             )
# end

# """
#     plot!(s::HPSolution)
# plots `s` using `Makie`.
# """
# function Makie.plot!(p::PlotSolHP)
#     (;mesh,u) = p
#     lift(p[1]) do mesh
#         (;points,trilist,edgelist) = mesh
#         lift(p[2]) do u
#             minu = minimum(u)
#             maxu = maximum(u)
#             cr = 0:1;#range(start=minu,stop=maxu,length=3length(u))
#             tris = hcat([[t...] for t in triangles(trilist)]...)'
#             poly!(p,points,tris,color=(u .-minu)/(maxu-minu),colormap=:coolwarm,colorrange=cr)
#         end
#     end
#     return p
# end

# """
#     plot_degs(mesh::HPMesh)
# plots `mesh` coloring the degrees of the edges, using `Makie`.
# """
# function plot_degs(mesh::HPMesh)
#     (;points,edgelist) = mesh
#     degs = degree.(edgelist)
#     nc  = Int(maximum(degs))
#     mc  = Int(minimum(degs))
#     pal = palette(:blues,nc-mc+1)
#     f = Figure()
#     Axis(f[1,1])
#     for e in edgelist
#         x = getindex.(points[e],1)
#         y = getindex.(points[e],2)
#         lines!(x,y,overdraw=true,linewidth=1,color=pal[degree(e)])
#     end
#     f
# end

# """
#     animate_refinement(meshes,path)
# Creates an animation, stored in `path` from a list of meshes.
# """
# function animate_refinement(meshes,path)
#     k = Observable(1)
#     msh = @lift(meshes[$k])
#     fig = plotmeshhp(msh,linewidth=0.25)
#     record(fig,path,1:length(meshes);framerate=2) do t
#         k[] = t
#     end
# end
