using RayTracingEngine
using FileIO
using JLD2
using GLMakie

function plotBBox!(ax, b)
    # Extract corners
    x1, y1, z1 = b.minx, b.miny, b.minz
    x2, y2, z2 = b.maxx, b.maxy, b.maxz

    segs = Point3f[
        (x1, y1, z1),
        (x2, y1, z1),
        (x2, y1, z1),
        (x2, y2, z1),
        (x2, y2, z1),
        (x1, y2, z1),
        (x1, y2, z1),
        (x1, y1, z1),
        (x1, y1, z2),
        (x2, y1, z2),
        (x2, y1, z2),
        (x2, y2, z2),
        (x2, y2, z2),
        (x1, y2, z2),
        (x1, y2, z2),
        (x1, y1, z2),
        (x1, y1, z1),
        (x1, y1, z2),
        (x2, y1, z1),
        (x2, y1, z2),
        (x2, y2, z1),
        (x2, y2, z2),
        (x1, y2, z1),
        (x1, y2, z2),
    ]

    linesegments!(ax, segs)
    return nothing
end

function plotBvh()
    model = load("C:/Users/francesco capolupo/OneDrive - ESA/Desktop/Work/Tools/rayTracer/examples/SL2_final_adj_5mpp_surf_crop.obj")
    bvh = buildBvh(model; maxLeafSize=500)
    @show length(bvh.nodes)

    fig = Figure();
    display(fig)
    ax = GLMakie.Axis3(fig[1, 1])
    for i in eachindex(bvh.nodes)
        if bvh.nodes[i].count > 0
            plotBBox!(ax, bvh.nodes[i].bbox)
        end
    end
    mesh!(ax, model; color=(:grey, 0.4))
    return nothing
end
plotBvh()
