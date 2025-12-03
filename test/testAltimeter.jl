using RayTracingEngine
using FileIO
using LinearAlgebra
using BenchmarkTools
using JLD2

function mainAlt()
    model = load("C:/Users/francesco capolupo/OneDrive - ESA/Desktop/Work/Tools/rayTracer/examples/SL2_final_adj_5mpp_surf_crop.obj")
    bbox = RayTracingEngine.BBox(model)
    argonautLat = -88.35006
    argonautLon = 300.14565
    pos3 = Float32.([
        cosd(argonautLon)*cosd(argonautLat);
        sind(argonautLon)*cosd(argonautLat);
        sind(argonautLat)
    ] .* (1.7374e6+300.0))
    ray = Ray(pos3, Float32[0.0; 0.0; 1.0])
    @time rayTracingAltimeter(model, ray, bbox)
    @show ray.t, ray.idxFace
    # @time rayTracingAltimeter(model, ray, bbox)
    # @btime rayTracingAltimeter($model, $ray, $bbox)
    # @profview rayTracingAltimeter(model, ray, bbox)
end

mainAlt()
# 12.775387

function mainBvh()
    model = load("C:/Users/francesco capolupo/OneDrive - ESA/Desktop/Work/Tools/rayTracer/examples/SL2_final_adj_5mpp_surf_crop.obj")

    @time bvh = buildBvh(model)

    argonautLat = -88.35006
    argonautLon = 300.1456
    pos3 = Float32.([
        cosd(argonautLon)*cosd(argonautLat);
        sind(argonautLon)*cosd(argonautLat);
        sind(argonautLat)
    ] .* (1.7374e6+300.0))
    ray = Ray(pos3, Float32[0.0; 0.0; 1.0])

    @time RayTracingEngine.intersect!(ray, bvh)
    # if ray.idxFace > 0
    #     println("Hit triangle $(ray.idxFace) at t=$(ray.t)")
    # end
    @show ray.t, ray.idxFace

    # FileIO.save("SL2_final_adj_5mpp_surf_crop_bvh.jld2", "bvh", bvh)

end

mainBvh()
