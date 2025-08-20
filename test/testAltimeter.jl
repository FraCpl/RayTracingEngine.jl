using RayTracingEngine
using FileIO
using LinearAlgebra
using BenchmarkTools

function mainAlt()
    model = load("C:/Users/francesco capolupo/OneDrive - ESA/Desktop/Work/Tools/rayTracer/examples/SL2_final_adj_5mpp_surf_crop.obj")
    bbox = BBox(model)
    argonautLat = -88.35006
    argonauLon = 300.14565
    pos3 = Float32.([cosd(argonauLon)*cosd(argonautLat); sind(argonauLon)*cosd(argonautLat); sind(argonautLat)].*(1.7374e6+300.0))
    ray = Ray(pos3, Float32[0.0; 0.0; 1.0])
    @time rayTracingAltimeter(model, ray, bbox)
    @show ray.t
    # @time rayTracingAltimeter(model, ray, bbox)
    @btime rayTracingAltimeter($model, $ray, $bbox)
    # @profview rayTracingAltimeter(model, ray, bbox)
end

mainAlt()
# 12.775387
