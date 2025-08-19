using RayTracingEngine
using FileIO
using LinearAlgebra
using BenchmarkTools

function mainAero2()
    model = load("C:/fc/data/3Dmodels/dragonLikeCapsule.obj")
    bbox = BBox(model)
    dirObs = convert.(Float32, [0; 0; 1.0])
    Nrays = 100*100
    # @btime rayTracingSurface($model, $dirObs; Nrays=$Nrays)

    # rayTracingSurface(model, dirObs; Nrays=Nrays)
    # @btime rayTracingSurface($model, $dirObs; Nrays=$Nrays)
    # @profview rayTracingSurface(model, dirObs; Nrays=Nrays)

    ray = Ray(Float32[-2.0; 0.0; 50.0], Float32[0.0; 0.0; -1.0])
    tmp1 = zero(ray.dir); tmp2 = zero(ray.dir); tmp3 = zero(ray.dir); tmp4 = zero(ray.dir);
    rayTracingAltimeter(model, ray, bbox, tmp1, tmp2, tmp3, tmp4)
    @show ray.t
    @time rayTracingAltimeter(model, ray, bbox, tmp1, tmp2, tmp3, tmp4)
    @btime rayTracingAltimeter($model, $ray, $bbox, $tmp1, $tmp2, $tmp3, $tmp4)
end

mainAero2()
# 12.775387
