using RayTracingEngine
using FileIO
using LinearAlgebra
using BenchmarkTools

function mainAero2()
    model = load("C:/fc/data/3Dmodels/dragonLikeCapsule.obj")
    dirObs = convert.(Float32, [0; 0; 1.0])
    Nrays = 100*100
    @btime rayTracingSurface($model, $dirObs; Nrays=$Nrays)

    rayTracingSurface(model, dirObs; Nrays=Nrays)
end

mainAero2()
# 12.775387
