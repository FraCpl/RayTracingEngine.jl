module RayTracingEngine

using LinearAlgebra
using GeometryBasics

const EPS_PARALLEL = 1e-8
const ZERO_WITH_TOL = -1e-8
const ONE_WITH_TOL = 1.0 + 1e-8
const EPS_MIN_DIST = 1e-4
const RAY_SHIFT = 1e-2
const DEFAULT_NRAYS = 1_000_000

include("math.jl")

export Ray, BBox
include("rtCore.jl")

export buildBvh, traverse!
include("bvh.jl")

export rayTracingDrag, rayTracingSrp, rayTracingSurface, rayTracingHypersonicAero, rayTracingAltimeter
include("rtFunctions.jl")

end
