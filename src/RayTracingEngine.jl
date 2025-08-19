module RayTracingEngine

const EPS_PARALLEL = 1e-8
const ZERO_WITH_TOL = -1e-8
const ONE_WITH_TOL = 1.0 + 1e-8
const EPS_MIN_DIST = 1e-4
const RAY_SHIFT = 1e-2
const DEFAULT_NRAYS = 1_000_000

include("math.jl")
include("rtCore.jl")

end
