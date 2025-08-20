mutable struct Ray{T}
    origin::Vector{T}   # Ray origin position
    dir::Vector{T}      # Ray direction
    tmp1::Vector{T}     # Allocations
    tmp2::Vector{T}     # Allocations
    tmp3::Vector{T}     # Allocations
    tmp4::Vector{T}     # Allocations
    t::T                # Distance to hit
    idxFace::Int64      # ID face hit
end

function Ray(origin::Vector{T}, dir::Vector{T}) where T
    return Ray(origin, normalize(dir), similar(dir), similar(dir), similar(dir), similar(dir), T(Inf), 0)
end

@inline function resetRay!(ray::Ray{T}) where T
    ray.t = T(Inf)
    ray.idxFace = 0
    return
end

@inline function rayHitPosition!(posHit::Vector{T}, ray::Ray{T}) where T
    posHit[1] = ray.origin[1] + ray.t*ray.dir[1]
    posHit[2] = ray.origin[2] + ray.t*ray.dir[2]
    posHit[3] = ray.origin[3] + ray.t*ray.dir[3]
    return
end

@inline function intersect!(ray::Ray{X}, tri, idx::Int64) where X
    v1, v2, v3 = tri
    E12 = ray.tmp1; E13 = ray.tmp2; P = ray.tmp3
    @inbounds for i in 1:3
        E12[i] = v2[i] - v1[i]
        E13[i] = v3[i] - v1[i]
    end
    cross!(P, ray.dir, E13)
    detM = dot3(P, E12)
    if abs(detM) ≤ X(EPS_PARALLEL)
        return
    end

    T = ray.tmp4
    @inbounds for i in 1:3
        T[i] = ray.origin[i] - v1[i]
    end
    u = dot3(P, T)/detM
    if u < X(ZERO_WITH_TOL); return; end
    if u > X(ONE_WITH_TOL); return; end

    cross!(P, T, E12)       # Q <- P
    v = dot3(P, ray.dir)/detM
    if v < X(ZERO_WITH_TOL) || u + v > X(ONE_WITH_TOL)
        return
    end

    t = dot3(E13, P)/detM
    if t > X(EPS_MIN_DIST) && t < ray.t
        ray.t = t
        ray.idxFace = idx
    end
end

struct BBox{T}
    min::Vector{T}
    max::Vector{T}
end

BBox(model::GeometryBasics.Mesh) = BBox(model.position)

function BBox(vertices)
    bmin = Vector(vertices[1])
    bmax = Vector(vertices[1])
    for v in vertices
        @inbounds for i in 1:3
            bmin[i] = min(bmin[i], v[i])
            bmax[i] = max(bmax[i], v[i])
        end
    end
    return BBox(bmin, bmax)
end

@inline function getT1T2(ray::Ray, bbox::BBox, idx)
    t1 = (bbox.min[idx] - ray.origin[idx])/ray.dir[idx]
    t2 = (bbox.max[idx] - ray.origin[idx])/ray.dir[idx]
    return t1, t2
end

@inline function intersect!(ray::Ray, bbox::BBox)
    t1, t2 = getT1T2(ray, bbox, 1)
    tmin = min(t1, t2)
    tmax = max(t1, t2)

    t1, t2 = getT1T2(ray, bbox, 2)
    tmin = max(tmin, min(t1, t2))
    tmax = min(tmax, max(t1, t2))

    t1, t2 = getT1T2(ray, bbox, 3)
    tmin = max(0.0, max(tmin, min(t1, t2)))
    tmax = min(tmax, max(t1, t2))
    return tmax ≥ tmin
end
