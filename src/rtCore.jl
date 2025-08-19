mutable struct Ray{T}
    origin::Vector{T}
    dir::Vector{T}
    t::T
    idxFace::Int64
end

function Ray(origin::Vector{T}, dir::Vector{T}) where T
    return Ray(origin, normalize(dir), convert(T, Inf), 0)
end

@inline function resetRay!(ray::Ray)
    ray.t = oftype(ray.t, Inf)
    ray.idxFace = 0
    return
end

@inline function rayHitPosition!(posHit, ray::Ray)
    posHit[1] = ray.origin[1] + ray.t*ray.dir[1]
    posHit[2] = ray.origin[2] + ray.t*ray.dir[2]
    posHit[3] = ray.origin[3] + ray.t*ray.dir[3]
    return
end

@inline function intersect!(ray::Ray, tri, idx::Int64, E12::Vector{X}, E13::Vector{X}, P::Vector{X}, T::Vector{X}) where X
    v1, v2, v3 = tri
    @inbounds for i in 1:3
        E12[i] = v2[i] - v1[i]
        E13[i] = v3[i] - v1[i]
    end
    cross!(P, ray.dir, E13)
    detM = dot3(P, E12)
    if abs(detM) ≤ X(EPS_PARALLEL)
        return
    end

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
    bmin, bmax = vertices[1], vertices[1]
    for v in vertices
        bmin = min.(bmin, v)
        bmax = max.(bmax, v)
    end
    return BBox(Vector(bmin), Vector(bmax))
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
