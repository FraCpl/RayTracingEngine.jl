mutable struct Ray{T}
    origin::Vector{T}
    dir::Vector{T}
    t::T
    idxFace::Int64
end

Ray(origin, dir) = Ray(origin, normalize(dir), Inf, 0)

@inline function resetRay!(ray::Ray)
    ray.t = Inf
    ray.idxFace = 0
    return
end

@inline function rayHitPosition!(posHit, ray::Ray)
    posHit[1] = ray.origin[1] + ray.t*ray.dir[1]
    posHit[2] = ray.origin[2] + ray.t*ray.dir[2]
    posHit[3] = ray.origin[3] + ray.t*ray.dir[3]
    return
end

@inline function intersect!(ray::Ray, tri, idx::Int64, E12, E13, P, T)
    v1, v2, v3 = tri
    @inbounds for i in 1:3
        E12[i] = v2[i] - v1[i]
        E13[i] = v3[i] - v1[i]
    end
    cross!(P, ray.dir, E13)
    detM = dot3(P, E12)
    if abs(detM) ≤ EPS_PARALLEL
        return
    end

    @inbounds for i in 1:3
        T[i] = ray.origin[i] - v1[i]
    end
    u = dot3(P, T)/detM
    if u < ZERO_WITH_TOL || u > ONE_WITH_TOL
        return
    end

    cross!(P, T, E12)       # Q <- P
    v = dot3(P, ray.dir)/detM
    if v < ZERO_WITH_TOL || u + v > ONE_WITH_TOL
        return
    end

    t = dot3(E13, P)/detM
    if t > EPS_MIN_DIST && t < ray.t
        ray.t = t
        ray.idxFace = idx
    end
end

struct BBox{T}
    boxMin::Vector{T}
    boxMax::Vector{T}
end

BBox(model::GeometryBasics.Mesh) = BBox(model.position)

function BBox(vertices)
    boxMin, boxMax = vertices[1], vertices[1]
    for v in vertices
        boxMin = min.(boxMin, v)
        boxMax = max.(boxMax, v)
    end
    return BBox(Vector(boxMin), Vector(boxMax))
end

@inline function getT1T2(ray::Ray, bbox::BBox, idx)
    t1 = (bbox.boxMin[idx] - ray.origin[idx])/ray.dir[idx]
    t2 = (bbox.boxMax[idx] - ray.origin[idx])/ray.dir[idx]
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
