mutable struct Ray{T}
    origin::Vector{T}   # Ray origin position
    dir::Vector{T}      # Ray direction
    invdir::Vector{T}   # 1/direction
    tmp1::Vector{T}     # Allocations
    tmp2::Vector{T}     # Allocations
    tmp3::Vector{T}     # Allocations
    tmp4::Vector{T}     # Allocations
    t::T                # Distance to hit
    idxFace::Int64      # ID face hit
    idxSkip::Int64      # ID face to skip when checking for intersections
end

function Ray(origin::Vector{T}, dir::Vector{T}) where {T}
    return Ray(origin, normalize(dir), similar(dir), similar(dir), similar(dir), similar(dir), similar(dir), T(Inf), 0, -1)
end

Ray() = Ray(zeros(Float32, 3), zeros(Float32, 3))

@inline function resetRay!(ray::Ray{T}) where {T}
    ray.t = T(Inf)
    ray.idxFace = 0
    ray.idxSkip = -1
    return nothing
end

@inline function invdirRay!(ray::Ray{T}) where {T}
    ray.invdir[1] = 1/ray.dir[1]
    ray.invdir[2] = 1/ray.dir[2]
    ray.invdir[3] = 1/ray.dir[3]
    return nothing
end

@inline function rayHitPosition!(posHit::Vector{T}, ray::Ray{T}) where {T}
    posHit[1] = ray.origin[1] + ray.t*ray.dir[1]
    posHit[2] = ray.origin[2] + ray.t*ray.dir[2]
    posHit[3] = ray.origin[3] + ray.t*ray.dir[3]
    return nothing
end

@inline function intersect!(ray::Ray{X}, tri::GeometryBasics.Triangle{3,X}, idx::Int64) where {X}
    v1, v2, v3 = tri
    #     intersect!(ray, v1, v2, v3, idx)
    #     return
    # end

    # @inline function intersect!(ray::Ray{X}, v1, v2, v3, idx::Int64) where X
    E12 = ray.tmp1;
    E13 = ray.tmp2;
    P = ray.tmp3
    @inbounds for i in 1:3
        E12[i] = v2[i] - v1[i]
        E13[i] = v3[i] - v1[i]
    end
    cross!(P, ray.dir, E13)
    detM = dot3(P, E12)
    if abs(detM) ≤ X(EPS_PARALLEL)
        return nothing
    end

    T = ray.tmp4
    @inbounds for i in 1:3
        T[i] = ray.origin[i] - v1[i]
    end
    u = dot3(P, T)/detM
    if u < X(ZERO_WITH_TOL)
        ;
        return nothing;
    end
    if u > X(ONE_WITH_TOL)
        ;
        return nothing;
    end

    cross!(P, T, E12)       # Q <- P
    v = dot3(P, ray.dir)/detM
    if v < X(ZERO_WITH_TOL) || u + v > X(ONE_WITH_TOL)
        return nothing
    end

    t = dot3(E13, P)/detM
    if t > X(EPS_MIN_DIST) && t < ray.t
        ray.t = t
        ray.idxFace = idx
    end
end

struct BBox{T}
    # We use scalars instead of vectors to allow for a fast BVH computation
    minx::T
    miny::T
    minz::T
    maxx::T
    maxy::T
    maxz::T
end

@inline function BBox(bmin::Vector{T}, bmax::Vector{T}) where {T}
    return BBox(bmin[1], bmin[2], bmin[3], bmax[1], bmax[2], bmax[3])
end

@inline function BBox(model::M) where {M<:GeometryBasics.Mesh}
    return BBox(model.position)
end

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

@inline function intersect(ray::Ray{T}, bbox::BBox{T}) where {T}

    # X-axis
    t1 = (bbox.minx - ray.origin[1])*ray.invdir[1]
    t2 = (bbox.maxx - ray.origin[1])*ray.invdir[1]
    tmin = min(t1, t2)
    tmax = max(t1, t2)

    # Y-axis
    t1 = (bbox.miny - ray.origin[2])*ray.invdir[2]
    t2 = (bbox.maxy - ray.origin[2])*ray.invdir[2]
    tmin = max(tmin, min(t1, t2))
    tmax = min(tmax, max(t1, t2))

    # Z-axis
    t1 = (bbox.minz - ray.origin[3])*ray.invdir[3]
    t2 = (bbox.maxz - ray.origin[3])*ray.invdir[3]
    tmin = max(tmin, min(t1, t2))
    tmax = min(tmax, max(t1, t2))

    return tmax ≥ max(tmin, zero(T))
end

struct SphereModel{T}
    origin::Vector{T}
    R::T
end

@inline function intersect!(ray::Ray{T}, sphere::SphereModel{T}) where {T}
    c = zero(T);
    d = zero(T)
    @inbounds for i in 1:3
        # position of ray origin wrt sphere origin
        x = ray.origin[i] - sphere.origin[i]
        c += ray.dir[i]*x
        d += x*x
    end
    δ = c*c - d + sphere.R*sphere.R
    if δ ≤ 0
        return nothing
    end
    ρ = - c - √δ
    if 0 < ρ < ray.t
        ray.t = ρ
    end
    return nothing
end
