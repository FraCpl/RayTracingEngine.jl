# using SpacecraftBuilder
using FileIO
# using JTools
using StaticArrays
using LinearAlgebra

struct Ray{T}
    o::SVector{3,T}
    d::SVector{3,T}
    invd::SVector{3,T}
end

Ray(o::SVector{3,T}, d::SVector{3,T}) where {T} = Ray(o, normalize(d), 1 ./ normalize(d))

mutable struct Hit{T}
    t::T
    idx::Int32
end

Hit(::Type{T}) where {T} = Hit{T}(typemax(T), -1)

struct BBox{T}
    min::SVector{3,T}
    max::SVector{3,T}
end

function BBox(vertices)
    mn, mx = vertices[1], vertices[1]
    @inbounds for v in vertices
        mn = min.(mn, v)
        mx = max.(mx, v)
    end
    return BBox(SVector{3}(mn[1], mn[2], mn[3]), SVector{3}(mx[1], mx[2], mx[3]))
end

@inline function intersect(ray::Ray{T}, tri::Tuple{SVector{3,T},SVector{3,T},SVector{3,T}}, idx::Int32, hit::Hit{T}) where {T}
    v1, v2, v3 = tri
    e1 = v2 - v1
    e2 = v3 - v1

    p = cross(ray.d, e2)
    det = dot(e1, p)
    if abs(det) ≤ T(1e-8)
        return hit
    end

    invdet = 1 / det
    tvec = ray.o - v1
    u = dot(tvec, p) * invdet
    if (u < -1e-8) | (u > 1 + 1e-8)
        return hit
    end

    q = cross(tvec, e1)
    v = dot(ray.d, q) * invdet
    if (v < -1e-8) | (u + v > 1 + 1e-8)
        return hit
    end

    t = dot(e2, q) * invdet
    if (t > 1e-4) & (t < hit.t)
        hit.t = t
        hit.idx = idx
    end
    return nothing
end

@inline function intersect(ray::Ray{T}, box::BBox{T}, hit::Hit{T}) where {T}
    t1 = (box.min .- ray.o) .* ray.invd
    t2 = (box.max .- ray.o) .* ray.invd

    tmin = maximum(min.(t1, t2))
    tmax = minimum(max.(t1, t2))
    return (tmax ≥ max(tmin, zero(T))) & (tmin < hit.t)
end

function mainAero()
    anyIntersection = true
    model = load("C:/fc/data/3Dmodels/dragonLikeCapsule.obj")
    mdl = [(SVector{3}(model.position[f[1]]), SVector{3}(model.position[f[2]]), SVector{3}(model.position[f[3]])) for f in model.faces]
    bbox = BBox(model.position)

    dirSun = @SVector Float32[0.0, 0.0, 1.0]
    R_IS = @SMatrix Float32[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    X0 = (bbox.min + bbox.max)/2
    R = maximum(bbox.max - X0)

    coords = Float32.(1.1*R*range(-1, 1, round(Int, 100)))
    dirRay = -normalize(dirSun)
    Ap = (coords[2] - coords[1])^2

    hit = Hit(Float32)
    surf = 0.0f0
    idx = Int32.(1:lastindex(mdl))
    for x in coords, y in coords
        u = SVector(x, y, 3R)
        ray = Ray(X0 + R_IS*u, dirRay)
        hit.t = Inf
        if intersect(ray, bbox, hit)
            # Intersect ray with entire model (closest intersection)
            @inbounds for k in idx
                intersect(ray, mdl[k], k, hit)
                if anyIntersection && hit.t < Inf
                    ;
                    break
                end
            end
            if hit.t < Inf
                surf += Ap
            end
        end
    end
    return surf
end

mainAero()
