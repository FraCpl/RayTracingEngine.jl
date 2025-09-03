mutable struct RayTracingCamera{T}
    width::Int
    height::Int
    img::Matrix{T}
    depth::Matrix{T}
    fxinv::T
    fyinv::T
    u0inv::T
    v0inv::T

    # Allocations
    ray::Ray{T}
    rayShadow::Ray{T}
    dir_S::Vector{T}
end

function RayTracingCamera(width, height, Kcam::Matrix{T}) where T
    img = zeros(T, height, width)
    depth = zeros(T, height, width)
    ray = Ray(zeros(Float32, 3), [0.0f0, 0.0f0, 1.0f0])
    rayShadow = Ray(zeros(Float32, 3), [0.0f0, 0.0f0, 1.0f0])
    dir_S = zeros(Float32, 3)
    fx = Kcam[1, 1]; fy = Kcam[2, 2]
    u0 = Kcam[1, 3]; v0 = Kcam[2, 3]
    return RayTracingCamera(width, height, img, depth, T(1/fx), T(1/fy), -u0/fx, -v0/fy, ray, rayShadow, dir_S)
end

function rayTracingImage(
        cam::RayTracingCamera{T},       # Ray tracing camera model
        bvh::BVHModel,                  # BVH of the scene
        pos_W::Vector{T},               # Position of the sensor in World frame
        q_WS::Vector{T},                # Attitude of the sensor in World frame (zS: boresight)
        dirSun_W::Vector{T},            # Direction of the Sun in World frame (must be unit vector)
        flipNormal::Bool=false,
    ) where T

    # Initialize outputs and rays
    # zeroT = zero(T)
    cam.img .= 0
    cam.depth .= 0
    cam.ray.origin .= pos_W
    cam.rayShadow.dir .= dirSun_W
    dir_S = cam.dir_S

    e12 = cam.ray.tmp1
    e23 = cam.ray.tmp2
    N = cam.ray.tmp3
    t1 = T(1)
    SHDW_SHIFT = T(1.0 + 1e-6)

    for i in 1:cam.width, j in 1:cam.height
        # Compute pixel line of sight direction in Sensor frame
        dir_S[1] = cam.fxinv*T(i) + cam.u0inv
        dir_S[2] = cam.fyinv*T(j) + cam.v0inv
        dir_S[3] = t1
        normalize!(dir_S)

        # Compute pixel line of sight direction in World frame
        q_transformVector!(cam.ray.dir, q_WS, dir_S)

        # Intersect with the scene
        if RayTracingEngine.intersect!(cam.ray, bvh)

            # Update depth map
            cam.depth[j, i] = cam.ray.t

            # Compute normal and check for preliminary illumination condition
            tri = bvh.model[cam.ray.idxFace]
            v1, v2, v3 = tri
            @inbounds for k in 1:3
                e12[k] = v2[k] - v1[k]
                e23[k] = v3[k] - v2[k]
            end
            cross!(N, e12, e23)     # Caution, non-unitary norm here!
            if flipNormal; N .*= -t1; end
            if dot3(N, dirSun_W) < 0; continue; end     # Check: surface normal and light direction shall lie in the same hemisphere
			if dot3(N, cam.ray.dir) > 0; continue; end  # Check: surface normal and observer direction shall lie in the same hemisphere

            # Update shadow ray origin
            RayTracingEngine.rayHitPosition!(cam.rayShadow.origin, cam.ray)
            cam.rayShadow.origin .*= SHDW_SHIFT         # Shift to avoid self-intersection
            cam.rayShadow.idxSkip = cam.ray.idxFace     # Avoid self-intersection

            # Intersect with light
            if !RayTracingEngine.intersect!(cam.rayShadow, bvh, true)
                # Normalize surface normal
                normalize!(N)

                # Compute pixel value (Lambert BRDF)
                cam.img[j, i] = dot3(N, dirSun_W)

                # cosi = dot3(dirSun_W, N)
                # cosr = -dot3(cam.ray.dir, N)
                # α = acos(clamp(-dot3(cam.ray.dir, dirSun_W), -t1, t1))
                # β = exp(-α/1.0473)
                # cam.img[j, i] = T((1.0 - β + β/(cosi + cosr))*cosi)
            end
        end
    end
    return cam.img
end

@inline function dot3(a::Vector{T}, b::Vector{T}) where T
    return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
end
