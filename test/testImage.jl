using RayTracingEngine
using Images
using ImageProcessing
using GLMakie
using LinearAlgebra
using Quaternions
using Statistics
using DelimitedFiles
using SparseArrays
using JTools

mutable struct RayTracingCamera{T, K}
    Nx::Int
    Ny::Int
    img::Matrix{T}
    Kinv::K

    # Allocations
    ray::Ray{T}
    rayShadow::Ray{T}
    u::Vector{T}
    dir_S::Vector{T}
end

function RayTracingCamera(Nx, Ny, Kinv::Matrix{T}) where T
    img = zeros(T, Ny, Nx)
    ray = Ray(zeros(Float32, 3), [0.0f0, 0.0f0, 1.0f0])
    rayShadow = Ray(zeros(Float32, 3), [0.0f0, 0.0f0, 1.0f0])
    u = [0.0f0, 0.0f0, 1.0f0]
    dir_S = zeros(Float32, 3)
    return RayTracingCamera(Nx, Ny, img, sparse(Kinv), ray, rayShadow, u, dir_S)
end

function imgCore(cam::RayTracingCamera{T}, bvh::BVHModel, pos_W::Vector{T}, q_WS::Vector{T}, dirSun_W::Vector{T}) where T
    cam.img .= zero(T)
    cam.ray.origin .= pos_W
    cam.rayShadow.dir .= dirSun_W
    u = cam.u; dir_S = cam.dir_S

    e12 = zeros(T, 3)
    e23 = zeros(T, 3)
    N = zeros(T, 3)

    for i in 1:cam.Nx, j in 1:cam.Ny
        # Compute pixel line of sight direction in Sensor frame
        u[1] = T(i); u[2] = T(j)
        mul!(dir_S, cam.Kinv, u)

        # Compute pixel line of sight direction in World frame
        q_transformVector!(cam.ray.dir, q_WS, dir_S)

        # Intersect with the scene
        if RayTracingEngine.intersect!(cam.ray, bvh)

            # Update shadow ray origin
            RayTracingEngine.rayHitPosition!(cam.rayShadow.origin, cam.ray)
            cam.rayShadow.origin .*= (1.0f0 + 1f-6)
            cam.rayShadow.idxSkip = cam.ray.idxFace

            # Intersect with light
            if !RayTracingEngine.intersect!(cam.rayShadow, bvh, true)

                # Compute pixel value (Lambert BRDF)
                tri = bvh.model[cam.ray.idxFace]
                v1, v2, v3 = tri
                @inbounds for k in 1:3
                    e12[k] = v2[k] - v1[k]
                    e23[k] = v3[k] - v2[k]
                end
                cross!(N, e12, e23)
                normalize!(N)
                cam.img[j, i] = N[1]*dirSun_W[1] + N[2]*dirSun_W[2] + N[3]*dirSun_W[3]
            end
        end
    end
    return cam.img
end

function mainImg()
    objFile = readdirext("C:/Users/francesco capolupo/OneDrive - ESA/Desktop/Work/Tools/ToolsJulia/Work/Argonaut/AbsNav/change3data/", ".obj"; join=true)[3]
    data = readdlm("C:/Users/francesco capolupo/OneDrive - ESA/Desktop/Work/Tools/ToolsJulia/Work/Argonaut/AbsNav/change3out/change3pose.txt")
    idx0 = findall(Int.(data[:, 1]) .== 1)[1]
    posEst_W = Float32.(data[idx0, 2:4])
    qEst_WS = Float32.(data[idx0, 5:8])
    dirSun_W  = Float32.(normalize([54718120892.43014; 124338657748.59918; 62652598859.22545]))

    N0 = 1024
    Nsup = 3
    N = Nsup*N0
    bvh = buildBvh(load(objFile))
    Kinv = Float32.(camMatrixInv(N, N, 46.6*π/180))
    cam = RayTracingCamera(N, N, Kinv)

    imgCore(cam, bvh, posEst_W, qEst_WS, dirSun_W)
    # imshow(cam.img)

    imOut = zeros(Float32, N0, N0)
    @inbounds for i in 1:N, j in 1:N
        ii = fld1(i, Nsup); jj = fld1(j, Nsup)
        imOut[ii, jj] += cam.img[i, j]
    end
    imOut ./= Nsup*Nsup
    # i1 = cam.img[1:2:end, 1:2:end]
    # i2 = cam.img[1:2:end, 2:2:end]
    # i3 = cam.img[2:2:end, 1:2:end]
    # i4 = cam.img[2:2:end, 2:2:end]
    # imOut = (i1 .+ i2 .+ i3.+ i4).*0.25
    imshow(imOut)
    return
end

mainImg()
