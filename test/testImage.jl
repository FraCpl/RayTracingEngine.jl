using RayTracingEngine
using ImageProcessing
using LinearAlgebra
using FileIO
using DelimitedFiles
using JTools

function mainImg()
    objFile = readdirext("C:/Users/francesco capolupo/OneDrive - ESA/Desktop/Work/Tools/ToolsJulia/Work/Argonaut/AbsNav/change3data/", ".obj"; join=true)[3]
    data = readdlm("C:/Users/francesco capolupo/OneDrive - ESA/Desktop/Work/Tools/ToolsJulia/Work/Argonaut/AbsNav/change3out/change3pose.txt")
    idx0 = findall(Int.(data[:, 1]) .== 1)[1]
    posEst_W = Float32.(data[idx0, 2:4])
    qEst_WS = Float32.(data[idx0, 5:8])
    dirSun_W  = Float32.(normalize([54718120892.43014; 124338657748.59918; 62652598859.22545]))

    N0 = 1024
    Nsup = 5
    N = Nsup*N0
    bvh = buildBvh(load(objFile))
    Kcam = Float32.(camMatrix(N, N, 46.6*π/180))
    cam = RayTracingCamera(N, N, Kcam)

    rayTracingImage(cam, bvh, posEst_W, qEst_WS, dirSun_W)
    imOut = imageBinning(cam.img, Nsup)

    K = ImageProcessing.gaussianKernel1D(2.0f0, 5)
    imOut = imageFilter(imOut, K, K)
    imshow(imOut)
    return
end

mainImg()
